/*
  Copyright (c) 2014 - 2019 University of Bergen
  
  This file is part of the BROOMStyx project.

  BROOMStyx is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  BROOMStyx is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with BROOMStyx.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the AUTHORS file
  for the list of copyright holders.
*/

#include "AlternateNonlinearMinimization.hpp"
#include <chrono>
#include <cmath>
#include <stdexcept>
#include <string>
#include <tuple>
#include <omp.h>

#include "Core/AnalysisModel.hpp"
#include "Core/ObjectFactory.hpp"
#include "Core/Diagnostics.hpp"
#include "Core/DofManager.hpp"
#include "Core/DomainManager.hpp"
#include "Core/LoadStep.hpp"
#include "Core/NumericsManager.hpp"
#include "MeshReaders/MeshReader.hpp"
#include "Core/SolutionManager.hpp"
#include "LinearSolvers/LinearSolver.hpp"
#include "Numerics/Numerics.hpp"
#include "SparseMatrix/SparseMatrix.hpp"
#include "Util/linearAlgebra.hpp"
#include "Util/readOperations.hpp"

#define REPORT_SUBSYSTEM_STATUS false

#define WRITE_DEBUG_OUTPUT false
#define DEBUG_SUBSTEP 1
#define DEBUG_SUBSYS 1

using namespace broomstyx;

registerBroomstyxObject(SolutionMethod, AlternateNonlinearMinimization)

AlternateNonlinearMinimization::AlternateNonlinearMinimization()
{
    _name = "AlternateNonlinearMinimization";
}
// ---------------------------------------------------------------------------
AlternateNonlinearMinimization::~AlternateNonlinearMinimization() {}

// Public methods
// ---------------------------------------------------------------------------
int AlternateNonlinearMinimization::computeSolutionFor
    ( int stage
    , const std::vector<BoundaryCondition>& bndCond
    , const std::vector<FieldCondition>& fldCond
    , const TimeData& time )
{
    std::chrono::time_point<std::chrono::system_clock> tic, toc, innertic, innertoc;
    std::chrono::duration<double> tictoc;
    
    _substepCount += 1;
    
    // Get all active DOFs at current stage
    std::vector<Dof*> dof = analysisModel().dofManager().giveActiveDofsAtStage(stage);
    
    // Impose constraints on nodal DOFs
    std::printf("    %-40s", "Imposing constraints ...");
    std::fflush(stdout);
    tic = std::chrono::high_resolution_clock::now();
    this->imposeConstraintsAt(stage, bndCond, time);
    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    std::printf("done (time = %f sec.)\n", tictoc.count());

    // Initialize global vectors
    tic = std::chrono::high_resolution_clock::now();
    std::vector<RealVector> dU(_nSubsystems, RealVector());
    auto resid = dU;
    auto rhs = dU;
    auto lhs = dU;
    auto sum_dU = dU;
    
    for ( int i = 0; i < _nSubsystems; i++ )
    {
        dU[i].init(_nUnknowns[i]);
        sum_dU[i].init(_nUnknowns[i]);
    }

    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    diagnostics().addSetupTime(tictoc.count());
    
    // Start iterations
    tic = std::chrono::high_resolution_clock::now();
    bool converged = false;
    int iterCount = 0;
    while ( !converged && iterCount <= _maxIter )
    {
        // Perform pre-iteration operations, if any ...
        innertic = std::chrono::high_resolution_clock::now();
        std::vector<Numerics*> numVec = analysisModel().numericsManager().giveAllNumerics();
        for ( Numerics* curNumerics : numVec )
            curNumerics->performPreIterationOperationsAt(stage, iterCount);
        innertoc = std::chrono::high_resolution_clock::now();
        tictoc = innertoc - innertic;
        diagnostics().addPostprocessingTime(tictoc.count());
        
        // Reinitialize data for calculation of residual convergence
        for ( int i = 0; i < _nDofGroups; i++ )
        {
            _sumAbsFlux(i) = 0;
            _fluxCount(i) = 0;
        }
        
        analysisModel().dofManager().resetSecondaryVariablesAtStage(stage);
        
        // Calculate residual for each subsystem
        // Note: We calculate global external force vectors in each iteration 
        // to accommodate numerics classes that implement special methods for imposing constraints        
        for ( int i = 0; i < _nSubsystems; i++)
        {
            rhs[i] = assembleRightHandSide(stage, _subsysNum[i], bndCond, fldCond, time);
            lhs[i] = assembleLeftHandSide(stage, _subsysNum[i], time);
            resid[i] = rhs[i] - lhs[i];
        }

        // Store residuals in DOFs
        std::vector<Dof*> activeDof = analysisModel().dofManager().giveActiveDofsAtStage(stage);

#pragma omp parallel for
        for ( int i = 0; i < (int)activeDof.size(); i++ )
        {
            int ssNum = analysisModel().dofManager().giveSubsystemNumberFor(activeDof[i]);
            int eqNum = analysisModel().dofManager().giveEquationNumberAt(activeDof[i]);

            for ( int j = 0; j < _nSubsystems; j++ )
                if ( ssNum == _subsysNum[j] )
                    analysisModel().dofManager().updateResidualAt(activeDof[i], -resid[j](eqNum));
        }
        
        if ( iterCount == 0 )
        {
            std::printf("\n    Initial guess (LS # %d, SS # %d)", _loadStep->giveLoadStepNum(), _substepCount);
            std::printf("\n   -----------------");
        }
        
        // Compute convergence norms
        RealMatrix normDat;
        
        // Disallow convergence at initial guess
        if ( iterCount == 0 )
            converged = false;
        else
        {
            std::tie(converged,normDat) = this->computeConvergenceNormsFrom(resid, dof);
            this->reportConvergenceStatus(normDat);
        }
        
        // Additional convergence checks from numerics
        bool numericsConverged = this->checkConvergenceOfNumericsAt(stage);
        if ( !numericsConverged )
            converged = false;
        
        
        if ( iterCount > 0 )
            _loadStep->writeConvergenceDataForStage(stage, normDat);
        
        toc = std::chrono::high_resolution_clock::now();
        tictoc = toc - tic;
        std::printf("\n\n    Total time for iteration = %f sec.\n", tictoc.count());
        
        // Print post iteration messages from each numerics
        if ( iterCount >= 0 )
            for ( Numerics* curNumerics : numVec )
                curNumerics->printPostIterationMessage(stage);
        
        if ( converged )
        {
            // Clear memory for solvers
            for ( int i = 0; i < _nSubsystems; i++ )
                _solver[i]->clearInternalMemory();
            
            innertic = std::chrono::high_resolution_clock::now();
            _loadStep->writeIterationDataForStage(stage, time.target, iterCount);
            innertoc = std::chrono::high_resolution_clock::now();
            tictoc = innertoc - innertic;
            diagnostics().addOutputWriteTime(tictoc.count());
        }
        else if ( iterCount == _maxIter )
        {
            // Clear memory for solvers
            for ( int i = 0; i < _nSubsystems; i++ )
                _solver[i]->clearInternalMemory();
            
            innertic = std::chrono::high_resolution_clock::now();
            _loadStep->writeIterationDataForStage(stage, time.target, iterCount);
            innertoc = std::chrono::high_resolution_clock::now();
            tictoc = innertoc - innertic;
            diagnostics().addOutputWriteTime(tictoc.count());
            ++iterCount;
        }
        else
        {
            // Start of new iteration
            ++iterCount;
            
            tic = std::chrono::high_resolution_clock::now();
            
            std::printf("\n    Iter.No. %-4d (LS # %d, SS # %d)", iterCount, _loadStep->giveLoadStepNum(), _substepCount);
            std::printf("\n   -----------------\n");
            
            // Cycle over all subsystems
            for ( int curSubsys = 0; curSubsys < _nSubsystems; curSubsys++ )
            {
                std::printf("    ");
                // Cumulative values of corrections for each subsystem
                sum_dU[curSubsys].init(_nUnknowns[curSubsys]);
                
                bool subsysConverged = false;
                int localIterCount = 0;
                while ( !subsysConverged && localIterCount < _maxSubsysIter[curSubsys] )
                {
                    // Reinitialize data for calculation of residual convergence
                    for ( int i = 0; i < _nDofGroups; i++ )
                    {
                        _sumAbsFlux(i) = 0;
                        _fluxCount(i) = 0;
                    }
                    
                    // Perform pre-iteration operations, if any ...
                    std::vector<Numerics*> numVec = analysisModel().numericsManager().giveAllNumerics();
                    for ( Numerics* curNumerics : numVec )
                        curNumerics->performPreIterationOperationsAt(stage, iterCount);
                    
                    // Recalculate right and left hand side vectors
                    rhs[curSubsys] = this->assembleRightHandSide(stage, _subsysNum[curSubsys], bndCond, fldCond, time);
                    lhs[curSubsys] = this->assembleLeftHandSide(stage, _subsysNum[curSubsys], time);
                    resid[curSubsys] = rhs[curSubsys] - lhs[curSubsys];
                    
                    // Check subsystem convergence
                    subsysConverged = this->checkSubsystemConvergenceAt(stage, curSubsys, resid[curSubsys], iterCount, localIterCount);
                    
                    if ( localIterCount == 0 )
                        subsysConverged = false;
                    
                    if ( !subsysConverged )
                    {
                        ++localIterCount;
                        if ( !REPORT_SUBSYSTEM_STATUS )
                        {
                            std::printf(".");
                            std::fflush(stdout);
                        }
                        
                        // Assemble system Jacobian
                        _spMatrix[curSubsys]->initializeValues();
                        this->assembleJacobian(stage, _subsysNum[curSubsys], time);

                        // Solve system
                        innertic = std::chrono::high_resolution_clock::now();
                        
                        _solver[curSubsys]->allocateInternalMemoryFor(_spMatrix[curSubsys]);
                        dU[curSubsys] = _solver[curSubsys]->solve(_spMatrix[curSubsys], resid[curSubsys]);
                        
                        innertoc = std::chrono::high_resolution_clock::now();
                        tictoc = innertoc - innertic;
                        diagnostics().addSolveTime(tictoc.count());
                        
                        // Execute line search
                        if ( _enableLineSearch[curSubsys] )
                        {
                            int lsIter = 0;
                            
                            // Compute initial slope
                            double s0 = dU[curSubsys].dot(resid[curSubsys]);
                            if ( s0 < 0 )
                                std::printf("  WARNING: Detected non-descent direction!\n");
                            
                            double sL = s0;
                            double etaL = 0;
                            double etaH = _initMultiplier(curSubsys);
                            
                            // Compute slope at initial value of eta
                            innertic = std::chrono::high_resolution_clock::now();
#ifdef _OPENMP
#pragma omp parallel for
#endif
                            for ( int i = 0; i < (int)dof.size(); i++ )
                                if ( analysisModel().dofManager().giveSubsystemNumberFor(dof[i]) == _subsysNum[curSubsys] )
                                {
                                    int eqNo = analysisModel().dofManager().giveEquationNumberAt(dof[i]);
                                    if ( eqNo != UNASSIGNED )
                                        analysisModel().dofManager().updatePrimaryVariableAt(dof[i], etaH*dU[curSubsys](eqNo), correction );
                                }
                            
                            innertoc = std::chrono::high_resolution_clock::now();
                            tictoc = innertoc - innertic;
                            diagnostics().addUpdateTime(tictoc.count());
                            
                            // Calculate initial residual
                            // Perform pre-iteration operations, if any ...
                            std::vector<Numerics*> numVec = analysisModel().numericsManager().giveAllNumerics();
                            for ( Numerics* curNumerics : numVec )
                                curNumerics->performPreIterationOperationsAt(stage, iterCount);
                    
                            RealVector lsRhs = this->assembleRightHandSide(stage, _subsysNum[curSubsys], bndCond, fldCond, time);
                            RealVector lsLhs = this->assembleLeftHandSide(stage, _subsysNum[curSubsys], time);
                            RealVector lsResid = lsRhs - lsLhs;
                            double sH = dU[curSubsys].dot(lsResid);
                            std::printf("\n\n        s0 = %E, sH = %E, sH/s0 = %E\n", s0, sH, sH/s0);
                            
                            bool lineSearchConverged;
                            if ( std::fabs(sH/s0) < _slackTolerance(curSubsys) )
                                lineSearchConverged = true;
                            else
                                lineSearchConverged = false;
                            
                            double eta = etaH;
                            while ( !lineSearchConverged && lsIter < _maxLineSearchIter[curSubsys] )
                            {
                                ++lsIter;
                                
                                // Implement regula falsi line search
                                eta = etaH - sH*(etaL - etaH)/(sL - sH);
                                if ( eta < _minMultiplier(curSubsys) )
                                {
                                    eta = _minMultiplier(curSubsys);
                                    lineSearchConverged = true;
                                }
                                if ( eta > _maxMultiplier(curSubsys) )
                                {
                                    eta = _maxMultiplier(curSubsys);
                                    lineSearchConverged = true;
                                }
                                
                                // Recalculate residual for current value of eta
                                innertic = std::chrono::high_resolution_clock::now();
#ifdef _OPENMP
#pragma omp parallel for
#endif
                                for ( int i = 0; i < (int)dof.size(); i++ )
                                    if ( analysisModel().dofManager().giveSubsystemNumberFor(dof[i]) == _subsysNum[curSubsys] )
                                    {
                                        int eqNo = analysisModel().dofManager().giveEquationNumberAt(dof[i]);
                                        if ( eqNo != UNASSIGNED )
                                            analysisModel().dofManager().updatePrimaryVariableAt(dof[i], eta*dU[curSubsys](eqNo), replacement_correction );
                                    }

                                innertoc = std::chrono::high_resolution_clock::now();
                                tictoc = innertoc - innertic;
                                diagnostics().addUpdateTime(tictoc.count());
                                
                                // Perform pre-iteration operations, if any ...
                                innertic = std::chrono::high_resolution_clock::now();
                                std::vector<Numerics*> numVec = analysisModel().numericsManager().giveAllNumerics();
                                for ( Numerics* curNumerics : numVec )
                                    curNumerics->performPreIterationOperationsAt(stage, iterCount);

                                innertoc = std::chrono::high_resolution_clock::now();
                                tictoc = innertoc - innertic;
                                diagnostics().addPostprocessingTime(tictoc.count());

                                lsRhs = this->assembleRightHandSide(stage, _subsysNum[curSubsys], bndCond, fldCond, time);
                                lsLhs = this->assembleLeftHandSide(stage, _subsysNum[curSubsys], time);
                                lsResid = lsRhs - lsLhs;
                                double sEta = dU[curSubsys].dot(lsResid);
                                
                                if ( sEta*sL < 0 )
                                {
                                    etaH = eta;
                                    sH = sEta;
                                }
                                else
                                {
                                    etaL = eta;
                                    sL = sEta;
                                }
                                
                                // Output status
                                double sRatio = std::fabs(sEta/s0);
                                std::printf("        LS # %d: eta = %E, s/s0 = %E\n", lsIter, eta, sRatio);
                                
                                if ( sRatio < _slackTolerance(curSubsys) )
                                    lineSearchConverged = true;
                            }

                            innertic = std::chrono::high_resolution_clock::now();
#ifdef _OPENMP
#pragma omp parallel for
#endif
                            for ( int i = 0; i < (int)dof.size(); i++ )
                                if ( analysisModel().dofManager().giveSubsystemNumberFor(dof[i]) == _subsysNum[curSubsys] )
                                {
                                    int eqNo = analysisModel().dofManager().giveEquationNumberAt(dof[i]);
                                    if ( eqNo != UNASSIGNED )
                                    {
                                        sum_dU[curSubsys](eqNo) += eta*dU[curSubsys](eqNo);
                                    
                                        // Update secondary variable at DOF with residual
                                        analysisModel().dofManager().updateResidualAt(dof[i], resid[curSubsys](eqNo));
                                    }
                                }

                            innertoc = std::chrono::high_resolution_clock::now();
                            tictoc = innertoc - innertic;
                            diagnostics().addUpdateTime(tictoc.count());
                        }
                        else
                        {
                            innertic = std::chrono::high_resolution_clock::now();
#ifdef _OPENMP
#pragma omp parallel for
#endif
                            for ( int i = 0; i < (int)dof.size(); i++ )
                                if ( analysisModel().dofManager().giveSubsystemNumberFor(dof[i]) == _subsysNum[curSubsys] )
                                {
                                    int eqNo = analysisModel().dofManager().giveEquationNumberAt(dof[i]);
                                    if ( eqNo != UNASSIGNED )
                                    {
                                        // Update primary variable at DOF
                                        analysisModel().dofManager().updatePrimaryVariableAt(dof[i], _initMultiplier(curSubsys)*dU[curSubsys](eqNo), correction);
                                        sum_dU[curSubsys](eqNo) += _initMultiplier(curSubsys)*dU[curSubsys](eqNo);
                                        
                                        // Update secondary variable at DOF with residual
                                        analysisModel().dofManager().updateResidualAt(dof[i], resid[curSubsys](eqNo));
                                    }
                                }

                            innertoc = std::chrono::high_resolution_clock::now();
                            tictoc = innertoc - innertic;
                            diagnostics().addUpdateTime(tictoc.count());
                        }
                    }
                }

                // Apply over-relaxation (also ensures that the cumulative correction is accounted for and not
                // only the final correction of the subsystem
                innertic = std::chrono::high_resolution_clock::now();
#ifdef _OPENMP
#pragma omp parallel for
#endif
                for ( int i = 0; i < (int)dof.size(); i++ )
                    if ( analysisModel().dofManager().giveSubsystemNumberFor(dof[i]) == _subsysNum[curSubsys] )
                    {
                        int eqNo = analysisModel().dofManager().giveEquationNumberAt(dof[i]);
                        if ( eqNo != UNASSIGNED )
                        {
                            analysisModel().dofManager().updatePrimaryVariableAt(dof[i], -sum_dU[curSubsys](eqNo), correction);
                            analysisModel().dofManager().updatePrimaryVariableAt(dof[i], sum_dU[curSubsys](eqNo), correction);
                        }   
                    }

                innertoc = std::chrono::high_resolution_clock::now();
                tictoc = innertoc - innertic;
                diagnostics().addUpdateTime(tictoc.count());

                std::printf("\n\n    Subsys # %d: Inner iterations performed = %d\n", _subsysNum[curSubsys], localIterCount);
                std::fflush(stdout);
            }
        }
    }
    
    std::printf("\n");
    if ( !converged )
        return 1;
    else
        return 0;
}
// ---------------------------------------------------------------------------
void AlternateNonlinearMinimization::readDataFromFile( FILE* fp )
{    
    std::string key, errmsg, src = "AlternateNonlinearMinimization (SolutionMethod)";
    
    // Read number of DOF groups
    verifyKeyword(fp, "DofGroups", _name);
    _nDofGroups = getIntegerInputFrom(fp, "Failed to read number of DOF groups from input file!", _name);
    _dofGrpNum.assign(_nDofGroups, 0);
    
    // Initialize solution control vectors
    _relTolCor.init(_nDofGroups);
    _relTolRes.init(_nDofGroups);
    _absTolCor.init(_nDofGroups);
    _absTolRes.init(_nDofGroups);
    _sumAbsFlux.init(_nDofGroups);
    _fluxCount.init(_nDofGroups);
    
    for ( int i = 0; i < _nDofGroups; i++ )
    {
        // Read DOF group number
        _dofGrpNum[i] = getIntegerInputFrom(fp, "Failed to read DOF group number from input file!", _name);

        // Read DOF group convergence parameters
        verifyKeyword(fp, key = "Parameters", _name);

        // Correction tolerance
        _relTolCor(i) = getRealInputFrom(fp, "Failed to read correction tolerance from input file!", _name);

        // Residual tolerance
        _relTolRes(i) = getRealInputFrom(fp, "Failed to read residual tolerance from input file!", _name);

        // Absolute tolerance for corrections
        _absTolCor(i) = getRealInputFrom(fp, "Failed to read absolute correction tolerance from input file!", _name);

        // Absolute tolerance for residuals
        _absTolRes(i) = getRealInputFrom(fp, "Failed to read absolute residual tolerance from input file!", _name);
    }
    
    // Read number of subsystems
    verifyKeyword(fp, "Subsystems", _name);
    _nSubsystems = getIntegerInputFrom(fp, "Failed reading number of subsystems from input file!", _name);

    _subsysNum.assign(_nSubsystems, 0);
    _subsysDofGroup.assign(_nSubsystems, std::vector<int>());
    _symmetry.assign(_nSubsystems, false);
    _solver.assign(_nSubsystems, nullptr);
    _spMatrix.assign(_nSubsystems, nullptr);
    _nUnknowns.assign(_nSubsystems, 0);
    _maxSubsysIter.assign(_nSubsystems, 0);
    _initMultiplier.init(_nSubsystems);
    _tolFactor.init(_nSubsystems, 2);
    _enableLineSearch.assign(_nSubsystems, false);
    _slackTolerance.init(_nSubsystems);
    _minMultiplier.init(_nSubsystems);
    _maxMultiplier.init(_nSubsystems);
    _maxLineSearchIter.assign(_nSubsystems, 0);
    
    for ( int i = 0; i < _nSubsystems; i++ )
    {
        // Read subsystem number
        _subsysNum[i] = getIntegerInputFrom(fp, "Failed reading subsystem number from input file!", _name);
        
        // Read which DOF groups are part of subsystem
        verifyKeyword(fp, "nDofGroups", _name);
        int nDofGrp = getIntegerInputFrom(fp, "Failed to read number of subsystem DOF groups from input file!", _name);
        
        _subsysDofGroup[i].assign(nDofGrp, UNASSIGNED);
        verifyKeyword(fp, "Label", _name);
        for ( int j = 0; j < nDofGrp; j++ )
            _subsysDofGroup[i][j] = getIntegerInputFrom(fp, "Failed to read subsystem DOF group from input file!", _name);
            
        // Read linear solver to be used
        verifyKeyword(fp, "LinearSolver", _name);
        key = getStringInputFrom(fp, "Failed to read linear solver from input file!", _name);
        _solver[i] = objectFactory().instantiateLinearSolver(key);
        _solver[i]->readDataFrom(fp);
        _symmetry[i] = _solver[i]->giveSymmetryOption();
                
        // Maximum number of iterations in each subsystem
        verifyKeyword(fp, key = "MaxSubsystemIterations", _name);
        _maxSubsysIter[i] = getIntegerInputFrom(fp, "Failed to read maximum number of subsystem iterations from input file!", _name);
        
        // Subsystem overrelaxation
        verifyKeyword(fp, key = "InitMultiplier", _name);
        _initMultiplier(i) = getRealInputFrom(fp, "Failed to read subsystem damping factor from input file!", _name);
        
        // Tolerance reduction factors
        verifyKeyword(fp, key = "TolReduction", _name);
        _tolFactor(i,0) = getRealInputFrom(fp, "Failed to read subsystem tolerance reduction factor for corrections from input file!", _name);
        _tolFactor(i,1) = getRealInputFrom(fp, "Failed to read subsystem tolerance reduction factor for residuals from input file!", _name);
        
        // Line search
        verifyKeyword(fp, "LineSearch", this->_name);
        key = getStringInputFrom(fp, "Failed to reach line search option from input file!", _name);
        if ( key == "On" )
        {
            _enableLineSearch[i] = true;
            
        }
        else if ( key == "Off" )
            _enableLineSearch[i] = false;
        else
        {
            std::printf("Error in line search option for solution method '%s': < On/Off > expected, '%s' found.", _name.c_str(), key.c_str());
            throw std::runtime_error("");
        }

        if ( _enableLineSearch[i] )
        {
            verifyKeyword(fp, "SlackTolerance", _name);
            _slackTolerance(i) = getRealInputFrom(fp, "Failed to read line search slack tolerance from input file!", _name);

            verifyKeyword(fp, "MinMultiplier", _name);
            _minMultiplier(i) = getRealInputFrom(fp, "Failed to read minimum line search multiplier from input file", _name);
            
            verifyKeyword(fp, "MaxMultiplier", _name);
            _maxMultiplier(i) = getRealInputFrom(fp, "Failed to read maximum line search multiplier from input file", _name);

            verifyKeyword(fp, "MaxLineSearchIterations", _name);
            _maxLineSearchIter[i] = getRealInputFrom(fp, "Failed to read maximum number of line search iterations from input file!", _name);
        }
        
        // Create sparse matrices
        _spMatrix[i] = objectFactory().instantiateSparseMatrix(_solver[i]->giveRequiredMatrixFormat());
        _spMatrix[i]->setSymmetryTo(_symmetry[i]);
    }

    // Maximum number of iterations
    verifyKeyword(fp, key = "MaxIterations", _name);
    _maxIter = getIntegerInputFrom(fp, "Failed to read maximum number of iterations from input file!", _name);
}

// Private methods
// --------------------------------------------------------------------------
bool AlternateNonlinearMinimization::checkSubsystemConvergenceAt( int stage
                                                                , int subsysIdx
                                                                , const RealVector& resid
                                                                , int globalIterCount
                                                                , int localIterCount )
{
    std::chrono::time_point<std::chrono::system_clock> tic, toc;
    std::chrono::duration<double> tictoc;

    tic = std::chrono::high_resolution_clock::now();

    // Get all active DOFs at current stage
    std::vector<Dof*> dof = analysisModel().dofManager().giveActiveDofsAtStage(stage);
    
    // Build subsystem vectors of DOFs
    int nSubsysDof = _nUnknowns[subsysIdx];
    std::vector<Dof*> subsysDof(nSubsysDof, nullptr);
    
    for ( int i = 0; i < (int)dof.size(); i++)
    {
        int ssNum = analysisModel().dofManager().giveSubsystemNumberFor(dof[i]);
        int eqNum = analysisModel().dofManager().giveEquationNumberAt(dof[i]);
        if ( ssNum == _subsysNum[subsysIdx] )
            subsysDof[eqNum] = dof[i];
    }
    
    // Initialize matrix to store convergence data
    //   col 0: norm
    //   col 1: criterion
    // -----------------------
    
    int nSubsysDofGroups = (_subsysDofGroup[subsysIdx]).size();
    RealMatrix normDat(2*nSubsysDofGroups, 2);

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < (int)subsysDof.size(); i++ )
    {
        // Find group number of DOF
        int grpNum = analysisModel().dofManager().giveGroupNumberFor(subsysDof[i]);
        
        int ssGrpIdx = 0;
        for ( int j = 0; j < nSubsysDofGroups; j++ )
            if ( _subsysDofGroup[subsysIdx][j] == grpNum )
                ssGrpIdx = j;

        // Contribution to L2-norm of corrections
        double corrVal = analysisModel().dofManager().giveValueOfPrimaryVariableAt(subsysDof[i], correction);

#ifdef _OPENMP
#pragma omp atomic
#endif
        normDat(2*ssGrpIdx,0) += corrVal*corrVal;

        // Contribution to L2-norm of incremental solution
        double incVal = analysisModel().dofManager().giveValueOfPrimaryVariableAt(subsysDof[i], incremental_value);

#ifdef _OPENMP
#pragma omp atomic
#endif
        normDat(2*ssGrpIdx,1) += incVal*incVal;

        // Contribution to L2-norm of residual
        double rVal = resid(i);

#ifdef _OPENMP
#pragma omp atomic
#endif
        normDat(2*ssGrpIdx+1,0) += rVal*rVal;
    }
        
    for ( int i = 0; i < nSubsysDofGroups; i++ )
    {
        // Normalize entries
        int grpIdx = this->giveIndexForDofGroup(_subsysDofGroup[subsysIdx][i]);
        normDat(2*i,0) = std::sqrt(normDat(2*i,0)/_dofGrpCount(grpIdx));
        normDat(2*i+1,0) = std::sqrt(normDat(2*i+1,0)/_dofGrpCount(grpIdx));
        normDat(2*i,1) = std::sqrt(normDat(2*i,1)/_dofGrpCount(grpIdx));
        
        // Calculate criteria using relative tolerances
        normDat(2*i,1) *= _relTolCor(grpIdx)*_tolFactor(subsysIdx,0);
        normDat(2*i+1,1) = _relTolRes(grpIdx)*_sumAbsFlux(grpIdx)/_fluxCount(grpIdx)*_tolFactor(subsysIdx,1);
        
        // Check against absolute tolerances
        if ( normDat(2*i,1) < _absTolCor(grpIdx) )
            normDat(2*i,1) = _absTolCor(grpIdx);
        if ( normDat(2*i+1,1) < _absTolRes(grpIdx) )
            normDat(2*i+1,1) = _absTolRes(grpIdx);
    }
    
    bool isConverged = true;
    for ( int i = 0; i < nSubsysDofGroups; i++ )
    {
        if ( normDat(2*i,0) > normDat(2*i,1) )
            isConverged = false;
        if ( normDat(2*i+1,0) > normDat(2*i+1,1) )
            isConverged = false;
    }
    
    if ( REPORT_SUBSYSTEM_STATUS )
    {
        std::printf("\n\n        Global Iter. %d, Subsystem %d, Local Iter. %d", globalIterCount, _subsysNum[subsysIdx], localIterCount);    
        std::printf("\n        DOF Grp   L2-Norm         Criterion");
        std::printf("\n       -------------------------------------------------------");

        for (int j = 0; j < nSubsysDofGroups; j++)
        {
            int grpNum = _subsysDofGroup[subsysIdx][j];
            // Convergence of corrections
            std::printf("\n         C  %-6d%e   %e ", grpNum, normDat(2*j,0), normDat(2*j,1));

            if ( normDat(2*j,0) <= normDat(2*j,1) )
                std::printf(" << CONVERGED");

            // Convergence of residuals
            std::printf("\n         R  %-6d%e   %e ", grpNum, normDat(2*j+1,0), normDat(2*j+1,1));

            if ( normDat(2*j+1,0) <= normDat(2*j+1,1) )
                std::printf(" << CONVERGED");
        }
        std::fflush(stdout);
    }

    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    diagnostics().addConvergenceCheckTime(tictoc.count());

    return isConverged;
}
