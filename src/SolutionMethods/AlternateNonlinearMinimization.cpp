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
#include "SolutionMethods/ConvergenceCriteria/ConvergenceCriterion.hpp"
#include "Util/linearAlgebra.hpp"
#include "Util/readOperations.hpp"

#define REPORT_SUBSYSTEM_STATUS false

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
        std::vector<Numerics*> numerics = analysisModel().numericsManager().giveAllNumerics();
        for ( Numerics* curNumerics : numerics )
            curNumerics->performPreIterationOperationsAt(stage, iterCount);
        innertoc = std::chrono::high_resolution_clock::now();
        tictoc = innertoc - innertic;
        diagnostics().addPostprocessingTime(tictoc.count());
        
        analysisModel().dofManager().resetSecondaryVariablesAtStage(stage);
        for ( int i = 0; i < _nDofGroups; i++ )
            _convergenceCriterion[i]->resetResidualCriteria();
        
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
        RealMatrix normDat(2*_nDofGroups, 2);
        
        // Disallow convergence at initial guess
        if ( iterCount == 0 )
            converged = false;
        else
        {
            converged = true;
            for ( int i = 0; i < _nSubsystems; i++ )
                for ( int j = 0; j < (int)_subsysDofGroup[i].size(); j++ )
                {
                    int idx = this->giveIndexForDofGroup(_subsysDofGroup[i][j]);
                    bool dofGrpConverged = _convergenceCriterion[idx]->checkConvergenceOf(resid[i], dof);
                    if ( !dofGrpConverged )
                        converged = false;
                }

            for ( int i = 0; i < _nDofGroups; i++ )
            {
                RealMatrix dofGrpNormDat = _convergenceCriterion[i]->giveConvergenceData();
                normDat(2*i, 0) = dofGrpNormDat(0, 0);
                normDat(2*i, 1) = dofGrpNormDat(0, 1);
                normDat(2*i+1, 0) = dofGrpNormDat(1, 0);
                normDat(2*i+1, 1) = dofGrpNormDat(1, 1);
            }
            
            // Report convergence status
            std::printf("\n\n    DOF Grp      Norm          Criterion");
            std::printf("\n   -------------------------------------------------------");

            for ( int i = 0; i < _nDofGroups; i++ )
                _convergenceCriterion[i]->reportConvergenceStatus();
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
            for ( Numerics* curNumerics : numerics )
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
                        _convergenceCriterion[i]->resetResidualCriteria();
                    
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
                        if ( _solver[curSubsys]->takesInitialGuess() )
                            _solver[curSubsys]->setInitialGuessTo(dU[curSubsys]);
                        dU[curSubsys] = _solver[curSubsys]->solve(_spMatrix[curSubsys], resid[curSubsys]);
                        
                        innertoc = std::chrono::high_resolution_clock::now();
                        tictoc = innertoc - innertic;
                        diagnostics().addSolveTime(tictoc.count());
                        

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
                                    analysisModel().dofManager().updatePrimaryVariableAt(dof[i], dU[curSubsys](eqNo), correction);
                                    sum_dU[curSubsys](eqNo) += dU[curSubsys](eqNo);
                                    
                                    // Update secondary variable at DOF with residual
                                    analysisModel().dofManager().updateResidualAt(dof[i], resid[curSubsys](eqNo));
                                }
                            }

                        innertoc = std::chrono::high_resolution_clock::now();
                        tictoc = innertoc - innertic;
                        diagnostics().addUpdateTime(tictoc.count());
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
    if ( !converged && _abortAtMaxIter )
        return 1;
    else
        return 0;
}
// ---------------------------------------------------------------------------
void AlternateNonlinearMinimization::readDataFromFile( FILE* fp )
{    
    // Read number of DOF groups
    verifyKeyword(fp, "DofGroups", _name);
    _nDofGroups = getIntegerInputFrom(fp, "Failed to read number of DOF groups from input file!", _name);
    _dofGrpNum.assign(_nDofGroups, -1);
    
    // Read convergence criterion and associated data
    _convergenceCriterion.assign(_nDofGroups, nullptr);
    for ( int i = 0; i < _nDofGroups; i++ )
    {
        _dofGrpNum[i] = getIntegerInputFrom(fp, "Failed to read DOF group number from input file!", _name);
        std::string convCritName = getStringInputFrom(fp, "Failed to read convergence criterion type from input file!", _name);
        _convergenceCriterion[i] = objectFactory().instantiateConvergenceCriterion(convCritName);
        _convergenceCriterion[i]->initialize(_dofGrpNum[i]);
        _convergenceCriterion[i]->readDataFromFile(fp);
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
        std::string linSolver = getStringInputFrom(fp, "Failed to read linear solver from input file!", _name);
        _solver[i] = objectFactory().instantiateLinearSolver(linSolver);
        _solver[i]->readDataFrom(fp);
        _symmetry[i] = _solver[i]->giveSymmetryOption();
                
        // Maximum number of iterations in each subsystem
        verifyKeyword(fp, "MaxSubsystemIterations", _name);
        _maxSubsysIter[i] = getIntegerInputFrom(fp, "Failed to read maximum number of subsystem iterations from input file!", _name);
        
        // Create sparse matrices
        _spMatrix[i] = objectFactory().instantiateSparseMatrix(_solver[i]->giveRequiredMatrixFormat());
        _spMatrix[i]->setSymmetryTo(_symmetry[i]);
    }

    // Maximum number of iterations
    verifyKeyword(fp, "MaxIterations", _name);
    _maxIter = getIntegerInputFrom(fp, "Failed to read maximum number of iterations from input file!", _name);

    // Directive upon reaching maximum iterations
    std::string directive = getStringInputFrom(fp, "Failed to read solution method directive from input file!", _name);
    if ( directive == "Abort" )
        _abortAtMaxIter = true;
    else if ( directive == "Continue" )
        _abortAtMaxIter = false;
    else
    {
        throw std::runtime_error("ERROR: Invalid directive to solution method encountered! Valid options are \"Abort\" or \"Continue\"\n");
    }
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
    
#pragma omp parallel for
    for ( int i = 0; i < (int)dof.size(); i++)
    {
        int ssNum = analysisModel().dofManager().giveSubsystemNumberFor(dof[i]);
        int eqNum = analysisModel().dofManager().giveEquationNumberAt(dof[i]);
        if ( ssNum == _subsysNum[subsysIdx] )
            subsysDof[eqNum] = dof[i];
    }
    
    bool converged = true;
    for ( int j = 0; j < (int)_subsysDofGroup[subsysIdx].size(); j++ )
    {
        int idx = this->giveIndexForDofGroup(_subsysDofGroup[subsysIdx][j]);
        bool dofGrpConverged = _convergenceCriterion[idx]->checkConvergenceOf(resid, subsysDof);
        if ( !dofGrpConverged )
            converged = false;
    }
    
    if ( REPORT_SUBSYSTEM_STATUS )
    {
        std::printf("\n\n        Global Iter. %d, Subsystem %d, Local Iter. %d", globalIterCount, _subsysNum[subsysIdx], localIterCount);    
        std::printf("\n        DOF Grp      Norm         Criterion");
        std::printf("\n       -------------------------------------------------------");

        for (int j = 0; j < (int)_subsysDofGroup[subsysIdx].size(); j++)
        {
            int grpNum = _subsysDofGroup[subsysIdx][j];
            int grpIdx = this->giveIndexForDofGroup(grpNum);
            _convergenceCriterion[grpIdx]->reportConvergenceStatus();
        }
        std::fflush(stdout);
    }

    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    diagnostics().addConvergenceCheckTime(tictoc.count());

    return converged;
}