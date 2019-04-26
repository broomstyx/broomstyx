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

#include "AlternateMinimization.hpp"
#include <chrono>
#include <cmath>
#include <stdexcept>
#include <string>
#include <tuple>
#include <omp.h>
#include <list>
#include <iterator>

#include "../Core/AnalysisModel.hpp"
#include "../Core/ObjectFactory.hpp"
#include "../Core/DofManager.hpp"
#include "../Core/DomainManager.hpp"
#include "../Core/LoadStep.hpp"
#include "../Core/NumericsManager.hpp"
#include "../MeshReaders/MeshReader.hpp"
#include "../Core/SolutionManager.hpp"
#include "../LinearSolvers/LinearSolver.hpp"
#include "../Numerics/Numerics.hpp"
#include "../Util/linearAlgebra.hpp"
#include "../Util/readOperations.hpp"
#include "../Util/reductions.hpp"

using namespace broomstyx;

registerBroomstyxObject(SolutionMethod, AlternateMinimization)

AlternateMinimization::AlternateMinimization()
{
    _name = "AlternateMinimization";
}

AlternateMinimization::~AlternateMinimization()
{
    for ( int i = 0; i < _nSubsystems; i++ )
    {
        delete _spMatrix[i];
        delete _solver[i];
    }
}

// Public methods
// ---------------------------------------------------------------------------
int AlternateMinimization::computeSolutionFor( int stage
                                             , const std::vector<BoundaryCondition>& bndCond
                                             , const std::vector<FieldCondition>& fldCond
                                             , const TimeData& time )
{
    std::chrono::time_point<std::chrono::system_clock> tic, toc;
    std::chrono::duration<double> tictoc;
    
    _substepCount += 1;
    
    // Get all active DOFs at current stage
    std::vector<Dof*> dof = analysisModel().dofManager().giveActiveDofsAtStage(stage);
    
    // Impose constraints on nodal DOFs
    std::printf("    %-40s", "Imposing constraints ...");
    tic = std::chrono::high_resolution_clock::now();
    this->imposeConstraintsAt(stage, bndCond, time);
    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    std::printf("done (time = %f sec.)\n", tictoc.count());
    
    // Initialize global vectors
    std::vector<RealVector> dU(_nSubsystems, RealVector());
    auto resid = dU;
    auto initResid = dU;
    auto rhs = dU;
    auto lhs = dU;

    for ( int i = 0; i < _nSubsystems; i++)
        dU[i].init(_nUnknowns[i]);
    
    // Start iterations
    tic = std::chrono::high_resolution_clock::now();
    bool converged = false;
    int iterCount = 0;
    while ( !converged && iterCount <= _maxIter )
    {
        // Perform pre-iteration operations, if any ...
        std::vector<Numerics*> numerics = analysisModel().numericsManager().giveAllNumerics();
        for ( Numerics* curNumerics : numerics )
            curNumerics->performPreIterationOperationsAt(stage, iterCount);
        
        // Reinitialize data for calculation of residual convergence
        for ( int i = 0; i < _nDofGroups; i++ )
        {
            _sumAbsFlux(i) = 0;
            _fluxCount(i) = 0;
        }
        
        // Calculate residual for each subsystem
        // Note: We calculate global right hand sides in each
        // iteration to accommodate numerics classes that implement special
        // methods for imposing constraints
        for ( int i = 0; i < _nSubsystems; i++)
        {
            rhs[i] = assembleRightHandSide(stage, _subsysNum[i], bndCond, fldCond, time);
            lhs[i] = assembleLeftHandSide(stage, _subsysNum[i], time);
            resid[i] = rhs[i] - lhs[i];
        }
        
        if ( iterCount == 0 )
            std::printf("\n    Initial guess (LS # %d, SS # %d)", _loadStep->giveLoadStepNum(), _substepCount);
        else
            std::printf("\n    Iter.No. %-4d (LS # %d, SS # %d)", iterCount, _loadStep->giveLoadStepNum(), _substepCount);
        std::printf("\n   -----------------");
        
        // Compute convergence norms
        RealMatrix normDat;
        
        // Disallow convergence at initial guess
        if ( iterCount == 0 )
            converged = false;
        else
        {
//            if ( _enableLineSearch )
//                this->performLineSearchAt(stage, dof, rhs, initResid, time, dU, resid);

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
            for ( Numerics* curNumerics : numerics )
            {
                curNumerics->printPostIterationMessage(stage);
            }
        
        if ( converged )
        {
            // Clear memory for solvers
            for ( int i = 0; i < _nSubsystems; i++ )
                _solver[i]->clearInternalMemory();
            
            _loadStep->writeIterationDataForStage(stage, time.target, iterCount);
        }
        else
        {
            // Start of new iteration
            ++iterCount;
            
            tic = std::chrono::high_resolution_clock::now();
            
            for ( int i = 0; i < _nSubsystems; i++ )
            {
                // Recalculate subsystem residual vector
                if ( i > 0 )
                {
                    rhs[i] = this->assembleRightHandSide(stage, _subsysNum[i], bndCond, fldCond, time);
                    lhs[i] = this->assembleLeftHandSide(stage, _subsysNum[i], time);
                    resid[i] = rhs[i] - lhs[i];
                }
                
                // Assemble system Jacobian
                _spMatrix[i]->initializeValues();
                this->assembleJacobian(stage, _subsysNum[i], time);
                
                // Solve system and apply over-relaxation
                _solver[i]->allocateInternalMemoryFor(_spMatrix[i]);
                if ( _solver[i]->takesInitialGuess() )
                    _solver[i]->setInitialGuessTo(dU[i]);
                dU[i] = _solver[i]->solve(_spMatrix[i], resid[i]);
                
#ifdef _OPENMP
#pragma omp parallel for
#endif
                for ( int j = 0; j < (int)dof.size(); j++ )
                    if ( analysisModel().dofManager().giveSubsystemNumberFor(dof[j]) == _subsysNum[i] )
                    {
                        int eqNo = analysisModel().dofManager().giveEquationNumberAt(dof[j]);
                        if ( eqNo != UNASSIGNED )
                        {
                            // Apply over-relaxation and update primary variable at DOF
                            dU[i](eqNo) *= _overRelaxation(i);
                            analysisModel().dofManager().updatePrimaryVariableAt(dof[j], dU[i](eqNo), correction);
                        }   
                    }
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
void AlternateMinimization::formSparsityProfileForStage( int stage )
{
    std::chrono::time_point<std::chrono::system_clock> tic, toc;
    std::chrono::duration<double> tictoc;
    
    // Count number of DOFs belonging to each DOF group
    _dofGrpCount.init(_nDofGroups);
    
    std::vector<Dof*> dof = analysisModel().dofManager().giveActiveDofsAtStage(stage);
    int nActiveDof = (int)dof.size();
    for ( int i = 0; i < nActiveDof; i++ )
    {
        int grpNo = analysisModel().dofManager().giveGroupNumberFor(dof[i]);
        int idx = this->giveIndexForDofGroup(grpNo);
        if ( idx >= 0 )
            _dofGrpCount(idx) += 1.;
    }
    // Assigning DOFs to subsystems
    std::printf("\n  %-40s", "Assigning DOFs to subsystems ...");
    tic = std::chrono::high_resolution_clock::now();
    
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < nActiveDof; i++ )
    {
        int grpNo = analysisModel().dofManager().giveGroupNumberFor(dof[i]);
        for ( int j = 0; j < _nSubsystems; j++ )
            for ( int k = 0; k < (int)_subsysDofGroup[j].size(); k++ )
                if ( grpNo == _subsysDofGroup[j][k] )
                {
                    analysisModel().dofManager().setSubsystemFor(dof[i], _subsysNum[j]);
#ifdef _OPENMP
#pragma omp atomic
#endif
                    _nUnknowns[j] += 1;
                }
    }
    
    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    std::printf("done (time = %f sec.)", tictoc.count());
    
    // Assigning DOF equation numbers
    std::printf("\n  %-40s", "Assigning DOF equation numbers ...");
    tic = std::chrono::high_resolution_clock::now();
    std::vector<int> subsysEqNo;
    subsysEqNo.assign(_nSubsystems, 0);
    
    for ( int i = 0; i < (int)dof.size(); i++)
    {
        // Check that DOF group is included in solution
        int grp = analysisModel().dofManager().giveGroupNumberFor(dof[i]);
        int idx = this->giveIndexForDofGroup(grp);
        
        if ( idx >= 0 )
        {
            int ssn = analysisModel().dofManager().giveSubsystemNumberFor(dof[i]);
            for ( int j = 0; j < _nSubsystems; j++ )
                if ( _subsysNum[j] == ssn )
                    analysisModel().dofManager().setEquationNumberFor(dof[i], subsysEqNo[j]++);
        }
    }
    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    std::printf("done (time = %f sec.)", tictoc.count());
    
    // Determine sparsity profile for each subsystem
    for ( int i = 0; i < _nSubsystems; i++ )
    {
        std::printf("\n  %-40s", "Determining sparsity pattern ...");
        std::fflush(stdout);
        tic = std::chrono::high_resolution_clock::now();

        // Initialize sparsity profile
        int nUnknowns = subsysEqNo[i];
        _spMatrix[i]->initializeProfile(nUnknowns, nUnknowns);
        _spMatrix[i]->setSymmetryTo(_symmetry[i]);
        
        int nCells = analysisModel().domainManager().giveNumberOfDomainCells();

        // Matrix assembly for determining sparsity profile cannot be 
        // parallelized because std::set is not thread-safe.
        for ( int j = 0; j < nCells; j++ )
        {
            Cell* curCell = analysisModel().domainManager().giveDomainCell(j);
            int label = analysisModel().domainManager().giveLabelOf(curCell);
            Numerics* numerics = analysisModel().domainManager().giveNumericsForDomain(label);

            TimeData time;
            std::vector<Dof*> rowDof, colDof;
            std::tie(rowDof,colDof,std::ignore) = numerics->giveStaticCoefficientMatrixAt(curCell, stage, _subsysNum[i], time);
            
            int lnnz = rowDof.size();
            for ( int k = 0; k < lnnz; k++)
            {
                if ( rowDof[k] && colDof[k] )
                {
                    int rowNum = analysisModel().dofManager().giveEquationNumberAt(rowDof[k]);
                    int rssNum = analysisModel().dofManager().giveSubsystemNumberFor(rowDof[k]);
                    int colNum = analysisModel().dofManager().giveEquationNumberAt(colDof[k]);
                    int cssNum = analysisModel().dofManager().giveSubsystemNumberFor(colDof[k]);
                    if ( rowNum != UNASSIGNED && rssNum == _subsysNum[i] && colNum != UNASSIGNED && cssNum == _subsysNum[i] )
                        _spMatrix[i]->insertNonzeroComponentAt(rowNum, colNum);
                }
            }
        }

        _spMatrix[i]->finalizeProfile();

        toc = std::chrono::high_resolution_clock::now();
        tictoc = toc - tic;
        std::printf("done (time = %f sec.)\n", tictoc.count());

        // Compute storage size in bytes
        int nnz = _spMatrix[i]->giveNumberOfNonzeros();
        double spratio = (double)nnz/((double)nUnknowns*(double)nUnknowns);
        std::printf("\n    Stage %-2d, Subsystem %-2d: nUnknowns = %d, nnz = %d", stage, _subsysNum[i], nUnknowns, nnz);
        std::printf("\n         Sparsity ratio = %.4e\n", spratio);
    }
}
// ---------------------------------------------------------------------------
void AlternateMinimization::initializeSolvers()
{
    for ( int i = 0; i < _nSubsystems; i++ )
        _solver[i]->initialize();
}
// ---------------------------------------------------------------------------
void AlternateMinimization::readDataFromFile( FILE* fp )
{    
    std::string key;
    
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
//    
//    _ctrlParam.init(_nDofGroups,6);
//    
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
    _overRelaxation.init(_nSubsystems);
    
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
                
        // Read over-relaxation parameter
        verifyKeyword(fp, key = "OverRelaxation", _name);
        _overRelaxation(i) = getRealInputFrom(fp, "Failed to read over-relaxation parameter from input file!", _name);
        
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
RealVector AlternateMinimization::assembleLeftHandSide( int stage
                                                      , int subsys
                                                      , const TimeData& time )
{
    // Find subsystem index
    int idx = this->giveIndexForSubsystem(subsys);
    
    RealVector lhs(_nUnknowns[idx]);
    
    // Assembly of global internal force vector for current subsystem
    int nCells = analysisModel().domainManager().giveNumberOfDomainCells();
    
    // ...because OpenMP doesn't allow class members to be shared by parallel threads
    RealVector sumAbsFlux = _sumAbsFlux;
    RealVector fluxCount = _fluxCount;
    
#ifdef _OPENMP
#pragma omp parallel for reduction ( + : sumAbsFlux, fluxCount )
#endif
    for ( int i = 0; i < nCells; i++ )
    {
        Cell* curCell = analysisModel().domainManager().giveDomainCell(i);
        Numerics* numerics = analysisModel().domainManager().giveNumericsFor(curCell);
        
        std::vector<Dof*> rowDof;
        RealVector localLhs;
        
        // Calculate cell internal forces
        std::tie(rowDof,localLhs) = numerics->giveStaticLeftHandSideAt(curCell, stage, subsys, time);
        
        // Assembly
        for ( int j = 0; j < (int)rowDof.size(); j++ )
        {
            if ( rowDof[j] )
            {
                int rowNum = analysisModel().dofManager().giveEquationNumberAt(rowDof[j]);
                int ssNum = analysisModel().dofManager().giveSubsystemNumberFor(rowDof[j]);

                if ( rowNum != UNASSIGNED && ssNum == subsys )
                {
                    int dgNum = analysisModel().dofManager().giveGroupNumberFor(rowDof[j]);
                    int dgIdx = this->giveIndexForDofGroup(dgNum);
                    
                    if ( std::fabs(localLhs(j)) > _absTolRes(dgIdx) )
                    {
                        sumAbsFlux(dgIdx) += std::fabs(localLhs(j));
                        fluxCount(dgIdx) += 1.;
                    }
#ifdef _OPENMP
#pragma omp atomic
#endif
                    lhs(rowNum) += localLhs(j);
                }
            }
        }
    }
    
    _sumAbsFlux = sumAbsFlux;
    _fluxCount = fluxCount;
    
    return lhs;
}
// ---------------------------------------------------------------------------
void AlternateMinimization::assembleJacobian( int stage
                                            , int subsys
                                            , const TimeData& time )
{
    int nCells = analysisModel().domainManager().giveNumberOfDomainCells();

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < nCells; i++)
    {
        Cell* curCell = analysisModel().domainManager().giveDomainCell(i);
        Numerics* numerics = analysisModel().domainManager().giveNumericsFor(curCell);
        
        std::vector<Dof*> rowDof, colDof;
        RealVector coefVal;

        std::tie(rowDof,colDof,coefVal) = numerics->giveStaticCoefficientMatrixAt(curCell, stage, subsys, time);

        for ( int j = 0; j < (int)rowDof.size(); j++)
        {
            if ( rowDof[j] && colDof[j] )
            {
                int rowNum = analysisModel().dofManager().giveEquationNumberAt(rowDof[j]);
                int rssNum = analysisModel().dofManager().giveSubsystemNumberFor(rowDof[j]);
                int colNum = analysisModel().dofManager().giveEquationNumberAt(colDof[j]);
                int cssNum = analysisModel().dofManager().giveSubsystemNumberFor(colDof[j]);
                
                if ( rowNum != UNASSIGNED && rssNum == subsys && colNum != UNASSIGNED && cssNum == subsys )
                {
                    int idx = this->giveIndexForSubsystem(subsys);
                    _spMatrix[idx]->atomicAddToComponent(rowNum, colNum, coefVal(j));
                }
            }
        }
    }
}
// ---------------------------------------------------------------------------
RealVector AlternateMinimization::assembleRightHandSide( int stage
                                                       , int subsys
                                                       , const std::vector<BoundaryCondition>& bndCond
                                                       , const std::vector<FieldCondition>& fldCond
                                                       , const TimeData& time )
{   
    // Get subsystem index
    int idx = this->giveIndexForSubsystem(subsys);
    RealVector rhs(_nUnknowns[idx]);
    
    // Loop through all field conditions
    int nCells = analysisModel().domainManager().giveNumberOfDomainCells();
    
    // ...because OpenMP doesn't allow class members to be shared by parallel threads
    RealVector sumAbsFlux = _sumAbsFlux;
    RealVector fluxCount = _fluxCount;
    
#ifdef _OPENMP
#pragma omp parallel for reduction ( + : sumAbsFlux, fluxCount )
#endif
    for ( int i = 0; i < nCells; i++ )
    {
        Cell* curCell = analysisModel().domainManager().giveDomainCell(i);
        int label = analysisModel().domainManager().giveLabelOf(curCell);
            
        for ( int ifc = 0; ifc < (int)fldCond.size(); ifc++ )
        {
            int fcLabel = analysisModel().domainManager().givePhysicalEntityNumberFor(fldCond[ifc].domainLabel());
            if ( label == fcLabel )
            {
                RealVector localRhs;
                std::vector<Dof*> rowDof;

                Numerics* numerics = analysisModel().domainManager().giveNumericsForDomain(label);
                std::tie(rowDof,localRhs) = numerics->giveStaticRightHandSideAt(curCell, stage, subsys, fldCond[ifc], time);

                for ( int j = 0; j < localRhs.dim(); j++)
                {
                    if ( rowDof[j] )
                    {
                        int rowNum = analysisModel().dofManager().giveEquationNumberAt(rowDof[j]);
                        int ssNum = analysisModel().dofManager().giveSubsystemNumberFor(rowDof[j]);

                        if ( rowNum != UNASSIGNED && ssNum == subsys )
                        {
                            int dgNum = analysisModel().dofManager().giveGroupNumberFor(rowDof[j]);
                            int dgIdx = this->giveIndexForDofGroup(dgNum);
                            
                            if ( std::fabs(localRhs(j)) > _absTolRes(dgIdx) )
                            {
                                sumAbsFlux(dgIdx) += std::fabs(localRhs(j));
                                fluxCount(dgIdx) += 1.;
                            }
#ifdef _OPENMP
#pragma omp atomic
#endif
                            rhs(rowNum) += localRhs(j);
                        }
                    }
                }
            }
        }
    }
    
    // Loop through all natural boundary conditions
    for (int ibc = 0; ibc < (int)bndCond.size(); ibc++)
    {
        int boundaryId = analysisModel().domainManager().givePhysicalEntityNumberFor(bndCond[ibc].boundaryName());
        Numerics* numerics = analysisModel().numericsManager().giveNumerics(bndCond[ibc].targetNumerics());
        
        int nBCells = analysisModel().domainManager().giveNumberOfBoundaryCells();
#ifdef _OPENMP
#pragma omp parallel for reduction ( + : sumAbsFlux, fluxCount )
#endif
        for ( int i = 0; i < nBCells; i++ ) 
        {
            Cell* curCell = analysisModel().domainManager().giveBoundaryCell(i);
            int label = analysisModel().domainManager().giveLabelOf(curCell);

            if ( label == boundaryId )
            {
                RealVector localRhs;
                std::vector<Dof*> rowDof;

                // Specifics of BC imposition are handled by numerics
                std::tie(rowDof,localRhs) = numerics->giveStaticRightHandSideAt(curCell, stage, subsys, bndCond[ibc], time);
                
                for ( int j = 0; j < localRhs.dim(); j++)
                {
                    if ( rowDof[j] )
                    {
                        int rowNum = analysisModel().dofManager().giveEquationNumberAt(rowDof[j]);
                        int ssNum = analysisModel().dofManager().giveSubsystemNumberFor(rowDof[j]);

                        if ( rowNum != UNASSIGNED && ssNum == subsys )
                        {
                            int dgNum = analysisModel().dofManager().giveGroupNumberFor(rowDof[j]);
                            int dgIdx = this->giveIndexForDofGroup(dgNum);
                            
                            if ( std::fabs(localRhs(j)) > _absTolRes(dgIdx) )
                            {
                                sumAbsFlux(dgIdx) += std::fabs(localRhs(j));
                                fluxCount(dgIdx) += 1.;
                            }
#ifdef _OPENMP
#pragma omp atomic
#endif
                            rhs(rowNum) += localRhs(j);
                        }
                    }
                }
            }
        }
    }
    
    _sumAbsFlux = sumAbsFlux;
    _fluxCount = fluxCount;
    
    return rhs;
}
// ---------------------------------------------------------------------------
std::tuple<bool,RealMatrix>
AlternateMinimization::computeConvergenceNormsFrom( const std::vector<RealVector>& resid
                                                  , const std::vector<Dof*>& dof )
{    
    /*********************************************** 
     * Initialize matrix to store convergence data
     *   col 0: norm
     *   col 1: criterion
     ***********************************************/
    
    RealMatrix normDat(2*_nDofGroups, 2);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < (int)dof.size(); i++ )
    {
        // Find group number of DOF
        int grpNum = analysisModel().dofManager().giveGroupNumberFor(dof[i]);
        int grpIdx = this->giveIndexForDofGroup(grpNum);
        
        // Find subsystem of DOF
        int ssNum = analysisModel().dofManager().giveSubsystemNumberFor(dof[i]);
        int ssIdx = this->giveIndexForSubsystem(ssNum);
        
        if ( grpIdx >= 0 && ssIdx >= 0 )
        {
            // Find equation number of DOF
            int eqNo = analysisModel().dofManager().giveEquationNumberAt(dof[i]);

            // Contribution to L2-norm of corrections
            double corrVal = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof[i], correction);
            
#ifdef _OPENMP
#pragma omp atomic
#endif
            normDat(2*grpIdx,0) += corrVal*corrVal;

            // Contribution to L2-norm of incremental solution
            double incVal = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof[i], incremental_value);

#ifdef _OPENMP
#pragma omp atomic
#endif
            normDat(2*grpIdx,1) += incVal*incVal;

            // Contribution to L2-norm of residual
            double rVal = resid[ssIdx](eqNo);

#ifdef _OPENMP
#pragma omp atomic
#endif
            normDat(2*grpIdx+1,0) += rVal*rVal;
        }
    }
    
    for ( int i = 0; i < _nDofGroups; i++ )
    {
        // Normalize entries
        normDat(2*i,0) = std::sqrt(normDat(2*i,0)/_dofGrpCount(i));
        normDat(2*i+1,0) = std::sqrt(normDat(2*i+1,0)/_dofGrpCount(i));
        normDat(2*i,1) = std::sqrt(normDat(2*i,1)/_dofGrpCount(i));
        
        // Calculate criteria using relative tolerances
        normDat(2*i,1) *= _relTolCor(i);
        normDat(2*i+1,1) = _relTolRes(i)*_sumAbsFlux(i)/_fluxCount(i);
        
        // Check against absolute tolerances
        if ( normDat(2*i,1) < _absTolCor(i) )
            normDat(2*i,1) = _absTolCor(i);
        if ( normDat(2*i+1,1) < _absTolRes(i) )
            normDat(2*i+1,1) = _absTolRes(i);
    }
    
    bool isConverged = true;
    for ( int i = 0; i < _nDofGroups; i++ )
    {
//        // Consider converged if either correction or residual is below absolute tolerance
//        if ( normDat(2*i,0) > _ctrlParam(i,2) && normDat(2*i+1,0) > _ctrlParam(i,3) )
//        {
            if ( normDat(2*i,0) > normDat(2*i,1) )
                isConverged = false;
            if ( normDat(2*i+1,0) > normDat(2*i+1,1) )
                isConverged = false;
//        }
    }
    
    return std::make_tuple(isConverged,normDat);
}
// ---------------------------------------------------------------------------
int AlternateMinimization::giveIndexForDofGroup( int dofGroupNum )
{
    int idx = -1;
    for ( int i = 0; i < _nDofGroups; i++ )
        if ( dofGroupNum == _dofGrpNum[i] )
            idx = i;
    
    return idx;
}
// ---------------------------------------------------------------------------
int AlternateMinimization::giveIndexForSubsystem( int subsysNum )
{
    int idx = -1;
    for ( int i = 0; i < _nSubsystems; i++ )
        if ( subsysNum == _subsysNum[i] )
            idx = i;
    
    return idx;
}
// ---------------------------------------------------------------------------
void AlternateMinimization::reportConvergenceStatus( const RealMatrix& normDat )
{
    // Report status
    std::printf("\n\n    DOF Grp   L2-Norm         Criterion");
    std::printf("\n   -------------------------------------------------------");
    
    for ( int i = 0; i < _nDofGroups; i++ )
    {
        // Convergence of corrections
        std::printf("\n     C  %-6d%e   %e ", _dofGrpNum[i], normDat(2*i, 0), normDat(2*i, 1));
        
        if ( normDat(2*i,0) <= normDat(2*i,1) )
            std::printf(" << CONVERGED");
        
        // Convergence of residuals
        std::printf("\n     R  %-6d%e   %e ", _dofGrpNum[i], normDat(2*i+1, 0), normDat(2*i+1, 1));
        
        if ( normDat(2*i+1,0) <= normDat(2*i+1,1) )
            std::printf(" << CONVERGED");
    }
}
