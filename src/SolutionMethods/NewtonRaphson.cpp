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

#include "NewtonRaphson.hpp"
#include <chrono>
#include <cmath>
#include <stdexcept>
#include <cstring>
#include <tuple>
#include "omp.h"

#include "ConvergenceCriteria/ConvergenceCriterion.hpp"
#include "Core/AnalysisModel.hpp"
#include "Core/ObjectFactory.hpp"
#include "Core/Diagnostics.hpp"
#include "Core/DofManager.hpp"
#include "Core/DomainManager.hpp"
#include "Core/LoadStep.hpp"
#include "Core/NumericsManager.hpp"
#include "Core/SolutionManager.hpp"
#include "LinearSolvers/LinearSolver.hpp"
#include "Numerics/Numerics.hpp"
#include "SparseMatrix/SparseMatrix.hpp"
#include "Util/linearAlgebra.hpp"
#include "Util/readOperations.hpp"
#include "Util/reductions.hpp"

using namespace broomstyx;

registerBroomstyxObject(SolutionMethod, NewtonRaphson)

// Constructor
NewtonRaphson::NewtonRaphson()
{
    _name = "NewtonRaphson";
}

// Destructor
NewtonRaphson::~NewtonRaphson()
{
    delete _spMatrix;
    delete _solver;
}

// Public methods
// ----------------------------------------------------------------------------
int NewtonRaphson::computeSolutionFor( int stage
                                     , const std::vector<BoundaryCondition>& bndCond
                                     , const std::vector<FieldCondition>& fldCond
                                     , TimeData& time )
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
    
    RealVector dU(_nUnknowns), rhs, lhs, resid;
        
    // Start iterations
    tic = std::chrono::high_resolution_clock::now();
    bool converged = false;
    int iterCount = 0;
    while ( !converged && iterCount <= _maxIter )
    {
        // Pre-iteration operations, if any ...
        auto innertic = std::chrono::high_resolution_clock::now();
        std::vector<Numerics*> numVec = analysisModel().numericsManager().giveAllNumerics();
        for ( Numerics* curNumerics : numVec )
            curNumerics->performPreIterationOperationsAt(stage, iterCount);
        auto innertoc = std::chrono::high_resolution_clock::now();
        tictoc = innertoc - innertic;
        diagnostics().addPostprocessingTime(tictoc.count());

        // Calculate global external force vectors
        // Note: this is normally redudant, but we recalculate them in each
        // iteration to accommodate numerics classes that implement special
        // methods for imposing constraints

        analysisModel().dofManager().resetSecondaryVariablesAtStage(stage);
        
        // Calculate residual for each subsystem
        rhs = assembleRightHandSide(stage, bndCond, fldCond, time);
        lhs = assembleLeftHandSide(stage, time);
        resid = rhs - lhs;
        
        if ( iterCount == 0 )
            std::printf("\n    Initial guess");
        else
            std::printf("\n    Iter.No. %-4d", iterCount);
        std::printf("\n   -----------------");
        
        // Compute convergence norms
        RealMatrix normDat(2*_nDofGroups, 2);
        
        // Disallow convergence at initial guess
        if ( iterCount == 0 )
            converged = false;
        else
        {
            converged = true;
            for ( int i = 0; i < _nDofGroups; i++ )
            {
                int idx = this->giveIndexForDofGroup(_dofGrpNum[i]);
                bool dofGrpConverged = _convergenceCriterion[idx]->checkConvergenceOf(resid, dof);
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
            for ( Numerics* curNumerics : numVec )
            {
                curNumerics->printPostIterationMessage(stage);
            }
        
        if ( converged )
        {
            // Clear memory for solvers
            _solver->clearInternalMemory();
            innertic = std::chrono::high_resolution_clock::now();
            _loadStep->writeIterationDataForStage(stage, time.giveTargetTime(), iterCount);
            innertoc = std::chrono::high_resolution_clock::now();
            tictoc = innertoc - innertic;
            diagnostics().addOutputWriteTime(tictoc.count());
        }
        else if ( iterCount == _maxIter )
        {
            // Clear memory for solver
            _solver->clearInternalMemory();
            
            innertic = std::chrono::high_resolution_clock::now();
            _loadStep->writeIterationDataForStage(stage, time.giveTargetTime(), iterCount);
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
            
            // Recalculate subsystem residual vector
            rhs = this->assembleRightHandSide(stage, bndCond, fldCond, time);
            lhs = this->assembleLeftHandSide(stage, time);
            resid = rhs - lhs;
                
            // Assemble system Jacobian
            _spMatrix->initializeValues();
            this->assembleJacobian(stage, time);
                
            // Solve system and apply over-relaxation
            innertic = std::chrono::high_resolution_clock::now();

            _solver->allocateInternalMemoryFor(_spMatrix);
            if ( _solver->takesInitialGuess() )
                _solver->setInitialGuessTo(dU);
            dU = _overRelaxation*_solver->solve(_spMatrix, resid);

            innertoc = std::chrono::high_resolution_clock::now();
            tictoc = innertoc - innertic;
            diagnostics().addSolveTime(tictoc.count());

            innertic = std::chrono::high_resolution_clock::now();
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for ( int j = 0; j < (int)dof.size(); j++ )
            {
                int eqNo = analysisModel().dofManager().giveEquationNumberAt(dof[j]);
                if ( eqNo != UNASSIGNED )
                {
                    // Update value primary variable at DOF
                    analysisModel().dofManager().updatePrimaryVariableAt(dof[j], dU(eqNo), correction);
                }   
            }
            innertoc = std::chrono::high_resolution_clock::now();
            tictoc = innertoc - innertic;
            diagnostics().addUpdateTime(tictoc.count());
        }
    }
    
    std::printf("\n");
    if ( !converged )
        return 1;
    else
        return 0;
}
// ----------------------------------------------------------------------------
void NewtonRaphson::formSparsityProfileForStage( int stage )
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
    
    // Assigning DOF equation numbers
    std::printf("\n  %-40s", "Assigning DOF equation numbers ...");
    tic = std::chrono::high_resolution_clock::now();
    
    int eqNo = 0;
    
    for ( int i = 0; i < (int)dof.size(); i++)
    {
        // Check that DOF group is included in solution
        int grp = analysisModel().dofManager().giveGroupNumberFor(dof[i]);
        int idx = this->giveIndexForDofGroup(grp);
        
        if ( idx >= 0 )
            analysisModel().dofManager().setEquationNumberFor(dof[i], eqNo++);
    }
    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    std::printf("done (time = %f sec.)", tictoc.count());
    
    // Determing sparsity profile
    std::printf("\n  %-40s", "Determining sparsity pattern ...");
    std::fflush(stdout);
    tic = std::chrono::high_resolution_clock::now();

    _nUnknowns = eqNo;
    _spMatrix->initializeProfile(_nUnknowns, _nUnknowns);
    _spMatrix->setSymmetryTo(_symmetry);

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
        std::tie(rowDof,colDof,std::ignore) = numerics->giveStaticCoefficientMatrixAt(curCell, stage, UNASSIGNED, time);

        int lnnz = rowDof.size();
        for ( int k = 0; k < lnnz; k++)
        {
            if ( rowDof[k] && colDof[k] )
            {
                int rowNum = analysisModel().dofManager().giveEquationNumberAt(rowDof[k]);
                int colNum = analysisModel().dofManager().giveEquationNumberAt(colDof[k]);
                if ( rowNum != UNASSIGNED && colNum != UNASSIGNED )
                    _spMatrix->insertNonzeroComponentAt(rowNum, colNum);
            }
        }
    }

    _spMatrix->finalizeProfile();

    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    std::printf("done (time = %f sec.)\n", tictoc.count());

    // Compute storage size in bytes
    int nnz = _spMatrix->giveNumberOfNonzeros();
    double spratio = (double)nnz/((double)_nUnknowns*(double)_nUnknowns);
    std::printf("\n    Stage %-2d: nUnknowns = %d, nnz = %d", stage, _nUnknowns, nnz);
    std::printf("\n         Sparsity ratio = %.4e\n", spratio);
}
// ---------------------------------------------------------------------------
void NewtonRaphson::initializeSolvers()
{
    _solver->initialize();
}
// ----------------------------------------------------------------------------
void NewtonRaphson::readDataFromFile( FILE* fp )
{    
    // Read number of DOF groups
    verifyKeyword(fp, "DofGroups", _name);
    _nDofGroups = getIntegerInputFrom(fp, "Failed to read number of DOF groups from input file!", _name);
    _dofGrpNum.assign(_nDofGroups, -1);

    // Read convergence criterion and associated data
    _convergenceCriterion.assign(_nDofGroups, nullptr);    
    for ( int i = 0; i < _nDofGroups; i++ )
    {
        _dofGrpNum[i] = getIntegerInputFrom(fp, "Failed to read DOF groupnumber from input file!", _name);
        std::string convCritName = getStringInputFrom(fp, "Failed to read convergence criterion type from input file!", _name);
        _convergenceCriterion[i] = objectFactory().instantiateConvergenceCriterion(convCritName);
        _convergenceCriterion[i]->initialize(_dofGrpNum[i]);
        _convergenceCriterion[i]->readDataFromFile(fp);
    }
    
    // Read linear solver to be used
    verifyKeyword(fp, "LinearSolver", _name);
    std::string key = getStringInputFrom(fp, "Failed to read linear solver from input file!", _name);
    _solver = objectFactory().instantiateLinearSolver(key);
    _solver->readDataFrom(fp);
    _symmetry = _solver->giveSymmetryOption();

    // Read over-relaxation parameter
    verifyKeyword(fp, "OverRelaxation", _name);
    _overRelaxation = getRealInputFrom(fp, "Failed to read over-relaxation parameter from input file!", _name);

    // Create sparse matrices
    _spMatrix = objectFactory().instantiateSparseMatrix(_solver->giveRequiredMatrixFormat());
    _spMatrix->setSymmetryTo(_symmetry);

    // Maximum number of iterations
    verifyKeyword(fp, key = "MaxIterations", _name);
    _maxIter = getIntegerInputFrom(fp, "Failed to read maximum number of iterations from input file!", _name);
}

// Private methods
// ----------------------------------------------------------------------------
RealVector NewtonRaphson::assembleLeftHandSide( int stage, const TimeData& time )
{
    std::chrono::time_point<std::chrono::system_clock> tic, toc;
    std::chrono::duration<double> tictoc;
    tic = std::chrono::high_resolution_clock::now();

    RealVector lhs(_nUnknowns);
    
    // Assembly of global internal force vector for current subsystem
    int nCells = analysisModel().domainManager().giveNumberOfDomainCells();
    
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        int threadNum = 0;
#ifdef _OPENMP
        threadNum = omp_get_thread_num();
#endif
#ifdef _OPENMP
#pragma omp for
#endif
        for (int i = 0; i < nCells; i++)
        {
            Cell* curCell = analysisModel().domainManager().giveDomainCell(i);
            Numerics* numerics = analysisModel().domainManager().giveNumericsFor(curCell);
            
            std::vector<Dof*> rowDof;
            RealVector localLhs;
            
            // Calculate cell internal forces
            std::tie(rowDof,localLhs) = numerics->giveStaticLeftHandSideAt(curCell, stage, UNASSIGNED, time);
            
            // Assembly
            for ( int j = 0; j < (int)rowDof.size(); j++ )
            {
                if ( rowDof[j] )
                {
                    int rowNum = analysisModel().dofManager().giveEquationNumberAt(rowDof[j]);
                    int dofGrp = analysisModel().dofManager().giveGroupNumberFor(rowDof[j]);
                    int idx = this->giveIndexForDofGroup(dofGrp);
                    _convergenceCriterion[idx]->processLocalResidualContribution(localLhs(j), threadNum);

                    if ( rowNum != UNASSIGNED )
                    {
#ifdef _OPENMP
#pragma omp atomic
#endif
                        lhs(rowNum) += localLhs(j);
                    }

                    analysisModel().dofManager().addToSecondaryVariableAt(rowDof[j], localLhs(j));
                }
            }
        }
    }

    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    diagnostics().addLhsAssemblyTime(tictoc.count());
    
    return lhs;
}
// ---------------------------------------------------------------------------
void NewtonRaphson::assembleJacobian( int stage, const TimeData& time )
{
    std::chrono::time_point<std::chrono::system_clock> tic, toc;
    std::chrono::duration<double> tictoc;
    tic = std::chrono::high_resolution_clock::now();

    int nCells = analysisModel().domainManager().giveNumberOfDomainCells();

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < nCells; i++ )
    {
        Cell* curCell = analysisModel().domainManager().giveDomainCell(i);
        Numerics* numerics = analysisModel().domainManager().giveNumericsFor(curCell);
        
        std::vector<Dof*> rowDof, colDof;
        RealVector coefVal;

        std::tie(rowDof,colDof,coefVal) = numerics->giveStaticCoefficientMatrixAt(curCell, stage, UNASSIGNED, time);

        for ( int j = 0; j < (int)rowDof.size(); j++ )
        {
            if ( rowDof[j] && colDof[j] )
            {
                int rowNum = analysisModel().dofManager().giveEquationNumberAt(rowDof[j]);
                int colNum = analysisModel().dofManager().giveEquationNumberAt(colDof[j]);

                if ( rowNum != UNASSIGNED && colNum != UNASSIGNED )
                    _spMatrix->atomicAddToComponent(rowNum, colNum, coefVal(j));
            }
        }
    }

    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    diagnostics().addCoefMatAssemblyTime(tictoc.count());
}
// ---------------------------------------------------------------------------
RealVector NewtonRaphson::assembleRightHandSide( int stage
                                               , const std::vector<BoundaryCondition>& bndCond
                                               , const std::vector<FieldCondition>& fldCond
                                               , const TimeData& time )
{
    std::chrono::time_point<std::chrono::system_clock> tic, toc;
    std::chrono::duration<double> tictoc;
    tic = std::chrono::high_resolution_clock::now();

    RealVector rhs(_nUnknowns);
    
    // Loop through all field conditions
    int nCells = analysisModel().domainManager().giveNumberOfDomainCells();
    
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        int threadNum = 0;
#ifdef _OPENMP
        threadNum = omp_get_thread_num();
#endif
#ifdef _OPENMP
#pragma omp for
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
                    std::tie(rowDof,localRhs) = numerics->giveStaticRightHandSideAt(curCell, stage, UNASSIGNED, fldCond[ifc], time);

                    for ( int j = 0; j < localRhs.dim(); j++)
                    {
                        if ( rowDof[j] )
                        {
                            int rowNum = analysisModel().dofManager().giveEquationNumberAt(rowDof[j]);
                            int dofGrp = analysisModel().dofManager().giveGroupNumberFor(rowDof[j]);
                            int idx = this->giveIndexForDofGroup(dofGrp);
                            _convergenceCriterion[idx]->processLocalResidualContribution(localRhs(j), threadNum);

                            if ( rowNum != UNASSIGNED )
                            {
#ifdef _OPENMP
#pragma omp atomic
#endif
                                rhs(rowNum) += localRhs(j);
                            }

                            analysisModel().dofManager().addToSecondaryVariableAt(rowDof[j], localRhs(j));
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
#pragma omp parallel
#endif
        {
            int threadNum = 0;
#ifdef _OPENMP
            threadNum = omp_get_thread_num();
#endif
#ifdef _OPENMP
#pragma omp for
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
                    std::tie(rowDof,localRhs) = numerics->giveStaticRightHandSideAt(curCell, stage, UNASSIGNED, bndCond[ibc], time);
                    
                    for ( int j = 0; j < localRhs.dim(); j++)
                    {
                        if ( rowDof[j] )
                        {
                            int rowNum = analysisModel().dofManager().giveEquationNumberAt(rowDof[j]);
                            int dofGrp = analysisModel().dofManager().giveGroupNumberFor(rowDof[j]);
                            int idx = this->giveIndexForDofGroup(dofGrp);
                            _convergenceCriterion[idx]->processLocalResidualContribution(localRhs(j), threadNum);
                            
                            if ( rowNum != UNASSIGNED )
                            {
#ifdef _OPENMP
#pragma omp atomic
#endif
                                rhs(rowNum) += localRhs(j);
                            }

                            analysisModel().dofManager().addToSecondaryVariableAt(rowDof[j], localRhs(j));
                        }
                    }
                }
            }
        }
    }
    
    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    diagnostics().addRhsAssemblyTime(tictoc.count());
    
    return rhs;
}
// ---------------------------------------------------------------------------
int NewtonRaphson::giveIndexForDofGroup( int dofGroupNum )
{
    int idx = -1;
    for ( int i = 0; i < _nDofGroups; i++ )
        if ( dofGroupNum == _dofGrpNum[i] )
            idx = i;
    
    return idx;
}
