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

#include "LinearStatic.hpp"
#include <chrono>
#include <omp.h>
#include <stdexcept>
#include <tuple>

#include "Core/AnalysisModel.hpp"
#include "Core/ObjectFactory.hpp"
#include "Core/Diagnostics.hpp"
#include "Core/DofManager.hpp"
#include "Core/DomainManager.hpp"
#include "Core/NumericsManager.hpp"
#include "Core/SolutionManager.hpp"
#include "LinearSolvers/LinearSolver.hpp"
#include "MeshReaders/MeshReader.hpp"
#include "Numerics/Numerics.hpp"
#include "Util/linearAlgebra.hpp"
#include "Util/readOperations.hpp"

using namespace broomstyx;

registerBroomstyxObject(SolutionMethod, LinearStatic)

// Constructor
LinearStatic::LinearStatic()
{
    _name = "LinearStatic";
}

// Destructor
LinearStatic::~LinearStatic() 
{
    delete _spMatrix;
    delete _solver;
}

// Public methods
// ---------------------------------------------------------------------------
int LinearStatic::computeSolutionFor( int stage
                                    , const std::vector<BoundaryCondition>& bndCond
                                    , const std::vector<FieldCondition>& fldCond
                                    , const TimeData& time )
{
    std::chrono::time_point<std::chrono::system_clock> tic, toc;
    std::chrono::duration<double> tictoc;

    // Impose constraints on nodal DOFs
    std::printf("    %-40s", "Imposing constraints ...");
    tic = std::chrono::high_resolution_clock::now();
    this->imposeConstraintsAt(stage, bndCond, time);
    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    std::printf("done (time = %f sec.)\n", tictoc.count());
    diagnostics().addSetupTime(tictoc.count());
    
    // Initialize global coefficient matrix
    _spMatrix->initializeValues();
    
    // Initialize global right hand side
    RealVector sysU, sysRHS;
    int nUnknowns;
    std::tie(nUnknowns, std::ignore) = _spMatrix->giveMatrixDimensions();
    sysRHS.init(nUnknowns);
    sysU.init(nUnknowns);
    
    // Assemble global equations
    std::printf("    %-40s", "Assembling equations ...");
    tic = std::chrono::high_resolution_clock::now();
    this->assembleEquations(stage, bndCond, fldCond, time, sysRHS);

    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    std::printf("done (time = %f sec.)\n", tictoc.count());
    
    // Solve system
    std::printf("    %-40s", "Solving system ...");
    tic = std::chrono::high_resolution_clock::now();
    _solver->allocateInternalMemoryFor(_spMatrix);
    if ( _solver->takesInitialGuess() )
        _solver->setInitialGuessTo(sysU);
    sysU = _solver->solve(_spMatrix, sysRHS);
    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    std::printf("done (time = %f sec.)\n", tictoc.count());
    diagnostics().addSolveTime(tictoc.count());

    // Update DOF values
    std::printf("    %-40s", "Updating DOF values ...");
    tic = std::chrono::high_resolution_clock::now();

    std::vector<Dof*> dof = analysisModel().dofManager().giveActiveDofsAtStage(stage);

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < (int)dof.size(); i++ )
    {
        int eqNo = analysisModel().dofManager().giveEquationNumberAt(dof[i]);
        analysisModel().dofManager().updatePrimaryVariableAt(dof[i], sysU(eqNo), current_value);
    }

    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    std::printf("done (time = %f sec.)\n", tictoc.count());
    diagnostics().addUpdateTime(tictoc.count());

    return 0;
}
// ---------------------------------------------------------------------------
void LinearStatic::formSparsityProfileForStage( int stage )
{
    std::chrono::time_point<std::chrono::system_clock> tic, toc;
    std::chrono::duration<double> tictoc;

    // Assign DOF equation numbers
    std::printf("\n  %-40s", "Assigning DOF equation numbers ...");
    std::fflush(stdout);
    tic = std::chrono::high_resolution_clock::now();

    std::vector<Dof*> dof = analysisModel().dofManager().giveActiveDofsAtStage(stage);

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < (int)dof.size(); i++ )
        analysisModel().dofManager().setEquationNumberFor(dof[i], i);

    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    std::printf("done (time = %f sec.)", tictoc.count());

    std::printf("\n  %-40s", "Determining sparsity pattern ...");
    tic = std::chrono::high_resolution_clock::now();

    // Get sparse matrix size
    int nvar = analysisModel().dofManager().giveNumberOfActiveDofsAtStage(stage);

    _spMatrix->initializeProfile(nvar, nvar);
    _spMatrix->setSymmetryTo(_solver->giveSymmetryOption());

    int nCells = analysisModel().domainManager().giveNumberOfDomainCells();

    for ( int i = 0; i < nCells; i++ )
    {
        Cell* curCell = analysisModel().domainManager().giveDomainCell(i);
        int label = analysisModel().domainManager().giveLabelOf(curCell);
        Numerics* numerics = analysisModel().domainManager().giveNumericsForDomain(label);

        // Dummy variables
        TimeData dummyTime;

        std::vector<Dof*> rowDof, colDof;
        std::tie(rowDof,colDof,std::ignore) = numerics->giveStaticCoefficientMatrixAt(curCell, stage,UNASSIGNED,dummyTime);

        int lnnz = rowDof.size();

        for ( int j = 0; j < lnnz; j++ )
        {
            if ( rowDof[j] && colDof[j] )
            {
                int rowNum, colNum;
                rowNum = analysisModel().dofManager().giveEquationNumberAt(rowDof[j]);
                colNum = analysisModel().dofManager().giveEquationNumberAt(colDof[j]);
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
    double spratio = (double)nnz/((double)nvar*(double)nvar);
    std::printf("\n    Stage %-2d: nUnknowns = %d, nnz = %d", stage, nvar, nnz);
    std::printf("\n         Sparsity ratio = %.4e\n", spratio);
}
// ---------------------------------------------------------------------------
void LinearStatic::initializeSolvers()
{
    _solver->initialize();
}
// ---------------------------------------------------------------------------
void LinearStatic::readDataFromFile( FILE* fp )
{
    std::string key;

    verifyKeyword(fp, key = "LinearSolver", _name);
    key = getStringInputFrom(fp, "Failed to read linear solver from input file!", _name);
    _solver = objectFactory().instantiateLinearSolver(key);
    _solver->readDataFrom(fp);
    
    // Initialize sparse matrix
    _spMatrix = objectFactory().instantiateSparseMatrix(_solver->giveRequiredMatrixFormat());
}

// Private methods
// -----------------------------------------------------------------------------
void LinearStatic::assembleEquations( int stage
                                    , const std::vector<BoundaryCondition>& bndCond
                                    , const std::vector<FieldCondition>& fldCond
                                    , const TimeData& time
                                    , RealVector& rhs )
{    
    std::chrono::time_point<std::chrono::system_clock> tic, toc;
    std::chrono::duration<double> tictoc;
    
    int nCells = analysisModel().domainManager().giveNumberOfDomainCells();
    tic = std::chrono::high_resolution_clock::now();

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int iCell = 0; iCell < nCells; iCell++ )
    {
        Cell* curCell = analysisModel().domainManager().giveDomainCell(iCell);
        int label = analysisModel().domainManager().giveLabelOf(curCell);
        Numerics* numerics = analysisModel().domainManager().giveNumericsForDomain(label);

        // Calculate local coefficient matrix
        std::vector<Dof*> rowDof, colDof;
        RealVector coefVal;
        std::tie(rowDof,colDof,coefVal) = numerics->giveStaticCoefficientMatrixAt(curCell, stage, UNASSIGNED, time);

        int lnnz = rowDof.size();

        for ( int j = 0; j < lnnz; j++)
        {
            if ( rowDof[j] && colDof[j] )
            {
                int rowNum, colNum;
                rowNum = analysisModel().dofManager().giveEquationNumberAt(rowDof[j]);
                colNum = analysisModel().dofManager().giveEquationNumberAt(colDof[j]);

                // Assembly to global coefficient matrix
                if ( rowNum != UNASSIGNED && colNum != UNASSIGNED )
                    _spMatrix->atomicAddToComponent(rowNum, colNum, coefVal(j));

                // Right hand side contribution arising from constraints
                if ( colNum == UNASSIGNED && rowNum != UNASSIGNED )
                {
                    double fmatVal = coefVal(j)*analysisModel().dofManager().giveValueOfConstraintAt(colDof[j], current_value);
#ifdef _OPENMP
#pragma omp atomic
#endif
                    rhs(rowNum) -= fmatVal;
                }
            }
        }
    }
    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    diagnostics().addCoefMatAssemblyTime(tictoc.count());

    // Assemble right hand side contribution of natural boundary conditions
    // and source/sink terms
    tic = std::chrono::high_resolution_clock::now();
    this->assembleRightHandSide(stage, bndCond, fldCond, time, rhs);
    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    diagnostics().addRhsAssemblyTime(tictoc.count());
}
// ---------------------------------------------------------------------------
void LinearStatic::assembleRightHandSide( int stage
                                        , const std::vector<BoundaryCondition>& bndCond
                                        , const std::vector<FieldCondition>& fldCond
                                        , const TimeData& time
                                        , RealVector& rhs )
{
    // Loop through all field conditions
    int nCells = analysisModel().domainManager().giveNumberOfDomainCells();

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int iCell = 0; iCell < nCells; iCell++)
    {
        Cell* curCell = analysisModel().domainManager().giveDomainCell(iCell);
        int label = analysisModel().domainManager().giveLabelOf(curCell);

        for (int ifc = 0; ifc < (int)fldCond.size(); ifc++)
        {
            int fcLabel = analysisModel().domainManager().givePhysicalEntityNumberFor(fldCond[ifc].domainLabel());
            if ( label == fcLabel )
            {
                RealVector localRhs;
                std::vector<Dof*> rowDof;

                Numerics* numerics = analysisModel().domainManager().giveNumericsForDomain(label);
                std::tie(rowDof,localRhs) = numerics->giveStaticRightHandSideAt(curCell, stage, UNASSIGNED, fldCond[ifc], time);

                for ( int i = 0; i < localRhs.dim(); i++)
                {
                    if ( rowDof[i] )
                    {
                        int rowNum = analysisModel().dofManager().giveEquationNumberAt(rowDof[i]);
                        if ( rowNum != UNASSIGNED )
                        {
#ifdef _OPENMP
#pragma omp atomic
#endif
                            rhs(rowNum) += localRhs(i);
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
#pragma omp parallel for
#endif
        for ( int iCell = 0; iCell < nBCells; iCell++ )
        {
            Cell* curCell = analysisModel().domainManager().giveBoundaryCell(iCell);
            int label = analysisModel().domainManager().giveLabelOf(curCell);

            if ( label == boundaryId )
            {
                RealVector localRhs;
                std::vector<Dof*> rowDof;

                // Specifics of boundary condition imposition are handled by numerics
                std::tie(rowDof,localRhs) = numerics->giveStaticRightHandSideAt(curCell, stage, UNASSIGNED, bndCond[ibc], time);

                for ( int i = 0; i < localRhs.dim(); i++)
                {
                    if ( rowDof[i])
                    {
                        int rowNum = analysisModel().dofManager().giveEquationNumberAt(rowDof[i]);
                        if ( rowNum != UNASSIGNED )
                        {
#ifdef _OPENMP
#pragma omp atomic
#endif
                            rhs(rowNum) += localRhs(i);
                        }
                    }
                }
            }                
        }
    }
}