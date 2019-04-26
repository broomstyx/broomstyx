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

#include "LinearTransient.hpp"
#include <chrono>
#include <omp.h>
#include <stdexcept>
#include <tuple>

#include "../Core/AnalysisModel.hpp"
#include "../Core/ObjectFactory.hpp"
#include "../Core/DofManager.hpp"
#include "../Core/DomainManager.hpp"
#include "../Core/NumericsManager.hpp"
#include "../Core/SolutionManager.hpp"
#include "../LinearSolvers/LinearSolver.hpp"
#include "../MeshReaders/MeshReader.hpp"
#include "../Numerics/Numerics.hpp"
#include "../Util/linearAlgebra.hpp"
#include "../Util/readOperations.hpp"

using namespace broomstyx;

registerBroomstyxObject(SolutionMethod, LinearTransient)

// Constructor
LinearTransient::LinearTransient()
{
    _name = "LinearTransient";
}

// Destructor
LinearTransient::~LinearTransient() {}

// Public methods
// ---------------------------------------------------------------------------
void LinearTransient::formSparsityProfileForStage( int stage )
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
    for (int i = 0; i < (int)dof.size(); i++)
        analysisModel().dofManager().setEquationNumberFor(dof[i], i);

    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    std::printf("done (time = %f sec.)", tictoc.count());

    std::printf("\n  %-40s", "Determining sparsity pattern ...");
    std::fflush(stdout);
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
        std::tie(rowDof,colDof,std::ignore) = numerics->giveStaticCoefficientMatrixAt(curCell, stage, UNASSIGNED, dummyTime);

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

        std::tie(rowDof,colDof,std::ignore) = numerics->giveTransientCoefficientMatrixAt(curCell, stage, UNASSIGNED, dummyTime);

        lnnz = rowDof.size();

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

// Private methods
// -----------------------------------------------------------------------------
void LinearTransient::assembleEquations( int stage
                                       , const std::vector<BoundaryCondition>& bndCond
                                       , const std::vector<FieldCondition>& fldCond
                                       , const TimeData& time
                                       , RealVector& rhs )
{    
    int nCells = analysisModel().domainManager().giveNumberOfDomainCells();

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int iCell = 0; iCell < nCells; iCell++ )
    {
        Cell* curCell = analysisModel().domainManager().giveDomainCell(iCell);
        int label = analysisModel().domainManager().giveLabelOf(curCell);
        Numerics* numerics = analysisModel().domainManager().giveNumericsForDomain(label);

        // Calculate local static coefficient matrix
        std::vector<Dof*> rowDof, colDof;
        RealVector coefVal;
        std::tie(rowDof,colDof,coefVal) = numerics->giveStaticCoefficientMatrixAt(curCell, stage, UNASSIGNED, time);

        for ( int i = 0; i < (int)rowDof.size(); i++ )
        {
            if ( rowDof[i] && colDof[i] )
            {
                int rowNum, colNum;
                rowNum = analysisModel().dofManager().giveEquationNumberAt(rowDof[i]);
                colNum = analysisModel().dofManager().giveEquationNumberAt(colDof[i]);

                // Assembly to global coefficient matrix
                if ( rowNum != UNASSIGNED && colNum != UNASSIGNED )
                    _spMatrix->atomicAddToComponent(rowNum, colNum, coefVal(i)*time.increment);

                // Right hand side contribution arising from constraints
                if ( colNum == UNASSIGNED && rowNum != UNASSIGNED )
                {
                    double constraintVal = analysisModel().dofManager().giveValueOfConstraintAt(colDof[i], current_value);
                    double fmatVal = coefVal(i)*constraintVal*time.increment;
#ifdef _OPENMP
#pragma omp atomic
#endif
                    rhs(rowNum) -= fmatVal;
                }
            }
        }

        // Calculate local transient coefficient matrix and rhs terms
        std::tie(rowDof,colDof,coefVal) = numerics->giveTransientCoefficientMatrixAt(curCell, stage, UNASSIGNED, time);

        for ( int i = 0; i < (int)rowDof.size(); i++ )
        {
            if ( rowDof[i] && colDof[i] )
            {
                int rowNum, colNum;
                rowNum = analysisModel().dofManager().giveEquationNumberAt(rowDof[i]);
                colNum = analysisModel().dofManager().giveEquationNumberAt(colDof[i]);

                // Assembly to global coefficient matrix
                if ( rowNum != UNASSIGNED && colNum != UNASSIGNED )
                    _spMatrix->atomicAddToComponent(rowNum, colNum, coefVal(i));

                // Right hand side contribution arising from constraints
                if ( colNum == UNASSIGNED && rowNum != UNASSIGNED )
                {
                    double constraintVal = analysisModel().dofManager().giveValueOfConstraintAt(colDof[i], current_value);
                    double fmatVal = coefVal(i)*constraintVal;
#ifdef _OPENMP
#pragma omp atomic
#endif
                    rhs(rowNum) -= fmatVal;
                }
            }
        }

        // Local right hand side arising from transient terms
        RealVector localRhs;
        std::tie(rowDof,localRhs) = numerics->giveTransientLeftHandSideAt(curCell, stage, UNASSIGNED, time, converged_value);

        for ( int i = 0; i < (int)rowDof.size(); i++ )
        {
            if ( rowDof[i] )
            {
                int rowNum = analysisModel().dofManager().giveEquationNumberAt(rowDof[i]);

                // Assemble to global right hand side
                if ( rowNum != UNASSIGNED )
    #ifdef _OPENMP
    #pragma omp atomic
    #endif
                    rhs(rowNum) += localRhs(i);
            }
        }
    }

    // Assemble right hand side contribution of natural boundary conditions
    // and source/sink terms
    this->assembleRightHandSide(stage, bndCond, fldCond, time, rhs);
}
// ----------------------------------------------------------------------------
void LinearTransient::assembleRightHandSide( int stage
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
                            rhs(rowNum) += localRhs(i)*time.increment;
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
                    int rowNum = analysisModel().dofManager().giveEquationNumberAt(rowDof[i]);
                    if ( rowNum != UNASSIGNED )
                    {
#ifdef _OPENMP
#pragma omp atomic
#endif
                        rhs(rowNum) += localRhs(i)*time.increment;
                    }
                }
            }                
        }
    }
}