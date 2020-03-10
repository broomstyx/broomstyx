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

#include "TransientAlternateMinimization.hpp"
#include <chrono>
#include <cmath>
#include <omp.h>
#include <stdexcept>
#include <tuple>

#include "Core/AnalysisModel.hpp"
#include "Core/ObjectFactory.hpp"
#include "Core/DofManager.hpp"
#include "Core/DomainManager.hpp"
#include "Core/LoadStep.hpp"
#include "Core/NumericsManager.hpp"
#include "Core/SolutionManager.hpp"
#include "LinearSolvers/LinearSolver.hpp"
#include "MeshReaders/MeshReader.hpp"
#include "Numerics/Numerics.hpp"
#include "SparseMatrix/SparseMatrix.hpp"
#include "Util/linearAlgebra.hpp"
#include "Util/readOperations.hpp"
#include "Util/reductions.hpp"

using namespace broomstyx;

registerBroomstyxObject(SolutionMethod, TransientAlternateMinimization)

// Constructor
TransientAlternateMinimization::TransientAlternateMinimization()
{
    _name = "TransientAlternateMinimization";
}

// Destructor
TransientAlternateMinimization::~TransientAlternateMinimization()
{}

// Public methods
// ----------------------------------------------------------------------------
void TransientAlternateMinimization::formSparsityProfileForStage( int stage )
{
    std::chrono::time_point<std::chrono::system_clock> tic, toc;
    std::chrono::duration<double> tictoc;

    // Count number of unknowns for each DOF group
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
        for ( int j = 0; j < nCells; j++)
        {
            Cell* curCell = analysisModel().domainManager().giveDomainCell(j);
            int label = analysisModel().domainManager().giveLabelOf(curCell);
            Numerics* numerics = analysisModel().domainManager().giveNumericsForDomain(label);
            
            TimeData time;
            std::vector<Dof*> rowDof, colDof;
            
            // A. Static part of coefficient matrix
            std::tie(rowDof,colDof,std::ignore) = numerics->giveStaticCoefficientMatrixAt(curCell, stage, _subsysNum[i], time);

            int lnnz = rowDof.size();
            for ( int k = 0; k < lnnz; k++ )
            {
                if ( rowDof[k] && colDof[k] )
                {
                    // Check that DOF group is included in solution
                    int rgrpNum = analysisModel().dofManager().giveGroupNumberFor(rowDof[k]);
                    int rgrpIdx = this->giveIndexForDofGroup(rgrpNum);
                    int cgrpNum = analysisModel().dofManager().giveGroupNumberFor(colDof[k]);
                    int cgrpIdx = this->giveIndexForDofGroup(cgrpNum);

                    if ( rgrpIdx >= 0 && cgrpIdx >= 0 )
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

            // B. Transient part of coefficient matrix
            std::tie(rowDof,colDof,std::ignore) = numerics->giveTransientCoefficientMatrixAt(curCell, stage, _subsysNum[i], time);

            lnnz = rowDof.size();
            for ( int k = 0; k < lnnz; k++ )
            {
                if ( rowDof[k] && colDof[k] )
                {
                    // Check that DOF group is included in solution
                    int rgrpNum = analysisModel().dofManager().giveGroupNumberFor(rowDof[k]);
                    int rgrpIdx = this->giveIndexForDofGroup(rgrpNum);
                    int cgrpNum = analysisModel().dofManager().giveGroupNumberFor(colDof[k]);
                    int cgrpIdx = this->giveIndexForDofGroup(cgrpNum);

                    if ( rgrpIdx >= 0 && cgrpIdx >= 0 )
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

// Private methods
// --------------------------------------------------------------------------
RealVector TransientAlternateMinimization::assembleLeftHandSide( int stage
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
        // A. Static part
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
                    sumAbsFlux(dgIdx) += std::fabs(localLhs(j)*time.increment);
                    fluxCount(dgIdx) += 1.;
#ifdef _OPENMP
#pragma omp atomic
#endif
                    lhs(rowNum) += localLhs(j)*time.increment;
                }
            }
        }
        
        /***************************************/

        // B. Transient part
        RealVector tLhsOld, tLhsNew;
        std::tie(rowDof,tLhsNew) = numerics->giveTransientLeftHandSideAt(curCell, stage, subsys, time, current_value);
        std::tie(rowDof,tLhsOld) = numerics->giveTransientLeftHandSideAt(curCell, stage, subsys, time, converged_value);
        localLhs = tLhsNew - tLhsOld;

        // Assembly
        for ( int j = 0; j < (int)rowDof.size(); j++)
        {
            if ( rowDof[j] )
            {
                int rowNum = analysisModel().dofManager().giveEquationNumberAt(rowDof[j]);
                int ssNum = analysisModel().dofManager().giveSubsystemNumberFor(rowDof[j]);

                if ( rowNum != UNASSIGNED && ssNum == subsys )
                {
                    int dgNum = analysisModel().dofManager().giveGroupNumberFor(rowDof[j]);
                    int dgIdx = this->giveIndexForDofGroup(dgNum);
                    sumAbsFlux(dgIdx) += std::fabs(localLhs(j));
                    fluxCount(dgIdx) += 1.;
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
void TransientAlternateMinimization::assembleJacobian( int stage
                                                     , int subsys
                                                     , const TimeData& time )
{
    int nCells = analysisModel().domainManager().giveNumberOfDomainCells();

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < nCells; i++ )
    {
        Cell* curCell = analysisModel().domainManager().giveDomainCell(i);
        int label = analysisModel().domainManager().giveLabelOf(curCell);
        Numerics* numerics = analysisModel().domainManager().giveNumericsForDomain(label);

        std::vector<Dof*> rowDof, colDof;
        RealVector coefVal;

        // A. Static terms
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
                    _spMatrix[idx]->atomicAddToComponent(rowNum, colNum, coefVal(j)*time.increment);
                }
            }
        }

        // B. Transient terms
        std::tie(rowDof,colDof,coefVal) = numerics->giveTransientCoefficientMatrixAt(curCell,stage,subsys,time);

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
RealVector TransientAlternateMinimization::assembleRightHandSide( int stage
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

                for ( int j = 0; j < localRhs.dim(); j++ )
                {
                    if ( rowDof[j] )
                    {
                        int rowNum = analysisModel().dofManager().giveEquationNumberAt(rowDof[j]);
                        int ssNum = analysisModel().dofManager().giveSubsystemNumberFor(rowDof[j]);

                        if ( rowNum != UNASSIGNED && ssNum == subsys )
                        {
                            int dgNum = analysisModel().dofManager().giveGroupNumberFor(rowDof[j]);
                            int dgIdx = this->giveIndexForDofGroup(dgNum);
                            sumAbsFlux(dgIdx) += std::fabs(localRhs(j)*time.increment);
                            fluxCount(dgIdx) += 1.;
#ifdef _OPENMP
#pragma omp atomic
#endif
                            rhs(rowNum) += localRhs(j)*time.increment;
                        }
                    }
                }
            }
        }
    }

    // Loop through all natural boundary conditions
    for ( int ibc = 0; ibc < (int)bndCond.size(); ibc++ )
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
                            sumAbsFlux(dgIdx) += std::fabs(localRhs(j)*time.increment);
                            fluxCount(dgIdx) += 1.;
#ifdef _OPENMP
#pragma omp atomic
#endif
                            rhs(rowNum) += localRhs(j)*time.increment;
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