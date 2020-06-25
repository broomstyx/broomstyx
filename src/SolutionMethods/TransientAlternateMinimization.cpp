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

#include "ConvergenceCriteria/ConvergenceCriterion.hpp"
#include "Core/AnalysisModel.hpp"
#include "Core/Diagnostics.hpp"
#include "Core/DofManager.hpp"
#include "Core/DomainManager.hpp"
#include "Core/LoadStep.hpp"
#include "Core/NumericsManager.hpp"
#include "Core/ObjectFactory.hpp"
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
    std::chrono::time_point<std::chrono::system_clock> tic, toc;
    std::chrono::duration<double> tictoc;
    tic = std::chrono::high_resolution_clock::now();

    // Find subsystem index
    int idx = this->giveIndexForSubsystem(subsys);
    
    RealVector lhs(_nUnknowns[idx]);
    
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
        for ( int i = 0; i < nCells; i++ )
        {
            Cell* curCell = analysisModel().domainManager().giveDomainCell(i);
            Numerics* numerics = analysisModel().domainManager().giveNumericsFor(curCell);
            
            std::vector<Dof*> rowDof;
            RealVector localLhs;
            
            // Calculate cell internal forces
            std::tie(rowDof,localLhs) = numerics->giveStaticLeftHandSideAt(curCell, stage, subsys, time);
            std::vector<int> localDofGrp(localLhs.dim(), -1);

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
                        localDofGrp[j] = analysisModel().dofManager().giveGroupNumberFor(rowDof[j]);
#ifdef _OPENMP
#pragma omp atomic
#endif
                        lhs(rowNum) += localLhs(j);
                    }

                    analysisModel().dofManager().addToSecondaryVariableAt(rowDof[j], localLhs(j));
                }
            }
            _convergenceCriterion->processLocalResidualContribution(localLhs, localDofGrp, threadNum);
            
            /***************************************/

            // B. Transient part
            RealVector tLhsOld, tLhsNew;
            std::tie(rowDof,tLhsNew) = numerics->giveTransientLeftHandSideAt(curCell, stage, subsys, time, current_value);
            std::tie(rowDof,tLhsOld) = numerics->giveTransientLeftHandSideAt(curCell, stage, subsys, time, converged_value);
            localLhs = (tLhsNew - tLhsOld);
            localDofGrp.assign(localLhs.dim(), -1);

            // Assembly
            for ( int j = 0; j < (int)rowDof.size(); j++)
            {
                if ( rowDof[j] )
                {
                    int rowNum = analysisModel().dofManager().giveEquationNumberAt(rowDof[j]);
                    int ssNum = analysisModel().dofManager().giveSubsystemNumberFor(rowDof[j]);

                    if ( rowNum != UNASSIGNED && ssNum == subsys )
                    {
                        localDofGrp[j] = analysisModel().dofManager().giveGroupNumberFor(rowDof[j]);
#ifdef _OPENMP
#pragma omp atomic
#endif
                        lhs(rowNum) += localLhs(j)/time.increment;
                    }

                    analysisModel().dofManager().addToSecondaryVariableAt(rowDof[j], localLhs(j)/time.increment);
                }
            }
            _convergenceCriterion->processLocalResidualContribution(localLhs, localDofGrp, threadNum);
        }
    }

    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    diagnostics().addLhsAssemblyTime(tictoc.count());
    
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
                    _spMatrix[idx]->atomicAddToComponent(rowNum, colNum, coefVal(j));
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
                    _spMatrix[idx]->atomicAddToComponent(rowNum, colNum, coefVal(j)/time.increment);
                }
            }
        }
    }
}