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
#include "MeshReaders/MeshReader.hpp"
#include "Numerics/Numerics.hpp"
#include "SparseMatrix/SparseMatrix.hpp"
#include "Util/linearAlgebra.hpp"
#include "Util/readOperations.hpp"
#include "Util/reductions.hpp"

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

    for ( int i = 0; i < _nDofGroups; i++ )
        delete _convergenceCriterion[i];
}

// Public methods
// ---------------------------------------------------------------------------
int AlternateMinimization::computeSolutionFor( int stage
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
    tic = std::chrono::high_resolution_clock::now();
    this->imposeConstraintsAt(stage, bndCond, time);
    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    std::printf("done (time = %f sec.)\n", tictoc.count());
    
    // Initialize global vectors
    tic = std::chrono::high_resolution_clock::now();
    std::vector<RealVector> dU(_nSubsystems, RealVector());
    auto resid = dU;
    auto initResid = dU;
    auto rhs = dU;
    auto lhs = dU;

    for ( int i = 0; i < _nSubsystems; i++)
        dU[i].init(_nUnknowns[i]);
    
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
        // Note: We calculate global right hand sides in each
        // iteration to accommodate numerics classes that implement special
        // methods for imposing constraints
        for ( int i = 0; i < _nSubsystems; i++)
        {
            rhs[i] = this->assembleRightHandSide(stage, _subsysNum[i], bndCond, fldCond, time);
            lhs[i] = this->assembleLeftHandSide(stage, _subsysNum[i], time);
            resid[i] = rhs[i] - lhs[i];
        }
        
        if ( iterCount == 0 )
            std::printf("\n    Initial guess (LS # %d, SS # %d)", _loadStep->giveLoadStepNum(), _substepCount);
        else
            std::printf("\n    Iter.No. %-4d (LS # %d, SS # %d)", iterCount, _loadStep->giveLoadStepNum(), _substepCount);
        std::printf("\n   -----------------");
        
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
            {
                curNumerics->printPostIterationMessage(stage);
            }
        
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
                innertic = std::chrono::high_resolution_clock::now();
                
                _solver[i]->allocateInternalMemoryFor(_spMatrix[i]);
                if ( _solver[i]->takesInitialGuess() )
                    _solver[i]->setInitialGuessTo(dU[i]);
                dU[i] = _solver[i]->solve(_spMatrix[i], resid[i]);
                
                innertoc = std::chrono::high_resolution_clock::now();
                tictoc = innertoc - innertic;
                diagnostics().addSolveTime(tictoc.count());
                
                innertic = std::chrono::high_resolution_clock::now();
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

                innertoc = std::chrono::high_resolution_clock::now();
                tictoc = innertoc - innertic;
                diagnostics().addUpdateTime(tictoc.count());
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
        std::string key = getStringInputFrom(fp, "Failed to read linear solver from input file!", _name);
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
RealVector AlternateMinimization::assembleLeftHandSide( int stage
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
            
            // Assembly
            for ( int j = 0; j < (int)rowDof.size(); j++ )
            {
                if ( rowDof[j] )
                {
                    int rowNum = analysisModel().dofManager().giveEquationNumberAt(rowDof[j]);
                    int ssNum = analysisModel().dofManager().giveSubsystemNumberFor(rowDof[j]);

                    int dofGrp = analysisModel().dofManager().giveGroupNumberFor(rowDof[j]);
                    int idx = this->giveIndexForDofGroup(dofGrp);
                    _convergenceCriterion[idx]->processLocalResidualContribution(localLhs(j), threadNum);

                    if ( rowNum != UNASSIGNED && ssNum == subsys )
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
void AlternateMinimization::assembleJacobian( int stage
                                            , int subsys
                                            , const TimeData& time )
{
    std::chrono::time_point<std::chrono::system_clock> tic, toc;
    std::chrono::duration<double> tictoc;
    tic = std::chrono::high_resolution_clock::now();

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

    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    diagnostics().addCoefMatAssemblyTime(tictoc.count());
}
// ---------------------------------------------------------------------------
RealVector AlternateMinimization::assembleRightHandSide( int stage
                                                       , int subsys
                                                       , const std::vector<BoundaryCondition>& bndCond
                                                       , const std::vector<FieldCondition>& fldCond
                                                       , const TimeData& time )
{
    std::chrono::time_point<std::chrono::system_clock> tic, toc;
    std::chrono::duration<double> tictoc;
    tic = std::chrono::high_resolution_clock::now();

    // Get subsystem index
    int idx = this->giveIndexForSubsystem(subsys);
    RealVector rhs(_nUnknowns[idx]);
    
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
                    std::tie(rowDof,localRhs) = numerics->giveStaticRightHandSideAt(curCell, stage, subsys, fldCond[ifc], time);
                    
                    for ( int j = 0; j < localRhs.dim(); j++)
                    {
                        if ( rowDof[j] )
                        {
                            int rowNum = analysisModel().dofManager().giveEquationNumberAt(rowDof[j]);
                            int ssNum = analysisModel().dofManager().giveSubsystemNumberFor(rowDof[j]);
                            
                            int dofGrp = analysisModel().dofManager().giveGroupNumberFor(rowDof[j]);
                            int idx = this->giveIndexForDofGroup(dofGrp);
                            _convergenceCriterion[idx]->processLocalResidualContribution(localRhs(j), threadNum);

                            if ( rowNum != UNASSIGNED && ssNum == subsys )
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
                    std::tie(rowDof,localRhs) = numerics->giveStaticRightHandSideAt(curCell, stage, subsys, bndCond[ibc], time);

                    for ( int j = 0; j < localRhs.dim(); j++)
                    {
                        if ( rowDof[j] )
                        {
                            int rowNum = analysisModel().dofManager().giveEquationNumberAt(rowDof[j]);
                            int ssNum = analysisModel().dofManager().giveSubsystemNumberFor(rowDof[j]);

                            int dofGrp = analysisModel().dofManager().giveGroupNumberFor(rowDof[j]);
                            int idx = this->giveIndexForDofGroup(dofGrp);
                            _convergenceCriterion[idx]->processLocalResidualContribution(localRhs(j), threadNum);

                            if ( rowNum != UNASSIGNED && ssNum == subsys )
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