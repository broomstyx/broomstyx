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

#include "LoadStep.hpp"
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <omp.h>
#include <stdexcept>
#include <string>

#include "AnalysisModel.hpp"
#include "Cell.hpp"
#include "DomainManager.hpp"
#include "DofManager.hpp"
#include "DomainManager.hpp"
#include "DomainManager.hpp"
#include "NumericsManager.hpp"
#include "OutputManager.hpp"
#include "SolutionManager.hpp"
#include "../MeshReaders/MeshReader.hpp"
#include "../Numerics/Numerics.hpp"
#include "../SolutionMethods/SolutionMethod.hpp"
#include "../Util/readOperations.hpp"
#include "ObjectFactory.hpp"

using namespace broomstyx;

// Constructor
LoadStep::LoadStep( int lsNum, int nStg )
    : _loadStepNum(lsNum)
    , _nStages(nStg)
{
    _name = "LoadStep";
    
    _convDatFile.assign(nStg, nullptr);
    _iterDatFile.assign(nStg, nullptr);
    _convDatCount.assign(nStg, 0);
    _iterDatCount.assign(nStg, 0);
    
    // Create directories
    const int dir_err = system("mkdir -p LoadStepData");
    if ( dir_err == -1 )
        throw std::runtime_error("Error creating directory 'LoadStepData'!\n");
    
    for ( int i = 0; i < nStg; i++)
    {
        std::string iterFilename = "./LoadStepData/LoadStep_" + std::to_string(_loadStepNum) + "_Stage_"
                + std::to_string(i+1) + "_IterCount.csv";
        
        std::string convFilename = "./LoadStepData/LoadStep_" + std::to_string(_loadStepNum) + "_Stage_"
                + std::to_string(i+1) + "_ConvDat.csv";
        
        _convDatFile[i] = std::fopen(convFilename.c_str(), "w");
        _iterDatFile[i] = std::fopen(iterFilename.c_str(), "w");
    }
}

// Destructor
LoadStep::~LoadStep()
{
    for ( int i = 0; i < _nStages; i++)
    {
        if ( _convDatFile[i] )
            std::fclose(_convDatFile[i]);
        if ( _iterDatFile[i] )
            std::fclose(_iterDatFile[i]);
    }
    
    for ( int i = 0; i < (int)_solutionMethod.size(); i++)
        if ( _solutionMethod[i] )
            delete _solutionMethod[i];
}

// Public methods
// ----------------------------------------------------------------------------
int LoadStep::giveLoadStepNum()
{
    return _loadStepNum;
}
// ----------------------------------------------------------------------------
void LoadStep::readDataFrom( FILE *fp )
{   
    // Special pre-processing
    verifyDeclaration(fp, "PREPROCESSING", _name);
    int nPreprocessing = getIntegerInputFrom(fp, "Failed reading number of special pre-processing procedures for load step # " + std::to_string(_loadStepNum), _name);
    
    _preProcess.assign(nPreprocessing, ProcessData());
    for ( int i = 0; i < nPreprocessing; i++)
    {
        _preProcess[i].domainTag = getStringInputFrom(fp, "Failed to read domain label for pre-processing in load step # " + std::to_string(_loadStepNum), _name);
        _preProcess[i].directive = getStringInputFrom(fp, "Failed to read directive for pre-processing in load step # " + std::to_string(_loadStepNum), _name);
    }
    
    // Time at start of load step
    verifyDeclaration(fp, "START_TIME", _name);
    _startTime = getRealInputFrom(fp, "Failed to read start time for load step # " + std::to_string(_loadStepNum) + " from input file!", _name);
    
    // Time at end of load step
    verifyDeclaration(fp, "END_TIME", _name);
    _endTime = getRealInputFrom(fp, "Failed to read end time for load step # " + std::to_string(_loadStepNum) + " from input file!", _name);
    
    // Initial time increment
    verifyDeclaration(fp, "INITIAL_TIME_INCREMENT", _name);
    _dtime = getRealInputFrom(fp, "Failed to read initial time increment for load step # " + std::to_string(_loadStepNum) + " from input file!", _name);
    
    // Maximum number of substeps
    verifyDeclaration(fp, "MAX_SUBSTEPS", _name);
    _maxSubsteps = getIntegerInputFrom(fp, "Failed to read maximum number of substeps for load step # " + std::to_string(_loadStepNum) + " from input file!", _name);
    
    // Read boundary conditions
    verifyDeclaration(fp, "BOUNDARY_CONDITIONS", _name);
    int nBC = getIntegerInputFrom(fp, "Failed to read number of boundary conditions for load step # " + std::to_string(_loadStepNum) + " from input file!", _name);
    
    _boundaryCondition.assign(nBC, BoundaryCondition());
    for ( int i = 0; i < nBC; i++)
        _boundaryCondition[i].readDataFrom(fp);
    
    // Read field conditions
    verifyDeclaration(fp, "FIELD_CONDITIONS", _name);
    int nFC = getIntegerInputFrom(fp, "Failed to read number of field conditions for load step # " + std::to_string(_loadStepNum) + " from input file!", _name);
    
    _fieldCondition.assign(nFC, FieldCondition());
    for ( int i = 0; i < nFC; i++)
        _fieldCondition[i].readDataFrom(fp);
    
    // Read solution methods for each stage
    verifyDeclaration(fp, "SOLUTION_METHODS", _name);
    _solutionMethod.assign(_nStages + 1, nullptr);
    for ( int i = 1; i <= _nStages; i++)
    {
        verifyKeyword(fp, "Stage", _name);
        int stg = getIntegerInputFrom(fp, "Failed to read stage number in solution method input for " + std::string("load step # ") + std::to_string(_loadStepNum), _name);
        std::string methStr = getStringInputFrom(fp, "Failed to solution method name for stage " + std::to_string(stg) + " of load step # " + std::to_string(_loadStepNum), _name);
        
        _solutionMethod[stg] = objectFactory().instantiateSolutionMethod(methStr);
        _solutionMethod[stg]->readDataFromFile(fp);
    }
    
    // Output frequency
    verifyDeclaration(fp, "WRITE_INTERVAL", _name);
    _writeInterval = getIntegerInputFrom(fp, "Failed to read number of output skips for load step # " + std::to_string(_loadStepNum), _name);
    if ( _writeInterval < 1 )
        _writeInterval = 1;
    
    // Special post-processing
    verifyDeclaration(fp, "POSTPROCESSING", _name);
    int nPostprocessing = getIntegerInputFrom(fp, "Failed reading number of special post-processing "
            + std::string("procedures for load step # ") + std::to_string(_loadStepNum), _name);
    
    _postProcess.assign(nPostprocessing, ProcessData());
    for ( int i = 0; i < nPostprocessing; i++)
    {
        _postProcess[i].domainTag = getStringInputFrom(fp, "Failed to domain label for post-processing in load step # " + std::to_string(_loadStepNum), _name);
        _postProcess[i].directive = getStringInputFrom(fp, "Failed to read directive for post-processing in load step # " + std::to_string(_loadStepNum), _name);
    }
}
// -----------------------------------------------------------------------------
void LoadStep::solveYourself()
{
    std::chrono::time_point<std::chrono::system_clock> tic, toc;
    std::chrono::duration<double> tictoc;
    
//    // HACK MESSAGES
//    FILE* fp = fopen("./HackMessages.txt","a");
    
    // Time data for immediate substep
    TimeData time;
    time.start = _startTime;
    time.end = _endTime;
    time.current = _startTime;
    time.target = _startTime + _dtime;
    time.increment = _dtime;
    // Note: time.increment can be used to implement adaptive time-stepping
    
    // Initialize solvers
    std::printf("\n  %-40s\n", "Initializing solvers ...");
    std::fflush(stdout);
    tic = std::chrono::high_resolution_clock::now();
    
    for ( int i = 0; i < (int)_solutionMethod.size(); i++ )
        if ( _solutionMethod[i] )
        {
            _solutionMethod[i]->getCurrentLoadStep();
            _solutionMethod[i]->initializeSolvers();
        }
    
    toc  = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    std::printf("  Done initializing solvers                    (time = %f sec.)\n", tictoc.count());
    
    std::printf("\n  ------------------------");
    std::printf("\n    LOADSTEP # %d", _loadStepNum);
    std::printf("\n  ------------------------\n");
    
    // Reset all constraints from previous load steps
    analysisModel().dofManager().removeAllDofConstraints();
    analysisModel().domainManager().removeAllCellConstraints();
    
    // Carry out special pre-processing procedures
    std::printf("\n  %-40s", "Running preprocesing routines ...");
    std::fflush(stdout);
    tic = std::chrono::high_resolution_clock::now();
    
    int nCells = analysisModel().domainManager().giveNumberOfDomainCells();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int iCell = 0; iCell < nCells; iCell++ )
    {
        Cell* curCell = analysisModel().domainManager().giveDomainCell(iCell);
        int cellLabel = analysisModel().domainManager().giveLabelOf(curCell);
        
        for ( int i = 0; i < (int)_preProcess.size(); i++ )
        {
            int domainLabel = analysisModel().domainManager().givePhysicalEntityNumberFor(_preProcess[i].domainTag);
            if ( cellLabel == domainLabel )
            {
                Numerics* numerics = analysisModel().domainManager().giveNumericsForDomain(domainLabel);
                numerics->performPreprocessingAt(curCell, _preProcess[i].directive);
            }
        }
    }
    
    toc  = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    std::printf("done (time = %f sec.)\n", tictoc.count());
    
    // Find constrained degrees of freedom
    std::printf("\n  %-40s", "Marking constrained DOFs ...");
    std::fflush(stdout);
    tic = std::chrono::high_resolution_clock::now();
    this->findConstrainedDofs();
    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    std::printf("done (time = %f sec.)\n", tictoc.count());
    
    // Find active DOFs
    std::printf("  %-40s", "Finding active DOFs ...");
    std::fflush(stdout);
    tic = std::chrono::high_resolution_clock::now();
    analysisModel().dofManager().findActiveDofs();
    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    std::printf("done (time = %f sec.)\n", tictoc.count());
    
    analysisModel().dofManager().reportNumberOfActiveDofs();
    
    // Construct sparse matrix profiles for each stage
    for (int i = 1; i <= _nStages; i++ )
    {
        // First impose constraints to allow numerics to set proper flags for cells where needed
        _solutionMethod[i]->imposeConstraintsAt(i, _boundaryCondition, time);
        _solutionMethod[i]->formSparsityProfileForStage(i);    
    }
    
    // ------------------------------------------------------------------
    //                       Solution proper
    
    // Cycle through substeps
    // ----------------------    
    bool endOfLoadStep = false;
    
    int curSubstep = 0;
    int skipCount = 0;
    bool forceBreak = false;
    
    while (!endOfLoadStep)
    {
//        if ( forceBreak )
//        {
//            analysisModel->dofManager->resetDofCurrentPrimaryValues();
//            forceBreak = false;
//        }
//        else
            curSubstep++;
        
        if ( curSubstep > _maxSubsteps )
            throw std::runtime_error("Maximum number of substeps exceeded!");
        
        std::printf("\n  -----------------------------------------");
        std::printf("\n    LOADSTEP # %d, Substep # %d", _loadStepNum, curSubstep);
        std::printf("\n  -----------------------------------------");
        
        tic = std::chrono::high_resolution_clock::now();
        std::printf("\n    Target time: %.14E\n", time.target);        
        
        int nStage = analysisModel().solutionManager().giveNumberOfSolutionStages();
        
        std::vector<bool> stageConverged(nStage, false);
        
        int  curSubstepIter = 0;
        bool substepConverged = false;
        forceBreak = false;
        
        do
        {
            ++curSubstepIter;
            
            std::printf("\n    Substep Iter. # %d", curSubstepIter);
            std::printf("\n  -----------------------------------------\n");
            
            for (int curStage = 1; curStage <= nStage; curStage++) {
            
                std::printf("\n    Stage # %d", curStage);
                std::printf("\n  -------------------\n");

                int error = _solutionMethod[curStage]->computeSolutionFor(curStage, _boundaryCondition, _fieldCondition, time);
                
                if ( error  == 0 )
                    stageConverged[curStage] = true;
                else if ( error == 1 )
                {
                    std::printf("\n**************************************\n");
                    std::printf("Maximum number of iterations exceeded!\n");
                    std::printf("(LoadStep # %d, Substep # %d\n\n", _loadStepNum, curSubstep);
                    forceBreak = true;
                }
                else
                {
                    std::printf("\n**************************************\n");
                    std::printf("Solution method terminated with error flag");
                    std::printf(" '%d'!\n\n", error);
                    throw std::runtime_error("Terminating LoadStep!");
                }
            }
            
            substepConverged = true;
            for (int curStage = 1; curStage <= nStage; curStage++) {
                if ( !stageConverged[curStage] )
                    substepConverged = false;
            }
            
        } while ( !substepConverged && !forceBreak );
        
//        if ( forceBreak )
//        {
//            std::fprintf(fp, "Substep %d did not converge.\n", curSubstep);
//            std::fflush(fp);
//        }
        
        if ( forceBreak )
        {
            // Perform any needed computations at cells before finalizing data
            this->performPrefinalCalculationsAtCells();
            
            // Finalize data (unconverged results)
            analysisModel().dofManager().finalizeDofPrimaryValues();
            analysisModel().domainManager().finalizeCellData();

            // Perform post-processing for nodal field values
            analysisModel().domainManager().performNodalPostProcessing();
            
            // Write unconverged results and then terminate program
            analysisModel().outputManager().writeOutput(time.target);
            std::printf("IMPORTANT: Above output contains non-converged results!!!\n\n");
            throw std::runtime_error("");
        }
        else
        {
            // Perform any needed computations at cells before finalizing data
            this->performPrefinalCalculationsAtCells();
            
            // Finalize data
            analysisModel().dofManager().finalizeDofPrimaryValues();
            analysisModel().domainManager().finalizeCellData();

            // Perform post-processing for nodal field values
            analysisModel().domainManager().performNodalPostProcessing();
            
            // Update time and target time for next substep
            time.current = time.target;
            time.target += time.increment;
            if (time.target > time.end)
                time.target = time.end;

            if ( std::fabs(time.current - time.end ) < 1.0e-13)
                endOfLoadStep = true;

            toc = std::chrono::high_resolution_clock::now();
            tictoc = toc - tic;
            std::printf("\n    Substep completed in %f sec.\n", tictoc.count());

            ++skipCount;
            if ( skipCount == _writeInterval )
            {
                analysisModel().outputManager().writeOutput(time.current);
                skipCount = 0;
            }
            else if ( endOfLoadStep )
                analysisModel().outputManager().writeOutput(time.current);

            analysisModel().outputManager().writeOutputQuantities(time.current); 
        }
    }
    
    // Carry out special post-processing procedures
    std::printf("  %-40s", "Running postprocessing routines ...");
    std::fflush(stdout);
    tic = std::chrono::high_resolution_clock::now();
    
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int iCell = 0; iCell < nCells; iCell++ )
    {
        Cell* curCell = analysisModel().domainManager().giveDomainCell(iCell);
        int cellLabel = analysisModel().domainManager().giveLabelOf(curCell);
        
        for ( int i = 0; i < (int)_postProcess.size(); i++ )
        {
            int domainLabel = analysisModel().domainManager().givePhysicalEntityNumberFor(_postProcess[i].domainTag);
            if ( cellLabel == domainLabel )
            {
                Numerics* numerics = analysisModel().domainManager().giveNumericsForDomain(domainLabel);
                numerics->performPostprocessingAt(curCell, _postProcess[i].directive);
            }
        }
    }
    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    std::printf("done (time = %f sec.)\n", tictoc.count());
    
//    // HACK MESSAGES
//    std::fclose(fp);
}
// -----------------------------------------------------------------------------
void LoadStep::writeConvergenceDataForStage( int stg, RealMatrix& convDat)
{
    int nRows = convDat.dim1();
    
    std::fprintf(_convDatFile[stg-1], "%d, ", ++_convDatCount[stg-1]);
    for (int i = 0; i < nRows; i++)
    {
        std::fprintf(_convDatFile[stg-1], "%.15e, %.15e", convDat(i,0), convDat(i,1));
        if (i == nRows - 1)
            std::fprintf(_convDatFile[stg-1], "\n");
        else
            std::fprintf(_convDatFile[stg-1], ", ");
    }
    std::fflush(_convDatFile[stg-1]);
}
// -----------------------------------------------------------------------------
void LoadStep::writeIterationDataForStage( int    stg
                                         , double time
                                         , int    nIter )
{
    std::fprintf(_iterDatFile[stg-1], "%d, ", ++_iterDatCount[stg-1]);
    std::fprintf(_iterDatFile[stg-1], "%.15e, %d\n", time, nIter);
    std::fflush(_iterDatFile[stg-1]);
}

// Private methods
// -----------------------------------------------------------------------------
void LoadStep::findConstrainedDofs()
{
    for ( int bcIdx = 0; bcIdx < (int)_boundaryCondition.size(); bcIdx++ )
    {
        BoundaryCondition curBC = _boundaryCondition[bcIdx];
        Numerics* targetNumerics = analysisModel().numericsManager().giveNumerics(curBC.targetNumerics());
        
        int boundaryId = analysisModel().domainManager().givePhysicalEntityNumberFor(curBC.boundaryName());
        
        // Check for essential boundary condition on nodes
        if ( curBC.conditionType() == "NodalConstraint" )
        {
            int nBCells = analysisModel().domainManager().giveNumberOfBoundaryCells();
            for ( int i = 0; i < nBCells; i++ )
            {
                Cell* curCell = analysisModel().domainManager().giveBoundaryCell(i);
                int label = analysisModel().domainManager().giveLabelOf(curCell);
                
                if ( label == boundaryId )
                {       
                    // Retrieve nodes of boundary element
                    std::vector<Node*> node;
                    node = analysisModel().domainManager().giveNodesOf(curCell);
                    
                    std::vector<Cell*> domCell = analysisModel().domainManager().giveDomainCellsAssociatedWith(curCell);
                    for ( Cell* curDomCell : domCell )
                    {
                        Numerics* domCellNumerics = analysisModel().domainManager().giveNumericsFor(curDomCell);
                        if ( domCellNumerics == targetNumerics )
                        {
                            int targetDofNum = targetNumerics->giveIndexOfNodalDof(curBC.targetDof());
                            
                            for ( Node* curNode : node )
                            {
                                Dof* targetDof = analysisModel().domainManager().giveNodalDof(targetDofNum, curNode);
                                analysisModel().dofManager().putDirichletConstraintOn(targetDof);
                            }
                        }
                    }
                }
            }
        }
        
        // Check for essential boundary condition on cells
        if ( curBC.conditionType() == "CellConstraint" )
        {
            int nDCells = analysisModel().domainManager().giveNumberOfDomainCells();
            for ( int i = 0; i < nDCells; i++ )
            {
                Cell* curCell = analysisModel().domainManager().giveDomainCell(i);
                int label = analysisModel().domainManager().giveLabelOf(curCell);
                
                if ( label == boundaryId )
                {
                    Numerics* numerics = analysisModel().domainManager().giveNumericsForDomain(label);
                    int targetDofNum = numerics->giveIndexOfCellDof(curBC.targetDof());
                    
                    Dof* targetDof = analysisModel().domainManager().giveCellDof(targetDofNum, curCell);
                    analysisModel().dofManager().putDirichletConstraintOn(targetDof);
                }
            }
        }
    }
}

void LoadStep::performPrefinalCalculationsAtCells()
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
        
        numerics->performPrefinalizationCalculationsAt(curCell);
    }
}