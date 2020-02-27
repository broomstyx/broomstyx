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

#include "AnalysisModel.hpp"
#include <chrono>
#include <stdexcept>
#include <string>
#include "Diagnostics.hpp"
#include "DofManager.hpp"
#include "DomainManager.hpp"
#include "MaterialManager.hpp"
#include "NumericsManager.hpp"
#include "ObjectFactory.hpp"
#include "OutputManager.hpp"
#include "SolutionManager.hpp"
#include "../MeshReaders/MeshReader.hpp"
#include "../Util/readOperations.hpp"

using namespace broomstyx;

typedef std::chrono::time_point<std::chrono::system_clock> TimePoint;
typedef std::chrono::duration<double> TimeDuration;
typedef std::chrono::high_resolution_clock Timer;

AnalysisModel::AnalysisModel()
    : _meshReader(nullptr)
{}

AnalysisModel::~AnalysisModel()
{
#ifdef VERBOSE_DESTRUCTION
    std::printf("\nDestroying AnalysisModel... ");
    std::fflush(stdout);
#endif
    
    if ( _meshReader )
        delete _meshReader;
    
    delete _outputManager;
    delete _solutionManager;
    delete _domainManager;
    delete _numericsManager;
    delete _materialManager;
    delete _dofManager;
    
#ifdef VERBOSE_DESTRUCTION
    std::printf("\nDone destroying AnalysisModel.\n\n");
    std::fflush(stdout);
#endif
}

// Public methods
void AnalysisModel::initializeYourself( std::string filename )
{
    TimePoint tic, toc;
    TimeDuration tictoc;
    
    tic = Timer::now();

    std::printf("\n----------------------------------------------------------------------");
    std::printf("\n               I N S T A N T I A T I O N    P H A S E");
    std::printf("\n----------------------------------------------------------------------\n\n");
    
    // Construct input filename
    _inputFilename = filename + ".inp";
    
    // Instantiate managers
    _dofManager = new DofManager();
    _domainManager = new DomainManager();
    _materialManager = new MaterialManager();
    _numericsManager = new NumericsManager();
    _outputManager = new OutputManager();
    _solutionManager = new SolutionManager();
    
    // Read input file
    try
    {
        readInputFile();
    }
    catch (std::exception& e)
    {
        std::printf("\n\nERROR reading input file: %s\n", e.what());
        throw std::runtime_error("Program run aborted!");
    }
    
    // Initialize materials
    _materialManager->initializeMaterials();
    
    // Read mesh file and create node and cell objects
    if ( !_meshReader )
        throw std::runtime_error("Error: MeshReader not defined!\n");
    else
        _meshReader->readMeshFile( _meshFilename );
    
    // Find domain cell neighbors
    _domainManager->findDomainCellNeighbors();
    
    // Construct cell faces
    _domainManager->constructCellFaces();
    
    // Find boundary cell associations
    _domainManager->findBoundaryAssociations();
    
    // Initialize numerics at cells
    _domainManager->initializeNumericsAtCells();
    
    // Initialize material data at cells
    _domainManager->initializeMaterialsAtCells();
    
    // Initialize output writer
    _outputManager->initializeOutputWriter();
    
    // Initialize CSV output quantities
    _outputManager->initializeCSVOutput();
    
    // Impose multi-freedom constraints
    _dofManager->imposeMultiFreedomConstraints();

    toc = Timer::now();
    tictoc = toc - tic;
    diagnostics().addSetupTime(tictoc.count());
}

void AnalysisModel::solveYourself()
{
    std::printf("\n----------------------------------------------------------------------");
    std::printf("\n                    S O L U T I O N    P H A S E");
    std::printf("\n----------------------------------------------------------------------\n\n");
    
    _solutionManager->commenceSolution();
}

// -----------------------------------------------------------------------------
DofManager& AnalysisModel::dofManager() { return *_dofManager; }
// -----------------------------------------------------------------------------
DomainManager& AnalysisModel::domainManager() { return *_domainManager; }
// -----------------------------------------------------------------------------
MaterialManager& AnalysisModel::materialManager() { return *_materialManager; }
// -----------------------------------------------------------------------------
NumericsManager& AnalysisModel::numericsManager() { return *_numericsManager; }
// -----------------------------------------------------------------------------
OutputManager& AnalysisModel::outputManager() { return *_outputManager; }
// -----------------------------------------------------------------------------
SolutionManager& AnalysisModel::solutionManager() { return *_solutionManager; }
// -----------------------------------------------------------------------------
MeshReader& AnalysisModel::meshReader() { return *_meshReader; }

// Private methods
// -----------------------------------------------------------------------------
void AnalysisModel::readMeshReaderFrom( FILE* fp )
{
    std::string name, src = "AnalysisModel";
    name = getStringInputFrom(fp, "Failed to read mesh format from input file.", src);
    
    _meshReader = objectFactory().instantiateMeshReader(name);
    
    if ( !_meshReader )
        throw std::runtime_error("Error: MeshReader not defined!\n");
}

void AnalysisModel::readInputFile()
{
    FILE* fp = nullptr;
    TimePoint tic, toc;
    TimeDuration tictoc;
    std::string src = "AnalysisModel";
    
    std::printf("\n\n  %-40s", "Reading input file ...");
    std::fflush(stdout);
    
    tic = Timer::now();
    
    fp = std::fopen( _inputFilename.c_str(), "r" );
    if ( !fp )
        throw std::runtime_error("\nSpecified input file '" + _inputFilename + "' not found!");

    std::string decl;
    do {
        decl = getDeclarationFrom(fp);
        
        if ( decl == "CSV_OUTPUT" )
            _outputManager->readDataForCSVOutputFrom(fp);
        else if ( decl == "DOF_PER_CELL" )
            _dofManager->readCellDofsFrom(fp);
        else if ( decl == "DOF_PER_FACE" )
            _dofManager->readFaceDofsFrom(fp);
        else if ( decl == "DOF_PER_NODE" )
            _dofManager->readNodalDofsFrom(fp);
        else if ( decl == "DOMAIN_ASSIGNMENTS" )
            _domainManager->readDomainAssignmentsFrom(fp);
        else if ( decl == "FIELDS_PER_CELL" )
            _domainManager->readNumberOfFieldsPerCellFrom(fp);
        else if ( decl == "FIELDS_PER_FACE" )
            _domainManager->readNumberOfFieldsPerFaceFrom(fp);
        else if ( decl == "FIELDS_PER_NODE" )
            _domainManager->readNumberOfFieldsPerNodeFrom(fp);
        else if ( decl == "INITIAL_CONDITIONS" )
            _solutionManager->readInitialConditionsFrom(fp);
        else if ( decl == "LOADSTEPS" )
            _solutionManager->readLoadStepsFrom(fp);
        else if ( decl == "MATERIALS" )
            _materialManager->readMaterialsFrom(fp);
        else if ( decl == "MESH_FILE" )
            _meshFilename = getStringInputFrom( fp, "\nFailed to read mesh filename from input file!", src);
        else if ( decl == "MESH_READER" )
            this->readMeshReaderFrom(fp);
        else if ( decl == "MULTIFREEDOM_CONSTRAINTS" )
            _dofManager->readMultiFreedomConstraintsFrom(fp);
        else if ( decl == "NUMERICS" )
            _numericsManager->readNumericsFrom(fp);
        else if ( decl == "OUTPUT_FORMAT" )
            _outputManager->readOutputWriterFromFile(fp);
        else if ( decl == "SOLUTION_STAGES" )
            _solutionManager->readNumberOfStagesFrom(fp);
        else if ( decl != "END" )
            throw std::runtime_error("Error: Unrecognized declaration '" + decl + "' encountered in input file!\n");
    }
    while ( decl != "END" );
    
    std::fclose(fp);
    toc = Timer::now();
    tictoc = toc - tic;
    std::printf("done (time = %f sec.)\n", tictoc.count());
}

// ----------------------------------------------------------------------------
AnalysisModel& broomstyx::analysisModel()
{
    static AnalysisModel enggModel;
    return enggModel;
}