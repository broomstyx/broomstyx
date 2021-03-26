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

#include "SolutionManager.hpp"
#include <chrono>
#include <stdexcept>

#include "AnalysisModel.hpp"
#include "Diagnostics.hpp"
#include "ObjectFactory.hpp"
#include "DomainManager.hpp"
#include "DomainManager.hpp"
#include "LoadStep.hpp"
#include "Node.hpp"
#include "OutputManager.hpp"
#include "Numerics/Numerics.hpp"
#include "User/UserFunction.hpp"
#include "Util/readOperations.hpp"

using namespace broomstyx;

// Constructor
SolutionManager::SolutionManager() {}

// Destructor
SolutionManager::~SolutionManager() 
{
#ifdef VERBOSE_DESTRUCTION
    std::printf("\n  Destroying SolutionManager... ");
    std::fflush(stdout);
#endif

    for ( int i = 0; i < (int)_loadStep.size(); i++ )
        if ( _loadStep[i] )
            delete _loadStep[i];
    
    for ( int i = 0; i < (int)_userFunction.size(); i++ )
        delete _userFunction[i];
    
#ifdef VERBOSE_DESTRUCTION
    std::printf("done.");
    std::fflush(stdout);
#endif
}

// Public methods
// ----------------------------------------------------------------------------
void SolutionManager::commenceSolution()
{
    std::chrono::time_point<std::chrono::system_clock> tic, toc;
    std::chrono::duration<double> tictoc;
    
    // Set stages for nodal and elemental degrees of freedom
    tic = std::chrono::high_resolution_clock::now();
    int nDomCells = analysisModel().domainManager().giveNumberOfDomainCells();
    for ( int i = 0; i < nDomCells; i++ )
    {
        Cell* curCell = analysisModel().domainManager().giveDomainCell(i);
        int label = analysisModel().domainManager().giveLabelOf(curCell);
        Numerics* numerics = analysisModel().domainManager().giveNumericsForDomain(label);
        numerics->setDofStagesAt(curCell);
    }
    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    diagnostics().addSetupTime(tictoc.count());

    // Impose initial conditions
    std::printf("\n  %-40s", "Imposing initial conditions ...");
    std::fflush(stdout);
    tic = std::chrono::high_resolution_clock::now();
    this->imposeInitialConditions();
    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    std::printf("done (time = %f sec.)\n", tictoc.count());
    diagnostics().addSetupTime(tictoc.count());
    
    // Write output corresponding to initial state
    TimeData time;
    time.start = 0.0;
    
    // analysisModel().outputManager().writeOutput(time.start);
    analysisModel().outputManager().writeOutputQuantities(time.start);
    
    // Solve load steps
    for (int i = 0; i < (int)_loadStep.size(); i++)
    {
        _curLoadStep = _loadStep[i];
        _curLoadStep->solveYourself();
    }
}
// ----------------------------------------------------------------------------
LoadStep* SolutionManager::giveCurrentLoadStep()
{
    return _curLoadStep;
}
// ----------------------------------------------------------------------------
int SolutionManager::giveNumberOfSolutionStages()
{ 
    return (int)_stage.size();
}
// ----------------------------------------------------------------------------
void SolutionManager::imposeInitialConditions()
{
    for ( int i = 0; i < (int)_initCond.size(); i++ )
    {
        std::string domainLabel = _initCond[i].domainLabel();
        int domainId = analysisModel().domainManager().givePhysicalEntityNumberFor(domainLabel);
        int nCells = analysisModel().domainManager().giveNumberOfDomainCells();

        std::string condType = _initCond[i].conditionType();

        if ( condType == "nodalDof" )
        {
            int nNodes = analysisModel().domainManager().giveNumberOfNodes();
            std::vector<bool> nodeIsInitialized(nNodes, false);

            for ( int j = 0; j < nCells; i++ )
            {
                Cell* curCell = analysisModel().domainManager().giveDomainCell(j);
                int cellLabel = analysisModel().domainManager().giveLabelOf(curCell);

                if ( cellLabel == domainId )
                {
                    std::vector<Node*> cellNode = analysisModel().domainManager().giveNodesOf(curCell);
                    for ( auto it = cellNode.begin(); it != cellNode.end(); it++ )
                    {
                        int nodeId = analysisModel().domainManager().giveIdOf(*it);
                        if ( !nodeIsInitialized[nodeId] )
                        {
                            RealVector coor = analysisModel().domainManager().giveCoordinatesOf(*it);
                            int dofNum = _initCond[i].targetDofNumber();
                            double val = _initCond[i].valueAt(coor);
                            
                            Dof* targetDof = analysisModel().domainManager().giveNodalDof(dofNum, *it);
                            analysisModel().dofManager().updatePrimaryVariableAt(targetDof, val, converged_value);
                        }
                    }
                }
            }
        }
        else if ( condType == "CellDof" )
        {
            Numerics* numerics = analysisModel().domainManager().giveNumericsForDomain(domainId);

            for ( int j = 0; j < nCells; j++ )
            {
                Cell* curCell = analysisModel().domainManager().giveDomainCell(j);
                int cellLabel = analysisModel().domainManager().giveLabelOf(curCell);

                if ( cellLabel == domainId )
                {
                    numerics->imposeInitialConditionAt(curCell, _initCond[i]);
                }
            }
        }
    }
}
// ----------------------------------------------------------------------------
std::vector<int> SolutionManager::giveRegisteredSolutionStages()
{
    std::vector<int> stage;
    stage.assign((int)_stage.size(), 0);
    int counter = 0;
    for ( auto i : _stage )
        stage[counter++] = i;
    
    return stage;
}
// ----------------------------------------------------------------------------
UserFunction* SolutionManager::makeNewUserFunction( std::string name )
{
    UserFunction* usrFcn = objectFactory().instantiateUserFunction(name);
    
    _userFunction.push_back(usrFcn);
    return usrFcn;
}
// ----------------------------------------------------------------------------
void SolutionManager::readInitialConditionsFrom( FILE* fp )
{
    std::string errmsg, src = "SolutionManager";
    
    int nInitCond = getIntegerInputFrom(fp, "Failed to read number of initial conditions from input file!", src);
    _initCond.assign(nInitCond, InitialCondition());
    
    for ( int i = 0; i < nInitCond; i++ )
        _initCond[i].readDataFrom(fp);
}
// ----------------------------------------------------------------------------
void SolutionManager::readLoadStepsFrom( FILE *fp )
{
    std::string src = "SolutionManager";
    
    int nLoadSteps = getIntegerInputFrom(fp, "Failed to read number of load steps from input file!", src);
    
    _loadStep.assign(nLoadSteps, nullptr);
    
    for ( int i = 0; i < nLoadSteps; i++)
    {
        int lsNum = getIntegerInputFrom(fp, "Failed to read load step number from input file!", src);
        _loadStep[i] = new LoadStep(lsNum, _stage.size());
        _loadStep[i]->readDataFrom(fp);
    }
}
// ----------------------------------------------------------------------------
void SolutionManager::readNumberOfStagesFrom( FILE* fp )
{
    std::string src = "SolutionManager";
    _nStages = getIntegerInputFrom(fp, "Failed to read number of solution stages from input file.", src);
}
// ----------------------------------------------------------------------------
void SolutionManager::registerStage( int stage, std::string tag )
{
    _stage.insert(stage);
}
// ----------------------------------------------------------------------------
void SolutionManager::reportRegisteredStages()
{
    std::printf("\nNumber of registered solution stages = %d.\n", (int)_stage.size());
    for ( auto curStage : _stage )
        std::printf("%d ", curStage);
    std::printf("\n");
}