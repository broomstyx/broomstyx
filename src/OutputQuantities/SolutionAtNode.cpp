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

#include "SolutionAtNode.hpp"
#include <cstring>
#include <stdexcept>
#include "../Core/AnalysisModel.hpp"
#include "../Core/ObjectFactory.hpp"
#include "../Core/Dof.hpp"
#include "../Core/DofManager.hpp"
#include "../Core/DomainManager.hpp"
#include "../MeshReaders/MeshReader.hpp"
#include "../Util/readOperations.hpp"

using namespace broomstyx;

registerBroomstyxObject(OutputQuantity, SolutionAtNode)

// Constructor
SolutionAtNode::SolutionAtNode()
{
    _name = "SolutionAtNode";
}

// Destructor
SolutionAtNode::~SolutionAtNode() {}

// Public methods
double SolutionAtNode::computeOutput()
{
    // Form set of all relevant DOFs associated with boundary
    int physNum = analysisModel().domainManager().givePhysicalEntityNumberFor(_physTag);
    std::set<Dof*> dof;
    
    int nBCells = analysisModel().domainManager().giveNumberOfBoundaryCells();
    for ( int i = 0; i < nBCells; i++)
    {
        Cell* curCell = analysisModel().domainManager().giveBoundaryCell(i);
        if ( analysisModel().domainManager().giveLabelOf(curCell) == physNum )
        {
            std::vector<Node*> cellNode = analysisModel().domainManager().giveNodesOf(curCell);

            for ( int j = 0; j < (int)cellNode.size(); j++)
            {
                Dof* curDof = analysisModel().domainManager().giveNodalDof(_dofNum, cellNode[j]);
                dof.insert(curDof);
            }
        }
    }
    
    double result;
    if ( dof.size() > 1 )
        throw std::runtime_error("Error: multiple nodes selected for nodal solution output!\nSource: " + _name);
    else
        result = analysisModel().dofManager().giveValueOfPrimaryVariableAt(*(dof.begin()), converged_value);
    
    return result;
}

void SolutionAtNode::initialize() {}

void SolutionAtNode::readDataFrom( FILE *fp )
{
    std::string key, errmsg;
    
    // Read physical tag of boundary
    _physTag = getStringInputFrom(fp, "Failed to read point label from input file!", _name);
    
    // Read nodal DOF associated with reaction
    verifyKeyword(fp, key = "NodalDof", _name);
    _dofNum = getIntegerInputFrom(fp, "Failed to read DOF number from input file!", _name);
}