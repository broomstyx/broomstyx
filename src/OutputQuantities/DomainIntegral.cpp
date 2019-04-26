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

#include "DomainIntegral.hpp"
#include <cstring>
#include <tuple>
#include "../Core/AnalysisModel.hpp"
#include "../Core/ObjectFactory.hpp"
#include "../Core/DomainManager.hpp"
#include "../MeshReaders/MeshReader.hpp"
#include "../Numerics/Numerics.hpp"
#include "../Util/RealVector.hpp"
#include "../Util/linearAlgebra.hpp"
#include "../Util/readOperations.hpp"

using namespace broomstyx;

registerBroomstyxObject(OutputQuantity, DomainIntegral)

// Constructor
DomainIntegral::DomainIntegral()
{
    _name = "DomainIntegral";
}

// Destructor
DomainIntegral::~DomainIntegral() {}

// Public methods
double DomainIntegral::computeOutput()
{
    double result = 0.0;
    
    // Retrieve relevant numerics
    int physNum = analysisModel().domainManager().givePhysicalEntityNumberFor(_physTag);
    Numerics* numerics = analysisModel().domainManager().giveNumericsForDomain(physNum);
    
    // Cycle through all domain cells
    int nCells = analysisModel().domainManager().giveNumberOfDomainCells();
    for ( int i = 0; i < nCells; i++ )
    {
        Cell* curCell = analysisModel().domainManager().giveDomainCell(i);
        if ( analysisModel().domainManager().giveLabelOf(curCell) == physNum )
        {
            RealVector val, wt;
            std::tie(val,wt) = numerics->giveCellFieldOutputAtEvaluationPointsOf(curCell, _cellFieldNum);
            result += val.dot(wt);
        }
    }
    
    return result;
}

void DomainIntegral::initialize() {}

void DomainIntegral::readDataFrom( FILE* fp )
{
    // Read domain tag
    _physTag = getStringInputFrom(fp, "Failed reading domain label from input file", _name);
            
    // Read cell field number
    verifyKeyword(fp, "CellField", _name);
    _cellFieldNum = getIntegerInputFrom(fp, "Failed reading cell field number from input file.", _name);
}