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

#include "CellFieldOutput.hpp"
#include <cstring>
#include <stdexcept>
#include "Core/AnalysisModel.hpp"
#include "Core/ObjectFactory.hpp"
#include "Core/DomainManager.hpp"
#include "Numerics/Numerics.hpp"
#include "Util/readOperations.hpp"

using namespace broomstyx;

registerBroomstyxObject(OutputQuantity, CellFieldOutput)

// Constructor
CellFieldOutput::CellFieldOutput()
{
    _name = "CellFieldOutput";
}

// Destructor
CellFieldOutput::~CellFieldOutput() {}

// Public methods
double CellFieldOutput::computeOutput()
{
    double value = 0.;
    // Find relevant cell
    int physNum = analysisModel().domainManager().givePhysicalEntityNumberFor(_physTag);
    Numerics* numerics = analysisModel().domainManager().giveNumericsForDomain(physNum);
    
    int nCells = analysisModel().domainManager().giveNumberOfDomainCells();
    int count = 0;
    for ( int i = 0; i < nCells; i++ )
    {
        Cell* curCell = analysisModel().domainManager().giveDomainCell(i);

        if ( analysisModel().domainManager().giveLabelOf(curCell) == physNum )
        {
            ++count;
            if ( count > 1 )
                throw std::runtime_error("ERROR: More than one cell detected under specified physical tag!\nSource = " + _name);
            
            value = numerics->giveCellFieldValueAt(curCell, _cellField);
        }
    }

    return value;
}

void CellFieldOutput::initialize() {}

void CellFieldOutput::readDataFrom( FILE *fp )
{
    std::string key, errmsg;
    
    // Read physical tag of boundary
    _physTag = getStringInputFrom(fp, "Failed to read cell label from input file!", _name);
    _cellField = getIntegerInputFrom(fp, "Failed to read cell field number from input file!", _name);
}