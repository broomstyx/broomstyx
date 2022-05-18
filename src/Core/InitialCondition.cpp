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

#include "InitialCondition.hpp"

#include <stdexcept>
#include "SolutionManager.hpp"
#include "User/UserFunction.hpp"
#include "Util/readOperations.hpp"
#include "AnalysisModel.hpp"

using namespace broomstyx;

InitialCondition::InitialCondition() {}

InitialCondition::~InitialCondition() {}

// ----------------------------------------------------------------------------
std::string InitialCondition::domainLabel() const { return _domainLabel; }
// ----------------------------------------------------------------------------
void InitialCondition::readDataFrom( FILE* fp )
{
    std::string str;
    std::string src = "class InitialCondition";
    
    // Domain on which to apply initial condition
    _domainLabel = getStringInputFrom(fp, "Failed to read domain label for initial condition from input file!", src);
    
    // Target DOF type
    _conditionType = getStringInputFrom(fp, "Failed to read target DOF type for initial condition from input file!", src);
    
    // Target DOF label
    if ( _conditionType == "NodalDof" || _conditionType == "CellDof" )
    {
        str = getStringInputFrom(fp, "Failed to read target DOF label for initial condition from input file!", src);
        _dofNumber = analysisModel().dofManager().giveIndexForNodalDof(str);
    }
    
    // Specification Type
    _specType = getStringInputFrom(fp, "Failed to read specification type for initial condition from input file!", src);
    
    if ( _specType == "Constant" )
        _val = getRealInputFrom(fp, "Failed to value for initial condition from input file!", src);
    else if ( _specType == "UserFunction" )
    {
        std::string fcnName = getStringInputFrom(fp, "Failed to read user function for initial condition from input file!", src);
        _usrFcn = analysisModel().solutionManager().makeNewUserFunction(fcnName);
        _usrFcn->readDataFrom(fp);
    }
    else
        throw std::runtime_error("Invalid specification type '" + str + "' encountered for initial condition in input file!");
}
// ----------------------------------------------------------------------------
int InitialCondition::targetDofNumber() const { return _dofNumber; }
// ----------------------------------------------------------------------------
std::string InitialCondition::conditionType() const { return _conditionType; }
// ----------------------------------------------------------------------------
double InitialCondition::valueAt( const RealVector& coor ) const
{
    if ( _specType == "Constant" )
        return _val;
    else // has to be user function
    {
        TimeData time; // Is automatically initialized to zero
        return _usrFcn->at(coor, time);
    }
}
