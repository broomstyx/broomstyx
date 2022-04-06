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

#include "BoundaryCondition.hpp"
#include <stdexcept>
#include "SolutionManager.hpp"
#include "User/UserFunction.hpp"
#include "Util/readOperations.hpp"
#include "AnalysisModel.hpp"
#include "DofManager.hpp"

using namespace broomstyx;

// Constructor
BoundaryCondition::BoundaryCondition() {}

// Destructor
BoundaryCondition::~BoundaryCondition() {}

// Public methods
// ----------------------------------------------------------------------------
std::string BoundaryCondition::boundaryName() const
{ 
    return _bndName; 
}
// ----------------------------------------------------------------------------
std::string BoundaryCondition::conditionType() const
{ 
    return _cndType; 
}
// ----------------------------------------------------------------------------
void BoundaryCondition::readDataFrom( FILE* fp )
{
    std::string src = "BoundaryCondition";
    
    // Boundary on which condition is acting
    _bndName = getStringInputFrom(fp, "Failed to read boundary name from input file!", src);
    
    // Target numerics
    _numericsTag = getIntegerInputFrom(fp, "Failed to read target numerics tag from input file!", src);
    
    // Condition type
    _cndType = getStringInputFrom(fp, "Failed to read boundary condition type from input file!", src);
    
    // Target DOF name
    _dofTag = getStringInputFrom(fp, "Failed to read DOF name from input file!", src);
    
    // Specification Type
    _specType = getStringInputFrom(fp, "Failed to read specification type for boundary condition from input file!", src);
    
    if ( _specType == "None" )
    {
        // Do nothing.
    }
    else if ( _specType == "Constant" )
    {
        // Read starting and ending values
        _startVal = broomstyx::getRealInputFrom(fp, "Failed to read value for boundary condition from input file!", src);
        _endVal = _startVal;
    }
    else if ( _specType == "Linear" )
    {
        // Read starting and ending values
        _startVal = broomstyx::getRealInputFrom(fp, "Failed to read starting value for boundary condition from input file!", src);
        _endVal = broomstyx::getRealInputFrom(fp, "Failed to read end value for boundary condition from input file!", src);
    }
    else if ( _specType ==  "UserFunction" )
    {
        // Read user function number
        std::string fcnName = getStringInputFrom(fp, "Failed to read user function for boundary condition from input file!", src);
        _usrFcn = analysisModel().solutionManager().makeNewUserFunction(fcnName);
        _usrFcn->readDataFrom(fp);
    }
    else
        throw std::runtime_error("Invalid specification type '" + _specType + "' for boundary condition encountered in input file!");
}
// ----------------------------------------------------------------------------
std::string BoundaryCondition::targetDof() const
{
    return _dofTag;
}
// ----------------------------------------------------------------------------
int BoundaryCondition::targetNumerics() const
{
    return _numericsTag;
}
// ----------------------------------------------------------------------------
double BoundaryCondition::valueAt( const RealVector& coor, const TimeData& time ) const
{
    if ( _specType == "Constant" )
        return _startVal;
    else if ( _specType == "Linear" )
        return _startVal + (_endVal - _startVal)*
        (time.giveTargetTime() - time.giveStartTime())/(time.giveEndTime() - time.giveStartTime());
    else if ( _specType == "UserFunction" )
        return _usrFcn->at(coor,time);
    
    return 0.0;
}
