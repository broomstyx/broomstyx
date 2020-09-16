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

#include "FieldCondition.hpp"

#include <stdexcept>
#include "SolutionManager.hpp"
#include "User/UserFunction.hpp"
#include "Util/readOperations.hpp"
#include "AnalysisModel.hpp"

using namespace broomstyx;

FieldCondition::FieldCondition() {}

FieldCondition::~FieldCondition() {}

// ----------------------------------------------------------------------------
std::string FieldCondition::conditionType() const 
{ 
    return _cndType;
}
// ----------------------------------------------------------------------------
std::string FieldCondition::domainLabel() const 
{ 
    return _domainLabel; 
}
// ----------------------------------------------------------------------------
void FieldCondition::readDataFrom( FILE* fp )
{
    std::string str;
    std::string src = "class FieldCondition";
    
    // Domain on which to apply field condition
    _domainLabel = getStringInputFrom(fp, "Failed to read domain label for field condition from input file!", src);
    
    // Name of field
    _cndType = getStringInputFrom(fp, "Failed to read type of field condition from input file!", src);
    
    // Specification Type
    _specType = getStringInputFrom(fp, "Failed to read specification type for field condition from input file!", src);
    
    if ( _specType == "None" )
    {
        // Do nothing.
    }
    else if ( _specType == "Constant" )
    {
        // Read starting and ending values
        _startVal = getRealInputFrom(fp, "Failed to read value for field condition from input file!", src);
        
        _endVal = _startVal;
    }
    else if ( _specType == "Linear" )
    {
        // Read starting and ending values
        _startVal = getRealInputFrom(fp, "Failed to read starting value for field condition from input file!", src);
        
        _endVal = getRealInputFrom(fp, "Failed to read end value for field condition from input file!", src);
    }
    else if ( _specType == "UserFunction" )
    {
        std::string fcnName = getStringInputFrom(fp, "Failed to read user function name for field condition from input file!", src);
        _usrFcn = analysisModel().solutionManager().makeNewUserFunction(fcnName);
        _usrFcn->readDataFrom(fp);
    }
    else
        throw std::runtime_error("Invalid specification type '" + str + "' for field condition in input file!");
}
// ----------------------------------------------------------------------------
double FieldCondition::valueAt( const RealVector& coor, const TimeData& time ) const
{
    if ( _specType == "Constant" )
        return _startVal;
    else if ( _specType == "Linear" )
        return _startVal + (_endVal - _startVal)*(time.target - time.start)/(time.end - time.start);
    else // has to be user function
        return _usrFcn->at(coor,time);
}