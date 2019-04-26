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

#include "ThresholdDegradation.hpp"
#include "../../../Core/ObjectFactory.hpp"
#include "../../../Util/readOperations.hpp"

using namespace broomstyx;

registerBroomstyxObject(Material, ThresholdDegradation)

// Constructor
ThresholdDegradation::ThresholdDegradation()
{
    _name = "ThresholdDegradation";
}

// Destructor
ThresholdDegradation::~ThresholdDegradation() {}

// Public methods
// ----------------------------------------------------------------------------
double ThresholdDegradation::givePotentialFrom( const RealVector& conState, const MaterialStatus* matStatus )
{
    double phi = conState(0);
    double val;
    
    if ( phi < _threshold )
        val = 1.;
    else
        val = 0.;
    
    return val;
}
// ----------------------------------------------------------------------------
RealVector ThresholdDegradation::giveForceFrom( const RealVector& conState, const MaterialStatus* matStatus )
{
    RealVector force(1);
    return force;
}
// ----------------------------------------------------------------------------
RealMatrix ThresholdDegradation::giveModulusFrom( const RealVector& conState, const MaterialStatus* matStatus )
{
    RealMatrix modulus(1,1);
    return modulus;
}
// ----------------------------------------------------------------------------
void ThresholdDegradation::readParamatersFrom( FILE* fp )
{
    _threshold = getRealInputFrom(fp, "Failed to read threshold value from input file!", _name);
}