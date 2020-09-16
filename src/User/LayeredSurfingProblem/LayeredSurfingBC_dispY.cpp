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

#include "LayeredSurfingBC_dispY.hpp"

#include <cmath>
#include "Core/ObjectFactory.hpp"
#include "Util/readOperations.hpp"

using namespace broomstyx;
        
registerBroomstyxObject(UserFunction, LayeredSurfingBC_dispY)

LayeredSurfingBC_dispY::LayeredSurfingBC_dispY()
{
    _name = "LayeredSurfingBC_dispY";
}

LayeredSurfingBC_dispY::~LayeredSurfingBC_dispY() {}

double LayeredSurfingBC_dispY::at( const RealVector& coor, const TimeData& time )
{
    double PI = 4.*std::atan(1.);

    double x = coor(0) - time.target;
    double y = coor(1);

    // convert to polar coordinates
    double r = std::sqrt(x*x + y*y);
    double theta = std::acos(x/r);
    if ( y < 0 )
        theta = -theta;

    // Kolosov constant for plane strain
    double kappa = 3. - 4.*_nu;

    // Calculate y-displacement
    return _KI * (1. + _nu)/_E * std::sqrt(r/(2.*PI)) * (kappa - x/r) * std::sin(theta/2.);
}

void LayeredSurfingBC_dispY::readDataFrom( FILE* fp )
{
    _E = getRealInputFrom(fp, "Failed to read value of Young's modulus from input file!", _name);
    _nu = getRealInputFrom(fp, "Failed to read value of Poisson ratio from input file!", _name);
    _KI = getRealInputFrom(fp, "Failed to read value of KI stress intensity factor from input file!", _name);
}
