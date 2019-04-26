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

#include "SurfingBC_dispY.hpp"

#include <cmath>
#include "../../Core/ObjectFactory.hpp"

using namespace broomstyx;
        
registerBroomstyxObject(UserFunction, SurfingBC_dispY)

SurfingBC_dispY::SurfingBC_dispY() {}

SurfingBC_dispY::~SurfingBC_dispY() {}

double SurfingBC_dispY::at( const RealVector& coor, const TimeData& time )
{
    double PI = 4.*std::atan(1.);

    // Hard-coded parameters
    double Gc = 2.7;
    double E  = 210.e3;
    double nu = 0.3;

    double x = coor(0) - time.target;
    double y = coor(1);

    // convert to polar coordinates
    double r = std::sqrt(x*x + y*y);
    double theta = std::acos(x/r);
    if ( y < 0 )
        theta = -theta;

    // Stress intensity factor
    double KI = std::sqrt(E*Gc);

    // Kolosov constant for plane strain
    double kappa = 3. - 4.*nu;

    // Calculate y-displacement
    return KI * (1. + nu)/E * std::sqrt(r/(2.*PI)) * (kappa - x/r) * std::sin(theta/2.);
}