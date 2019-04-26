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

#include "Line_P2.hpp"

using namespace broomstyx;

Line_P2::Line_P2() {}

Line_P2::~Line_P2() {}

RealVector Line_P2::giveBasisFunctionsAt( const RealVector& coor )
{
    double xi = coor(0);
    RealVector psi({0.5*xi*(xi - 1.), 0.5*xi*(xi + 1.), 1. - xi*xi});
    return psi;
}

std::vector<RealVector>
Line_P2::giveBasisFunctionDerivativesAt( const RealVector& coor )
{
    std::vector<RealVector> dpsi;
    dpsi.assign(1, RealVector());
    
    double xi = coor(0);
    dpsi[0] = {xi - 0.5, xi + 0.5, -2.*xi};
    return dpsi;
}