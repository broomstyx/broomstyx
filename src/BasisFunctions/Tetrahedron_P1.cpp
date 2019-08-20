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

#include "Tetrahedron_P1.hpp"

#include "../Core/AnalysisModel.hpp"
#include "../Core/DomainManager.hpp"

using namespace broomstyx;

Tetrahedron_P1::Tetrahedron_P1() {}

Tetrahedron_P1::~Tetrahedron_P1() {}

RealVector Tetrahedron_P1::giveBasisFunctionsAt( const RealVector& coor )
{
    double xi = coor(0);
    double eta = coor(1);
    double zeta = coor(2);
    
    RealVector psi({1. - xi - eta - zeta,
                    xi,
                    eta,
                    zeta});
                    
    return psi;
}

std::vector<RealVector> Tetrahedron_P1::giveBasisFunctionDerivativesAt( const RealVector& coor )
{
    std::vector<RealVector> dpsi;
    dpsi.assign(3, RealVector());
    
    dpsi[0] = {-1., 1., 0., 0.};
    dpsi[1] = {-1., 0., 1., 0.};
    dpsi[2] = {-1., 0., 0., 1.};
    
    return dpsi;
}