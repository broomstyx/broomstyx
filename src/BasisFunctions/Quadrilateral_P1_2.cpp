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

#include "Quadrilateral_P1_2.hpp"
#include <cmath>
#include "../Core/AnalysisModel.hpp"
#include "../Core/DomainManager.hpp"
#include "../Util/linearAlgebra.hpp"

using namespace broomstyx;

Quadrilateral_P1_2::Quadrilateral_P1_2() {}

Quadrilateral_P1_2::~Quadrilateral_P1_2() {}

RealVector Quadrilateral_P1_2::giveBasisFunctionsAt( const RealVector& coor )
{
    double xi = coor(0);
    double eta = coor(1);
    
    RealVector psi({0.25*(1. - xi)*(1. - eta),
                    0.25*(1. + xi)*(1. - eta),
                    0.25*(1. + xi)*(1. + eta),
                    0.25*(1. - xi)*(1. + eta)});
    
    return psi;
}
// ----------------------------------------------------------------------------
std::vector<RealVector> Quadrilateral_P1_2::giveBasisFunctionDerivativesAt( const RealVector& coor )
{
    std::vector<RealVector> dpsi;
    dpsi.assign(2, RealVector());
    
    double xi = coor(0);
    double eta = coor(1);
    
    dpsi[0] = {-0.25*(1. - eta),
               0.25*(1. - eta),
               0.25*(1. + eta),
               -0.25*(1. + eta)};
    
    dpsi[1] = {-0.25*(1. - xi),
               -0.25*(1. + xi),
               0.25*(1. + xi),
               0.25*(1. - xi)};
    
    return dpsi;
}