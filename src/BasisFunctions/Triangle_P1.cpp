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

#include "Triangle_P1.hpp"

#include "../Core/AnalysisModel.hpp"
#include "../Core/DomainManager.hpp"

using namespace broomstyx;

Triangle_P1::Triangle_P1() {}

Triangle_P1::~Triangle_P1() {}

double Triangle_P1::giveAreaOf( const std::vector<Node*>& node )
{
    RealVector coor_0, coor_1, coor_2;
    
    coor_0 = analysisModel().domainManager().giveCoordinatesOf(node[0]);
    coor_1 = analysisModel().domainManager().giveCoordinatesOf(node[1]);
    coor_2 = analysisModel().domainManager().giveCoordinatesOf(node[2]);
    
    double area = 0.5*(coor_0(0)*(coor_1(1) - coor_2(1)) + 
                       coor_1(0)*(coor_2(1) - coor_0(1)) +
                       coor_2(0)*(coor_0(1) - coor_1(1)));
    
    if ( area <= 0 )
        throw std::runtime_error("Calculation of negative area detected!\nSource: Triangle_P1 (ShapeFunction)");
    
    return area;
}

RealVector Triangle_P1::giveBasisFunctionsAt( const RealVector& coor, const std::vector<Node*>& node )
{
    RealVector coor_0, coor_1, coor_2;
    
    coor_0 = analysisModel().domainManager().giveCoordinatesOf(node[0]);
    coor_1 = analysisModel().domainManager().giveCoordinatesOf(node[1]);
    coor_2 = analysisModel().domainManager().giveCoordinatesOf(node[2]);
    
    double twiceArea = 2.*this->giveAreaOf(node);
    
    RealVector beta(3), gamma(3);
    beta(0) = (coor_1(1) - coor_2(1))/twiceArea;
    beta(1) = (coor_2(1) - coor_0(1))/twiceArea;
    beta(2) = (coor_0(1) - coor_1(1))/twiceArea;
    
    gamma(0) = (coor_2(0) - coor_1(0))/twiceArea;
    gamma(1) = (coor_0(0) - coor_2(0))/twiceArea;
    gamma(2) = (coor_1(0) - coor_0(0))/twiceArea;
    
    RealVector psi(3);
    
    psi(0) = (coor_1(0)*coor_2(1) - coor_2(0)*coor_1(1))/twiceArea + coor(0)*beta(0) + coor(1)*gamma(0);
    psi(1) = (coor_2(0)*coor_0(1) - coor_0(0)*coor_2(1))/twiceArea + coor(0)*beta(1) + coor(1)*gamma(1);
    psi(2) = (coor_0(0)*coor_1(1) - coor_1(0)*coor_0(1))/twiceArea + coor(0)*beta(2) + coor(1)*gamma(2);
    
    return psi;
}

std::vector<RealVector> Triangle_P1::giveBasisFunctionDerivativesAt( const std::vector<Node*>& node )
{
    RealVector coor_0, coor_1, coor_2;
    
    coor_0 = analysisModel().domainManager().giveCoordinatesOf(node[0]);
    coor_1 = analysisModel().domainManager().giveCoordinatesOf(node[1]);
    coor_2 = analysisModel().domainManager().giveCoordinatesOf(node[2]);
    
    double twiceArea = 2.*this->giveAreaOf(node);
    
    std::vector<RealVector> dpsi;
    dpsi.assign(2, RealVector(3));
    
    dpsi[0](0) = (coor_1(1) - coor_2(1))/twiceArea;
    dpsi[0](1) = (coor_2(1) - coor_0(1))/twiceArea;
    dpsi[0](2) = (coor_0(1) - coor_1(1))/twiceArea;
    
    dpsi[1](0) = (coor_2(0) - coor_1(0))/twiceArea;
    dpsi[1](1) = (coor_0(0) - coor_2(0))/twiceArea;
    dpsi[1](2) = (coor_1(0) - coor_0(0))/twiceArea;
    
    return dpsi;
}