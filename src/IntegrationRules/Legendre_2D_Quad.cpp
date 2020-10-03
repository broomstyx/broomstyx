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

#include "Legendre_2D_Quad.hpp"
#include "Legendre_1D.hpp"
#include <cmath>

using namespace broomstyx;

Legendre_2D_Quad::Legendre_2D_Quad( int nPoints ) : IntegrationRule(nPoints) {}

Legendre_2D_Quad::~Legendre_2D_Quad() {}

std::tuple< std::vector<RealVector>, RealVector > 
Legendre_2D_Quad::giveIntegrationPointsAndWeights()
{
    std::vector<RealVector> gpLoc(_nIntegrationPoints, RealVector());
    RealVector gpWt(_nIntegrationPoints);

    // 2D integration rules are derived from appropriate lower-dimensional (1D) rules
    int _nIntegrationPoints_1d;

    switch (_nIntegrationPoints)
    {
        case 4:
            _nIntegrationPoints_1d = 2;
            break;
        case 9:
            _nIntegrationPoints_1d = 3;
            break;
        case 16:
            _nIntegrationPoints_1d = 4;
            break;
        case 25:
            _nIntegrationPoints_1d = 5;
            break;
        default:
            throw std::runtime_error("Invalid number of Gauss points specified for 'Legendre_2D_Quad'!");
    }

    Legendre_1D integRule_1d(_nIntegrationPoints_1d);
    std::vector<RealVector> gpLoc_1d;
    RealVector gpWt_1d;

    std::tie(gpLoc_1d, gpWt_1d) = integRule_1d.giveIntegrationPointsAndWeights();
    for ( int i = 0; i < _nIntegrationPoints_1d; i++ )
        for ( int j = 0; j < _nIntegrationPoints_1d; j++ )
        {
            gpLoc[2*i + j] = {gpLoc_1d[i](0), gpLoc_1d[j](0)};
            gpWt(2*i + j) = gpWt_1d(i)*gpWt_1d(j);
        }
    
    return std::make_tuple(gpLoc, gpWt);
}