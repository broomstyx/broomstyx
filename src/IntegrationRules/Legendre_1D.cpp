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

#include "Legendre_1D.hpp"
#include <cmath>

using namespace broomstyx;

Legendre_1D::Legendre_1D( int nPoints ) : IntegrationRule(nPoints) {}

Legendre_1D::~Legendre_1D() {}

std::tuple< std::vector<RealVector>, RealVector > 
Legendre_1D::giveIntegrationPointsAndWeights()
{
    const double loc_2 = std::sqrt(3.)/3;
    
    const double loc_3_1 = std::sqrt(3./5.);
    const double wt_3_0 = 8./9.;
    const double wt_3_1 = 5./9.;
    
    const double loc_4_1 = std::sqrt(3./7. - 2./7.*std::sqrt(6./5.));
    const double loc_4_2 = std::sqrt(3./7. + 2./7.*std::sqrt(6./5.));
    const double wt_4_1 = (18. + std::sqrt(30.))/36.;
    const double wt_4_2 = (18. - std::sqrt(30.))/36.;
    
    const double loc_5_1 = std::sqrt(5. - 2.*std::sqrt(10./7.))/3.;
    const double loc_5_2 = std::sqrt(5. + 2.*std::sqrt(10./7.))/3.;
    const double wt_5_0 = 128./225.;
    const double wt_5_1 = (322. + 13.*std::sqrt(70.))/900.;
    const double wt_5_2 = (322. - 13.*std::sqrt(70.))/900.;
    
    std::vector<RealVector> gpLoc(_nIntegrationPoints, RealVector(1));
    RealVector gpWt;
    
    switch (_nIntegrationPoints)
    {
    	case 1:
    		gpLoc[0](0) = 0.;
    		gpWt = {2.};
    		break;
        case 2:
            gpLoc[0](0) = -loc_2;
            gpLoc[1](0) = loc_2;
            gpWt = {1., 1.};
            break;
        case 3:
            gpLoc[0](0) = 0.;
            gpLoc[1](0) = -loc_3_1;
            gpLoc[2](0) = loc_3_1;
            gpWt = {wt_3_0, wt_3_1, wt_3_1};
            break;
        case 4:
            gpLoc[0](0) = -loc_4_1;
            gpLoc[1](0) = loc_4_1;
            gpLoc[2](0) = -loc_4_2;
            gpLoc[3](0) = loc_4_2;
            gpWt = {wt_4_1, wt_4_1, wt_4_2, wt_4_2};
            break;
        case 5:
            gpLoc[0](0) = 0.;
            gpLoc[1](0) = -loc_5_1;
            gpLoc[2](0) = loc_5_1;
            gpLoc[3](0) = -loc_5_2;
            gpLoc[4](0) = loc_5_2;
            gpWt = {wt_5_0, wt_5_1, wt_5_1, wt_5_2, wt_5_2};
            break;
        default:
            throw std::runtime_error("Invalid number of Gauss points specified for 'Legendre_1D'!");
    }
    
    return std::make_tuple(gpLoc, gpWt);
}
