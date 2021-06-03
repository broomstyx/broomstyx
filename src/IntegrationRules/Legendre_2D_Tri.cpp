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

#include "Legendre_2D_Tri.hpp"
#include <cmath>

using namespace broomstyx;

Legendre_2D_Tri::Legendre_2D_Tri( int nPoints ) : IntegrationRule(nPoints) {}

Legendre_2D_Tri::~Legendre_2D_Tri() {}

std::tuple< std::vector<RealVector>, RealVector > 
Legendre_2D_Tri::giveIntegrationPointsAndWeights()
{
    std::vector<RealVector> gpLoc(_nIntegrationPoints, RealVector());
    RealVector gpWt;
    
    switch (_nIntegrationPoints)
    {
    	case 1: // Degree of precision = 1
    		gpLoc[0] = {0.333333333333333, 0.333333333333333};
    		gpWt = {0.5};
    		break;
        case 3: // Degree of precision = 2
            gpLoc[0] = {0.166666666666667, 0.166666666666667};
            gpLoc[1] = {0.666666666666667, 0.166666666666667};
            gpLoc[2] = {0.166666666666667, 0.666666666666667};
            gpWt = {0.166666666666667,
                    0.166666666666667,
                    0.166666666666667};
            break;
        case 4: // Degree of precision = 3
            gpLoc[0] = {0.333333333333333, 0.333333333333333};
            gpLoc[1] = {0.2, 0.2};
            gpLoc[2] = {0.6, 0.2};
            gpLoc[3] = {0.2, 0.6};
            gpWt = {-0.28125, 
                    0.260416666666667,
                    0.260416666666667,
                    0.260416666666667};
            break;
        case 6: // Degree of precision = 4
            gpLoc[0] = {0.816847572980459, 0.091576213509771};
            gpLoc[1] = {0.091576213509771, 0.816847572980459};
            gpLoc[2] = {0.091576213509771, 0.091576213509771};
            gpLoc[3] = {0.108103018168070, 0.445948490915965};
            gpLoc[4] = {0.445948490915965, 0.108103018168070};
            gpLoc[5] = {0.445948490915965, 0.445948490915965};
            gpWt = {0.054975871827661,
                    0.054975871827661,
                    0.054975871827661,
                    0.111690794839006,
                    0.111690794839006,
                    0.111690794839006};
            break;
        default:
            throw std::runtime_error("Invalid number of Gauss points specified for 'Legendre_2D_Tri'!");
    }
    
    return std::make_tuple(gpLoc, gpWt);
}
