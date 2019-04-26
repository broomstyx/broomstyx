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

#include "SneddonProblem.hpp"
#include <cmath>

double sneddonPressure( const TimeData& time )
{
    double p_start = 2.25e5;
    double p_end = 2.25e5;

    double p = p_start + (p_end - p_start)*(time.target - time.start)/(time.end - time.start);

    return p;
}

void sneddonStresses( RealVector coor, const TimeData& time, double sigma[3] )
{
    const double PI = 3.14159265358979323846;

    // crack radius
    double cRad = 2.15;
    double c[3] = {0.0, cRad, -cRad};

    double y = coor(1);
    double x[3], r[3], theta[3];

    // Get radial coordinates
    for ( int i = 0; i < 3; i++ )
    {
        x[i] = coor(0) - c[i];
        r[i] = std::sqrt(x[i]*x[i] + y*y);
        theta[i] = std::atan(y/x[i]);
            if ( x[i] < 0 )
                theta[i] += PI;
    }

    double p = sneddonPressure(time);

    // Calculate sigma_xx and sigma_xy
    double halfSum = p*(r[0]/(std::sqrt(r[1]*r[2]))*std::cos(theta[0] - 0.5*theta[1] - 0.5*theta[2]) - 1.0);
    double halfDiff = p*r[0]*std::sin(theta[0])/cRad*std::pow(cRad*cRad/(r[1]*r[2]),1.5)*std::sin(1.5*(theta[1] + theta[2]));

    sigma[0] = halfSum - halfDiff;
    sigma[1] = halfSum + halfDiff;
    sigma[2] = p*r[0]*std::sin(theta[0])/cRad*std::pow(cRad*cRad/(r[1]*r[2]),1.5)*std::cos(1.5*(theta[1] + theta[2]));
}