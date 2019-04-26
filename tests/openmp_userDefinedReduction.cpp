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

#include <cstdio>
#include <algorithm>
#include <functional>
#include "omp.h"
#include "../source/Util/RealVector.hpp"
#include "../source/Util/reductions.hpp"
#include <iostream>
#include <cmath>

using namespace broomstyx;

void omp_sumRealVector()
{
    std::cout << "********************************" << std::endl;
    std::cout << "Parellel reduction on RealVector" << std::endl;
    std::cout << "********************************" << std::endl;
    
    RealVector S(2);
    S(0) = 0;
    S(1) = -10;
    
#pragma omp parallel for reduction ( + : S )
    for ( int i = 1; i <= 100; i++ )
    {
        S(0) += (double)i;
        S(1) += -(double)i;
    }
    
    S.print("S", 3);
    if ( std::fabs(S(0) - 5050.) > 1.e-15 || std::fabs(S(1) + 5060) > 1.e-15 )
        std::cout << "%TEST_FAILED% time=0 testname=omp_sumRealVector (OpenMP_userDefinedRecution) message=Wrong result." << std::endl;
}

int main( int argc, char** argv )
{
    std::cout << "%SUITE_STARTING% OpenMP_userDefinedRecution" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;

    std::cout << "%TEST_STARTED% test1 (omp_sumRealVector)" << std::endl;
    try
    {
        omp_sumRealVector();
    }
    catch (std::exception& e)
    {
        std::cout << "\n" << e.what() << "\n" << std::endl;
        std::cout << "%TEST_FAILED% time=0 testname=omp_sumRealVector (OpenMP_userDefinedRecution) message=Exception occurred." << std::endl;
    }
    std::cout << "%TEST_FINISHED% time=0 test1 (omp_sumRealVector)" << std::endl;

    std::cout << "%SUITE_FINISHED% time=0" << std::endl;
    
    return (EXIT_SUCCESS);
}