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

#ifndef TRIANGLE_P2_HPP
#define TRIANGLE_P2_HPP

/***********************************************************************
 * Basis functions computed according to the following nodal ordering:
 *
 *                               3
 *                              / \
 *                             /   \
 *                            6     5
 *                           /       \
 *                          /         \                      
 *                         1 --- 4 --- 2
 * 
 ***********************************************************************/

#include "ScalarBasisFunction.hpp"

namespace broomstyx
{
    class Triangle_P2 final : public ScalarBasisFunction
    {
    public:
        Triangle_P2();
        virtual ~Triangle_P2();
        
        RealVector giveBasisFunctionsAt( const RealVector& coor ) override;
        std::vector<RealVector> giveBasisFunctionDerivativesAt( const RealVector& coor ) override;
    };
}

#endif /* TRIANGLE_P2_HPP */