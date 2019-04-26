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

#include "ScalarBasisFunction.hpp"

using namespace broomstyx;

ScalarBasisFunction::ScalarBasisFunction() {}

ScalarBasisFunction::~ScalarBasisFunction() {}

void ScalarBasisFunction::error_unimplemented( std::string method )
{
    throw std::runtime_error("\nError: Call to unimplemented method '"
            + _name + "::" + method + "' encountered!\n");
}

RealVector ScalarBasisFunction::giveBasisFunctionsAt( const RealVector& coor )
{
    this->error_unimplemented("giveShapeFunctionsAt(...)");
    
    // Return statement is only meant to suppress warnings during compilation.
    RealVector dummy;
    return dummy;
}

std::vector<RealVector>
ScalarBasisFunction::giveBasisFunctionDerivativesAt( const RealVector& coor )
{
    this->error_unimplemented("giveShapeFunctionDerivativesAt(...)");
    
    // Return statement is only meant to suppress warnings during compilation.
    std::vector<RealVector> dummy;
    return dummy;
}

std::vector< std::vector<RealVector> >
ScalarBasisFunction::giveBasisFunction2ndDerivativesAt( const RealVector& coor )
{
    this->error_unimplemented("giveShapeFunction2ndDerivativesAt(...)");
    
    // Return statement is only meant to suppress warnings during compilation.
    std::vector< std::vector<RealVector> > dummy;
    return dummy;
}