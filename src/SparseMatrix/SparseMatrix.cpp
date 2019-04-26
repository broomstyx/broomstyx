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

#include "SparseMatrix.hpp"
#include <cstdio>
#include <cmath>
#include <stdexcept>
#include <string>
#include <tuple>
#include "omp.h"

using namespace broomstyx;

// Constructor
SparseMatrix::SparseMatrix()
    : _dim1(0)
    , _dim2(0)
    , _symFlag(false)
{}

// Destructor
SparseMatrix::~SparseMatrix() {}

// Public methods
// ----------------------------------------------------------------------------
std::tuple<int,int> SparseMatrix::giveMatrixDimensions()
{
    return std::make_tuple( _dim1, _dim2);
}
// ----------------------------------------------------------------------------
int SparseMatrix::giveNumberOfNonzeros() 
{
    return _nnz; 
}
// ----------------------------------------------------------------------------
bool SparseMatrix::isSymmetric()
{
    return _symFlag;
}
// ----------------------------------------------------------------------------
void SparseMatrix::setSymmetryTo(bool true_or_false )
{ 
    _symFlag = true_or_false;
}