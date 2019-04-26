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

#ifndef REDUCTIONS_HPP
#define REDUCTIONS_HPP

#include <omp.h>
#include <algorithm>
#include <functional>
#include "RealVector.hpp"

namespace broomstyx
{
#pragma omp declare reduction( + : RealVector : \
    std::transform(omp_in.ptr(), omp_in.ptr() + omp_in.dim(), \
    omp_out.ptr(), omp_out.ptr(), std::plus<double>())) initializer(omp_priv = RealVector(omp_orig.dim()))

#pragma omp declare reduction( + : RealMatrix : \
    std::transform(omp_in.ptr(), omp_in.ptr() + omp_in.dim1()*omp_in.dim2(), \
    omp_out.ptr(), omp_out.ptr(), std::plus<double>())) initializer(omp_priv = RealMatrix(omp_orig.dim1(), omp_orig.dim2()))
}
#endif /* REDUCTIONS_HPP */