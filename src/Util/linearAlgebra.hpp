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

#ifndef LINEARALGEBRA_HPP
#define	LINEARALGEBRA_HPP

#include "RealMatrix.hpp"
#include "RealVector.hpp"

namespace broomstyx
{
    // Matrix addition
    RealMatrix operator+( RealMatrix&& A, RealMatrix&& B );
    RealMatrix operator+( RealMatrix&& A, const RealMatrix& B );
    RealMatrix operator+( const RealMatrix& A, RealMatrix&& B );
    RealMatrix operator+( const RealMatrix& A, const RealMatrix& B );
    
    // Vector addition
    RealVector operator+( RealVector&& A, RealVector&& B );
    RealVector operator+( RealVector&& A, const RealVector& B );
    RealVector operator+( const RealVector& A, RealVector&& B );
    RealVector operator+( const RealVector& A, const RealVector& B );
    
    // Matrix subtraction
    RealMatrix operator-( RealMatrix&& A, RealMatrix&& B );
    RealMatrix operator-( RealMatrix&& A, const RealMatrix& B );
    RealMatrix operator-( const RealMatrix& A, RealMatrix&& B );
    RealMatrix operator-( const RealMatrix& A, const RealMatrix& B );
    
    // Vector subtraction
    RealVector operator-( RealVector&& A, RealVector&& B );
    RealVector operator-( RealVector&& A, const RealVector& B );
    RealVector operator-( const RealVector& A, RealVector&& B );
    RealVector operator-( const RealVector& A, const RealVector& B );
    
    // Scalar multiplication for matrices
    RealMatrix operator*( RealMatrix&& A, double b );
    RealMatrix operator*( const RealMatrix& A, double b );
    RealMatrix operator*( double a, RealMatrix&& B );
    RealMatrix operator*( double a, const RealMatrix& B );
    
    // Scalar division for matrices
    RealMatrix operator/( RealMatrix&& A, double b );
    RealMatrix operator/( const RealMatrix& A, double b );
    
    // Scalar multiplication for vectors
    RealVector operator*( RealVector&& A, double b );
    RealVector operator*( const RealVector& A, double b );
    RealVector operator*( double a, RealVector&& B );
    RealVector operator*( double a, const RealVector& B );
    
    // Scalar division for vectors
    RealVector operator/( RealVector&& A, double b );
    RealVector operator/( const RealVector& A, double b );
    
    // Matrix multiplicaton
    RealMatrix operator*( RealMatrix&& A, RealMatrix&& B );
    RealMatrix operator*( RealMatrix&& A, const RealMatrix& B );
    RealMatrix operator*( const RealMatrix& A, RealMatrix&& B );
    RealMatrix operator*( const RealMatrix& A, const RealMatrix& B );
    
    // Matrix-vector multiplication
    RealVector operator*( RealMatrix&& A, RealVector&& B );
    RealVector operator*( RealMatrix&& A, const RealVector& B);
    RealVector operator*( const RealMatrix& A, RealVector&& B );
    RealVector operator*( const RealMatrix& A, const RealVector& B);
    
    // Vector-matrix multiplication
    RealVector operator*( RealVector&& A, RealMatrix&& B );
    RealVector operator*( RealVector&& A, const RealMatrix& B );
    RealVector operator*( const RealVector& A, RealMatrix&& B );
    RealVector operator*( const RealVector& A, const RealMatrix& B );
    
    // Matrix transpose
    RealMatrix trp( RealMatrix&& A );
    RealMatrix trp( const RealMatrix& A );
    
    // Simplification (state resolution)
    RealMatrix simplify( RealMatrix&& A );
    RealVector simplify( RealVector&& A );
    
    // Matrix inverse
    RealMatrix inv( RealMatrix&& A );
    RealMatrix inv( const RealMatrix& B);
}

#endif	/* LINEARALGEBRA_HPP */