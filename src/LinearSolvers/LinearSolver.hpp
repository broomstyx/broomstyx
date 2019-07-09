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

#ifndef LINEARSOLVER_HPP
#define	LINEARSOLVER_HPP

#include <cstdio>
#include <string>

namespace broomstyx
{
    class RealVector;
    class SparseMatrix;

    class LinearSolver
    {
    public:
        LinearSolver();
        virtual ~LinearSolver();
        
        // Disable copy constructor and assignment operator
        LinearSolver( const LinearSolver& ) = delete;
        LinearSolver& operator=( const LinearSolver& ) = delete;

        virtual void allocateInternalMemoryFor( SparseMatrix* coefMat );
        virtual RealVector backSubstitute( SparseMatrix* coefMat, RealVector& rhs );
        virtual void factorize( SparseMatrix* coefMat );
        virtual bool giveSymmetryOption();
        virtual void initialize();
        virtual void clearInternalMemory();
        virtual void setInitialGuessTo( RealVector& initGuess );
        virtual bool takesInitialGuess();
        
        virtual std::string giveRequiredMatrixFormat() = 0;
        virtual void        readDataFrom( FILE* fp ) = 0;
        virtual RealVector  solve( SparseMatrix* coefMat, RealVector& rhs ) = 0;
    };
}

#endif	/* LINEARSOLVER_HPP */