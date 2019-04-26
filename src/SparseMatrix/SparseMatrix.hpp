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

#ifndef SPARSEMATRIX_HPP
#define	SPARSEMATRIX_HPP

#include <cstdio>
#include <set>
#include <tuple>
#include <vector>
#include "../Util/RealVector.hpp"

namespace broomstyx
{
    class SparseMatrix
    {
    public:
        SparseMatrix();
        virtual ~SparseMatrix();

        std::tuple< int,int > giveMatrixDimensions();
        int     giveNumberOfNonzeros();
        bool    isSymmetric();
        void    setSymmetryTo( bool true_or_false );
        
        virtual void addToComponent( int rowNum, int colNum, double val ) = 0;
        virtual void atomicAddToComponent( int rowNum, int colNum, double val ) = 0;
        virtual void finalizeProfile() = 0;
        virtual std::tuple< int*,int* > giveProfileArrays() = 0;
        virtual double* giveValArray() = 0;
        virtual void initializeProfile( int dim1, int dim2 ) = 0;
        virtual void initializeValues() = 0;
        virtual void insertNonzeroComponentAt( int rowIdx, int colIdx ) = 0;
        virtual RealVector lumpRows() = 0;
        virtual void printTo( FILE* fp, int n ) = 0;
        virtual RealVector times( const RealVector& x ) = 0;

    protected:
        int _dim1;
        int _dim2;
        int _nnz;
        bool _symFlag;
    };
}

#endif	/* SPARSEMATRIX_HPP */