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

#ifndef ISTLMAT_CPP
#define	ISTLMAT_CPP

#include "ISTLMat.hpp"
#include "../Core/ObjectFactory.hpp"


#if HAVE_DUNE_ISTL

using namespace broomstyx;

registerBroomstyxObject(SparseMatrix, ISTLMat)


ISTLMat::ISTLMat() : _matrix()
{
}

void
ISTLMat::addToComponent( int rowNum, int colNum, double val )
{
    // Do nothing if component is in lower triangular portion of a symmetric
    // sparse matrix
    //if ( _symFlag && colNum < rowNum)
    //    return;

    (*_matrix)[ rowNum ][ colNum ][ 0 ][ 0 ] += val;
}

void
ISTLMat::atomicAddToComponent( int rowNum, int colNum, double val )
{
    // Do nothing if component is in lower triangular portion of a symmetric
    // sparse matrix
    //if ( _symFlag && colNum < rowNum)
    //    return;

    (*_matrix)[ rowNum ][ colNum ][ 0 ][ 0 ] += val;
}

void
ISTLMat::finalizeProfile()
{
    _matrix->compress();
    _nnz = _matrix->nonzeroes();
}

std::tuple< int*,int* >
ISTLMat::giveProfileArrays() { return std::tuple< int*,int* > (nullptr, nullptr ) ; }

double*
ISTLMat::giveValArray() { std::abort(); return nullptr; }

void
ISTLMat::initializeProfile( int dim1, int dim2 )
{
    _symFlag = false;
    _dim1 = dim1;
    _dim2 = dim2;
    // TODO: estimate nnz, here 30
    _matrix.reset( new MatrixType( dim1, dim2, 15, 1.4, MatrixType::implicit ) );
}

void
ISTLMat::initializeValues()
{
    (*_matrix) = 0.0;
}

void
ISTLMat::insertNonzeroComponentAt( int rowIdx, int colIdx)
{
    _matrix->entry( rowIdx, colIdx ) = 0.0;
}

RealVector
ISTLMat::lumpRows()
{
    RealVector y( _matrix->N() );
    const auto endrow = _matrix->end();
    for( auto row = _matrix->begin(); row != endrow; ++row )
    {
        const auto colend = row->end();
        double sum = 0.0;
        for( auto col = row->begin(); col != colend; ++ col )
        {
            sum += (*col)[ 0 ][ 0 ];
        }
        y( row.index() ) = sum;
    }
    return y;
}


void
ISTLMat::printTo( FILE* fp, int n )
{
    int width1 = (int)std::log10((double)_dim2);
    int width2 = (int)std::log10((double)_dim2);

    std::fprintf(fp, "\nnRows = %d", _dim1);
    std::fprintf(fp, "\nnCols = %d", _dim2);
    std::fprintf(fp, "\nnNonzeros = %d", _nnz);
    if ( _symFlag )
        std::fprintf(fp, "\nSymmetric\n");
    else
        std::fprintf(fp, "\nNonsymmetric\n");

    int colCount = -1;

    for ( int i = 0; i < _dim1; i++)
    {
        const auto end = (*_matrix)[ i ].end();
        for( auto j = (*_matrix)[ i ].begin(); j!=end; ++j )
        {
            colCount = j.index();
            double _val = (*j)[0][0];
            std::fprintf(fp, "\n%*d  %*d  %.*e", width1, i, width2, colCount, n, _val);
        }
    }
    fprintf(fp, "\n");

}

RealVector
ISTLMat::times( const RealVector& x )
{
    assert( x.dim() == _matrix->M() );
    RealVector y( _matrix->N() );
    const auto endrow = _matrix->end();
    for( auto row = _matrix->begin(); row != endrow; ++row )
    {
        const auto colend = row->end();
        double sum = 0.0;
        for( auto col = row->begin(); col != colend; ++ col )
        {
            sum += (*col)[ 0 ][ 0 ] * x( col.index() );
        }
        y( row.index() ) = sum;
    }
    return y;
}

#endif // HAVE_DUNE_ISTL

#endif	/* ISTLMAT_HPP */
