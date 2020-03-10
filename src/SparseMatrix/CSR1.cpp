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

#include "CSR1.hpp"
#include <cmath>
#include <cstdio>
#include "Core/ObjectFactory.hpp"

using namespace broomstyx;

registerBroomstyxObject(SparseMatrix, CSR1)

// Constructor
CSR1::CSR1() {}

// Destructor
CSR1::~CSR1() {}

// Public methods
// ----------------------------------------------------------------------------
void CSR1::addToComponent( int rowNum, int colNum, double val )
{
    
    // Do nothing if component is in lower triangular portion of a symmetric
    // sparse matrix
    if ( _symFlag && colNum < rowNum)
        return;
    
    int rowStart = _prfVec1[rowNum] - 1;
    int nextRow = _prfVec1[rowNum + 1] - 1;

    bool success = false;
    int jdx = rowStart;
    while ( !success )
    {
        if ( _prfVec2[jdx] - 1 == colNum )
        {
            _val(jdx) += val;
            success = true;
        }
        else
            jdx++;

        // Check for instance where specified location cannot be found
        if ( jdx == nextRow )
            throw std::runtime_error("\nAttempted access to non-existing component location in sparse matrix:\n\trow = "
                    + std::to_string(rowNum) + ", col = " + std::to_string(colNum));
    }
}
// ----------------------------------------------------------------------------
void CSR1::atomicAddToComponent( int rowNum, int colNum, double val )
{
    
    // Do nothing if component is in lower triangular portion of a symmetric
    // sparse matrix
    if ( _symFlag && colNum < rowNum)
        return;
    
    int rowStart = _prfVec1[rowNum] - 1;
    int nextRow = _prfVec1[rowNum + 1] - 1;

    bool success = false;
    int jdx = rowStart;
    while ( !success )
    {
        if ( _prfVec2[jdx] - 1 == colNum )
        {
#ifdef _OPENMP
#pragma omp atomic
#endif
            _val(jdx) += val;
            success = true;
        }
        else
            jdx++;

        // Check for instance where specified location cannot be found
        if ( jdx == nextRow )
            throw std::runtime_error("\nAttempted access to non-existing component location in sparse matrix:\n\trow = "
                    + std::to_string(rowNum) + ", col = " + std::to_string(colNum));
    }
}
// ----------------------------------------------------------------------------
void CSR1::finalizeProfile()
{
    // Count number of nonzeros
    _nnz = 0;
    for ( int i = 0; i < _dim1; i++ )
        _nnz += _nz[i].size();

    // Construct arrays
    _prfVec1.assign( _dim1 + 1, 0 );
    _prfVec2.assign( _nnz, 0 );

    int curIdx = 0;
    _prfVec1[0] = 1;
    for ( int i = 0; i < _dim1; i++ )
    {
        _prfVec1[i+1] = _prfVec1[i] + _nz[i].size();
        for ( auto it = _nz[i].begin(); it != _nz[i].end(); ++it )
            _prfVec2[curIdx++] = (*it) + 1;
    }
}
// ----------------------------------------------------------------------------
void CSR1::initializeProfile( int dim1, int dim2 )
{
    _dim1 = dim1;
    _dim2 = dim2;
    _nz.assign( _dim1, std::set<int>() );
}
// ----------------------------------------------------------------------------
void CSR1::insertNonzeroComponentAt( int rowIdx, int colIdx )
{
    if ( colIdx >= rowIdx || !_symFlag )
        _nz[rowIdx].insert(colIdx);
}
// ----------------------------------------------------------------------------
RealVector CSR1::lumpRows()
{
    RealVector b(_dim1);
    
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < _dim1; i++ )
    {
        for ( int j = _prfVec1[i]-1; j < _prfVec1[i+1]-1; j++ )
        {
            b(i) += std::fabs(_val(j));
            if ( _symFlag && _prfVec2[j] != i+1 )
                b(_prfVec2[j]-1) += std::fabs(_val(j));
        }
    }
    
    return b;
}
// ----------------------------------------------------------------------------
std::tuple<int*,int*> CSR1::giveProfileArrays()
{
    return std::make_tuple(_prfVec1.data(),_prfVec2.data());
}
// ----------------------------------------------------------------------------
double* CSR1::giveValArray() 
{ 
    return _val.ptr(); 
}
// ----------------------------------------------------------------------------
void CSR1::initializeValues()
{
    _val.init(_nnz);
}
// ----------------------------------------------------------------------------
void CSR1::printTo( FILE* fp, int n )
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
        int startColIdx = _prfVec1[i];
        int endColIdx = _prfVec1[i+1];
        int diff = endColIdx - startColIdx;
        for ( int j = 0; j < diff; j++)
        {
            colCount++;
            std::fprintf(fp, "\n%*d  %*d  %.*e", width1, i+1, width2, _prfVec2[colCount], n, _val(colCount));
        }
    }
    fprintf(fp, "\n");
}
// ----------------------------------------------------------------------------
RealVector CSR1::times(const RealVector& x)
{
    if ( x.dim() != _dim2 )
        throw std::runtime_error("\nSize mismatch in sparse matrix - vector multiplication.\n\tdim(A) = [ "
                + std::to_string(_dim1) + " x " + std::to_string(_dim2) + " ], dim(B) = " 
                + std::to_string(x.dim()));
    
    RealVector b(_dim1);
    
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < _dim1; i++ )
    {
        for ( int j = _prfVec1[i]-1; j < _prfVec1[i+1]-1; j++ )
        {
            b(i) += _val(j) * x(_prfVec2[j]-1);            
            if ( _symFlag && _prfVec2[j] != i+1 )
                b(_prfVec2[j]-1) += _val(j) * x(i);
        }
    }
    
    return b;
}