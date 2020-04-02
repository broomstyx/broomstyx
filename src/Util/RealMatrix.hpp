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

#ifndef REALMATRIX_HPP
#define	REALMATRIX_HPP

#include <config.h>

#include <cstdio>
#include <stdexcept>
#include <initializer_list>

#ifdef HAVE_MKL
    #include "mkl_cblas.h"
#else
    #include "cblas.h"
#endif

namespace broomstyx
{        
    class RealMatrix final
    {
        friend RealMatrix operator*( RealMatrix&& A, double b );
        friend RealMatrix operator*( double a, RealMatrix&& B );
        friend RealMatrix operator/( RealMatrix&& A, double b );
        friend RealMatrix trp( RealMatrix&& A );
        
    public:
        // Default constructor
        RealMatrix()
            : _dim1(0)
            , _dim2(0)
            , _ptr(nullptr)
            , _ownsItsPointer(true)            
            , _isScaled(false)
            , _scaling(1.)
            , _isTransposed(false)
        {}

        // Constructor with initial memory allocation
        RealMatrix( int dim1, int dim2 )\
            : _ownsItsPointer(true)
            , _isScaled(false)
            , _scaling(1.)
            , _isTransposed(false)

        {
#ifndef NDEBUG
            if ( dim1 < 1 )
                throw std::runtime_error("Cannot construct RealMatrix with dim1 = " + std::to_string(dim1));
            if ( dim2 < 1 )
                throw std::runtime_error("Cannot construct RealMatrix with dim1 = " + std::to_string(dim2));
#endif
            _dim1 = dim1;
            _dim2 = dim2;
            _ptr = new double[ _dim1 * _dim2 ]();
        }
        
        // Constructor with initializer list
        RealMatrix( std::initializer_list< std::initializer_list<double>> initList )
            : _ownsItsPointer(true)
            , _isScaled(false)
            , _scaling(1.)
            , _isTransposed(false)

        {
            _dim1 = (int)initList.size();
            _dim2 = (int)(*(initList.begin())).size();
#ifndef NDEBUG
            for ( int i = 1; i < _dim1; i++ )
                if ( _dim2 != (int)(*(initList.begin() + i)).size() )
                    throw std::runtime_error("RealMatrix initializer list has inconsistent dimensions!");
#endif
            _ptr = new double[ _dim1 * _dim2 ];
            
            for ( int i = 0; i < _dim1; i++ )
            {
#pragma GCC ivdep
                for ( int j = 0; j < _dim2; j++ )
                    _ptr[ j*_dim1 + i ] = *((*(initList.begin() + i)).begin() + j);
            }
        }
        
        // Copy constructor
        RealMatrix( const RealMatrix& source )
            : _dim1(source._dim1)
            , _dim2(source._dim2)
            , _ptr(new double[ _dim1 * _dim2 ])
            , _ownsItsPointer(true)
            , _isScaled(source._isScaled)
            , _scaling(source._scaling)
            , _isTransposed(source._isTransposed)
        {
#ifdef VERBOSE_REALMATRIX_CONSTRUCTION
            std::printf("...RealMatrix Copy constructor called.\n");
            std::fflush(stdout);
#endif
            std::copy(source._ptr, source._ptr + _dim1 * _dim2, _ptr);
            this->simplify();
        }

        // Move constructor
        RealMatrix( RealMatrix&& source )
            : _dim1(source._dim1)
            , _dim2(source._dim2)
            , _ptr(source._ptr)
            , _ownsItsPointer(source._ownsItsPointer)
            , _isScaled(source._isScaled)
            , _scaling(source._scaling)
            , _isTransposed(source._isTransposed)
        {
#ifdef VERBOSE_REALMATRIX_CONSTRUCTION
            std::printf("...RealMatrix Move constructor called.\n");
            std::fflush(stdout);
#endif
            source._dim1 = 0;
            source._dim2 = 0;
            source._ptr = nullptr;
        }

        // Destructor
        virtual ~RealMatrix()
        {
            if ( _ptr && _ownsItsPointer ) 
                delete[] _ptr;
        }

        // Assignment with initializer list
        RealMatrix& operator= ( std::initializer_list< std::initializer_list<double>> initList )
        {
            _dim1 = (int)initList.size();
            _dim2 = (int)(*(initList.begin())).size();
#ifndef NDEBUG
            for ( int i = 1; i < _dim1; i++ )
                if ( _dim2 != (int)(*(initList.begin() + i)).size() )
                    throw std::runtime_error("RealMatrix initializer list has inconsistent dimensions!");
#endif
            if ( _ptr && _ownsItsPointer )
                delete[] _ptr;
                
            _ptr = new double[ _dim1 * _dim2 ];
            
            for ( int i = 0; i < _dim1; i++ )
            {
#pragma GCC ivdep
                for ( int j = 0; j < _dim2; j++ )
                    _ptr[ j*_dim1 + i ] = *((*(initList.begin() + i)).begin() + j);
            }
            
            _ownsItsPointer = true;
            _isScaled = false;
            _scaling = 1.;
            _isTransposed = false;
            
            return *this;
        }
        
        // Copy assignment operator
        RealMatrix& operator= ( const RealMatrix& source )
        {
#ifdef VERBOSE_REALMATRIX_CONSTRUCTION
            std::printf("...RealMatrix Copy assignment called.\n");
            std::fflush(stdout);
#endif
            // Check for self-assignment
            if ( this == &source )
            {
                this->simplify();
                return *this;
            }

            if ( _ptr && _ownsItsPointer )
                delete[] _ptr;

            _dim1 = source._dim1;
            _dim2 = source._dim2;
            _ptr = new double[ _dim1 * _dim2 ];
            std::copy(source._ptr, source._ptr + _dim1 * _dim2, _ptr);
            
            _ownsItsPointer = true;
            _isScaled = source._isScaled;
            _scaling = source._scaling;
            _isTransposed = source._isTransposed;
            
            this->simplify();
                
            return *this;
        }

        // Move assignment operator
        RealMatrix& operator= ( RealMatrix&& source )
        {
#ifdef VERBOSE_REALMATRIX_CONSTRUCTION
            std::printf("...RealMatrix Move assignment called.\n");
            std::fflush(stdout);
#endif
            if ( _ptr && _ownsItsPointer )
                delete[] _ptr;
            
            _dim1 = source._dim1;
            _dim2 = source._dim2;
            _ptr = source.ptr();
            _ownsItsPointer = source._ownsItsPointer;
            _isScaled = source._isScaled;
            _scaling = source._scaling;
            _isTransposed = source._isTransposed;
            
            source._dim1 = 0;
            source._dim2 = 0;
            source._ptr = nullptr;    
            
            this->simplify();
            
            return *this;
        }
        
        // Addition to self
        RealMatrix& operator+= ( const RealMatrix& source )
        {
#ifndef NDEBUG
            if ( source._dim1 != _dim1 || source._dim2 != _dim2 )
                throw std::runtime_error("\nSize mismatch in operands for operator '+='!\n\tdim(A) = [ " + std::to_string(_dim1) + " x " 
                    + std::to_string(_dim2) + " ], dim(B) = [ " + std::to_string(source._dim1) + " x " + std::to_string(source._dim2) + " ]");
#endif
            this->simplify();
            cblas_daxpy(_dim1*_dim2, source._scaling, source._ptr, 1, _ptr, 1);

            return *this;
        }
        
        RealMatrix& operator+= ( RealMatrix&& source )
        {
            if ( source._isTransposed )
                source.simplify();
#ifndef NDEBUG
            if ( source._dim1 != _dim1 || source._dim2 != _dim2 )
                throw std::runtime_error("\nSize mismatch in operands for operator '+='!\n\tdim(A) = [ " + std::to_string(_dim1) + " x " 
                    + std::to_string(_dim2) + " ], dim(B) = [ " + std::to_string(source._dim1) + " x " + std::to_string(source._dim2) + " ]");
#endif
            this->simplify();
            cblas_daxpy(_dim1*_dim2, source._scaling, source._ptr, 1, _ptr, 1);

            return *this;
        }
        
        // Subtraction from self
        RealMatrix& operator-= ( const RealMatrix& source )
        {
#ifndef NDEBUG
            if ( source._dim1 != _dim1 || source._dim2 != _dim2 )
                throw std::runtime_error("\nSize mismatch in operands for operator '+='!\n\tdim(A) = [ " + std::to_string(_dim1) + " x " 
                    + std::to_string(_dim2) + " ], dim(B) = [ " + std::to_string(source._dim1) + " x " + std::to_string(source._dim2) + " ]");
#endif
            this->simplify();
            cblas_daxpy(_dim1*_dim2, -source._scaling, source._ptr, 1, _ptr, 1);

            return *this;
        }
        
        RealMatrix& operator-= ( RealMatrix&& source )
        {
            if ( source._isTransposed )
                source.simplify();
#ifndef NDEBUG
            if ( source._dim1 != _dim1 || source._dim2 != _dim2 )
                throw std::runtime_error("\nSize mismatch in operands for operator '+='!\n\tdim(A) = [ " + std::to_string(_dim1) + " x " 
                    + std::to_string(_dim2) + " ], dim(B) = [ " + std::to_string(source._dim1) + " x " + std::to_string(source._dim2) + " ]");
#endif
            this->simplify();
            cblas_daxpy(_dim1*_dim2, -source._scaling, source._ptr, 1, _ptr, 1);

            return *this;
        }

        // NOTE: Storage of matrix components uses column-major format
        // Matric component access as variableName(i,j)
        double& operator()( int idx1, int idx2 )
        {
#ifndef NDEBUG
            if ( _isScaled )
                throw std::runtime_error("\nDenied access to address of RealMatrix component due to non-unity scaling!");
            
            if ( _isTransposed )
                throw std::runtime_error("\nDenied access to address of RealMatrix component due to unresolved transposition!");
            
            if ( !_ptr )
                throw std::runtime_error("\nCannot access RealMatrix component (" + std::to_string(idx1) + "," + std::to_string(idx2) + ") -- "
                        + "\nmatrix is not initialized.");

            if ( idx1 < 0 || idx1 >= _dim1 || idx2 < 0 || idx2 >= _dim2 )
                throw std::runtime_error("\nCannot access RealMatrix component (" + std::to_string(idx1) + "," + std::to_string(idx2) + ")! "
                        + "Valid range is (0-" + std::to_string(_dim1 - 1) + ",0-" + std::to_string(_dim2 - 1) + ").");
#endif
            return _ptr[ idx2 * _dim1 + idx1 ];
        }

        double operator()( int idx1, int idx2 ) const
        {
#ifndef NDEBUG
            if ( _isScaled )
                throw std::runtime_error("\nDenied access to address of RealMatrix component due to non-unity scaling!");
            
            if ( _isTransposed )
                throw std::runtime_error("\nDenied access to address of RealMatrix component due to unresolved transposition!");

            if ( !_ptr )
                throw std::runtime_error("\nCannot access RealMatrix component (" + std::to_string(idx1) + "," + std::to_string(idx2) + ") -- "
                        + "\nmatrix is not initialized.");
            
            if ( idx1 < 0 || idx1 >= _dim1 || idx2 < 0 || idx2 >= _dim2 )
                throw std::runtime_error("\nCannot access RealMatrix component (" + std::to_string(idx1) + "," + std::to_string(idx2) + ")! "
                        + "Valid range is (0-" + std::to_string(_dim1 - 1) + ",0-" + std::to_string(_dim2 - 1) + ").");
#endif
            return _ptr[ idx2 * _dim1 + idx1 ];
        }

        // Dimensions of matrix
        int dim1() const { return _dim1; }
        int dim2() const { return _dim2; }

        // Erase contents
        void erase()
        {
            if ( _ptr && _ownsItsPointer )
                delete[] _ptr;
                
            _ptr = nullptr;
            _dim1 = 0;
            _dim2 = 0;
            _scaling = 1.;
            _isScaled = false;
            _ownsItsPointer = true;
        }

        // Allocate/reallocate space for matrix of size (dim1 x dim2), with all
        // values set to zero
        void init( int dim1, int dim2 )
        {
#ifndef NDEBUG
            if ( dim1 < 1 )
                throw std::runtime_error("\nCannot initialize RealMatrix with dim1 = " + std::to_string(dim1));
            if ( dim2 < 1 )
                throw std::runtime_error("\nCannot initialize RealMatrix with dim2 = " + std::to_string(dim2));
#endif
            if ( _ptr)
                delete[] _ptr;

            _dim1 = dim1;
            _dim2 = dim2;
            _ptr = new double[ _dim1 * _dim2 ]();
            _isScaled = false;
            _scaling = 1.;
            _ownsItsPointer = true;
        }

        // Query scaling status
        bool isScaled() const { return _isScaled; }
        
        // Query transposition status
        bool isTransposed() const { return _isTransposed; }
        
        // Query pointer ownership
        bool ownsItsPointer() const { return _ownsItsPointer; }
        
        // Show contents of matrix in scientific precision
        void print( const char *s, int n ) const 
        {
            std::printf("\nRealMatrix %s:\n\n", s);
            
            if ( _isScaled )
                std::printf("...is scaled with factor %e\n", _scaling);
            else
                std::printf("...is not scaled (scaling = %e)\n", _scaling);
            
            if ( _ownsItsPointer )
                std::printf("...owns its pointer\n");
            else
                std::printf("...does not own its pointer\n");
            
            if ( _isTransposed )
                std::printf("...is transposed\n");
            else
                std::printf("...is not transposed\n");

            if ( !_ptr )
                std::printf("...is empty\n");
            else
            {
                std::printf("...size = %d x %d\n\n", _dim1, _dim2);
                for ( int i = 0; i < _dim1; i++ ) 
                {
                    for ( int j = 0; j < _dim2; j++ )
                        std::printf("%*.*e", n+10, n, _ptr[ j * _dim1 + i ]);
                    std::printf("\n");
                }
            }
            std::printf("\n");
        }

        // Give pointer to matrix components
        double* ptr() const { return _ptr; }
        
        // Scale by real number
        RealMatrix scaleBy( double val ) const
        {
            RealMatrix A;
            
            A._dim1 = _dim1;
            A._dim2 = _dim2;
            A._isScaled = true;
            A._scaling = _scaling*val;
            A._ptr = _ptr;
            A._ownsItsPointer = false;
            A._isTransposed = _isTransposed;
            
            return A;
        }
        
        // Scaling factor
        double scaling() const { return _scaling; }
        
        // Resolve transposition, pointer ownership and scaling
        void simplify()
        {
            if ( _isTransposed )
            {
#ifdef VERBOSE_REALMATRIX_STATE_RESOLUTION
                std::printf("...Resolving transposition of RealMatrix.\n");
                std::fflush(stdout);
#endif
                double* t_ptr = new double[ _dim1*_dim2 ];
                for ( int i = 0; i < _dim1; i++ )
                {
#pragma GCC ivdep
                    for ( int j = 0; j < _dim2; j++ )
                        t_ptr[ i*_dim2 + j ] = _ptr[ j * _dim1 + i ];
                }
                
                double dum = _dim1;
                _dim1 = _dim2;
                _dim2 = dum;
                
                if ( _ownsItsPointer )
                    delete[] _ptr;
                    
                _ptr = t_ptr;
                _isTransposed = false;
                _ownsItsPointer = true;
            }
            else if ( !_ownsItsPointer )
            {
#ifdef VERBOSE_REALMATRIX_STATE_RESOLUTION
                std::printf("...Resolving pointer ownership of RealMatrix.\n");
                std::fflush(stdout);
#endif
                double* temp = _ptr;
                _ptr = new double[ _dim1*_dim2 ];
                std::copy(temp, temp + _dim1*_dim2, _ptr);
                _ownsItsPointer = true;
            }
            
            if ( _isScaled )
            {
#ifdef VERBOSE_REALMATRIX_STATE_RESOLUTION
                std::printf("...Resolving scaling of RealMatrix.\n");
                std::fflush(stdout);
#endif
                cblas_dscal(_dim1*_dim2, _scaling, _ptr, 1);
                _isScaled = false;
                _scaling = 1.;
            }
        }
        
        // Transpose (only marks matrix as transposed)
        RealMatrix transpose() const
        {
            RealMatrix A;
            
            A._dim1 = _dim1;
            A._dim2 = _dim2;
            A._isScaled = _isScaled;
            A._scaling = _scaling;
            A._ptr = _ptr;
            A._ownsItsPointer = false;
            
            if ( _isTransposed )
                A._isTransposed = false;
            else
                A._isTransposed = true;
            
            return A;
        }
        
    private:
        int     _dim1;
        int     _dim2;
        double* _ptr;
        bool    _ownsItsPointer;
        bool    _isScaled;
        double  _scaling;
        bool    _isTransposed;
    };
}

#endif	/* REALMATRIX_HPP */