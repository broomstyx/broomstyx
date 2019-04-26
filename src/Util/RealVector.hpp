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

#ifndef REALVECTOR_HPP
#define	REALVECTOR_HPP

#include "../../config.h"
#include <cstdio>
#include <stdexcept>
#include <initializer_list>
#include "RealMatrix.hpp"

#ifdef HAVE_MKL
    #include "mkl_cblas.h"
#else
    #include "cblas.h"
#endif

namespace broomstyx
{
    class RealVector final
    {
        friend RealVector operator*( RealVector&& A, double b );
        friend RealVector operator*( double a, RealVector&& B );
        friend RealVector operator/( RealVector&& A, double b );
        
    public:
        // Default constructor
        RealVector()
            : _dim(0)
            , _ptr(nullptr)
            , _ownsItsPointer(true)            
            , _isScaled(false)
            , _scaling(1.)
        {}

        // Constructor with initial memory allocation
        RealVector( int dim )
            : _ownsItsPointer(true)
            , _isScaled(false)
            , _scaling(1.)
        {
#ifndef NDEBUG
            if ( dim < 1 )
                throw std::runtime_error("\nCannot construct RealVector with dim = " + std::to_string(dim));
#endif

            _dim = dim;
            _ptr = new double[ _dim ]();
            
        }
        
        // Constructor with initializer list
        RealVector( std::initializer_list<double> initList )
            : _ownsItsPointer(true)
            , _isScaled(false)
            , _scaling(1.)
        {
            _dim = (int)initList.size();
            _ptr = new double[ _dim ];
            
            std::copy(initList.begin(), initList.begin() + _dim, _ptr);
        }

        // Copy constructor
        RealVector( const RealVector& source )
            : _dim(source._dim)
            , _ownsItsPointer(true)
            , _isScaled(source._isScaled)
            , _scaling(source._scaling)
        {
#ifdef VERBOSE_REALVECTOR_CONSTRUCTION
            std::printf("...RealVector Copy constructor called.\n");
            std::fflush(stdout);
#endif
            _ptr = new double[ _dim ];
            std::copy(source._ptr, source._ptr + _dim, _ptr);
            
            this->simplify();
        }

        // Move constructor
        RealVector( RealVector&& source )
            : _dim(source._dim)
            , _ptr(source._ptr)
            , _ownsItsPointer(source._ownsItsPointer)
            , _isScaled(source._isScaled)
            , _scaling(source._scaling)
        {
#ifdef VERBOSE_REALVECTOR_CONSTRUCTION
            std::printf("...RealVector Move constructor called.\n");
            std::fflush(stdout);
#endif
            source._ptr = nullptr;
            source._dim = 0;
        }

        // Destructor
        virtual ~RealVector()
        {
            if ( _ptr && _ownsItsPointer )
                delete[] _ptr;
        }

        // Assignment with initializer list
        RealVector& operator= ( std::initializer_list<double> initList )
        {
            if ( _ptr && _ownsItsPointer && _dim != (int)initList.size() )
            {
                delete[] _ptr;
                _ptr = nullptr;
            }
            
            _dim = (int)initList.size();
            if ( !_ptr || !_ownsItsPointer )
                _ptr = new double[ _dim ];
            
            _ownsItsPointer = true;
            _isScaled = false;
            _scaling = 1.;
            
            std::copy(initList.begin(), initList.begin() + _dim, _ptr);
            
            return *this;
        }
        
        // Copy assignment operator
        RealVector& operator= ( const RealVector& source )
        {
#ifdef VERBOSE_REALVECTOR_CONSTRUCTION
            std::printf("...RealVector Copy assignment called.\n");
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

            _dim = source._dim;
            _ptr = new double[ _dim ];
            std::copy(source._ptr, source._ptr + _dim, _ptr);
            
            _ownsItsPointer = true;
            _isScaled = source._isScaled;
            _scaling = source._scaling;
            
            this->simplify();
            
            return *this;
        }

        // Move assignment operator
        RealVector& operator= ( RealVector&& source )
        {
#ifdef VERBOSE_REALVECTOR_CONSTRUCTION
            std::printf("...RealVector Move assignment called.\n");
            std::fflush(stdout);
#endif
            if ( _ptr )
                delete[] _ptr;
                
            _dim = source._dim;
            _ptr = source._ptr;
            _ownsItsPointer = source._ownsItsPointer;
            _isScaled = source._isScaled;
            _scaling = source._scaling;
            
            source._dim = 0;
            source._ptr = nullptr;
            this->simplify();
            
            return *this;
        }

        // Addition to self
        RealVector& operator+= ( const RealVector& source )
        {
#ifndef NDEBUG
            if ( source._dim != _dim )
                throw std::runtime_error("\nSize mismatch in operands for operator '+='!\n\tdim(A) = " + std::to_string(_dim) + ", dim(B) = " + std::to_string(source._dim));
#endif
            this->simplify();
            cblas_daxpy(_dim, source._scaling, source._ptr, 1, _ptr, 1);

            return *this;
        }
        
        RealVector& operator+= ( RealVector&& source )
        {
#ifndef NDEBUG
            if ( source._dim != _dim )
                throw std::runtime_error("\nSize mismatch in operands for operator '+='!\n\tdim(A) = " + std::to_string(_dim) + ", dim(B) = " + std::to_string(source._dim));
#endif
            this->simplify();
            cblas_daxpy(_dim, source._scaling, source._ptr, 1, _ptr, 1);

            return *this;
        }
        
        // Subtraction from self
        RealVector& operator-= ( const RealVector& source )
        {
#ifndef NDEBUG
            if ( source._dim != _dim )
                throw std::runtime_error("\nSize mismatch in operands for operator '-='!\n\tdim(A) = " + std::to_string(_dim) + ", dim(B) = " + std::to_string(source._dim));
#endif
            this->simplify();
            cblas_daxpy(_dim, -source._scaling, source._ptr, 1, _ptr, 1);

            return *this;
        }
        
        RealVector& operator-= ( RealVector&& source )
        {
#ifndef NDEBUG
            if ( source._dim != _dim )
                throw std::runtime_error("\nSize mismatch in operands for operator '-='!\n\tdim(A) = " + std::to_string(_dim) + ", dim(B) = " + std::to_string(source._dim));
#endif
            this->simplify();
            cblas_daxpy(_dim, -source._scaling, source._ptr, 1, _ptr, 1);

            return *this;
        }
        
        // Vector component access as variableName(i)
        double& operator()( int idx )
        {
#ifndef NDEBUG
            if ( _isScaled )
                throw std::runtime_error("\nDenied access to address of RealVector component due to non-unity scaling!");

            if ( !_ptr )
                throw std::runtime_error("\nCannot access RealVector component (" + std::to_string(idx) + ") -- vector is not initialized.");

            if ( idx < 0 || idx >= _dim )
                throw std::runtime_error("\nCannot access RealVector component (" + std::to_string(idx) + ")! Valid range is (0-" 
                        + std::to_string(_dim - 1) + ").");
#endif
            return _ptr[ idx ];
        }

        double operator()( int idx ) const
        {
#ifndef NDEBUG
            if ( _isScaled )
                throw std::runtime_error("\nDenied access to address of RealVector component due to non-unity scaling!");

            if ( !_ptr )
                throw std::runtime_error("\nCannot access RealVector component (" + std::to_string(idx) + ") -- vector is not initialized.");

            if ( idx < 0 || idx >= _dim )
                throw std::runtime_error("\nCannot access RealVector component (" + std::to_string(idx) + ")! Valid range is (0-" 
                        + std::to_string(_dim - 1) + ").");
#endif
            return _ptr[ idx ];
        }

        // Dimensions of vector
        int dim() const { return _dim; }

        // Vector dot product
        double dot( const RealVector& B )
        {
#ifndef NDEBUG
            if ( B._dim != _dim )
                throw std::runtime_error("\nSize mismatch in vector dot product!\ndim(A) = " + std::to_string(_dim) + ", dim(B) = " 
                        + std::to_string(B._dim));
#endif
            return _scaling*B._scaling*cblas_ddot(_dim, _ptr, 1, B._ptr, 1);
        }

        double dot( RealVector&& B )
        {
#ifndef NDEBUG
            if ( B._dim != _dim )
                throw std::runtime_error("\nSize mismatch in vector dot product!\ndim(A) = " + std::to_string(_dim) + ", dim(B) = "
                        + std::to_string(B._dim));
#endif
            return _scaling*B._scaling*cblas_ddot(_dim, _ptr, 1, B._ptr, 1);
        }

        // Erase contents
        void erase()
        {
            if ( _ptr && _ownsItsPointer )
                delete[] _ptr;
                
            _ptr = nullptr;
            _dim = 0;
            _scaling = 1.;
            _isScaled = false;
            _ownsItsPointer = true;
        }

        // Allocate/reallocate space for vector, initializing all values to zero
        void init( int dim )
        {
#ifndef NDEBUG
            if ( dim < 1 )
                throw std::runtime_error("\nCannot initialize RealVector with dim = " + std::to_string(dim));
#endif
            if ( _ptr && _ownsItsPointer )
                delete[] _ptr;

            _dim = dim;
            _ptr = new double[ _dim ]();
            _isScaled = false;
            _scaling = 1.;
            _ownsItsPointer = true;
        }
        
        // Query scaling status
        bool isScaled() const { return _isScaled; }
        
        // Query pointer ownership
        bool ownsItsPointer() const { return _ownsItsPointer; }
        
        // Show contents of matrix in scientific precision
        void print( const char *s, int n ) const
        {
            std::printf("\nRealVector %s:\n\n", s);
            
            if ( _isScaled )
                std::printf("...is scaled with factor %e\n", _scaling);
            else
                std::printf("...is not scaled (scaling = %e)\n", _scaling);
            
            if ( _ownsItsPointer )
                std::printf("...owns its pointer\n");
            else
                std::printf("...does not own its pointer\n");

            if ( !_ptr )
                std::printf("...is empty\n");
            else
            {
                std::printf("\n");
                for ( int i = 0; i < _dim; i++ )
                    std::printf("%*.*e\n", n+10, n, _ptr[i]);
            }
            std::printf("\n");
        }

        // Print contents to file
        void printTo( FILE* fp, int n ) const
        {
            if ( !_ptr )
                std::fprintf(fp, "...is empty\n");
            else
                for ( int i = 0; i < _dim; i++ )
                    std::fprintf(fp, "%*.*e\n", n+10, n, _scaling*_ptr[i]);
        }
        
        // Give pointer to vector components
        double* ptr() const { return _ptr; }
        
        // Scale by real number
        RealVector scaleBy( double val ) const
        {
            RealVector A;
            
            A._dim = _dim;
            A._isScaled = true;
            A._scaling = _scaling*val;
            A._ptr = _ptr;
            A._ownsItsPointer = false;
            
            return A;
        }
        
        // Scaling factor
        double scaling() const { return _scaling; }
        
        // Resolve pointer ownership and scaling
        void simplify()
        {
            if ( !_ownsItsPointer )
            {
#ifdef VERBOSE_REALVECTOR_STATE_RESOLUTION
                std::printf("...Resolving pointer ownership of RealVector.\n");
                std::fflush(stdout);
#endif
                double* temp = _ptr;
                _ptr = new double[ _dim ];
                std::copy(temp, temp + _dim, _ptr);
                _ownsItsPointer = true;
            }
            
            if ( _isScaled )
            {
#ifdef VERBOSE_REALVECTOR_STATE_RESOLUTION
                std::printf("...Resolving scaling of RealVector.\n");
                std::fflush(stdout);
#endif
                cblas_dscal(_dim, _scaling, _ptr, 1);
                _isScaled = false;
                _scaling = 1.;
            }
        }
        
        // Vector tensor product
        RealMatrix xMen( const RealVector& B ) const
        {
            int dimB = B.dim();
            RealMatrix C( _dim, dimB );
            
            cblas_dger(CblasColMajor, _dim, dimB, _scaling*B._scaling, _ptr, 1, B._ptr, 1, C.ptr(), _dim);

            return C;
        }

        RealMatrix xMen( RealVector&& B ) const
        {
            int dimB = B.dim();
            RealMatrix C( _dim, dimB );
            
            cblas_dger(CblasColMajor, _dim, dimB, _scaling*B._scaling, _ptr, 1, B._ptr, 1, C.ptr(), _dim);

            return C;
        }
        
    private:
        int     _dim;
        double* _ptr;
        bool    _ownsItsPointer;
        bool    _isScaled;
        double  _scaling;
    };
}

#endif	/* REALVECTOR_HPP */