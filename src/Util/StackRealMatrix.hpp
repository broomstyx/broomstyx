#ifndef STACKREALMATRIX_HPP
#define	STACKREALMATRIX_HPP

#include <config.h>

#include <cstdio>
#include <algorithm>
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
    template< int nrows, int ncols>
    class StackRealMatrix final
    {
    public:
        // Default constructor
        StackRealMatrix() {}

        // Constructor with initializer list
        StackRealMatrix( std::initializer_list< std::initializer_list<double>> initList )
        {
            int dim1 = (int)initList.size();
            int dim2 = (int)(*(initList.begin())).size();

            if ( dim1 != nrows || dim2 != ncols )
            	throw std::runtime_error("Mismatch detected between StackRealMatrix dimensions and initializer list!");

#ifndef NDEBUG
            for ( int i = 1; i < dim1; i++ )
                if ( dim2 != (int)(*(initList.begin() + i)).size() )
                    throw std::runtime_error("RealMatrix initializer list has inconsistent dimensions!");
#endif

            for ( int i = 0; i < dim1; i++ )
            {
#pragma GCC ivdep
                for ( int j = 0; j < dim2; j++ )
                    _data[ j*dim1 + i ] = *((*(initList.begin() + i)).begin() + j);
            }
        }
        
        // Copy constructor
        StackRealMatrix( const StackRealMatrix<nrows,ncols>& source )
        {
            std::copy(source._data, source._data + nrows*ncols, _data);
        }

        // Constructor from RealVector
        StackRealMatrix( const RealMatrix& source )
        {
            int dim1 = source.dim1();
            int dim2 = source.dim2();

            if ( dim1 != nrows || dim2 != ncols )
            	throw std::runtime_error("Mismatch detected between StackRealMatrix and source dimensions!");

            std::copy(source.ptr(), source.ptr() + nrows*ncols, _data);
        }

        // Destructor
        virtual ~StackRealMatrix() {}

        // Assignment with initializer list
        StackRealMatrix& operator= ( std::initializer_list< std::initializer_list<double>> initList )
        {
            int dim1 = (int)initList.size();
            int dim2 = (int)(*(initList.begin())).size();

            if ( dim1 != nrows || dim2 != ncols )
				throw std::runtime_error("Mismatch detected between StackRealMatrix dimensions and initializer list!");

#ifndef NDEBUG
            for ( int i = 1; i < dim1; i++ )
                if ( dim2 != (int)(*(initList.begin() + i)).size() )
                    throw std::runtime_error("RealMatrix initializer list has inconsistent dimensions!");
#endif
            
            for ( int i = 0; i < nrows; i++ )
            {
#pragma GCC ivdep
                for ( int j = 0; j < ncols; j++ )
                    _data[ j*nrows + i ] = *((*(initList.begin() + i)).begin() + j);
            }
            
            return *this;
        }
        
        // Copy assignment operator
        StackRealMatrix& operator= ( const StackRealMatrix<nrows,ncols>& source )
        {
            // Check for self-assignment
            if ( this == &source )
                return *this;
            
            std::copy(source._data, source._data + nrows*ncols, _data);
            
            return *this;
        }
        
        // Assignment from RealMatrix
        StackRealMatrix& operator= ( const RealMatrix& source )
        {
			if ( source.dim1() != nrows || source.dim2() != ncols )
				throw std::runtime_error("Mismatch detected between StackRealMatrix and RealMatrix source dimensions!");

			std::copy(source.ptr(), source.ptr() + nrows*ncols, _data);

			return *this;
        }

        // Addition to self
        StackRealMatrix& operator+= ( const StackRealMatrix<nrows,ncols>& source )
        {
            cblas_daxpy(nrows*ncols, 1.0, source._data, 1, _data, 1);

            return *this;
        }
        
        // Subtraction from self
        StackRealMatrix& operator-= ( const StackRealMatrix<nrows,ncols>& source )
        {
            cblas_daxpy(nrows*ncols, -1.0, source._data, 1, _data, 1);

            return *this;
        }
        
        // NOTE: Storage of matrix components uses column-major format
        // Matrix component address access as variableName(i,j)
        double& operator()( int idx1, int idx2 )
        {
#ifndef NDEBUG
        	if ( idx1 < 0 || idx1 >= nrows || idx2 < 0 || idx2 >= ncols )
                throw std::runtime_error("\nCannot access RealMatrix component (" + std::to_string(idx1) + "," + std::to_string(idx2) + ")! "
                        + "Valid range is (0-" + std::to_string(nrows - 1) + ",0-" + std::to_string(ncols - 1) + ").");
#endif
            return _data[ idx2*nrows + idx1 ];
        }

        // Matrix component value access as variableName.at(i,j)
        double at( int idx1, int idx2 ) const
        {
#ifndef NDEBUG
        	if ( idx1 < 0 || idx1 >= nrows || idx2 < 0 || idx2 >= ncols )
                throw std::runtime_error("\nCannot access RealMatrix component (" + std::to_string(idx1) + "," + std::to_string(idx2) + ")! "
                        + "Valid range is (0-" + std::to_string(nrows - 1) + ",0-" + std::to_string(ncols - 1) + ").");
#endif
            return _data[ idx2*nrows + idx1 ];
        }

        // Dimensions of matrix
        int dim1() const { return nrows; }
        int dim2() const { return ncols; }

        // Allocate/reallocate space for matrix of size (dim1 x dim2), with all
        // values set to zero
        void init()
        {
        	std::fill(_data, _data + nrows*ncols, 0.);
        }

        // Show contents of matrix in scientific precision
        void print( const char *s, int n ) const 
        {
            std::printf("\nStackRealMatrix %s: size = %d x %d\n\n", s, nrows, ncols);
            
			for ( int i = 0; i < nrows; i++ )
			{
				for ( int j = 0; j < ncols; j++ )
					std::printf("%*.*e", n+10, n, _data[ j*nrows + i ]);
				std::printf("\n");
			}
            std::printf("\n");
        }

        // Give pointer to matrix components
        double* ptr()
        {
        	return _data;
        }
        
    private:
        double _data[nrows*ncols] = {0.};
    };
}

#endif	/* STACKREALMATRIX_HPP */
