#ifndef STACKREALVECTOR_HPP
#define	STACKREALVECTOR_HPP

#include <config.h>

#include <cstdio>
#include <stdexcept>
#include <initializer_list>
#include "StackRealMatrix.hpp"
#include "RealVector.hpp"

#ifdef HAVE_MKL
    #include "mkl_cblas.h"
#else
    #include "cblas.h"
#endif

namespace broomstyx
{
    template<int nrows>
    class StackRealVector final
    {
    public:
        // Default constructor
        StackRealVector() {}

        // Constructor with initializer list
        StackRealVector( std::initializer_list<double> initList )
        {
        	int dim = (int)initList.size();
			if ( dim != nrows )
				throw std::runtime_error("Mismatch detected between StackRealVector dimension and initializer list!");
            
            std::copy(initList.begin(), initList.begin() + nrows, _data);
        }

        // Copy constructor
        StackRealVector( const StackRealVector& source )
        {
        	if ( source.dim() != nrows )
				throw std::runtime_error("Mismatch detected between StackRealVector and source dimensions!");

            std::copy(source._data, source._data + nrows, _data);
        }

        // Constructor from RealVector
        StackRealVector( const RealVector& source )
        {
            if ( source.dim() != nrows )
				throw std::runtime_error("Mismatch detected between StackRealVector and source dimensions!");

            std::copy(source.ptr(), source.ptr() + nrows, _data);
        }

        // Destructor
        virtual ~StackRealVector() {}

        // Assignment with initializer list
        StackRealVector& operator= ( std::initializer_list<double> initList )
        {
        	int dim = (int)initList.size();
        	if ( dim != nrows )
				throw std::runtime_error("Mismatch detected between StackRealVector dimensions and initializer list!");

        	std::copy(initList.begin(), initList.begin() + nrows, _data);
            
            return *this;
        }
        
        // Copy assignment operator
        StackRealVector& operator= ( const StackRealVector<nrows>& source )
        {
            // Check for self-assignment
            if ( this == &source )
                return *this;
            
			std::copy(source._data, source._data + nrows, _data);

            return *this;
        }

        // Assignment from RealVector
        StackRealVector& operator= ( const RealVector& source )
		{
			if ( source.dim() != nrows )
				throw std::runtime_error("Mismatch detected between StackRealVector and RealVector source dimensions!");

			std::copy(source.ptr(), source.ptr() + nrows, _data);

			return *this;
		}

        // Addition to self
        StackRealVector<nrows>& operator+= ( StackRealVector<nrows>& source )
        {
            cblas_daxpy(nrows, 1.0, source._data, 1, _data, 1);
            return *this;
        }
        
        // Subtraction from self
        StackRealVector<nrows>& operator-= ( StackRealVector<nrows>& source )
        {
            cblas_daxpy(nrows, -1.0, source._data, 1, _data, 1);
            return *this;
        }

        // In-place scalar multiplication
        StackRealVector<nrows>& operator*= ( double factor )
        {
            cblas_dscal(nrows, factor, _data, 1);
            return *this;         
        }

        // In-place scalar division
        StackRealVector<nrows>& operator/= ( double factor )
        {
            cblas_dscal(nrows, 1./factor, _data, 1);
            return *this;         
        }
        
        // Vector component address access as variableName(i)
        double& operator()( int idx )
        {
#ifndef NDEBUG
            if ( idx < 0 || idx >= nrows )
                throw std::runtime_error("\nCannot access StackRealVector component (" + std::to_string(idx) + ")! Valid range is (0-"
                        + std::to_string(nrows - 1) + ").");
#endif
            return _data[idx];
        }

        // Vector component value access as variableName.at(i)
        double at( int idx ) const
        {
#ifndef NDEBUG
            if ( idx < 0 || idx >= nrows )
                throw std::runtime_error("\nCannot access StackRealVector component (" + std::to_string(idx) + ")! Valid range is (0-"
                        + std::to_string(nrows - 1) + ").");
#endif
            return _data[idx];
        }

        // Vector cross product
        StackRealVector<3> cross( const StackRealVector<3>& B )
        {
            double c0 = _data[1]*B._data[2] - _data[2]*B._data[1];
            double c1 = _data[2]*B._data[0] - _data[0]*B._data[2];
            double c2 = _data[0]*B._data[1] - _data[1]*B._data[0];

            return StackRealVector<3>({c0, c1, c2});
        }

        // Dimensions of vector
        int dim() const { return nrows; }

        // Vector dot product
		double dot( StackRealVector<nrows> B )
        {
            return cblas_ddot(nrows, _data, 1, B._data, 1);
        }

        // Allocate/reallocate space for vector, initializing all values to zero
        void init()
        {
        	std::fill(_data, _data + nrows, 0.);
        }
        
        // Show contents of matrix in scientific precision
        void print( const char *s, int n ) const
        {
            std::printf("\nStackRealVector %s: size = %d\n\n", s, nrows);
            
			for ( int i = 0; i < nrows; i++ )
				std::printf("%*.*e\n", n+10, n, _data[i]);

            std::printf("\n");
        }

        // Print contents to file
        void printTo( FILE* fp, int n ) const
        {
            for ( int i = 0; i < nrows; i++ )
				std::fprintf(fp, "%*.*e\n", n+10, n, _data[i]);
        }
        
        // Give pointer to vector components
        double* ptr()
        {
        	return _data;
        }
        
        // Vector tensor product
        template<int ncols>
        StackRealMatrix<nrows,ncols> xMen( StackRealVector<ncols>& B ) const
        {
        	StackRealMatrix<nrows,ncols> C;
        	for ( int j = 0; j < ncols; j++ )
#pragma GCC ivdep
        		for ( int i = 0; i < nrows; i++ )
        			C(i,j) = _data[i]*B.at(j);

//			cblas_dger(CblasColMajor, nrows, ncols, 1.0, const_cast<double*>(_data), 1, const_cast<double*>(B.ptr()), 1, const_cast<double*>(C.ptr()), nrows);

            return C;
        }
        
    private:
        double _data[nrows] = {0.};
    };
}

#endif	/* STACKREALVECTOR_HPP */
