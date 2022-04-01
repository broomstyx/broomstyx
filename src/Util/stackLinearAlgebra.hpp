#include <cstdio>
#include <stdexcept>

#include <config.h>
#include "StackRealVector.hpp"
#include "StackRealMatrix.hpp"
#include "RealVector.hpp"
#include "RealMatrix.hpp"

#ifdef HAVE_MKL
    #include "mkl_cblas.h"
    #include "mkl_lapacke.h"
#else
    #include "cblas.h"
    #include "lapacke.h"
#endif

namespace broomstyx
{
    // Matrix addition
    template<int nrows, int ncols>
    StackRealMatrix<nrows, ncols> operator+( StackRealMatrix<nrows, ncols> A, StackRealMatrix<nrows, ncols> B )
    {
        cblas_daxpy(nrows*ncols, 1.0, B.ptr(), 1, A.ptr(), 1);
        
        return A;
    }

    // Vector addition
    template<int nrows>
	StackRealVector<nrows> operator+( StackRealVector<nrows> A, StackRealVector<nrows> B )
    {
        cblas_daxpy(nrows, 1.0, B.ptr(), 1, A.ptr(), 1);

        return A;
    }

    // Matrix subtraction
    template<int nrows, int ncols>
    StackRealMatrix<nrows,ncols> operator-( StackRealMatrix<nrows,ncols> A, StackRealMatrix<nrows,ncols> B )
    {
    	cblas_daxpy(nrows*ncols, -1.0, B.ptr(), 1, A.ptr(), 1);

        return A;
    }

    // Vector subtraction
    template<int nrows>
    StackRealVector<nrows> operator-( StackRealVector<nrows> A, StackRealVector<nrows> B )
    {
        cblas_daxpy(nrows, -1.0, B.ptr(), 1, A.ptr(), 1);

        return A;
    }

    // Scalar multiplication for matrices
    template<int nrows, int ncols>
    StackRealMatrix<nrows,ncols> operator*( StackRealMatrix<nrows,ncols> A, double b )
    {
    	cblas_dscal(nrows*ncols, b, A.ptr(), 1);

        return A;
    }

    template<int nrows, int ncols>
	StackRealMatrix<nrows,ncols> operator*( double a, StackRealMatrix<nrows,ncols> B )
	{
		cblas_dscal(nrows*ncols, a, B.ptr(), 1);

		return B;
	}

    // Scalar division for matrices
    template<int nrows, int ncols>
	StackRealMatrix<nrows,ncols> operator/( StackRealMatrix<nrows,ncols> A, double b )
	{
		cblas_dscal(nrows*ncols, 1./b, A.ptr(), 1);

		return A;
	}

    // Scalar multiplication for vectors
    template<int nrows>
	StackRealVector<nrows> operator*( StackRealVector<nrows> A, double b )
	{
		cblas_dscal(nrows, b, A.ptr(), 1);

		return A;
	}

	template<int nrows>
	StackRealVector<nrows> operator*( double a, StackRealVector<nrows> B )
	{
		cblas_dscal(nrows, a, B.ptr(), 1);

		return B;
	}

    // Scalar division for vectors
	template<int nrows>
	StackRealVector<nrows> operator/( StackRealVector<nrows> A, double b )
	{
		cblas_dscal(nrows, 1./b, A.ptr(), 1);

		return A;
	}

    // Matrix multiplication
	template<int nrows, int innerdim, int ncols>
    StackRealMatrix<nrows,ncols> operator*( StackRealMatrix<nrows,innerdim> A, StackRealMatrix<innerdim,ncols> B )
    {
        StackRealMatrix<nrows,ncols> C;
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, nrows, ncols, innerdim, 1.0, A.ptr(), nrows, B.ptr(), innerdim, 0.0, C.ptr(), nrows);

        return C;
    }

    // Matrix-vector multiplication
	template<int nrows, int ncols>
    StackRealVector<nrows> operator*( StackRealMatrix<nrows,ncols> A, StackRealVector<ncols> B)
    {
        StackRealVector<nrows> C;
        cblas_dgemv(CblasColMajor, CblasNoTrans, nrows, ncols, 1.0, A.ptr(), nrows, B.ptr(), 1, 1.0, C.ptr(), 1);

        return C;
    }

    // Vector-matrix multiplication
	template<int nrows, int ncols>
    StackRealVector<ncols> operator*( StackRealVector<nrows> A, StackRealMatrix<nrows, ncols> B )
    {   
        StackRealVector<ncols> C;
        cblas_dgemv(CblasColMajor, CblasTrans, nrows, ncols, 1.0, B.ptr(), nrows, A.ptr(), 1, 1.0, C.ptr(), 1);

        return C;
    }

    // Matrix transpose
	template<int nrows, int ncols>
    StackRealMatrix<ncols,nrows> trp( StackRealMatrix<nrows,ncols> A )
    {
		StackRealMatrix<ncols,nrows> B;

        for ( int i = 0; i < nrows; i++ )
        {
#pragma GCC ivdep
        	for ( int j = 0; j < ncols; j++ )
        		B(j,i) = A.at(i,j);
        }
        
        return B;
    }

    // Matrix inverse
    template<int nrows>
    StackRealMatrix<nrows,nrows> inv( StackRealMatrix<nrows,nrows> B )
    {
        int info;
        lapack_int* ipiv;

        ipiv = new lapack_int[nrows]();
        double *locB = B.ptr();

        //dgetrf(&nrows, &nrows, locA, &nrows, ipiv, &info);
        info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, nrows, nrows, locB, nrows, ipiv);
        if ( info < 0 )
        {
            delete[] ipiv;
            throw std::runtime_error("Error in matrix inversion: illegal value encountered during factorization!");
        }
        else if (info > 0)
        {
            delete[] ipiv;
            throw std::runtime_error("Cannot invert singular matrix!");
        }

        //dgetri(&nrows, locA, &nrows, ipiv, work, &lwork, &info);
        info = LAPACKE_dgetri(LAPACK_COL_MAJOR, nrows, locB, nrows, ipiv);

        delete[] ipiv;

        if ( info < 0 )
            throw std::runtime_error("Error in matrix inversion: illegal value encountered during factorization!");
        else if (info > 0)
            throw std::runtime_error("Cannot invert singular matrix!");

        return B;
    }
}
