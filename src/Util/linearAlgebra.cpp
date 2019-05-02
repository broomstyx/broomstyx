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

#include "linearAlgebra.hpp"
#include <cstdio>
#include <stdexcept>
#include "../../config.h"

#ifdef HAVE_MKL
    #include "mkl_cblas.h"
    #include "mkl_lapacke.h"
#else
    #include "cblas.h"

    typedef int lapack_int;
//    #include "lapacke.h"

extern "C" {
// LU decomoposition of a general matrix
void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
//// generate inverse of a matrix given its LU decomposition
void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int*    lwork, int* INFO);
void dgetrs_(char* C, int* N, int* NRHS, double* A, int* LDA, int* IPIV, double* B, int* LDB, int* INFO);
void dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);
}

#endif

namespace broomstyx
{
    // Matrix addition
    // ---------------------------------------------------------------------------
    RealMatrix operator+( RealMatrix&& A, RealMatrix&& B )
    {
        if ( A.isTransposed() )
            A.simplify();
        if ( B.isTransposed() )
            B.simplify();

        int A_dim1 = A.dim1();
        int A_dim2 = A.dim2();
        int B_dim1 = B.dim1();
        int B_dim2 = B.dim2();

#ifndef NDEBUG
        if ( A_dim1 != B_dim1 || A_dim2 != B_dim2 )
            throw std::runtime_error("\nSize mismatch in operands for matrix addition!\n\tdim(A) = [ " + std::to_string(A_dim1) + " x "
                    + std::to_string(A_dim2) + " ], dim(B) = [ " + std::to_string(B_dim1) + " x " + std::to_string(B_dim2) + " ]");
#endif
        if ( A.isScaled() || !A.ownsItsPointer() )
        {
            B.simplify();
            cblas_daxpy(A_dim1*A_dim2, A.scaling(), A.ptr(), 1, B.ptr(), 1);
            return (RealMatrix&&)B;
        }
        else
        {
            A.simplify();
            cblas_daxpy(A_dim1*A_dim2, B.scaling(), B.ptr(), 1, A.ptr(), 1);
            return (RealMatrix&&)A;
        }
    }

    RealMatrix operator+( RealMatrix&& A, const RealMatrix& B )
    {
        A.simplify();
#ifndef NDEBUG
        if ( B.isTransposed() )
            throw std::runtime_error("\nUnresolved transposition of lvalue RealMatrix B detected during operation A + B");
#endif
        int A_dim1 = A.dim1();
        int A_dim2 = A.dim2();
        int B_dim1 = B.dim1();
        int B_dim2 = B.dim2();
#ifndef NDEBUG
        if ( A_dim1 != B_dim1 || A_dim2 != B_dim2 )
            throw std::runtime_error("\nSize mismatch in operands for matrix addition!\n\tdim(A) = [ " + std::to_string(A_dim1) + " x "
                    + std::to_string(A_dim2) + " ], dim(B) = [ " + std::to_string(B_dim1) + " x " + std::to_string(B_dim2) + " ]");
#endif
        cblas_daxpy(A_dim1*A_dim2, B.scaling(), B.ptr(), 1, A.ptr(), 1);
        return (RealMatrix&&)A;
    }

    RealMatrix operator+( const RealMatrix& A, RealMatrix&& B )
    {
        B.simplify();
#ifndef NDEBUG
        if ( A.isTransposed() )
            throw std::runtime_error("\nUnresolved transposition of lvalue RealMatrix A detected during operation A + B");
#endif
        int A_dim1 = A.dim1();
        int A_dim2 = A.dim2();
        int B_dim1 = B.dim1();
        int B_dim2 = B.dim2();
#ifndef NDEBUG
        if ( A_dim1 != B_dim1 || A_dim2 != B_dim2 )
            throw std::runtime_error("\nSize mismatch in operands for matrix addition!\n\tdim(A) = [ " + std::to_string(A_dim1) + " x "
                    + std::to_string(A_dim2) + " ], dim(B) = [ " + std::to_string(B_dim1) + " x " + std::to_string(B_dim2) + " ]");
#endif
        cblas_daxpy(B_dim1*B_dim2, A.scaling(), A.ptr(), 1, B.ptr(), 1);

        return (RealMatrix&&)B;
    }

    RealMatrix operator+( const RealMatrix& A, const RealMatrix& B )
    {
#ifndef NDEBUG
        if ( A.isTransposed() )
            throw std::runtime_error("\nUnresolved transposition of lvalue RealMatrix A detected during operation A + B");
        if ( B.isTransposed() )
            throw std::runtime_error("\nUnresolved transposition of lvalue RealMatrix B detected during operation A + B");
#endif
        int A_dim1 = A.dim1();
        int A_dim2 = A.dim2();
        int B_dim1 = B.dim1();
        int B_dim2 = B.dim2();
#ifndef NDEBUG
        if ( A_dim1 != B_dim1 || A_dim2 != B_dim2 )
            throw std::runtime_error("\nSize mismatch in operands for matrix addition!\n\tdim(A) = [ " + std::to_string(A_dim1) + " x "
                    + std::to_string(A_dim2) + " ], dim(B) = [ " + std::to_string(B_dim1) + " x " + std::to_string(B_dim2) + " ]");
#endif
        RealMatrix C = A;
        C.simplify();
        cblas_daxpy(A_dim1*A_dim2, B.scaling(), B.ptr(), 1, C.ptr(), 1);

        return C;
    }

    // Vector addition
    // ---------------------------------------------------------------------------
    RealVector operator+( RealVector&& A, RealVector&& B )
    {
        int A_dim = A.dim();
        int B_dim = B.dim();
#ifndef NDEBUG
        if ( A_dim != B_dim )
            throw std::runtime_error("\nSize mismatch in operands for vector addition!\n\tdim(A) = " + std::to_string(A_dim) + ", dim(B) = " + std::to_string(B_dim));
#endif
        if ( A.isScaled() || !A.ownsItsPointer() )
        {
            B.simplify();
            cblas_daxpy(A_dim, A.scaling(), A.ptr(), 1, B.ptr(), 1);
            return (RealVector&&)B;
        }
        else
        {
            A.simplify();
            cblas_daxpy(A_dim, B.scaling(), B.ptr(), 1, A.ptr(), 1);
            return (RealVector&&)A;
        }
    }

    RealVector operator+( RealVector&& A, const RealVector& B )
    {
        int A_dim = A.dim();
        int B_dim = B.dim();
#ifndef NDEBUG
        if ( A_dim != B_dim )
            throw std::runtime_error("\nSize mismatch in operands for vector addition!\n\tdim(A) = " + std::to_string(A_dim) + ", dim(B) = " + std::to_string(B_dim));
#endif
        A.simplify();
        cblas_daxpy(A_dim, B.scaling(), B.ptr(), 1, A.ptr(), 1);

        return (RealVector&&)A;
    }

    RealVector operator+( const RealVector& A, RealVector&& B )
    {
        int A_dim = A.dim();
        int B_dim = B.dim();
#ifndef NDEBUG
        if ( A_dim != B_dim )
            throw std::runtime_error("\nSize mismatch in operands for vector addition!\n\tdim(A) = " + std::to_string(A_dim) + ", dim(B) = " + std::to_string(B_dim));
#endif
        B.simplify();
        cblas_daxpy(B_dim, A.scaling(), A.ptr(), 1, B.ptr(), 1);

        return (RealVector&&)B;
    }

    RealVector operator+( const RealVector& A, const RealVector& B )
    {
        int A_dim = A.dim();
        int B_dim = B.dim();
#ifndef NDEBUG
        if ( A_dim != B_dim )
            throw std::runtime_error("\nSize mismatch in operands for vector addition!\n\tdim(A) = " + std::to_string(A_dim) + ", dim(B) = " + std::to_string(B_dim));
#endif
        RealVector C = A;
        C.simplify();
        cblas_daxpy(A_dim, B.scaling(), B.ptr(), 1, C.ptr(), 1);

        return C;
    }

    // Matrix subtraction
    // ---------------------------------------------------------------------------
    RealMatrix operator-( RealMatrix&& A, RealMatrix&& B )
    {
        A.simplify();
        if ( B.isTransposed() )
            B.simplify();

        int A_dim1 = A.dim1();
        int A_dim2 = A.dim2();
        int B_dim1 = B.dim1();
        int B_dim2 = B.dim2();
#ifndef NDEBUG
        if ( A_dim1 != B_dim1 || A_dim2 != B_dim2 )
            throw std::runtime_error("\nSize mismatch in operands for matrix subtraction!\n\tdim(A) = [ " + std::to_string(A_dim1) + " x "
                    + std::to_string(A_dim2) + " ], dim(B) = [ " + std::to_string(B_dim1) + " x " + std::to_string(B_dim2) + " ]");
#endif
        cblas_daxpy(A_dim1*A_dim2, -1.*B.scaling(), B.ptr(), 1, A.ptr(), 1);

        return (RealMatrix&&)A;
    }

    RealMatrix operator-( RealMatrix&& A, const RealMatrix& B )
    {
        A.simplify();
#ifndef NDEBUG
        if ( B.isTransposed() )
            throw std::runtime_error("\nUnresolved transposition of lvalue RealMatrix B detected during operation A - B");
#endif
        int A_dim1 = A.dim1();
        int A_dim2 = A.dim2();
        int B_dim1 = B.dim1();
        int B_dim2 = B.dim2();
#ifndef NDEBUG
        if ( A_dim1 != B_dim1 || A_dim2 != B_dim2 )
            throw std::runtime_error("\nSize mismatch in operands for matrix subtraction!\n\tdim(A) = [ " + std::to_string(A_dim1) + " x "
                    + std::to_string(A_dim2) + " ], dim(B) = [ " + std::to_string(B_dim1) + " x " + std::to_string(B_dim2) + " ]");
#endif
        cblas_daxpy(A_dim1*A_dim2, -1.*B.scaling(), B.ptr(), 1, A.ptr(), 1);

        return (RealMatrix&&)A;
    }

    RealMatrix operator-( const RealMatrix& A, RealMatrix&& B )
    {
        if ( B.isTransposed() )
            B.simplify();

        RealMatrix C = A;
        C.simplify();

        int A_dim1 = C.dim1();
        int A_dim2 = C.dim2();
        int B_dim1 = B.dim1();
        int B_dim2 = B.dim2();
#ifndef NDEBUG
        if ( A_dim1 != B_dim1 || A_dim2 != B_dim2 )
            throw std::runtime_error("\nSize mismatch in operands for matrix subtraction!\n\tdim(A) = [ " + std::to_string(A_dim1) + " x "
                    + std::to_string(A_dim2) + " ], dim(B) = [ " + std::to_string(B_dim1) + " x " + std::to_string(B_dim2) + " ]");
#endif
        cblas_daxpy(A_dim1*A_dim2, -1.*B.scaling(), B.ptr(), 1, C.ptr(), 1);

        return C;
    }

    RealMatrix operator-( const RealMatrix& A, const RealMatrix& B )
    {
#ifndef NDEBUG
        if ( B.isTransposed() )
            throw std::runtime_error("\nUnresolved transposition of lvalue RealMatrix B detected during operation A - B");
#endif
        RealMatrix C = A;
        C.simplify();

        int A_dim1 = C.dim1();
        int A_dim2 = C.dim2();
        int B_dim1 = B.dim1();
        int B_dim2 = B.dim2();
#ifndef NDEBUG
        if ( A_dim1 != B_dim1 || A_dim2 != B_dim2 )
            throw std::runtime_error("\nSize mismatch in operands for matrix subtraction!\n\tdim(A) = [ " + std::to_string(A_dim1) + " x "
                    + std::to_string(A_dim2) + " ], dim(B) = [ " + std::to_string(B_dim1) + " x " + std::to_string(B_dim2) + " ]");
#endif

        cblas_daxpy(A_dim1*A_dim2, -1.*B.scaling(), B.ptr(), 1, C.ptr(), 1);

        return C;
    }

    // Vector subtraction
    // ---------------------------------------------------------------------------
    RealVector operator-( RealVector&& A, RealVector&& B )
    {
        A.simplify();

        int A_dim = A.dim();
        int B_dim = B.dim();
#ifndef NDEBUG
        if ( A_dim != B_dim )
            throw std::runtime_error("\nSize mismatch in operands for vector subtraction!\n\tdim(A) = " + std::to_string(A_dim) + ", dim(B) = " + std::to_string(B_dim));
#endif
        cblas_daxpy(A_dim, -1.*B.scaling(), B.ptr(), 1, A.ptr(), 1);

        return (RealVector&&)A;
    }

    RealVector operator-( RealVector&& A, const RealVector& B )
    {
        A.simplify();

        int A_dim = A.dim();
        int B_dim = B.dim();
#ifndef NDEBUG
        if ( A_dim != B_dim )
            throw std::runtime_error("\nSize mismatch in operands for vector subtraction!\n\tdim(A) = " + std::to_string(A_dim) + ", dim(B) = " + std::to_string(B_dim));
#endif
        cblas_daxpy(A_dim, -1.*B.scaling(), B.ptr(), 1, A.ptr(), 1);

        return (RealVector&&)A;
    }

    RealVector operator-( const RealVector& A, RealVector&& B )
    {
        int A_dim = A.dim();
        int B_dim = B.dim();
#ifndef NDEBUG
        if ( A_dim != B_dim )
            throw std::runtime_error("\nSize mismatch in operands for vector subtraction!\n\tdim(A) = " + std::to_string(A_dim) + ", dim(B) = " + std::to_string(B_dim));
#endif
        RealVector C = A;
        C.simplify();
        cblas_daxpy(A_dim, -1.*B.scaling(), B.ptr(), 1, C.ptr(), 1);

        return C;
    }

    RealVector operator-( const RealVector& A, const RealVector& B )
    {
        int A_dim = A.dim();
        int B_dim = B.dim();
#ifndef NDEBUG
        if ( A_dim != B_dim )
            throw std::runtime_error("\nSize mismatch in operands for vector subtraction!\n\tdim(A) = " + std::to_string(A_dim) + ", dim(B) = " + std::to_string(B_dim));
#endif
        RealVector C = A;
        C.simplify();
        cblas_daxpy(A_dim, -1.*B.scaling(), B.ptr(), 1, C.ptr(), 1);

        return C;
    }

    // Scalar multiplication for matrices
    // ---------------------------------------------------------------------------
    RealMatrix operator*( RealMatrix&& A, double b )
    {
        A._isScaled = true;
        A._scaling *= b;
        return (RealMatrix&&)A;
    }

    RealMatrix operator*( const RealMatrix& A, double b )
    {
        return A.scaleBy(b);
    }

    RealMatrix operator*( double a, RealMatrix&& B )
    {
        B._isScaled = true;
        B._scaling *= a;
        return (RealMatrix&&)B;
    }

    RealMatrix operator*( double a, const RealMatrix& B )
    {
        return B.scaleBy(a);
    }

    // Scalar division for matrices
    // ---------------------------------------------------------------------------
    RealMatrix operator/( RealMatrix&& A, double b )
    {
        A._isScaled = true;
        A._scaling /= b;
        return (RealMatrix&&)A;
    }

    RealMatrix operator/( const RealMatrix& A, double b )
    {
        return A.scaleBy(1./b);
    }

    // Scalar multiplication for vectors
    // ---------------------------------------------------------------------------
    RealVector operator*( RealVector&& A, double b )
    {
        A._isScaled = true;
        A._scaling *= b;
        return (RealVector&&)A;
    }

    RealVector operator*( const RealVector& A, double b )
    {
        return A.scaleBy(b);
    }

    RealVector operator*( double a, RealVector&& B )
    {
        B._isScaled = true;
        B._scaling *= a;
        return (RealVector&&)B;
    }

    RealVector operator*( double a, const RealVector& B )
    {
        return B.scaleBy(a);
    }

    // Scalar division for vectors
    // ---------------------------------------------------------------------------
    RealVector operator/( RealVector&& A, double b )
    {
        A._isScaled = true;
        A._scaling /= b;
        return (RealVector&&)A;
    }

    RealVector operator/( const RealVector& A, double b )
    {
        return A.scaleBy(1./b);
    }

    // Matrix multiplication
    // ---------------------------------------------------------------------------
    RealMatrix operator*( RealMatrix&& A, RealMatrix&& B )
    {
        CBLAS_TRANSPOSE transA, transB;
        int opA_dim1, opA_dim2, opB_dim1, opB_dim2;

        int ldA = A.dim1();
        int ldB = B.dim1();

        if ( A.isTransposed() )
        {
            opA_dim1 = A.dim2();
            opA_dim2 = A.dim1();
            transA = CblasTrans;
        }
        else
        {
            opA_dim1 = A.dim1();
            opA_dim2 = A.dim2();
            transA = CblasNoTrans;
        }

        if ( B.isTransposed() )
        {
            opB_dim1 = B.dim2();
            opB_dim2 = B.dim1();
            transB = CblasTrans;
        }
        else
        {
            opB_dim1 = B.dim1();
            opB_dim2 = B.dim2();
            transB = CblasNoTrans;
        }
#ifndef NDEBUG
        if ( opA_dim2 != opB_dim1 )
            throw std::runtime_error("\nSize mismatch in operands for matrix multiplication!\n\tdim(opA) = [ " + std::to_string(opA_dim1) + " x "
                    + std::to_string(opA_dim2) + " ], dim(opB) = [ " + std::to_string(opB_dim1) + " x " + std::to_string(opB_dim2) + " ]");
#endif
        RealMatrix C(opA_dim1, opB_dim2);
        cblas_dgemm(CblasColMajor, transA, transB, opA_dim1, opB_dim2, opA_dim2, A.scaling()*B.scaling(), A.ptr(), ldA, B.ptr(), ldB, 0.0, C.ptr(), opA_dim1);

        return C;
    }

    RealMatrix operator*( RealMatrix&& A, const RealMatrix& B )
    {
        CBLAS_TRANSPOSE transA, transB;
        int opA_dim1, opA_dim2, opB_dim1, opB_dim2;

        int ldA = A.dim1();
        int ldB = B.dim1();

        if ( A.isTransposed() )
        {
            opA_dim1 = A.dim2();
            opA_dim2 = A.dim1();
            transA = CblasTrans;
        }
        else
        {
            opA_dim1 = A.dim1();
            opA_dim2 = A.dim2();
            transA = CblasNoTrans;
        }

        if ( B.isTransposed() )
        {
            opB_dim1 = B.dim2();
            opB_dim2 = B.dim1();
            transB = CblasTrans;
        }
        else
        {
            opB_dim1 = B.dim1();
            opB_dim2 = B.dim2();
            transB = CblasNoTrans;
        }
#ifndef NDEBUG
        if ( opA_dim2 != opB_dim1 )
            throw std::runtime_error("\nSize mismatch in operands for matrix multiplication!\n\tdim(opA) = [ " + std::to_string(opA_dim1) + " x "
                    + std::to_string(opA_dim2) + " ], dim(opB) = [ " + std::to_string(opB_dim1) + " x " + std::to_string(opB_dim2) + " ]");
#endif
        RealMatrix C(opA_dim1, opB_dim2);
        cblas_dgemm(CblasColMajor, transA, transB, opA_dim1, opB_dim2, opA_dim2, A.scaling()*B.scaling(), A.ptr(), ldA, B.ptr(), ldB, 0.0, C.ptr(), opA_dim1);

        return C;
    }

    RealMatrix operator*( const RealMatrix& A, RealMatrix&& B )
    {
        CBLAS_TRANSPOSE transA, transB;
        int opA_dim1, opA_dim2, opB_dim1, opB_dim2;

        int ldA = A.dim1();
        int ldB = B.dim1();

        if ( A.isTransposed() )
        {
            opA_dim1 = A.dim2();
            opA_dim2 = A.dim1();
            transA = CblasTrans;
        }
        else
        {
            opA_dim1 = A.dim1();
            opA_dim2 = A.dim2();
            transA = CblasNoTrans;
        }

        if ( B.isTransposed() )
        {
            opB_dim1 = B.dim2();
            opB_dim2 = B.dim1();
            transB = CblasTrans;
        }
        else
        {
            opB_dim1 = B.dim1();
            opB_dim2 = B.dim2();
            transB = CblasNoTrans;
        }
#ifndef NDEBUG
        if ( opA_dim2 != opB_dim1 )
            throw std::runtime_error("\nSize mismatch in operands for matrix multiplication!\n\tdim(opA) = [ " + std::to_string(opA_dim1) + " x "
                    + std::to_string(opA_dim2) + " ], dim(opB) = [ " + std::to_string(opB_dim1) + " x " + std::to_string(opB_dim2) + " ]");
#endif
        RealMatrix C(opA_dim1, opB_dim2);
        cblas_dgemm(CblasColMajor, transA, transB, opA_dim1, opB_dim2, opA_dim2, A.scaling()*B.scaling(), A.ptr(), ldA, B.ptr(), ldB, 0.0, C.ptr(), opA_dim1);

        return C;
    }

    RealMatrix operator*( const RealMatrix& A, const RealMatrix& B )
    {
        CBLAS_TRANSPOSE transA, transB;
        int opA_dim1, opA_dim2, opB_dim1, opB_dim2;

        int ldA = A.dim1();
        int ldB = B.dim1();

        if ( A.isTransposed() )
        {
            opA_dim1 = A.dim2();
            opA_dim2 = A.dim1();
            transA = CblasTrans;
        }
        else
        {
            opA_dim1 = A.dim1();
            opA_dim2 = A.dim2();
            transA = CblasNoTrans;
        }

        if ( B.isTransposed() )
        {
            opB_dim1 = B.dim2();
            opB_dim2 = B.dim1();
            transB = CblasTrans;
        }
        else
        {
            opB_dim1 = B.dim1();
            opB_dim2 = B.dim2();
            transB = CblasNoTrans;
        }
#ifndef NDEBUG
        if ( opA_dim2 != opB_dim1 )
            throw std::runtime_error("\nSize mismatch in operands for matrix multiplication!\n\tdim(opA) = [ " + std::to_string(opA_dim1) + " x "
                    + std::to_string(opA_dim2) + " ], dim(opB) = [ " + std::to_string(opB_dim1) + " x " + std::to_string(opB_dim2) + " ]");
#endif
        RealMatrix C(opA_dim1, opB_dim2);
        cblas_dgemm(CblasColMajor, transA, transB, opA_dim1, opB_dim2, opA_dim2, A.scaling()*B.scaling(), A.ptr(), ldA, B.ptr(), ldB, 0.0, C.ptr(), opA_dim1);

        return C;
    }

    // Matrix-vector multiplication
    // ---------------------------------------------------------------------------
    RealVector operator*( RealMatrix&& A, RealVector&& B)
    {
        CBLAS_TRANSPOSE transA;
        int A_dim1, A_dim2, opA_dim1, opA_dim2, B_dim;

        A_dim1 = A.dim1();
        A_dim2 = A.dim2();

        if ( A.isTransposed() )
        {
            opA_dim1 = A_dim2;
            opA_dim2 = A_dim1;
            transA = CblasTrans;
        }
        else
        {
            opA_dim1 = A_dim1;
            opA_dim2 = A_dim2;
            transA = CblasNoTrans;
        }

        int ldA = A_dim1;
        B_dim = B.dim();
#ifndef NDEBUG
        if ( opA_dim2 != B_dim )
            throw std::runtime_error("\nSize mismatch in operands for matrix-vector multiplication!\n\tdim(A) = [ "
                    + std::to_string(opA_dim1) + " x " + std::to_string(opA_dim2) + " ], dim(B) =" + std::to_string(B_dim));
#endif
        RealVector C(opA_dim1);
        cblas_dgemv(CblasColMajor, transA, A_dim1, A_dim2, A.scaling()*B.scaling(), A.ptr(), ldA, B.ptr(), 1, 1.0, C.ptr(), 1);

        return C;
    }

    RealVector operator*( RealMatrix&& A, const RealVector& B)
    {
        CBLAS_TRANSPOSE transA;
        int A_dim1, A_dim2, opA_dim1, opA_dim2, B_dim;

        A_dim1 = A.dim1();
        A_dim2 = A.dim2();

        if ( A.isTransposed() )
        {
            opA_dim1 = A_dim2;
            opA_dim2 = A_dim1;
            transA = CblasTrans;
        }
        else
        {
            opA_dim1 = A_dim1;
            opA_dim2 = A_dim2;
            transA = CblasNoTrans;
        }

        int ldA = A_dim1;
        B_dim = B.dim();
#ifndef NDEBUG
        if ( opA_dim2 != B_dim )
            throw std::runtime_error("\nSize mismatch in operands for matrix-vector multiplication!\n\tdim(A) = [ "
                    + std::to_string(opA_dim1) + " x " + std::to_string(opA_dim2) + " ], dim(B) =" + std::to_string(B_dim));
#endif
        RealVector C(opA_dim1);
        cblas_dgemv(CblasColMajor, transA, A_dim1, A_dim2, A.scaling()*B.scaling(), A.ptr(), ldA, B.ptr(), 1, 1.0, C.ptr(), 1);

        return C;
    }

    RealVector operator*( const RealMatrix& A, RealVector&& B)
    {
        CBLAS_TRANSPOSE transA;
        int A_dim1, A_dim2, opA_dim1, opA_dim2, B_dim;

        A_dim1 = A.dim1();
        A_dim2 = A.dim2();

        if ( A.isTransposed() )
        {
            opA_dim1 = A_dim2;
            opA_dim2 = A_dim1;
            transA = CblasTrans;
        }
        else
        {
            opA_dim1 = A_dim1;
            opA_dim2 = A_dim2;
            transA = CblasNoTrans;
        }

        int ldA = A_dim1;
        B_dim = B.dim();
#ifndef NDEBUG
        if ( opA_dim2 != B_dim )
            throw std::runtime_error("\nSize mismatch in operands for matrix-vector multiplication!\n\tdim(A) = [ "
                    + std::to_string(opA_dim1) + " x " + std::to_string(opA_dim2) + " ], dim(B) =" + std::to_string(B_dim));
#endif
        RealVector C(opA_dim1);
        cblas_dgemv(CblasColMajor, transA, A_dim1, A_dim2, A.scaling()*B.scaling(), A.ptr(), ldA, B.ptr(), 1, 1.0, C.ptr(), 1);

        return C;
    }

    RealVector operator*( const RealMatrix& A, const RealVector& B)
    {
        CBLAS_TRANSPOSE transA;
        int A_dim1, A_dim2, opA_dim1, opA_dim2, B_dim;

        A_dim1 = A.dim1();
        A_dim2 = A.dim2();

        if ( A.isTransposed() )
        {
            opA_dim1 = A_dim2;
            opA_dim2 = A_dim1;
            transA = CblasTrans;
        }
        else
        {
            opA_dim1 = A_dim1;
            opA_dim2 = A_dim2;
            transA = CblasNoTrans;
        }

        int ldA = A_dim1;
        B_dim = B.dim();
#ifndef NDEBUG
        if ( opA_dim2 != B_dim )
            throw std::runtime_error("\nSize mismatch in operands for matrix-vector multiplication!\n\tdim(A) = [ "
                    + std::to_string(opA_dim1) + " x " + std::to_string(opA_dim2) + " ], dim(B) =" + std::to_string(B_dim));
#endif
        RealVector C(opA_dim1);
        cblas_dgemv(CblasColMajor, transA, A_dim1, A_dim2, A.scaling()*B.scaling(), A.ptr(), ldA, B.ptr(), 1, 1.0, C.ptr(), 1);

        return C;
    }

    // Vector-matrix multiplication
    // ---------------------------------------------------------------------------
    RealVector operator*( RealVector&& A, RealMatrix&& B )
    {
        CBLAS_TRANSPOSE transB;
        int B_dim1, B_dim2, opB_dim1, opB_dim2, A_dim;

        B_dim1 = B.dim1();
        B_dim2 = B.dim2();

        if ( !B.isTransposed() )
        {
            opB_dim1 = B_dim2;
            opB_dim2 = B_dim1;
            transB = CblasTrans;
        }
        else
        {
            opB_dim1 = B_dim1;
            opB_dim2 = B_dim2;
            transB = CblasNoTrans;
        }

        int ldB = B_dim1;
        A_dim = A.dim();
#ifndef NDEBUG
        if ( opB_dim2 != A_dim )
            throw std::runtime_error("\nSize mismatch in operands for matrix-vector multiplication!\n\tdim(B) = [ "
                    + std::to_string(opB_dim1) + " x " + std::to_string(opB_dim2) + " ], dim(A) =" + std::to_string(A_dim));
#endif
        RealVector C(opB_dim1);
        cblas_dgemv(CblasColMajor, transB, B_dim1, B_dim2, A.scaling()*B.scaling(), B.ptr(), ldB, A.ptr(), 1, 1.0, C.ptr(), 1);

        return C;
    }

    RealVector operator*( RealVector&& A, const RealMatrix& B )
    {
        CBLAS_TRANSPOSE transB;
        int B_dim1, B_dim2, opB_dim1, opB_dim2, A_dim;

        B_dim1 = B.dim1();
        B_dim2 = B.dim2();

        if ( !B.isTransposed() )
        {
            opB_dim1 = B_dim2;
            opB_dim2 = B_dim1;
            transB = CblasTrans;
        }
        else
        {
            opB_dim1 = B_dim1;
            opB_dim2 = B_dim2;
            transB = CblasNoTrans;
        }

        int ldB = B_dim1;
        A_dim = A.dim();
#ifndef NDEBUG
        if ( opB_dim2 != A_dim )
            throw std::runtime_error("\nSize mismatch in operands for matrix-vector multiplication!\n\tdim(B) = [ "
                    + std::to_string(opB_dim1) + " x " + std::to_string(opB_dim2) + " ], dim(A) =" + std::to_string(A_dim));
#endif
        RealVector C(opB_dim1);
        cblas_dgemv(CblasColMajor, transB, B_dim1, B_dim2, A.scaling()*B.scaling(), B.ptr(), ldB, A.ptr(), 1, 1.0, C.ptr(), 1);

        return C;
    }

    RealVector operator*( const RealVector& A, RealMatrix&& B )
    {
        CBLAS_TRANSPOSE transB;
        int B_dim1, B_dim2, opB_dim1, opB_dim2, A_dim;

        B_dim1 = B.dim1();
        B_dim2 = B.dim2();

        if ( !B.isTransposed() )
        {
            opB_dim1 = B_dim2;
            opB_dim2 = B_dim1;
            transB = CblasTrans;
        }
        else
        {
            opB_dim1 = B_dim1;
            opB_dim2 = B_dim2;
            transB = CblasNoTrans;
        }

        int ldB = B_dim1;
        A_dim = A.dim();
#ifndef NDEBUG
        if ( opB_dim2 != A_dim )
            throw std::runtime_error("\nSize mismatch in operands for matrix-vector multiplication!\n\tdim(B) = [ "
                    + std::to_string(opB_dim1) + " x " + std::to_string(opB_dim2) + " ], dim(A) =" + std::to_string(A_dim));
#endif
        RealVector C(opB_dim1);
        cblas_dgemv(CblasColMajor, transB, B_dim1, B_dim2, A.scaling()*B.scaling(), B.ptr(), ldB, A.ptr(), 1, 1.0, C.ptr(), 1);

        return C;
    }

    RealVector operator*( const RealVector& A, const RealMatrix& B )
    {
        CBLAS_TRANSPOSE transB;
        int B_dim1, B_dim2, opB_dim1, opB_dim2, A_dim;

        B_dim1 = B.dim1();
        B_dim2 = B.dim2();

        if ( !B.isTransposed() )
        {
            opB_dim1 = B_dim2;
            opB_dim2 = B_dim1;
            transB = CblasTrans;
        }
        else
        {
            opB_dim1 = B_dim1;
            opB_dim2 = B_dim2;
            transB = CblasNoTrans;
        }

        int ldB = B_dim1;
        A_dim = A.dim();
#ifndef NDEBUG
        if ( opB_dim2 != A_dim )
            throw std::runtime_error("\nSize mismatch in operands for matrix-vector multiplication!\n\tdim(B) = [ "
                    + std::to_string(opB_dim1) + " x " + std::to_string(opB_dim2) + " ], dim(A) =" + std::to_string(A_dim));
#endif
        RealVector C(opB_dim1);
        cblas_dgemv(CblasColMajor, transB, B_dim1, B_dim2, A.scaling()*B.scaling(), B.ptr(), ldB, A.ptr(), 1, 1.0, C.ptr(), 1);

        return C;
    }

    // Matrix transpose
    // ---------------------------------------------------------------------------
    RealMatrix trp( RealMatrix&& A )
    {
        if ( A._isTransposed )
            A._isTransposed = false;
        else
            A._isTransposed = true;

        return (RealMatrix&&)A;
    }

    RealMatrix trp( const RealMatrix& A )
    {
        return A.transpose();
    }

    // State resolution
    // ------------------------------------------------------------------------
    RealMatrix simplify( RealMatrix&& A )
    {
        A.simplify();
        return (RealMatrix&&)A;
    }

    RealVector simplify( RealVector&& A )
    {
        A.simplify();
        return (RealVector&&)A;
    }

    // Matrix inverse
    // ------------------------------------------------------------------------
    RealMatrix inv( RealMatrix&& A )
    {
        A.simplify();

        int info;
        lapack_int* ipiv;

        int nrows = A.dim1();
        int ncols = A.dim2();

        if ( ncols != nrows )
            throw std::runtime_error("Cannot invert non-square matrix!\ndim(A) = [ " + std::to_string(nrows) + " x " + std::to_string(ncols) + " ]");

        ipiv = new lapack_int[nrows]();
        double *locA = A.ptr();

        dgetrf_(&nrows, &nrows, locA, &nrows, ipiv, &info);
        //info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, nrows, nrows, locA, nrows, ipiv);
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

        double* work;
        int lwork;
        dgetri_(&nrows, locA, &nrows, ipiv, work, &lwork, &info);
        //info = LAPACKE_dgetri(LAPACK_COL_MAJOR, nrows, locA, nrows, ipiv);

        delete[] ipiv;

        if ( info < 0 )
            throw std::runtime_error("Error in matrix inversion: illegal value encountered during factorization!");
        else if (info > 0)
            throw std::runtime_error("Cannot invert singular matrix!");

        return (RealMatrix&&)A;
    }

    RealMatrix inv( const RealMatrix& B )
    {
        RealMatrix A = B;
        A.simplify();

        int info;
        lapack_int* ipiv;
        int nrows = A.dim1();
        int ncols = A.dim2();

        if ( ncols != nrows )
            throw std::runtime_error("Cannot invert non-square matrix!\ndim(A) = [ " + std::to_string(nrows) + " x " + std::to_string(ncols) + " ]");

        ipiv = new lapack_int[nrows]();
        double *locA = A.ptr();

        dgetrf_(&nrows, &nrows, locA, &nrows, ipiv, &info);
        //info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, nrows, nrows, locA, nrows, ipiv);
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

        double * work;
        int lwork;
        dgetri_(&nrows, locA, &nrows, ipiv, work, &lwork, &info);
        //info = LAPACKE_dgetri(LAPACK_COL_MAJOR, nrows, locA, nrows, ipiv);

        delete[] ipiv;

        if ( info < 0 )
            throw std::runtime_error("Error in matrix inversion: illegal value encountered during factorization!");
        else if (info > 0)
            throw std::runtime_error("Cannot invert singular matrix!");

        return A;
    }
}
