#include <stdlib.h>
#include <cstdio>
#include <chrono>
#include <omp.h>

#ifdef USE_MKL_BLAS
#include "mkl.h"
#endif

#include "Util/StackRealVector.hpp"
#include "Util/StackRealMatrix.hpp"
#include "Util/stackLinearAlgebra.hpp"

using namespace broomstyx;

void test_stackmatrix_transposition()
{
	// Matrix transposition
	StackRealMatrix<2,3> A({{1,2,3},{4,5,6}});
	A.print("A", 1);

	(trp(A)).print("trp(A)", 1);

	auto B = trp(A);
	B.print("B = trp(A)", 1);

	(trp(trp(A))).print("trp(trp(A))", 1);

	auto C = trp(trp(A));
	C.print("C = trp(trp(A))", 1);

	auto D = trp(trp(trp(A)));
	D.print("D = trp(trp(trp(A)))", 1);
}

void test_stack_scalar_multiplication_division()
{
    StackRealMatrix<2,3> A({{1,2,3},{4,5,6}});
    A.print("A", 1);

    StackRealMatrix<2,3> B;
    B = 2*A;
    B.print("B = 2*A", 1);

    B = A*2;
    B.print("B = A*2", 1);

    StackRealMatrix<3,2> C;
    C = 2.0*trp(A);
    C.print("C = 2*trp(A)", 1);

    C = 2.*trp(A)*2.;
    C.print("C = 2*trp(A)*2", 1);

    C = 2.*trp(A)/2.;
    C.print("C = 2*trp(A)/2", 1);

    StackRealVector<3> D({1,2,3});
    StackRealVector<3> E;

    D.print("D", 1);
    E = 2*D;
    E.print("E = 2*D", 1);

    E = D*2.;
    E.print("E = D*2", 1);

    E = 2.*D/2.;
    E.print("E = 2*D/2", 1);
}

void test_stack_matrixAddition()
{
    StackRealMatrix<2,3> A({{1,2,3},{4,5,6}});
    StackRealMatrix<2,3> B({{1,2,3},{4,5,6}});
    StackRealMatrix<3,2> C({{1,4},{2,5},{3,6}});

    auto D = A + B;
    D.print("D = A + B", 1);

    D = A + trp(C);
    D.print("D = A + trp(C)", 1);

    auto E = trp(A) + C;
    E.print("E = trp(A) + C", 1);

    D = 2*A + 2*B;
    D.print("D = 2*A + 2*B", 1);
}

void test_stack_vectorAddition()
{
    StackRealVector<3> A({1,2,3});
    StackRealVector<3> B({1,2,3});

    A.print("A", 1);
    B.print("B", 1);

    auto C = A + B;
    C.print("C = A + B", 3);

    auto D = 2*A + B*2;
    D.print("D = 2*A + B*2", 3);
}

void test_stack_matrixMultiplication()
{
    StackRealMatrix<2,3> A({{1,2,3},{4,5,6}});
    auto B = trp(A);

    A.print("A", 1);
    B.print("B", 1);

    auto C = A*B;
    C.print("C = A*B", 1);

    auto D = trp(A)*trp(B);
    D.print("D = trp(A)*trp(B)", 1);

    (0.5*A).print("0.5*A",1);

    trp(0.5*A).print("trp(0.5*A)",1);

    auto E = trp(A/2)*(2*A);
    E.print("E = trp(A/2)*(2*A)", 1);
}

void test_stack_matrixVectorMultiplication()
{
    StackRealMatrix<4,3> A({{1,2,3},{2,4,6},{4,8,12},{8,16,24}});
    StackRealVector<3> B({1,1,1});

    A.print("A",1);
    B.print("B",1);

    auto C = trp(A);
    C.print("C",1);

    auto D = A*B;
    D.print("D = A*B", 1);

    D = trp(C)*B;
    D.print("D = trp(C)*B", 1);

    auto E = B*trp(A);
    E.print("E = B*trp(A)", 1);

    E = B*C;
    E.print("E = B*C", 1);

    double alpha = 1, beta = 1;
    StackRealVector<4> F({1,1,1,1});

    auto G = alpha*trp(C)*B + beta*F;
    G.print("G = alpha*trp(C)*B + beta*F",1);
}

void test_stack_linear_algebra()
{
	std::printf("\n========================");
	std::printf("\n  Matrix Transposition");
	std::printf("\n========================\n");
	test_stackmatrix_transposition();

	std::printf("\n=====================");
	std::printf("\n  Scalar operations");
	std::printf("\n=====================\n");
	test_stack_scalar_multiplication_division();

	std::printf("\n=========================");
	std::printf("\n  Stack Matrix Addition");
	std::printf("\n=========================\n");
	test_stack_matrixAddition();

	std::printf("\n=========================");
	std::printf("\n  Stack Vector Addition");
	std::printf("\n=========================\n");
	test_stack_vectorAddition();

	std::printf("\n===============================");
	std::printf("\n  Stack Matrix Multiplication");
	std::printf("\n===============================\n");
	test_stack_matrixMultiplication();

	std::printf("\n======================================");
	std::printf("\n  Stack Matrix-Vector Multiplication");
	std::printf("\n======================================");
	test_stack_matrixVectorMultiplication();
}
