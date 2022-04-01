#include "Util/StackRealMatrix.hpp"

#include <chrono>
#include <cstdio>

using namespace broomstyx;

void test_StackRealMatrix_Implementation()
{
	std::printf("---------------------------------------\n");
	std::printf("Test for StackRealMatrix implementation\n");
	std::printf("---------------------------------------\n");

	// Default constructor
	StackRealMatrix<2,3> A;
	A.print("A", 3);

	// Constructor with initializer list
	StackRealMatrix<3,2> B({{1., 2.}, {3., 4.}, {5., 6.}});
	B.print("B", 3);

	// Copy constructor
	StackRealMatrix<3,2> C(B);
	C.print("C", 3);

	// Assignment via list
	StackRealMatrix<3,2> D;
	D = {{7., 8.}, {9., 10.}, {11., 12.}};
	D.print("D", 3);

	// Copy assignment
	StackRealMatrix<3,2> E;
	E = D;
	E.print("E", 3);

	// Addition to self
	E += D;
	E.print("E", 3);

	// Subtraction from self
	E -= D;
	E.print("E", 3);

	// Component access
	A(0,0) = 1.;
	A(0,1) = 2.;
	A.print("A", 3);

	std::printf("Value of A(0,0) is %f\n", A.at(0,0));
	std::printf("Value of A(0,1) is %f\n\n", A.at(0,1));

	// Dimensions
	std::printf("StackRealMatrix A is %d x %d.\n", A.dim1(), A.dim2());

	// Initialize to zero
	A.init();
	A.print("A", 3);
}
