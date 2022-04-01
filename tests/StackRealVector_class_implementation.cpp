#include "Util/StackRealVector.hpp"

#include <chrono>
#include <cstdio>

using namespace broomstyx;

void test_StackRealVector_Implementation()
{
	std::printf("---------------------------------------\n");
	std::printf("Test for StackRealVector implementation\n");
	std::printf("---------------------------------------\n");

	// Default constructor
	StackRealVector<2> A;
	A.print("A", 3);

	// Constructor with initializer list
	StackRealVector<3> B({1., 2., 3.});
	B.print("B", 3);

	// Copy constructor
	StackRealVector<3> C(B);
	C.print("C", 3);

	// Assignment via list
	StackRealVector<4> D;
	D = {5., 6., 7., 8.};
	D.print("D", 3);

	// Copy assignment
	StackRealVector<4> E;
	E = D;
	E.print("E", 3);

	// Copy RealVector to StackRealVector
	RealVector RE({1., 2., 3., 4., 5.});
	StackRealVector<5> SE;
	SE = RE;
	SE.print("SE", 3);

	// Addition to self
	E += D;
	E.print("E", 3);

	// Subtraction from self
	E -= D;
	E.print("E", 3);

	// In-place scalar multiplication
	E *= 2.;
	E.print("E", 3);

	// In-place scalar division
	E /= 2.;
	E.print("E", 3);

	// Component access
	A(0) = 1.;
	A(1) = 2.;
	A.print("A", 3);

	std::printf("Value of A(0) is %f\n", A.at(0));
	std::printf("Value of A(1) is %f\n\n", A.at(1));

	// Dimensions
	std::printf("StackRealVector A has length %d.\n", A.dim());

	// Vector cross product
	StackRealVector<3> G({1.0, 2.0, 3.0});
	G.print("G", 3);
	StackRealVector<3> H({4.0, 5.0, 6.0});
	H.print("H", 3);

	auto I = G.cross(H);
	I.print("G x H", 3);

	// Vector dot product
	std::printf("G . H = %lf\n\n", G.dot(H));

	// Initialize to zero
	A.init();
	A.print("A", 3);

	// Tensor product
	auto F = D.xMen<3>(B);
	F.print("F", 3);
}
