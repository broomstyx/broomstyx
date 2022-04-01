#ifndef _TEST_HPP_
#define _TEST_HPP_

void test_StackRealMatrix_Implementation();
void test_StackRealVector_Implementation();
void test_stack_linear_algebra();

int perform_tests()
{
//	test_StackRealMatrix_Implementation();
//	test_StackRealVector_Implementation();
	test_stack_linear_algebra();

	return 0;
}

#endif /* _TEST_HPP_ */
