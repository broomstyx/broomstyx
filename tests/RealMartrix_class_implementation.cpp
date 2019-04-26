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

#include <stdlib.h>
#include <iostream>

#define VERBOSE_REALMATRIX_STATE_RESOLUTION

#include "../source/Util/RealVector.hpp"
#include "../source/Util/RealMatrix.hpp"

using namespace broomstyx;

void testRealMatrixConstruction()
{
    std::cout << "********************************" << std::endl;
    std::cout << "    RealMatrix construction" << std::endl;
    std::cout << "********************************" << std::endl;
    
    RealMatrix A, B(3,2), C({{1, 2, 3}, {4, 5, 6}});
    
    std::cout << "Constructing empty matrix A." << std::endl;
    A.print("A", 3);
    
    std::cout << "Constructing matrix B with dim 3 x 2" << std::endl;
    B.print("B", 3);
    
    std::cout << "Constructing matrix C from initializer list." << std::endl; 
    C.print("c", 3);
}

void testRealMatrixAssignment()
{
    std::cout << "********************************" << std::endl;
    std::cout << "    RealMatrix construction" << std::endl;
    std::cout << "********************************" << std::endl;
    
    std::cout << "Constructing matrix b from initializer list." << std::endl; 
    RealMatrix A, B({{1, 1, 1}, {2, 3, 4}});
    B.print("B", 3);
    
    std::cout << "Assigning A = B." << std::endl; 
    A = B;
    
    A.print("A", 3);
    
    A = {{1,2},{3,4},{5,6}};
    A.print("A = {{1,2},{3,4},{5,6}}",1);
}

void testRealMatrixScaling()
{
    std::cout << "********************************" << std::endl;
    std::cout << "      RealMatrix scaling" << std::endl;
    std::cout << "********************************" << std::endl;
    
    RealMatrix A({{1, 2}, {3, 4}, {5, 6}}), B;
    
    A.print("A", 3);
    std::cout << "Calculating 2.0*A ..." << std::endl;
    (A.scaleBy(2.)).print("2.0*A", 3);
    
    std::cout << "Calculating B = 3.0*(2.0*A) ..." << std::endl;
    B = (A.scaleBy(2.)).scaleBy(3.);
    B.print("B", 3);
    
    std::cout << "Calculating C = 2.0*a ..." << std::endl;
    RealMatrix C = A.scaleBy(2.);
    C.print("C", 3);
    
    std::cout << "Important: The compiler may perform return value optimization/copy elision," << std::endl;
    std::cout << "resulting in C having non-unity scaling and/or not having ownership of its pointer." << std::endl;
    std::cout << "In this case, component access by address is forbidden and will throw an exception." << std::endl;
    std::cout << std::endl;
    
    std::cout << "Resolving state of C ..." << std::endl;
    C.simplify();
    C.print("C", 3);
}

void testRealMatrixTransposition()
{
    std::cout << "********************************" << std::endl;
    std::cout << "      RealMatrix transposition" << std::endl;
    std::cout << "********************************" << std::endl;
    
    RealMatrix A({{1, 2}, {3, 4}, {5, 6}}), B;
    
    A.print("A", 3);
    std::cout << "Calculating A^T ..." << std::endl;
    A.transpose().print("A.transpose()", 3);
    
    std::cout << "Assigning B = A^T" << std::endl;
    B = A.transpose();
    B.print("B", 3);
    
    std::cout << "Constructing C = (2.0*A)^T" << std::endl;
    RealMatrix C = A.scaleBy(2.).transpose();
    C.print("C", 3);
    
    C.simplify();
    C.print("C", 3);
}

int main(int argc, char** argv)
{
    std::cout << "%SUITE_STARTING% RealMatrix_class_implementation" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;

    std::cout << "%TEST_STARTED% Construction (RealMatrix_class_implementation)\n" << std::endl;
    try 
    {
        testRealMatrixConstruction();
    }
    catch (std::exception& e)
    {
        std::cout << "\n" << e.what() << "\n" << std::endl;
        std::cout << "%TEST_FAILED% time=0 testname=Construction (RealMatrix_class_implementation) message=Exception occurred." << std::endl;
    }
    std::cout << "%TEST_FINISHED% time=0 Construction (RealMatrix_class_implementation)" << std::endl;
    
    std::cout << "%TEST_STARTED% Assignment (RealMatrix_class_implementation)\n" << std::endl;
    try 
    {
        testRealMatrixAssignment();
    }
    catch (std::exception& e)
    {
        std::cout << "\n" << e.what() << "\n" << std::endl;
        std::cout << "%TEST_FAILED% time=0 testname=Assignment (RealMatrix_class_implementation) message=Exception occurred." << std::endl;
    }
    std::cout << "%TEST_FINISHED% time=0 Assignment (RealMatrix_class_implementation)" << std::endl;
    
    std::cout << "%TEST_STARTED% Scaling (RealMatrix_class_implementation)\n" << std::endl;
    try
    {
        testRealMatrixScaling();
    }
    catch (std::exception& e)
    {
        std::cout << "\n" << e.what() << "\n" << std::endl;
        std::cout << "%TEST_FAILED% time=0 testname=Scaling (RealMatrix_class_implementation) message=Exception occurred." << std::endl;
    }
    std::cout << "%TEST_FINISHED% time=0 Scaling (RealMatrix_class_implementation)" << std::endl;
    
    std::cout << "%TEST_STARTED% Transposition (RealMatrix_class_implementation)\n" << std::endl;
    try
    {
        testRealMatrixTransposition();
    }
    catch (std::exception& e)
    {
        std::cout << "\n" << e.what() << "\n" << std::endl;
        std::cout << "%TEST_FAILED% time=0 testname=Transposition (RealMatrix_class_implementation) message=Exception occurred." << std::endl;
    }
    std::cout << "%TEST_FINISHED% time=0 Transposition (RealMatrix_class_implementation)" << std::endl;
    
    std::cout << "%SUITE_FINISHED% time=0" << std::endl;

    return (EXIT_SUCCESS);
}

