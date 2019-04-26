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

#define VERBOSE_REALVECTOR_STATE_RESOLUTION

#include "../source/Util/RealVector.hpp"
#include "../source/Util/RealMatrix.hpp"

using namespace broomstyx;

void testRealVectorConstruction()
{
    std::cout << "********************************" << std::endl;
    std::cout << "    RealVector construction" << std::endl;
    std::cout << "********************************" << std::endl;
    
    RealVector a, b(3), c({-1, 0, 1, 2, 3});
    
    std::cout << "Constructing empty vector a." << std::endl;
    a.print("a", 3);
    
    std::cout << "Constructing vector b, initialized to 3 entries." << std::endl;
    b.print("b", 3);
    
    std::cout << "Constructing vector c from initializer list." << std::endl; 
    c.print("c", 3);
}

void testRealVectorAssignment()
{
    std::cout << "********************************" << std::endl;
    std::cout << "    RealVector construction" << std::endl;
    std::cout << "********************************" << std::endl;
    
    std::cout << "Constructing vector b from initializer list." << std::endl; 
    RealVector a, b({-1, 0, 1, 2, 3});
    b.print("b", 3);
    
    std::cout << "Assigning a = b." << std::endl; 
    a = b;
    
    a.print("a", 3);
}

void testRealVectorDotProduct()
{
    std::cout << "********************************" << std::endl;
    std::cout << "RealVector dot product method" << std::endl;
    std::cout << "********************************" << std::endl;
    
    RealVector a({1,1,1,1,1});
    a.print("a", 1);
    
    double d = a.dot(a);
    std::cout << "a.dot(a) = " << d << std::endl;
    
    d = (a.scaleBy(2)).dot(a.scaleBy(2));
    std::cout << "(2*a).dot(2*a) = " << d << std::endl;
}

void testRealVectorOuterProduct()
{
    std::cout << "********************************" << std::endl;
    std::cout << "RealVector outer product method" << std::endl;
    std::cout << "********************************" << std::endl;
    
    RealVector a({1, 2, 3, 4, 5}), b({-1, -2});
    RealMatrix c = a.xMen(a);
    RealMatrix d = b.xMen(a);
    
    c.print("c", 3);
    d.print("d", 3);
}

void testRealVectorScaling()
{
    std::cout << "********************************" << std::endl;
    std::cout << "      RealVector scaling" << std::endl;
    std::cout << "********************************" << std::endl;
    
    RealVector a({1, 2, 3}), b;
    
    a.print("a", 3);
    std::cout << "Calculating 2.0*a ..." << std::endl;
    (a.scaleBy(2.)).print("2.0*a", 3);
    
    std::cout << "Calculating b = 3.0*(2.0*a) ..." << std::endl;
    b = (a.scaleBy(2.)).scaleBy(3.);
    b.print("b", 3);
    
    std::cout << "Calculating c = 2.0*a ..." << std::endl;
    RealVector c = a.scaleBy(2.);
    c.print("c", 3);
    
    std::cout << "Important: The compiler may perform return value optimization/copy elision," << std::endl;
    std::cout << "resulting in c having non-unity scaling and/or not having ownership of its pointer." << std::endl;
    std::cout << "In this case, component access by address is forbidden and will throw an exception." << std::endl;
    std::cout << std::endl;
    
    std::cout << "Resolving state of c ..." << std::endl;
    c.simplify();
    c.print("c", 3);
}

void testRealVectorOperator_plusEquals()
{
    std::cout << "********************************" << std::endl;
    std::cout << "      RealVector operator +=" << std::endl;
    std::cout << "********************************" << std::endl;
    
    RealVector a({1, 2, 3});
    RealVector b;
    b = a;
    
    a += b;
    a.print("a += b", 3);
}

int main(int argc, char** argv)
{
    std::cout << "%SUITE_STARTING% RealVector_class_implementation" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;

    std::cout << "%TEST_STARTED% Construction (RealVector_class_implementation)" << std::endl;
    try
    {
        testRealVectorConstruction();
    }
    catch (std::exception& e)
    {
        std::cout << "\n" << e.what() << "\n" << std::endl;
        std::cout << "%TEST_FAILED% time=0 testname=Construction (RealVector_class_implementation) message=Exception occurred." << std::endl;
    }
    std::cout << "%TEST_FINISHED% time=0 Construction (RealVector_class_implementation)" << std::endl;

    std::cout << "%TEST_STARTED% Assignment (RealVector_class_implementation)\n" << std::endl;
    try
    {
        testRealVectorAssignment();
    }
    catch (std::exception& e)
    {
        std::cout << "\n" << e.what() << "\n" << std::endl;
        std::cout << "%TEST_FAILED% time=0 testname=Assignment (RealVector_class_implementation) message=Exception occurred." << std::endl;
    }
    std::cout << "%TEST_FINISHED% time=0 Assignment (RealVector_class_implementation)" << std::endl;
    
    std::cout << "%TEST_STARTED% dotProduct (RealVector_class_implementation)\n" << std::endl;
    try
    {
        testRealVectorDotProduct();
    }
    catch (std::exception& e)
    {
        std::cout << "\n" << e.what() << "\n" << std::endl;
        std::cout << "%TEST_FAILED% time=0 testname=dotProduct (RealVector_class_implementation) message=Exception occurred." << std::endl;
    }
    std::cout << "%TEST_FINISHED% time=0 dotProduct (RealVector_class_implementation)" << std::endl;
    
    std::cout << "%TEST_STARTED% OuterProduct (RealVector_class_implementation)\n" << std::endl;
    try
    {
        testRealVectorOuterProduct();
    }
    catch (std::exception& e)
    {
        std::cout << "\n" << e.what() << "\n" << std::endl;
        std::cout << "%TEST_FAILED% time=0 testname=OuterProduct (RealVector_class_implementation) message=Exception occurred." << std::endl;
    }
    std::cout << "%TEST_FINISHED% time=0 OuterProduct (RealVector_class_implementation)" << std::endl;
    
    std::cout << "%TEST_STARTED% Scaling (RealVector_class_implementation)\n" << std::endl;
    try
    {
        testRealVectorScaling();
    }
    catch (std::exception& e)
    {
        std::cout << "\n" << e.what() << "\n" << std::endl;
        std::cout << "%TEST_FAILED% time=0 testname=Scaling (RealVector_class_implementation) message=Exception occurred." << std::endl;
    }
    std::cout << "%TEST_FINISHED% time=0 Scaling (RealVector_class_implementation)" << std::endl;
    
    std::cout << "%TEST_STARTED% Operator+= (RealVector_class_implementation)\n" << std::endl;
    try
    {
        testRealVectorOperator_plusEquals();
    }
    catch (std::exception& e)
    {
        std::cout << "\n" << e.what() << "\n" << std::endl;
        std::cout << "%TEST_FAILED% time=0 testname=Operator+= (RealVector_class_implementation) message=Exception occurred." << std::endl;
    }
    std::cout << "%TEST_FINISHED% time=0 Operator+= (RealVector_class_implementation)" << std::endl;
    
    std::cout << "%SUITE_FINISHED% time=0" << std::endl;

    return (EXIT_SUCCESS);
}