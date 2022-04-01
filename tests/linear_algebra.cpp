#if 0
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
#include <chrono>
#include <omp.h>

#ifdef USE_MKL_BLAS
#include "mkl.h"
#endif

#include "../source/Util/RealVector.hpp"
#include "../source/Util/RealMatrix.hpp"
#include "../source/Util/linearAlgebra.hpp"

using namespace broomstyx;

void matrixTransposition()
{
    std::cout << "********************************" << std::endl;
    std::cout << "      Matrix transposition" << std::endl;
    std::cout << "********************************" << std::endl;
    
    RealMatrix A({{1,2,3},{4,5,6}});
    A.print("A", 1);
    
    (trp(A)).print("trp(A)", 1);
    
    RealMatrix B;
    B = trp(A);
    B.print("B = trp(A)", 1);
    
    (trp(trp(A))).print("trp(trp(A))", 1);
    
    B = trp(trp(A));
    B.print("B = trp(trp(A))", 1);
    
    B = trp(trp(trp(A)));
    B.print("B = trp(trp(trp(A)))", 1);
}

void scalarOperations()
{
    std::cout << "*****************************************" << std::endl;
    std::cout << "   Scalar multiplication and division" << std::endl;
    std::cout << "*****************************************" << std::endl;
    
    RealMatrix A({{1,2,3},{4,5,6}});
    RealMatrix B;
    
    A.print("A", 1);
    
    B = 2*A;
    B.print("B = 2*A", 1);
    
    B = A*2;
    B.print("B = A*2", 1);
    
    B = 2*trp(A);
    B.print("B = 2*trp(A)", 1);
    
    B = 2*trp(A)*2;
    B.print("B = 2*trp(A)*2", 1);
    
    B = 2*trp(A)/2;
    B.print("B = 2*trp(A)/2", 1);
    
    RealVector C({1,2,3});
    RealVector D;
    
    C.print("C", 1);
    D = 2*C;
    D.print("D = 2*C", 1);
    
    D = C*2;
    D.print("D = C*2", 1);
    
    D = 2*C/2;
    D.print("D = 2*C/2", 1);
}

void matrixAddition()
{
    std::cout << "********************************" << std::endl;
    std::cout << "        Matrix addition" << std::endl;
    std::cout << "********************************" << std::endl;
    
    RealMatrix A({{1,2,3},{4,5,6}});
    RealMatrix B({{1,2,3},{4,5,6}});
    RealMatrix C({{1,4},{2,5},{3,6}});
    RealMatrix D;
    
    D = A + B;
    D.print("D = A + B", 1);
    
    D = A + trp(C);
    D.print("D = A + trp(C)", 1);
    
    D = trp(A) + C;
    D.print("D = trp(A) + C", 1);
    
    D = 2*A + 2*B;
    D.print("D = 2*A + 2*B", 1);
}

void vectorAddition()
{
    std::cout << "********************************" << std::endl;
    std::cout << "        Vector addition" << std::endl;
    std::cout << "********************************" << std::endl;
    
    RealVector A({1,2,3});
    RealVector B({1,2,3});
    RealVector C;
    
    A.print("A", 1);
    B.print("B", 1);
    
    C = A + B;
    C.print("C = A + B", 3);
    
    RealVector D = 2*A + B*2;
    D.print("D = 2*A + B*2", 3);
}

void matrixMultiplication()
{
    std::cout << "********************************" << std::endl;
    std::cout << "        Matrix multiplication" << std::endl;
    std::cout << "********************************" << std::endl;
    
    RealMatrix A({{1,2,3},{4,5,6}});
    RealMatrix B = trp(A);
    B.simplify();
    
    A.print("A", 1);
    B.print("B", 1);
    
    RealMatrix C = A*B;
    C.print("C = A*B", 1);
    
    C = trp(A)*trp(B);
    C.print("C = trp(A)*trp(B)", 1);
    
    //C = (2*A)*trp(0.5*A);
    //C.print("C = (2*A)*trp(0.5*A)", 1);
    (0.5*A).print("0.5*A",1);
    
    trp(0.5*A).print("trp(0.5*A)",1);
    
    C = trp(A/2)*(2*A);
    C.print("C = trp(A/2)*(2*A)", 1);
}

void matrixVectorMultiplication()
{
    std::cout << "********************************" << std::endl;
    std::cout << "  Matrix-vector multiplication" << std::endl;
    std::cout << "********************************" << std::endl;
    
    RealMatrix A({{1,2,3},{2,4,6},{4,8,12},{8,16,24}});
    RealVector B({1,1,1});
    
    A.print("A",1);
    B.print("B",1);
    
    RealMatrix C;
    C = trp(A);
    C.print("C",1);
    
    RealVector D = A*B;
    D.print("D = A*B", 1);
    
    D = trp(C)*B;
    D.print("D = trp(C)*B", 1);
    
    D = B*trp(A);
    D.print("D = B*trp(A)", 1);
    
    D = B*C;
    D.print("D = B*C", 1);
    
    double alpha = 1, beta = 1;
    RealVector E({1,1,1,1});
    
    RealVector F;
    F = alpha*trp(C)*B + beta*E;
    F.print("alpha*trp(C)*B + beta*E;",1);
    
}

void comparison()
{
    std::chrono::time_point<std::chrono::system_clock> tic, toc;
    std::chrono::duration<double> tictoc;
    
    std::cout << "********************************" << std::endl;
    std::cout << "           Comparison" << std::endl;
    std::cout << "********************************" << std::endl;
    
    RealMatrix A(4,6);
    
    int n = 10000000;
    
    double total_time_manual;
    double total_time_lazy;
    
    tic = std::chrono::high_resolution_clock::now();
#pragma omp parallel for
    for ( int i = 0; i < n; i++ )
    {
        RealMatrix B, C;

        B = trp(A);
        B.simplify();
        C = A*B;
    }
    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    total_time_manual = tictoc.count();
    
    tic = std::chrono::high_resolution_clock::now();
#pragma omp parallel for
    for ( int i = 0; i < n; i++ )
    {
        RealMatrix D;
        D = A*trp(A);
    }
    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    total_time_lazy = tictoc.count();
    
    std::cout << "Time using manual transposition: " << total_time_manual << std::endl;
    std::cout << "Time using lazy transposition: " << total_time_lazy << std::endl;
    std::cout << "Speed-up = " << (total_time_manual - total_time_lazy)/total_time_manual*100 << "%" << std::endl;
    std::cout << std::endl;
    
    tic = std::chrono::high_resolution_clock::now();        
#pragma omp parallel for
    for ( int i = 0; i < n; i++ )
    {
        RealMatrix B, C;

        B = trp(A);
        B.simplify();
        C = B*A;
    }
    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    total_time_manual = tictoc.count();
    
    tic = std::chrono::high_resolution_clock::now();   
#pragma omp parallel for
    for ( int i = 0; i < n; i++ )    
    {
        RealMatrix D;
        D = trp(A)*A;
    }
    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    total_time_lazy = tictoc.count();
    
    std::cout << "Time using manual transposition: " << total_time_manual << std::endl;
    std::cout << "Time using lazy transposition: " << total_time_lazy << std::endl;
    std::cout << "Speed-up = " << (total_time_manual - total_time_lazy)/total_time_manual*100 << "%" << std::endl;
    std::cout << std::endl;
}

void test2() {
    std::cout << "linear_algebra test 2" << std::endl;
    std::cout << "%TEST_FAILED% time=0 testname=test2 (linear_algebra) message=error message sample" << std::endl;
}

int main(int argc, char** argv) {
    std::cout << "%SUITE_STARTING% linear_algebra" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;

    std::cout << "%TEST_STARTED% matrixTransposition (linear_algebra)" << std::endl;
    try
    {
        matrixTransposition();
    }
    catch (std::exception& e)
    {
        std::cout << "\n" << e.what() << "\n" << std::endl;
        std::cout << "%TEST_FAILED% time=0 testname=matrixTransposition (linear_algebra) message=Exception occurred." << std::endl;
    }
    std::cout << "%TEST_FINISHED% time=0 matrixTransposition (linear_algebra)" << std::endl;
    
    std::cout << "%TEST_STARTED% scalarOperations (linear_algebra)" << std::endl;
    try
    {
        scalarOperations();
    }
    catch (std::exception& e)
    {
        std::cout << "\n" << e.what() << "\n" << std::endl;
        std::cout << "%TEST_FAILED% time=0 testname=scalarOperations (linear_algebra) message=Exception occurred." << std::endl;
    }
    std::cout << "%TEST_FINISHED% time=0 scalarOperations (linear_algebra)" << std::endl;
    
    std::cout << "%TEST_STARTED% matrixAddition (linear_algebra)" << std::endl;
    try
    {
        matrixAddition();
    }
    catch (std::exception& e)
    {
        std::cout << "\n" << e.what() << "\n" << std::endl;
        std::cout << "%TEST_FAILED% time=0 testname=matrixAddition (linear_algebra) message=Exception occurred." << std::endl;
    }
    std::cout << "%TEST_FINISHED% time=0 matrixAddition (linear_algebra)" << std::endl;
    
    std::cout << "%TEST_STARTED% vectorAddition (linear_algebra)" << std::endl;
    try
    {
        vectorAddition();
    }
    catch (std::exception& e)
    {
        std::cout << "\n" << e.what() << "\n" << std::endl;
        std::cout << "%TEST_FAILED% time=0 testname=vectorAddition (linear_algebra) message=Exception occurred." << std::endl;
    }
    std::cout << "%TEST_FINISHED% time=0 matrixAddition (linear_algebra)" << std::endl;
    
    std::cout << "%TEST_STARTED% matrixMultiplication (linear_algebra)" << std::endl;
    try
    {
        matrixMultiplication();
    }
    catch (std::exception& e)
    {
        std::cout << "\n" << e.what() << "\n" << std::endl;
        std::cout << "%TEST_FAILED% time=0 testname=matrixMultiplication (linear_algebra) message=Exception occurred." << std::endl;
    }
    std::cout << "%TEST_FINISHED% time=0 matrixMultiplication (linear_algebra)" << std::endl;
    
    std::cout << "%TEST_STARTED% matrixVectorMultiplication (linear_algebra)" << std::endl;
    try
    {
        matrixVectorMultiplication();
    }
    catch (std::exception& e)
    {
        std::cout << "\n" << e.what() << "\n" << std::endl;
        std::cout << "%TEST_FAILED% time=0 testname=matrixVectorMultiplication (linear_algebra) message=Exception occurred." << std::endl;
    }
    std::cout << "%TEST_FINISHED% time=0 matrixVectorMultiplication (linear_algebra)" << std::endl;
    
//    std::cout << "%TEST_STARTED% comparison (linear_algebra)" << std::endl;
//    try
//    {
//        comparison();
//    }
//    catch (std::exception& e)
//    {
//        std::cout << "\n" << e.what() << "\n" << std::endl;
//        std::cout << "%TEST_FAILED% time=0 testname=comparison (linear_algebra) message=Exception occurred." << std::endl;
//    }
//    std::cout << "%TEST_FINISHED% time=0 comparison (linear_algebra)" << std::endl;

    std::cout << "%SUITE_FINISHED% time=0" << std::endl;

    return (EXIT_SUCCESS);
}
#endif
