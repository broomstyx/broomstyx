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

#include <cstdio>
#include <stdexcept>
#include <chrono>
#include "omp.h"

#ifdef USE_OPENBLAS
    #include "cblas.h"
    #include "lapacke.h"
#endif

#ifdef USE_MKL_BLAS
    #include "mkl_cblas.h"
    #include "mkl_lapacke.h"
#endif

void test_dgemm()
{
    std::chrono::time_point<std::chrono::system_clock> tic, toc;
    std::chrono::duration<double> tictoc;
    
    int n = 4;
    int m = 6;
    
    double* A = new double[n*m];
    double* B = new double[n*m];
    double* C = new double[n*m];
    
    int nloop = 1000;
    
    tic = std::chrono::high_resolution_clock::now();
#pragma omp parallel for
    for ( int i = 0; i < nloop; i++ )
    {
        // A is n x m, B is m x n
        A[0] = i;
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, m, 1, A, n, B, m, 0, C, n);
    }
    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    double time1 = tictoc.count();
    
    tic = std::chrono::high_resolution_clock::now();
#pragma omp parallel for
    for ( int i = 0; i < nloop; i++ )
    {
        // A is m x n, B is m x n
        A[0] = i;
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, n, n, m, 1, A, m, B, m, 0, C, n);
    }
    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    double time2 = tictoc.count();
    
    std::cout << "Time without transpose = " << time1 << std::endl;
    std::cout << "Time with transpose = " << time2 << std::endl;
    
    delete[] A;
    delete[] B;
    delete[] C;
}

void test2() {
    std::cout << "directBlasCalls test 2" << std::endl;
    std::cout << "%TEST_FAILED% time=0 testname=test2 (directBlasCalls) message=error message sample" << std::endl;
}

int main(int argc, char** argv) {
    std::cout << "%SUITE_STARTING% directBlasCalls" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;

    std::cout << "%TEST_STARTED% test_dgemm (directBlasCalls)" << std::endl;
    test_dgemm();
    std::cout << "%TEST_FINISHED% time=0 test_dgemm (directBlasCalls)" << std::endl;

    std::cout << "%SUITE_FINISHED% time=0" << std::endl;

    return (EXIT_SUCCESS);
}

