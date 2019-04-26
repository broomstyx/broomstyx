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

#ifdef USE_UB_PARDISO

#ifndef UB_PARDISO_HPP
#define	UB_PARDISO_HPP

#include "LinearSolver.hpp"
#include <string>

// PARDISO function prototypes
extern "C" void pardisoinit(void*, int*, int*, int*, double*, int*);
extern "C" void pardiso(void*, int*, int*, int*, int*, int*, double*, int*, int*, int*, int*, int*, int*, double*, double*, int*, double*);
extern "C" void pardiso_chkmatrix(int*, int*, double*, int*, int*, int*);
extern "C" void pardiso_chkvec(int*, int*, double*, int*);
extern "C" void pardiso_printstats(int*, int*, double*, int*, int*, int*, double*, int*);

namespace broomstyx
{
    class AnalysisModel;

    class UB_Pardiso final : public LinearSolver
    {
    public:
        UB_Pardiso();
        ~UB_Pardiso();

        void        allocateInternalMemoryFor( SparseMatrix* coefMat ) override;
        RealVector  backSubstitute( SparseMatrix* coefMat, RealVector& rhs ) override;
        void        initialize() override;
        void        clearInternalMemory() override;
        void        factorize( SparseMatrix* coefMat ) override;
        std::string giveRequiredMatrixFormat() override;
        bool        giveSymmetryOption() override;
        void        readDataFrom( FILE* fp ) override;
        RealVector  solve( SparseMatrix* coefMat, RealVector& rhs ) override;

    private:
        int   _solver;
        void* _pt[64];
        int   _iparm[64];
        int   _nThreads;
        bool  _memoryIsAllocated;
        int   _mtype;
        bool  _symmetry;
    };
}

#endif	/* UB_PARDISO_HPP */

#endif