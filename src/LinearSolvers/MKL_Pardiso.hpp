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

#ifndef MKL_PARDISO_HPP
#define MKL_PARDISO_HPP

#include "../../config.h"

#ifdef HAVE_MKL

#include "LinearSolver.hpp"
#include "mkl_pardiso.h"
#include "mkl_types.h"

namespace broomstyx
{
    class AnalysisModel;
    
    class MKL_Pardiso final : public LinearSolver
    {
    public:
        MKL_Pardiso();
        ~MKL_Pardiso();

        void        allocateInternalMemoryFor( SparseMatrix* coefMat ) override;
        void        initialize() override;
        void        clearInternalMemory() override;
        std::string giveRequiredMatrixFormat() override;
        bool        giveSymmetryOption() override;
        void        readDataFrom( FILE* fp ) override;
        RealVector  solve( SparseMatrix* coefMat, RealVector& rhs ) override;

    private:
        int  _nThreads;
        bool _memoryIsAllocated;
        bool _symmetry;
        
        void* _pt[64];
        int   _iparm[64];
        int   _mtype;
        
        int _n;
        int _maxfct;
        int _mnum;
        int _msglvl;
        int _nrhs;
        
        void giveErrorMessage( int error );
    };
}
    
#endif /* HAVE_MKL */
#endif /* MKL_PARDISO_HPP */
