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

#ifndef VIENNACL_CUDA_HPP
#define VIENNACL_CUDA_HPP

#include "../../config.h"

#ifdef HAVE_VIENNACL
#include "LinearSolver.hpp"
#include <string>
#include "../Util/RealVector.hpp"

namespace broomstyx
{
    class ViennaCL_cuda : public LinearSolver
    {
    public:
        ViennaCL_cuda();
        virtual ~ViennaCL_cuda();
        
        std::string giveRequiredMatrixFormat() override;
        void        readDataFrom( FILE* fp ) override;
        void        setInitialGuessTo( RealVector& initGuess ) override;
        RealVector  solve ( SparseMatrix* coefMat, RealVector& rhs ) override;
        bool        takesInitialGuess() override;
    
//        struct PrecondParam
//        {
//            int    _chowPatel_sweep;
//            int    _chowPatel_nJacIter;
//            int    _ilut_entriesPerRow;
//            double _ilut_dropTol;
//            bool   _ilut_levelSched;
//            bool   _ilu0_levelSched;
//            int    _rowScal_norm;
//            double _amg_strConThresh;
//            double _amg_jacSmthWt;
//            int    _amg_nPreSmooth;
//            int    _amg_nPostSmooth;
//            int    _amg_maxCrsLvl;
//            int    _amg_crsLvlCutoff;
//            int    _amg_precSetup;
//            int    _amg_precApp;
//        };
        
    private:
        double _tol;
        int    _maxIter;
        int    _restart;
        
        std::string  _algorithm;
        std::string  _preconditioner;
        RealVector   _initGuess;
        
        // Preconditioner parameters
        int    _chowPatel_sweep;
        int    _chowPatel_nJacIter;
    };
}

#endif /* HAVE_VIENNACL */
#endif /* VIENNACL_CUDA_HPP */
