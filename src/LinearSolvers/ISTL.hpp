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

#include <config.h>
#include "LinearSolver.hpp"

#ifdef HAVE_DUNE_ISTL

#include <dune/istl/solver.hh>
#include <string>
#include "../Util/RealVector.hpp"

namespace broomstyx
{
    class ISTL : public LinearSolver
    {
        using BlockType = Dune::FieldVector< double, 1 >;
        using BlockVectorType = Dune::BlockVector< BlockType >;
        using InverseOperatorType = Dune::InverseOperator< BlockVectorType, BlockVectorType >;
        
    public:
        ISTL();
        ~ISTL();

        std::string giveRequiredMatrixFormat() override;
        bool        giveSymmetryOption() override;
        void        readDataFrom( FILE* fp ) override;
        void        setInitialGuessTo( RealVector& initGuess ) override;
        RealVector  solve ( SparseMatrix* coefMat, RealVector& rhs ) override;
        bool        takesInitialGuess() override { return true; }

    private:
        double _tol;
        double _abstol;
        int    _maxIter;
        int    _verbosity;

        std::string _preconditioner;
        RealVector  _initGuess;

        std::unique_ptr< InverseOperatorType > _solver;
    };
}

#endif /* HAVE_DUNE_ISTL */
