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

#ifndef ISTL_CPP
#define ISTL_CPP

#include "../../config.h"

#ifdef HAVE_DUNE_ISTL

#include "../Core/ObjectFactory.hpp"
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <string>
#include "../Util/RealVector.hpp"
#include "../SparseMatrix/ISTLMat.hpp"

#include "ISTL.hpp"

using namespace broomstyx;

registerBroomstyxObject(LinearSolver, ISTLSolver)


ISTLSolver::ISTLSolver() {}

std::string
ISTLSolver::giveRequiredMatrixFormat() { return std::string("ISTLMat"); }
void
ISTLSolver::readDataFrom( FILE* fp )
{
    _tol = 1e-8;
    _maxIter = 1000;

    /*
    std::string src = "ViennaCL_cuda (LinearSolver)";

    verifyKeyword(fp, "Algorithm", src);
    _algorithm = getStringInputFrom(fp, "Failed to read iterative algorithm for linear solver from input file!", src);

    _tol = getRealInputFrom(fp, "Failed to read relative tolerance for iterative linear solver from input file!", src);
    _maxIter = getIntegerInputFrom(fp, "Failed to read max. iterations for iterative linear solver from input file!", src);

    if ( _algorithm == "GMRES" )
        _restart = getIntegerInputFrom(fp, "Failed to read number of iterations before restarting GMRES solver from input file!", src);

    verifyKeyword(fp, "Preconditioner", src);
    _preconditioner = getStringInputFrom(fp, "Failed to read preconditioner for linear solver from input file!", src);

    if ( _preconditioner == "Chow_Patel_ILU0" )
    {
        _chowPatel_sweep = getIntegerInputFrom(fp, "Failed to read number of sweeps for preconditioner Chow_Patel_ILU0 from input file!", src);
        _chowPatel_nJacIter = getIntegerInputFrom(fp, "Failed to read number of Jacobi iterations for preconditioner Chow_Patel_ILU0 from input file!", src);
    }
    else if ( _preconditioner == "ILU0" || _preconditioner == "none" )
    {
        // Do nothing.
    }
    else
    {
        std::string errMsg = "Invalid preconditioner tag '" + _preconditioner + "' encountered while reading input file!\nSource: ViennaCL_cuda (LinearSolver)";
        throw std::runtime_error(errMsg);
    }
    */
}

void
ISTLSolver::setInitialGuessTo( RealVector& initGuess )
{
    _initGuess = initGuess;
}

RealVector
ISTLSolver::solve ( SparseMatrix* coefMat, RealVector& rhs )
{
    ISTLMat* matrix = dynamic_cast< ISTLMat* > (coefMat);
    if( ! matrix )
    {
        std::cerr << "ERROR" << std::endl;
        std::abort();
    }

    typedef typename ISTLMat :: MatrixType M;

    M& mat = matrix->exportMatrix();

    typedef Dune::MatrixAdapter< M, BlockVectorType, BlockVectorType >
        AssembledOperatorType;
    AssembledOperatorType op( mat );
    typedef Dune::SeqILU< M, BlockVectorType, BlockVectorType > PreconditionerType;
    PreconditionerType precon( mat, 1.0, true );

    Dune::BiCGSTABSolver< BlockVectorType > solver( op, precon, _tol, int(_maxIter), int(0));

    const int dim = rhs.dim();
    BlockVectorType B( rhs.dim() );
    BlockVectorType X( rhs.dim() );
    for( int i=0; i<dim; ++i )
    {
        B[ i ][0] = rhs( i );
        X[ i ][0] = _initGuess( i );
    }

    Dune::InverseOperatorResult info;
    solver.apply( X, B, _tol, info );

    RealVector x( dim );
    for( int i=0; i<dim; ++i )
    {
        x( i ) = X[ i ][ 0 ];
    }

    return x;
}

#endif /* HAVE_DUNE_ISTL */
#endif /* ISTL_CPP */
