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

#include <config.h>

#ifdef HAVE_DUNE_ISTL

#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/fastamg.hh>
#include <dune/istl/paamg/fastamgsmoother.hh>
#include <dune/istl/paamg/kamg.hh>
#include <string>
#include "../Core/ObjectFactory.hpp"
#include "../SparseMatrix/ISTLMat.hpp"
#include "../Util/readOperations.hpp"
#include "../Util/RealVector.hpp"

#include "ISTL.hpp"

using namespace broomstyx;

registerBroomstyxObject(LinearSolver, ISTL)

namespace detail
{
    template<class M, class X, class Y>
    class OpenMPMatrixAdapter : public Dune::AssembledLinearOperator<M,X,Y>
    {
    public:
        //! export types
        using matrix_type = M;
        using domain_type = X;
        using range_type = Y;
        using field_type = typename X::field_type;

        //! constructor: just store a reference to a matrix
        explicit OpenMPMatrixAdapter (const M& A) : _A_(A) {}

        //! apply operator to x:  \f$ y = A(x) \f$
        virtual void apply (const X& x, Y& y) const
        {
#ifdef _OPENMP
            const size_t nRows = _A_.N();

#pragma omp parallel for
            for( size_t i = 0; i < nRows; ++i )
            {
                const auto& row = _A_[ i ];

                y[ i ] = 0;
                const auto endj = row.end();
                for (auto j=row.begin(); j!=endj; ++j)
                {
                    (*j).umv( x[ j.index() ], y[ i ]);
                }
            }
#else
            _A_.mv(x,y);
#endif
        }

        //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
        virtual void applyscaleadd (field_type alpha, const X& x, Y& y) const
        {
#ifdef _OPENMP
            const size_t nRows = _A_.N();

#pragma omp parallel for
            for( size_t i = 0; i<nRows; ++i )
            {
              const auto& row = _A_[ i ];

              const auto endj = row.end();
              for (auto j=row.begin(); j!=endj; ++j)
              {
                (*j).usmv(alpha, x[ j.index() ], y[ i ]);
              }
            }
#else
            _A_.usmv(alpha,x,y);
#endif
        }

        //! get matrix via *
        virtual const M& getmat () const
        {
            return _A_;
        }

        //! Category of the solver (see SolverCategory::Category)
        virtual Dune::SolverCategory::Category category() const
        {
            return Dune::SolverCategory::sequential;
        }

    private:
        const M& _A_;
    };
}

ISTL::ISTL() {}

ISTL::~ISTL() {}

std::string ISTL::giveRequiredMatrixFormat()
{
    return std::string("ISTLMat");
}

bool ISTL::giveSymmetryOption()
{
    if ( _solverType == "BiCGSTAB" )
        return false;
    else if ( _solverType == "CG" )
        return false;
    else
        throw std::runtime_error("ERROR! Invalid ISTL solver type '" + _solverType + "' encountered during evaluation of symmetry option.\n");
}

void ISTL::readDataFrom( FILE* fp )
{
	std::string src = "ISTL (LinearSolver)";

    _solverType = getStringInputFrom(fp, "Failed to read ISTL solver type from input file!\n", src);
    _tol = getRealInputFrom(fp, "Failed to read relative tolerance for iterative linear solver from input file!\n", src);
    _abstol = getRealInputFrom(fp, "Failed to read absolute tolerance for iterative linear solver from input file!\n", src);
    _maxIter = getIntegerInputFrom(fp, "Failed to read max. iterations for iterative linear solver from input file!\n", src);
    verifyKeyword(fp, "Preconditioner", src);
    _preconditionerType = getStringInputFrom(fp, "Failed to read preconditioner type from input file!\n", src);
    if ( _preconditionerType == "ILU0" || _preconditionerType == "ILDL" )
    {
        // No specific parameters needed
    }
    else if ( _preconditionerType == "AMG" || _preconditionerType == "KAMG" || _preconditionerType == "FastAMG" )
    {
        _maxDepth = getIntegerInputFrom(fp, "Failed to read max. depth for AMG/KAMG preconditioner solver from input file!\n", src);
        _coarseningLimit = getIntegerInputFrom(fp, "Failed to coarsening limit for AMG/KAMG preconditioner solver from input file!\n", src);
        _cycleType = getStringInputFrom(fp, "Failed to read cycle type for AMG/KAMG preconditioner from input file!\n", src);
        if ( _cycleType != "V" && _cycleType != "W" )
            throw std::runtime_error("ERROR! Invalid cycle type '" + _cycleType + "' encountered in input file!\n");
        
        _addmultType = getStringInputFrom(fp, "Failed to read additive/multiplicative multigrid type from input file!\n", src);
        if ( _addmultType != "A" && _addmultType != "M" )
            throw std::runtime_error("ERROR! Invalid additive/multiplicative multigrid type '" + _cycleType + "' encountered in input file!\n");
        
        if ( _preconditionerType != "FastAMG" )
        {
            _smootherType = getStringInputFrom(fp, "Failed to read smoother type for AMG/KAMG preconditioner from input file!\n", src);
            if ( _smootherType != "SSOR" && 
                 _smootherType != "SOR" && 
                 _smootherType != "GS" && 
                 _smootherType != "Jac" && 
                 _smootherType != "ILU" )
                throw std::runtime_error("ERROR! Invalid smoother type '" + _smootherType + "' encountered in input file!\n");
        }
    }
    else
        throw std::runtime_error("ERROR! Invalid preconditioner type '" + _preconditionerType + "' encountered in input file!\n");
        
    _verbosity = getIntegerInputFrom(fp, "Failed to read verbosity level from input file!", src);
    if ( _verbosity < 0 || _verbosity > 2 )
        throw std::runtime_error("ERROR! Invalid verbosity level '" + std::to_string(_verbosity) + "' encountered in input file!\n");
}

void ISTL::setInitialGuessTo( RealVector& initGuess )
{
    _initGuess = initGuess;
}

RealVector ISTL::solve( SparseMatrix* coefMat, RealVector& rhs )
{
    // Proceed only when residual norm is larger than absolute tolerance
    if ( std::sqrt(rhs.dot(rhs)) <= _abstol )
        return _initGuess;
    else
    {
        ISTLMat* matrix = dynamic_cast< ISTLMat* > (coefMat);
        if( ! matrix )
        {
            std::cerr << "ERROR" << std::endl;
            std::abort();
        }

        using MatrixType = typename ISTLMat::MatrixType;

        MatrixType& mat = matrix->exportMatrix();

        using AssembledOperatorType = detail::OpenMPMatrixAdapter< MatrixType, BlockVectorType, BlockVectorType >;
        _op = std::make_shared< AssembledOperatorType >( mat );

        if ( _verbosity > 0 || _preconditionerType == "FastAMG" )
            std::printf("\n"); // For pretty printing

        if ( _preconditionerType == "ILU0" )
            _preconditioner = std::make_shared< Dune::SeqILU< MatrixType, BlockVectorType, BlockVectorType > >( mat, 1.0, true );
        else if ( _preconditionerType == "ILDL" )
            _preconditioner = std::make_shared< Dune::SeqILDL< MatrixType, BlockVectorType, BlockVectorType > >( mat, 1.0 );
        else if ( _preconditionerType == "AMG" || _preconditionerType == "KAMG" || _preconditionerType == "FastAMG" )
        {
            using Operator = Dune::MatrixAdapter< MatrixType, BlockVectorType, BlockVectorType >;
            std::shared_ptr< Operator > fop = std::make_shared< Operator >(mat);
            _amg_op = fop;

            Dune::Amg::Parameters params(_maxDepth, _coarseningLimit);
            
            using Criterion = Dune::Amg::CoarsenCriterion< Dune::Amg::SymmetricCriterion<MatrixType, Dune::Amg::RowSum > >;
            Criterion criterion(params);
            
            if ( _cycleType == "V" )
                criterion.setGamma(1);
            else
                criterion.setGamma(2);
            
            if ( _addmultType == "A" )
                criterion.setAdditive(true);
            else
                criterion.setAdditive(false);
            
            criterion.setDebugLevel(_verbosity);
            if ( _verbosity > 0 )
                std::printf("\n");

            if ( _preconditionerType == "FastAMG" )
            {
                _preconditioner = std::make_shared< Dune::Amg::FastAMG< Operator, BlockVectorType > >(*fop, criterion);
            }
            else
            {
                if ( _smootherType == "SSOR" )
                {
                    using Smoother = Dune::SeqSSOR< MatrixType, BlockVectorType, BlockVectorType >;
                    using SmootherArgs = Dune::Amg::SmootherTraits< Smoother >::Arguments;
                    SmootherArgs smootherArgs;
                    
                    if ( _preconditionerType == "AMG" )
                        _preconditioner = std::make_shared< Dune::Amg::AMG< Operator, BlockVectorType, Smoother > >(*fop, criterion, smootherArgs);
                    else
                        _preconditioner = std::make_shared< Dune::Amg::KAMG< Operator, BlockVectorType, Smoother > >(*fop, criterion, smootherArgs);
                }
                else if ( _smootherType == "SOR" )
                {
                    using Smoother = Dune::SeqSOR< MatrixType, BlockVectorType, BlockVectorType >;
                    using SmootherArgs = Dune::Amg::SmootherTraits< Smoother >::Arguments;
                    SmootherArgs smootherArgs;
                    
                    if ( _preconditionerType == "AMG" )
                        _preconditioner = std::make_shared< Dune::Amg::AMG< Operator, BlockVectorType, Smoother > >(*fop, criterion, smootherArgs);
                    else
                        _preconditioner = std::make_shared< Dune::Amg::KAMG< Operator, BlockVectorType, Smoother > >(*fop, criterion, smootherArgs);
                }
                else if ( _smootherType == "GS" )
                {
                    using Smoother = Dune::SeqGS< MatrixType, BlockVectorType, BlockVectorType >;
                    using SmootherArgs = Dune::Amg::SmootherTraits< Smoother >::Arguments;
                    SmootherArgs smootherArgs;
                    
                    if ( _preconditionerType == "AMG" )
                        _preconditioner = std::make_shared< Dune::Amg::AMG< Operator, BlockVectorType, Smoother > >(*fop, criterion, smootherArgs);
                    else
                        _preconditioner = std::make_shared< Dune::Amg::KAMG< Operator, BlockVectorType, Smoother > >(*fop, criterion, smootherArgs);
                }
                else if ( _smootherType == "Jac" )
                {
                    using Smoother = Dune::SeqJac< MatrixType, BlockVectorType, BlockVectorType >;
                    using SmootherArgs = Dune::Amg::SmootherTraits< Smoother >::Arguments;
                    SmootherArgs smootherArgs;
                    
                    if ( _preconditionerType == "AMG" )
                        _preconditioner = std::make_shared< Dune::Amg::AMG< Operator, BlockVectorType, Smoother > >(*fop, criterion, smootherArgs);
                    else
                        _preconditioner = std::make_shared< Dune::Amg::KAMG< Operator, BlockVectorType, Smoother > >(*fop, criterion, smootherArgs);
                }
                else if ( _smootherType == "ILU" )
                {
                    using Smoother = Dune::SeqILU< MatrixType, BlockVectorType, BlockVectorType >;
                    using SmootherArgs = Dune::Amg::SmootherTraits< Smoother >::Arguments;
                    SmootherArgs smootherArgs;
                    
                    if ( _preconditionerType == "AMG" )
                        _preconditioner = std::make_shared< Dune::Amg::AMG< Operator, BlockVectorType, Smoother > >(*fop, criterion, smootherArgs);
                    else
                        _preconditioner = std::make_shared< Dune::Amg::KAMG< Operator, BlockVectorType, Smoother > >(*fop, criterion, smootherArgs);
                }
            }
        }

        if ( _solverType == "CG" )
            _solver = std::make_shared< Dune::CGSolver< BlockVectorType > >(*_op, *_preconditioner, _tol, _maxIter, _verbosity);
        else if ( _solverType == "BiCGSTAB" )
            _solver = std::make_shared< Dune::BiCGSTABSolver< BlockVectorType > >(*_op, *_preconditioner, _tol, _maxIter, _verbosity);
        else
            throw std::runtime_error("ERROR! Invalid ISTL solver type '" + _solverType + "' encountered in input file!\n");
        
        const int dim = rhs.dim();
        BlockVectorType B( rhs.dim() );
        BlockVectorType X( rhs.dim() );
        for( int i = 0; i < dim; ++i )
        {
            B[ i ][0] = rhs( i );
            X[ i ][0] = _initGuess( i );
        }

        Dune::InverseOperatorResult info;
        _solver->apply( X, B, _tol, info );
        
        if ( _verbosity > 0 )
            std::printf("\n    %-40s"," "); // hack for pretty printing

        RealVector x( dim );
        for( int i = 0; i < dim; ++i )
        {
            x( i ) = X[ i ][ 0 ];
        }

        return x;
    }
}

#endif /* HAVE_DUNE_ISTL */
#endif /* ISTL_CPP */