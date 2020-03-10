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

#include "ViennaCL_cuda.hpp"

#ifdef HAVE_VIENNACL

#include "Core/ObjectFactory.hpp"
#include "Util/readOperations.hpp"
#include "Util/RealVector.hpp"
#include "SparseMatrix/SparseMatrix.hpp"

#include "viennacl/scalar.hpp"
#include "viennacl/vector.hpp"
#include "viennacl/compressed_matrix.hpp"
#include "viennacl/linalg/prod.hpp"
#include "viennacl/linalg/ilu_operations.hpp"
#include "viennacl/linalg/ilu.hpp"
#include "viennacl/linalg/bicgstab.hpp"
#include "viennacl/linalg/gmres.hpp"

void copy_csr0_matrix( broomstyx::SparseMatrix* coefMat, viennacl::compressed_matrix<double>& gpuMatrix )
{
    int* csr_rows;
    int* csr_cols;
    double * csr_elements;
    int num_rows, num_cols, num_nnz;
    
    std::tie(csr_rows, csr_cols) = coefMat->giveProfileArrays();
    csr_elements = coefMat->giveValArray();
    std::tie(num_rows, num_cols) = coefMat->giveMatrixDimensions();
    num_nnz = coefMat->giveNumberOfNonzeros();
    
    if ( num_rows > 0 && num_cols > 0 && num_nnz > 0)
    {
        viennacl::backend::typesafe_host_array<unsigned int> row_buffer(gpuMatrix.handle1(), num_rows + 1);

        if (sizeof(int) != row_buffer.element_size()) // check whether indices are of the same length (same number of bits)
        {
            viennacl::backend::typesafe_host_array<unsigned int> col_buffer(gpuMatrix.handle2(), num_nnz);

            for ( int i=0; i<=num_rows; ++i )
                row_buffer.set( i, csr_rows[i] );
            for ( int i=0; i<num_nnz; ++i )
                col_buffer.set(i, csr_cols[i]);

            gpuMatrix.set(row_buffer.get(), col_buffer.get(), csr_elements, num_rows, num_cols, num_nnz);
        }
        else
            gpuMatrix.set(static_cast<const void*>(csr_rows), static_cast<const void*>(csr_cols), csr_elements, num_rows, num_cols, num_nnz);
    }
}

using namespace broomstyx;

registerBroomstyxObject(LinearSolver, ViennaCL_cuda)

// Constructor
ViennaCL_cuda::ViennaCL_cuda() : LinearSolver() {}

// Destructor
ViennaCL_cuda::~ViennaCL_cuda() {}

// Public methods
// ---------------------------------------------------------------------------
std::string ViennaCL_cuda::giveRequiredMatrixFormat()
{
    return std::string("CSR0");
}
// ---------------------------------------------------------------------------
void ViennaCL_cuda::readDataFrom( FILE* fp )
{
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
}
// ---------------------------------------------------------------------------
void ViennaCL_cuda::setInitialGuessTo( RealVector& initGuess )
{
//    // Temporary hack: last solution becomes initial guess :)
//    if ( _initGuess.dim() == 0 )
        _initGuess = initGuess;
}
// ---------------------------------------------------------------------------
RealVector ViennaCL_cuda::solve( SparseMatrix* coefMat, RealVector& rhs )
{
    // Instantiate objects
    int nUnknowns = rhs.dim();
    
    RealVector soln(nUnknowns);
    viennacl::compressed_matrix<double> gpuMatrix;
    viennacl::vector<double> gpuRhs(nUnknowns);
    viennacl::vector<double> gpuResult(nUnknowns);
    viennacl::vector<double> gpuInitGuess(nUnknowns);
    
    // Copy data from CPU into GPU
    copy_csr0_matrix(coefMat, gpuMatrix);
    viennacl::copy(rhs.ptr(), rhs.ptr() + nUnknowns, gpuRhs.begin());
    viennacl::copy(_initGuess.ptr(), _initGuess.ptr() + nUnknowns, gpuInitGuess.begin());
    
    // Setup iterative solver
    if ( _algorithm == "BiCGStab" )
    {
        viennacl::linalg::bicgstab_tag solverTag(_tol, _maxIter);
        viennacl::linalg::bicgstab_solver<viennacl::vector<double> > iterSolver(solverTag);
        iterSolver.set_initial_guess(gpuInitGuess);
        
        // Setup preconditioner
        if ( _preconditioner == "Chow_Patel_ILU0" )
        {
            viennacl::linalg::chow_patel_tag pcConfig;
            pcConfig.sweeps(_chowPatel_sweep);
            pcConfig.jacobi_iters(_chowPatel_nJacIter);
            viennacl::linalg::chow_patel_ilu_precond< viennacl::compressed_matrix<double> > pcObject(gpuMatrix, pcConfig);
            std::printf("\n      Preconditioner setup completed.");
            std::fflush(stdout);
            gpuResult = iterSolver(gpuMatrix, gpuRhs, pcObject);
        }
        else if ( _preconditioner == "ILU0" )
        {
            viennacl::linalg::ilu0_tag pcConfig;
            viennacl::linalg::ilu0_precond< viennacl::compressed_matrix<double> > pcObject(gpuMatrix, pcConfig);
            std::printf("\n      Preconditioner setup completed.");
            std::fflush(stdout);
            gpuResult = iterSolver(gpuMatrix, gpuRhs, pcObject);
        }
        else if ( _preconditioner == "none" )
        {
            gpuResult = iterSolver(gpuMatrix, gpuRhs);
        }
        else
        {
            std::string errmsg = "Preconditioner '" + _preconditioner + "' is not yet programmed!\n";
            throw std::runtime_error(errmsg);
        }
        
        std::printf("\n      System solved.");
        std::printf("\n      Num iters = %d, est. error = %e\n", (int)iterSolver.tag().iters(), iterSolver.tag().error());
    }
    else if ( _algorithm == "GMRES" )
    {
        viennacl::linalg::gmres_tag solverTag(_tol, _maxIter, _restart);
        viennacl::linalg::gmres_solver<viennacl::vector<double> > iterSolver(solverTag);
        iterSolver.set_initial_guess(gpuInitGuess);
        
        // Setup preconditioner
        if ( _preconditioner == "Chow_Patel_ILU0" )
        {
            viennacl::linalg::chow_patel_tag pcConfig;
            pcConfig.sweeps(_chowPatel_sweep);
            pcConfig.jacobi_iters(_chowPatel_nJacIter);
            viennacl::linalg::chow_patel_ilu_precond< viennacl::compressed_matrix<double> > pcObject(gpuMatrix, pcConfig);
            std::printf("\n      Preconditioner setup completed.");
            std::fflush(stdout);
            gpuResult = iterSolver(gpuMatrix, gpuRhs, pcObject);
        }
        else if ( _preconditioner == "ILU0" )
        {
            viennacl::linalg::ilu0_tag pcConfig;
            viennacl::linalg::ilu0_precond< viennacl::compressed_matrix<double> > pcObject(gpuMatrix, pcConfig);
            std::printf("\n      Preconditioner setup completed.");
            std::fflush(stdout);
            gpuResult = iterSolver(gpuMatrix, gpuRhs, pcObject);
        }
        else if ( _preconditioner == "none" )
        {
            gpuResult = iterSolver(gpuMatrix, gpuRhs);
        }
        else
        {
            std::string errmsg = "Preconditioner '" + _preconditioner + "' is not yet programmed!\n";
            throw std::runtime_error(errmsg);
        }
        
        std::printf("\n      System solved.");
        std::printf("\n      Num iters = %d, est. error = %e\n", (int)iterSolver.tag().iters(), iterSolver.tag().error());
    }
    else
    {
        std::string errmsg = "Iterative algorithm '" + _algorithm + "' is not yet programmed!\n";
        throw std::runtime_error(errmsg);
    }
    
    // Copy result vector from GPU into CPU
    viennacl::copy(gpuResult.begin(), gpuResult.end(), soln.ptr());
    
//    // Copy to initial guess (just a hack at the moment)
//    _initGuess = soln;
    
//    std::printf("    %-40s", "");
    
    return soln;
}
// ---------------------------------------------------------------------------
bool ViennaCL_cuda::takesInitialGuess()
{
    return true;
}

#endif /* HAVE_VIENNACL */
