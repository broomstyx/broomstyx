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

#include "MKL_Pardiso.hpp"

#ifdef HAVE_MKL

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <stdexcept>
#include <string>
#include <tuple>

#include "Core/ObjectFactory.hpp"
#include "SparseMatrix/SparseMatrix.hpp"
#include "Util/RealVector.hpp"
#include "Util/readOperations.hpp"
#include "mkl.h"

#define SKIP_SYMBOLIC_FACTORIZATION true

using namespace broomstyx;

registerBroomstyxObject(LinearSolver, MKL_Pardiso)

// Constructor
MKL_Pardiso::MKL_Pardiso() : LinearSolver()
{
    _nThreads = 0;
    _memoryIsAllocated = false;
    _mtype = 0;
    _symmetry = false;
    
    for ( int i = 0; i < 64; i++ )
    {
        _pt[i] = 0;
        _iparm[i] = 0;
    }
}

// Destructor
MKL_Pardiso::~MKL_Pardiso()
{
    this->clearInternalMemory();
}

// Public Methods
// ----------------------------------------------------------------------------
void MKL_Pardiso::allocateInternalMemoryFor( SparseMatrix* coefMat )
{
    if ( !_memoryIsAllocated )
    {
        // Check that coefficient is square
        int dim1, dim2;
        std::tie(dim1,dim2) = coefMat->giveMatrixDimensions();
        
        if ( dim1 != dim2 )
        {
            std::printf("\n\nCannot solve specified linear system!\n\tCoefficient matrix is not square");
            std::printf("\n\t--> dim(A) = [ %d x %d ]\n", dim1, dim2);
            std::fflush(stdout);
            throw std::runtime_error("Source: Pardiso");
        }
        int idum;
        
        int* ia;
        int* ja;
        
        std::tie(ia, ja) = coefMat->giveProfileArrays();
        double* a = coefMat->giveValArray();
        
        double ddum;
        
        // ******************************************************************
        // Reordering phase and symbolic factorization. This also allocates
        // all memory necessary for the factorization
        // ******************************************************************
        
        int error = 0;
        int phase = 11;

        pardiso (_pt, &_maxfct, &_mnum, &_mtype, &phase, &dim1, a, ia, ja, &idum, &_nrhs, _iparm, &_msglvl, &ddum, &ddum, &error);
        if ( error != 0 )
        {
            std::printf("\nError during symbolic factorization.\n");
            this->giveErrorMessage(error);
        }
        
        _memoryIsAllocated = true;
    }
}
// ----------------------------------------------------------------------------
void MKL_Pardiso::clearInternalMemory()
{
    if ( _memoryIsAllocated )
    {
        double ddum;
        int idum;
        int error = 0;
        // Memory release
        int phase = -1;
        pardiso(_pt, &_maxfct, &_mnum, &_mtype, &phase, &_n, &ddum, &idum, &idum, &idum, &_nrhs, _iparm, &_msglvl, &ddum, &ddum, &error);
        this->giveErrorMessage(error);
        
        _memoryIsAllocated = false;
    }
}
// ----------------------------------------------------------------------------
void MKL_Pardiso::initialize()
{
    // *********************************************************************
    // Note that this class calls the DIRECT in-core version of MKL Pardiso
    // and only works with REAL-valued matrices.
    //
    // The iterative version will be implemented in a separate class.
    // *********************************************************************
    
    // Set number of threads for Pardiso
#ifdef _OPENMP
    int error = mkl_domain_set_num_threads(_nThreads, MKL_DOMAIN_PARDISO);
#else
    int error = mkl_domain_set_num_threads(1, MKL_DOMAIN_PARDISO);
#endif
    if ( error == 0 )
        throw std::runtime_error("ERROR: Failed to set specified number of threads for MKL Pardiso!\n");
    
    _maxfct = 1;
    _mnum = 1;
    _msglvl = 0;
    _nrhs = 1;
    
    _memoryIsAllocated = false;
    
    // Set IPARM values
    _iparm[0] = 1;   // Manually specify values for IPARM
    _iparm[1] = 2;   // Use nested dissection from METIS
    _iparm[3] = 0;   // No iterative-direct algorithm
    _iparm[4] = 0;   // No user fill-in reducing permutation
    _iparm[5] = 0;   // Write solution into x
    _iparm[7] = 2;   // Max number of iterative refinement steps
    _iparm[9] = 13;  // Pivoting perturbation for unsymmetric matrices (1E-13)
    _iparm[10] = 1;  // Use nonsymmetric permutation and scaling MPS
    _iparm[11] = 0;  // Conjugate transposed/transpose solve
    _iparm[12] = 1;  // Maximum weighted matching algorithm is switched on (default for non-symmetric)
    _iparm[20] = 1;  // Default pivoting for symmetric indefinite matrices

    if ( _nThreads >= 8 )
        _iparm[23] = 10; // Use improved two-level factorization
    else
        _iparm[23] = 0;
    
#ifdef _OPENMP
    _iparm[24] = 1; // Use parallel algorithm for forward/back substitution
#else
    _iparm[24] = 0; // Use sequential algorithm for forward/back substitution
#endif
    _iparm[26] = 0; // Do not perform matrix checking
    _iparm[34] = 0; // Use Fortran-style indexing in the sparse matrix profile
    
    pardisoinit(_pt, &_mtype, _iparm);
}
// ----------------------------------------------------------------------------
std::string MKL_Pardiso::giveRequiredMatrixFormat()
{
    return std::string("CSR1");
}
// ----------------------------------------------------------------------------
bool MKL_Pardiso::giveSymmetryOption()
{
    return _symmetry;
}
// ----------------------------------------------------------------------------
void MKL_Pardiso::readDataFrom( FILE* fp )
{
    std::string key, errmsg, src = "Pardiso (LinearSolver)";

    verifyKeyword(fp, key = "nThreads", src);
    _nThreads = getIntegerInputFrom(fp, errmsg = "Failed to read number of threads for Pardiso solver in input file!", src);
    
    key = getStringInputFrom(fp, "Failed to read symmetry option in input file!", src);
    if ( key == "Symmetric" )
        _mtype = -2;
    else if ( key == "SymmetricPositiveDefinite")
        _mtype = 2;
    else if ( key == "Unsymmetric" )
        _mtype = 11;
    else
        throw std::runtime_error("Unrecognized option '" + key + "' for Pardiso solver encountered in input file!");
    
    if ( _mtype == 11 )
        _symmetry = false;
    else
        _symmetry = true;
}
// ----------------------------------------------------------------------------
RealVector MKL_Pardiso::solve( SparseMatrix* coefMat, RealVector& rhs )
{
    if ( !_memoryIsAllocated )
        throw std::runtime_error("MKL Pardiso solve phase called without proper memory allocation!");
    
    int dim1, dim2;
    std::tie(dim1,dim2) = coefMat->giveMatrixDimensions();

    // Check that sparse matrix and RHS vector have compatible dimensions
    if ( dim1 != rhs.dim() )
    {
        std::printf("\nCannot solve specified linear system!");
        std::printf("\n\tSparse coefficient matrix has dimensions [ %d x %d ]", dim1, dim2);
        std::printf("\n\tRHS vector has dimension [ %d ]\n", rhs.dim());

        throw std::runtime_error("Source: Pardiso");
    }
    _n = dim1;
    
    // Initialize solution vector
    RealVector u(_n);

    // Matrix Data
    int* ia;
    int* ja;
    double* a;

    std::tie(ia,ja) = coefMat->giveProfileArrays();
    a  = coefMat->giveValArray();

    // RHS and solution vectors
    double *b, *x;
    b = rhs.ptr();
    x = u.ptr();

    int error = 0;
    
    // ***********************************************************************
    // Numerical factorization, back-substitution and iterative refinement.
    // 
    // Note that symbolic factorization and memory allocation is actually also
    // done making use of the permutation vector stored in the _perm array.
    // ***********************************************************************

    // Auxiliary (dummy) variables
    int idum;

    error = 0;
    int phase;
    if ( SKIP_SYMBOLIC_FACTORIZATION )
        phase = 23;
    else
        phase = 13;

    pardiso(_pt, &_maxfct, &_mnum, &_mtype, &phase, &_n, a, ia, ja, &idum, &_nrhs, _iparm, &_msglvl, b, x, &error);
    
    if ( error != 0 )
    {
        // In case of error, matrix is reanalyzed and one more attempt at solution is made
        // A. Deallocate memory
        error = 0;
        phase = -1;
        pardiso(_pt, &_maxfct, &_mnum, &_mtype, &phase, &_n, a, ia, ja, &idum, &_nrhs, _iparm, &_msglvl, b, x, &error);
        this->giveErrorMessage(error);
        
        // B. Symbolic factorization, numerical factorization and back-substitution
        error = 0;
        phase = 13;
        pardiso(_pt, &_maxfct, &_mnum, &_mtype, &phase, &_n, a, ia, ja, &idum, &_nrhs, _iparm, &_msglvl, b, x, &error);
        this->giveErrorMessage(error);
    }

    return u;
}

// Private methods
void MKL_Pardiso::giveErrorMessage( int error )
{
    std::string msg;
    if ( error != 0 )
    {        
        switch (error)
        {
            case -1:
                msg = "Input inconsistent."; break;
            case -2:
                msg = "Not enough memory."; break;
            case -3:
                msg = "Reordering problem."; break;
            case -4:
                msg = "Zero pivot, numerical factorization or iterative refinement problem."; break;
            case -5:
                msg = "Unclassified (internal) error."; break;
            case -6:
                msg = "Reordering failed (matrix types 11 and 13 only)."; break;
            case -7:
                msg = "Diagonal matrix is singular."; break;
            case -8:
                msg = "32-bit integer overflow problem."; break;
            case -9:
                msg = "Not enough memory for OOC."; break;
            case -10:
                msg = "Error opening OOC files."; break;
            case -11:
                msg = "Read/write error with OOC files."; break;
            case -12:
                msg = "Pardiso_64 called from 32-bit library."; break;
            default:
                msg = "Error code = " + std::to_string(error) + ".";
        }

        throw std::runtime_error("MKL Pardiso error: " + msg);
    }
}

#endif