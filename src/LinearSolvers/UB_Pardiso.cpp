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

#include "UB_Pardiso.hpp"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <omp.h>
#include <stdexcept>
#include <string>
#include <tuple>

#include "../Core/ObjectFactory.hpp"
#include "../SparseMatrix/SparseMatrix.hpp"
#include "../Util/RealVector.hpp"
#include "../Util/readOperations.hpp"

using namespace broomstyx;

registerBroomstyxObject(LinearSolver, UB_Pardiso)

// Constructor
UB_Pardiso::UB_Pardiso() : LinearSolver()
{
    _solver = 0;
    _nThreads = 0;
    _memoryIsAllocated = false;
    _mtype = 0;
    _symmetry = false;
    
    for ( int i = 0; i < 64; i++ )
    {
        _pt[i] = nullptr;
        _iparm[i] = 0;
    }
}

// Destructor
UB_Pardiso::~UB_Pardiso()
{
    this->clearInternalMemory();
}

// Public Methods
// ----------------------------------------------------------------------------
void UB_Pardiso::allocateInternalMemoryFor( SparseMatrix* coefMat )
{
    if ( !_memoryIsAllocated )
    {
        // Pardiso control parameters
        double dparm[64];
        int    maxfct, mnum, phase, error, msglvl;

        int dim1, dim2;
        std::tie(dim1,dim2) = coefMat->giveMatrixDimensions();

        // Check that coefficient is square
        if ( dim1 != dim2 )
        {
            std::printf("\n\nCannot solve specified linear system!\n\tCoefficient matrix is not square");
            std::printf("\n\t--> dim(A) = [ %d x %d ]\n", dim1, dim2);
            throw std::runtime_error("Source: Pardiso");
        }

        // Matrix Data
        int n = dim1;
        int* ia;
        int* ja;
        double* a;

        std::tie(ia,ja) = coefMat->giveProfileArrays();
        a  = coefMat->giveValArray();

        int nrhs = 1; // Number of right hand sides

        // Auxiliary variables
        double ddum;   // dummy variable
        int    idum;   // dummy variable

        maxfct = 1; /* Maximum number of numerical factorizations.  */
        mnum   = 1; /* Which factorization to use. */

        msglvl = 0; /* Print statistical information  */
        error  = 0; /* Initialize error flag */

    /* -------------------------------------------------------------------- *
     * Reordering and Symbolic Factorization. This step also allocates all memory
     * that is necessary for the factorization */
    /* -------------------------------------------------------------------- */

        phase = 11;
        pardiso (_pt, &maxfct, &mnum, &_mtype, &phase, &n, a, ia, ja, &idum, &nrhs, _iparm, &msglvl, &ddum, &ddum, &error, dparm);

        if ( error != 0 )
        {
            std::printf("\nPardiso error during symbolic factorization: %d", error);
            throw std::runtime_error("\n");
        }

        _memoryIsAllocated = true;
    }
}
// ----------------------------------------------------------------------------
RealVector UB_Pardiso::backSubstitute( SparseMatrix* coefMat, RealVector& rhs )
{
    // Pardiso control parameters
    double dparm[64];
    int    maxfct, mnum, phase, error, msglvl;

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

    // Initialize solution vector
    RealVector u(dim1);

    // Matrix Data
    int n = dim1;
    int* ia;
    int* ja;
    double* a;

    std::tie(ia,ja) = coefMat->giveProfileArrays();
    a  = coefMat->giveValArray();

    // RHS and solution vectors
    double *b, *x; //, *bs, res, res0;
    b = rhs.ptr();
    x = u.ptr();

    int nrhs = 1; // Number of right hand sides

    // Auxiliary variables
    int    idum;   // dummy variable

    maxfct = 1; /* Maximum number of numerical factorizations.  */
    mnum   = 1; /* Which factorization to use. */

    msglvl = 0; /* Print statistical information  */
    error  = 0; /* Initialize error flag */

/* -------------------------------------------------------------------- */
/* .. Back substitution and iterative refinement. */
/* -------------------------------------------------------------------- */
    phase = 33;
    pardiso (_pt, &maxfct, &mnum, &_mtype, &phase, &n, a, ia, ja, &idum, &nrhs, _iparm, &msglvl, b, x, &error, dparm);
    
    maxfct = 1; /* Maximum number of numerical factorizations.  */
    mnum   = 1; /* Which factorization to use. */

    msglvl = 0; /* Print statistical information  */
    error  = 0; /* Initialize error flag */

    if (error != 0)
    {
        std::printf("\nPardiso ERROR during solution: %d", error);
        throw std::runtime_error("\n");
    }

    return u;
}
// ----------------------------------------------------------------------------
void UB_Pardiso::clearInternalMemory()
{
    if ( _memoryIsAllocated )
    {
        // Release internal memory for Pardiso

        // Dummy pardiso control parameters
        double dparm[64] = {0.}, ddum = 0., a[1] = {0.};
        int    maxfct = 0, mnum = 1, error = 0, msglvl = 0, n = 0, idum = 0, nrhs = 0, ia[1] = {0}, ja[1] = {0};

        int phase = -1;
        pardiso (_pt, &maxfct, &mnum, &_mtype, &phase, &n, a, ia, ja, &idum, &nrhs, _iparm, &msglvl, &ddum, &ddum, &error, dparm);

        _memoryIsAllocated = false;
    }
}
// ----------------------------------------------------------------------------
void UB_Pardiso::factorize( SparseMatrix* coefMat )
{
    int    mtype;

    // Pardiso control parameters
    double dparm[64];
    int    maxfct, mnum, phase, error, msglvl;

    int dim1, dim2;
    std::tie(dim1,dim2) = coefMat->giveMatrixDimensions();

    // Matrix Data
    int n = dim1;
    int* ia;
    int* ja;
    double* a;

    std::tie(ia,ja) = coefMat->giveProfileArrays();
    a  = coefMat->giveValArray();

    int nrhs = 1; // Number of right hand sides

    // Auxiliary variables
    double ddum;   // dummy variable
    int    idum;   // dummy variable

    maxfct = 1; /* Maximum number of numerical factorizations.  */
    mnum   = 1; /* Which factorization to use. */

    msglvl = 0; /* Print statistical information  */
    error  = 0; /* Initialize error flag */

/* -------------------------------------------------------------------- */
/* .. Numerical factorization. */
/* -------------------------------------------------------------------- */

    phase = 22;
    pardiso (_pt, &maxfct, &mnum, &_mtype, &phase, &n, a, ia, ja, &idum, &nrhs, _iparm, &msglvl, &ddum, &ddum, &error, dparm);

    if (error != 0)
    {
        // Try again using full symbolic factorization
        phase = -1;
        pardiso (_pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, _iparm, &msglvl, &ddum, &ddum, &error, dparm);
        phase = 12;
        pardiso (_pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, _iparm, &msglvl, &ddum, &ddum, &error, dparm);

        if (error != 0)
        {
            std::printf("\nPardiso error during numerical factorization: %d",error);
            throw std::runtime_error("\n");
        }
    }
}
// ----------------------------------------------------------------------------
void UB_Pardiso::initialize()
{
    // Pardiso control parameters
    double dparm[64];
    int    error = 0;
    _solver = 0;
    _memoryIsAllocated = false;

    // Initialize internal solver memory pointer. This is only necessary for the FIRST call
    // of the PARDISO solver
    for (int i = 0; i < 64; i++)
        _pt[i] = 0;

    // Setup Pardiso control parameters
#ifdef _OPENMP
    _iparm[2] = _nThreads;
#else
    _iparm[2] = 1;
#endif

    _iparm[7] = 100;     // Max number of iterative refinement
    _iparm[23] = 1;      // Parallel factorization algorithm
    _iparm[27] = 1;      // Parallel METIS reordering

    pardisoinit(_pt, &_mtype, &_solver, _iparm, dparm, &error);
    if ( error != 0 )
    {
        if (error == -10 )
            throw std::runtime_error("Pardiso: No license file found!\n");
        if (error == -11 )
            throw std::runtime_error("Pardiso: License is expired\n");
        if (error == -12 )
            throw std::runtime_error("Pardiso: Wrong username or hostname\n");
    }
    std::printf("\n");
//    for ( int i = 0; i < 75; i++)
//        std::printf("*");
//    std::printf("\n\n");
}
// ----------------------------------------------------------------------------
void UB_Pardiso::readDataFrom( FILE* fp )
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
std::string UB_Pardiso::giveRequiredMatrixFormat()
{
    return std::string("CSR1");
}
// ----------------------------------------------------------------------------
bool UB_Pardiso::giveSymmetryOption()
{
    return _symmetry;
}
// ----------------------------------------------------------------------------
RealVector UB_Pardiso::solve( SparseMatrix* coefMat, RealVector&   rhs )
{
    // Pardiso control parameters
    double dparm[64];
    int    maxfct, mnum, phase, error, msglvl;

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

    // Initialize solution vector
    RealVector u(dim1);

    // Matrix Data
    int n = dim1;
    int* ia;
    int* ja;
    double* a;

    std::tie(ia,ja) = coefMat->giveProfileArrays();
    a  = coefMat->giveValArray();

    // RHS and solution vectors
    double *b, *x; //, *bs, res, res0;
    b = rhs.ptr();
    x = u.ptr();

    int nrhs = 1; // Number of right hand sides

    // Auxiliary variables
    double ddum;   // dummy variable
    int    idum;   // dummy variable

    maxfct = 1; /* Maximum number of numerical factorizations.  */
    mnum   = 1; /* Which factorization to use. */

    msglvl = 0; /* Print statistical information  */
    error  = 0; /* Initialize error flag */

/* -------------------------------------------------------------------- */
/* .. Numerical factorization. */
/* -------------------------------------------------------------------- */

    phase = 22;
    pardiso (_pt, &maxfct, &mnum, &_mtype, &phase, &n, a, ia, ja, &idum, &nrhs, _iparm, &msglvl, &ddum, &ddum, &error, dparm);

    if (error != 0)
    {
        // Try again using full symbolic factorization
        phase = -1;
        pardiso (_pt, &maxfct, &mnum, &_mtype, &phase, &n, a, ia, ja, &idum, &nrhs, _iparm, &msglvl, &ddum, &ddum, &error, dparm);
        phase = 12;
        pardiso (_pt, &maxfct, &mnum, &_mtype, &phase, &n, a, ia, ja, &idum, &nrhs, _iparm, &msglvl, &ddum, &ddum, &error, dparm);

        if (error != 0)
        {
            std::printf("\nPardiso error during numerical factorization: %d",error);
            throw std::runtime_error("\n");
        }
    }

/* -------------------------------------------------------------------- */
/* .. Back substitution and iterative refinement. */
/* -------------------------------------------------------------------- */
    phase = 33;
    pardiso (_pt, &maxfct, &mnum, &_mtype, &phase, &n, a, ia, ja, &idum, &nrhs, _iparm, &msglvl, b, x, &error, dparm);

    if (error != 0)
    {
        std::printf("\nPardiso ERROR during solution: %d", error);
        throw std::runtime_error("\n");
    }

    return u;
}

#endif