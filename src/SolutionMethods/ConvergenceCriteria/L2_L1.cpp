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

#include "L2_L1.hpp"
#include <omp.h>
#include <chrono>
#include <cmath>
#include "Core/AnalysisModel.hpp"
#include "Core/Diagnostics.hpp"
#include "Core/DofManager.hpp"
#include "Core/ObjectFactory.hpp"
#include "Util/readOperations.hpp"

using namespace broomstyx;

registerBroomstyxObject(ConvergenceCriterion, L2_L1)

L2_L1::L2_L1()
{
    _name = "L2_L1 (Convergence criterion)";
}

L2_L1::~L2_L1() {}

// Public methods
// ----------------------------------------------------------------------------------------
bool L2_L1::checkConvergenceOf( const RealVector& resid, const std::vector<Dof*>& dof )
{
    std::chrono::time_point<std::chrono::system_clock> tic, toc;
    std::chrono::duration<double> tictoc;
    tic = std::chrono::high_resolution_clock::now();

    _threadCorrCrit.init(_nThreads);
    _threadCorrNorm.init(_nThreads);
    _threadResidNorm.init(_nThreads);
    _threadDofCount.init(_nThreads);

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        int threadNum = 0;
#ifdef _OPENMP
        threadNum = omp_get_thread_num();
#endif

#ifdef _OPENMP
#pragma omp for
#endif
        for ( int i = 0; i < (int)dof.size(); i++ )
        {
            // Find group number of DOF
            int grpNum = analysisModel().dofManager().giveGroupNumberFor(dof[i]);
            
            if ( grpNum == _dofGrpNum )
            {
                // Find equation number of DOF
                int eqNo = analysisModel().dofManager().giveEquationNumberAt(dof[i]);

                // Contribution to L2-norm of corrections
                double corrVal = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof[i], correction);
                _threadCorrNorm(threadNum) += corrVal*corrVal;

                // Contribution to L2-norm of incremental solution
                double incVal = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof[i], incremental_value);
                _threadCorrCrit(threadNum) += incVal*incVal;
                
                // Contribution to L2-norm of residual
                double rVal = resid(eqNo);
                _threadResidNorm(threadNum) += rVal*rVal;

                // Increase DOF count
                _threadDofCount(threadNum) += 1.;
            }
        }
    }

    // Accumulate values from different threads
    _corrNorm = 0.;
    _corrCrit = 0.;
    _residNorm = 0.;
    _residCrit = 0.;
    _dofCount = 0.;
    _contribCount = 0.;

    for ( int i = 0; i < _nThreads; i++ )
    {
        _corrNorm += _threadCorrNorm(i);
        _corrCrit += _threadCorrCrit(i);
        _residNorm += _threadResidNorm(i);
        _residCrit += _threadResidCrit(i);
        _dofCount += _threadDofCount(i);
        _contribCount += _threadContribCount(i);
    }
    if ( _contribCount < 1.0 )
        _contribCount = 1.0;

    // Normalize values
    _corrNorm = std::sqrt(_corrNorm/_dofCount);
    _corrCrit = std::sqrt(_corrCrit/_dofCount);
    _residNorm = std::sqrt(_residNorm/_dofCount);
    _residCrit /= _contribCount;

    // Apply relative tolerances
    _corrCrit *= _relTolCorr;
    _residCrit *= _relTolRes;

    // Check against absolute tolerances
    if ( _corrCrit < _absTolCorr )
        _corrCrit = _absTolCorr;
    if ( _residCrit < _absTolRes )
        _residCrit = _absTolRes;
    
    bool convergenceStatus = true;
    if ( _corrNorm > _corrCrit )
        convergenceStatus = false;
    if ( _residNorm > _residCrit )
        convergenceStatus = false;
    
    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    diagnostics().addConvergenceCheckTime(tictoc.count());

    return convergenceStatus;
}
// ----------------------------------------------------------------------------------------
RealMatrix L2_L1::giveConvergenceData()
{
    RealMatrix convData(2, 2);

    convData(0, 0) = _corrNorm;
    convData(0, 1) = _corrCrit;
    convData(1, 0) = _residNorm;
    convData(1, 1) = _residCrit;
    
    return convData;
}
// ----------------------------------------------------------------------------------------
void L2_L1::initialize( int dofGrpNum )
{
    _dofGrpNum = dofGrpNum;

    #ifdef _OPENMP
    #pragma omp parallel
    {
        #pragma omp single
        {
            _nThreads = omp_get_num_threads();
        }
    }
    #else
    _nThreads = 1;
    #endif
    
    _relTolCorr = 0.;
    _relTolRes = 0.;
    _absTolCorr = 0.;
    _absTolRes = 0.;

    _corrNorm = 0.;
    _corrCrit = 0.;
    _residNorm = 0.;
    _residCrit = 0.;
}
// ----------------------------------------------------------------------------------------
void L2_L1::processLocalResidualContribution( double contrib, int threadNum )
{
    double val = std::fabs(contrib);
    if ( val > _absTolRes )
    {
        _threadResidCrit(threadNum) += val;
        _threadContribCount(threadNum) += 1.;
    }
}
// ----------------------------------------------------------------------------------------
void L2_L1::processLocalResidualContribution( const RealVector& contrib, const std::vector<int>& dofGrp, int threadNum )
{
    for ( int i = 0; i < contrib.dim(); i++ )
    {
        double val = std::fabs(contrib(i));
        if ( dofGrp[i] == _dofGrpNum && val > _absTolRes )
        {
            _threadResidCrit(threadNum) += val;
            _threadContribCount(threadNum) += 1.;
        }
    }
}
// ----------------------------------------------------------------------------------------
void L2_L1::readDataFromFile( FILE* fp )
{
    // Read DOF group convergence parameters
    _trackingOption = getStringInputFrom(fp, "Failed to read criterion option from input file!", _name);
    
    if ( _trackingOption != "C" && _trackingOption != "R" && _trackingOption != "CR" )
        throw std::runtime_error("ERROR: Invalid tracking option encountered in input file! Valid: (C/R/CR)");

    if ( _trackingOption == "C" || _trackingOption == "CR" )
    {
        // Correction tolerance
        _relTolCorr = getRealInputFrom(fp, "Failed to read correction tolerance from input file!", _name);

        // Residual tolerance
        _relTolRes = getRealInputFrom(fp, "Failed to read residual tolerance from input file!", _name);
    }
    if ( _trackingOption == "R" || _trackingOption == "CR" )
    {
        // Absolute tolerance for corrections
        _absTolCorr = getRealInputFrom(fp, "Failed to read absolute correction tolerance from input file!", _name);

        // Absolute tolerance for residuals
        _absTolRes = getRealInputFrom(fp, "Failed to read absolute residual tolerance from input file!", _name);
    }
}
// ----------------------------------------------------------------------------------------
void L2_L1::reportConvergenceStatus()
{
    // Convergence of corrections
    if ( _trackingOption == "C" || _trackingOption == "CR" )
    {
        std::printf("\n     C  %-6d%e   %e ", _dofGrpNum, _corrNorm, _corrCrit);
        if ( _corrNorm <= _corrCrit )
            std::printf(" << CONVERGED");
    }
        
    // Convergence of residuals
    if ( _trackingOption == "R" || _trackingOption == "CR" )
    {
        std::printf("\n     R  %-6d%e   %e ", _dofGrpNum, _residNorm, _residCrit);
        if ( _residNorm <= _residCrit )
            std::printf(" << CONVERGED");
    }
}
// ----------------------------------------------------------------------------------------
void L2_L1::resetResidualCriteria()
{
    _threadContribCount.init(_nThreads);
    _threadResidCrit.init(_nThreads);
}