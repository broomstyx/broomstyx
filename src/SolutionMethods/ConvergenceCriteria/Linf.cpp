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

#include "LInf.hpp"
#include <omp.h>
#include <chrono>
#include <cmath>
#include "Core/AnalysisModel.hpp"
#include "Core/Diagnostics.hpp"
#include "Core/DofManager.hpp"
#include "Core/ObjectFactory.hpp"
#include "Util/readOperations.hpp"

using namespace broomstyx;

registerBroomstyxObject(ConvergenceCriterion, LInf)

LInf::LInf()
{
    _name = "LInf (Convergence criterion)";
}

LInf::~LInf() {}

// Public methods
// ----------------------------------------------------------------------------------------
bool LInf::checkConvergenceOf( const RealVector& resid, const std::vector<Dof*>& dof )
{
    std::chrono::time_point<std::chrono::system_clock> tic, toc;
    std::chrono::duration<double> tictoc;
    tic = std::chrono::high_resolution_clock::now();

    _threadCorrCrit.init(_nThreads);
    _threadCorrNorm.init(_nThreads);
    _threadResidNorm.init(_nThreads);

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
            // int grpIdx = this->giveIndexForDofGroup(grpNum);
            
            if ( grpNum == _dofGrpNum )
            {
                // Find equation number of DOF
                int eqNo = analysisModel().dofManager().giveEquationNumberAt(dof[i]);

                // Contribution to L2-norm of corrections
                double corrVal = std::fabs(analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof[i], correction));
                if ( corrVal > _threadCorrNorm(threadNum) )
                    _threadCorrNorm(threadNum) = corrVal;

                // Contribution to L2-norm of incremental solution
                double incVal = std::fabs(analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof[i], incremental_value));
                if ( incVal > _threadCorrCrit(threadNum) )
                    _threadCorrCrit(threadNum) = incVal;
                
                // Contribution to L2-norm of residual
                double rVal = std::fabs(resid(eqNo));
                if ( rVal > _threadResidNorm(threadNum) )
                    _threadResidNorm(threadNum) = rVal;
            }
        }
    }

    // Accumulate values from different threads
    _corrNorm = 0.;
    _corrCrit = 0.;
    _residNorm = 0.;
    _residCrit = 0.;

    for ( int i = 0; i < _nThreads; i++ )
    {
        if ( _corrNorm < _threadCorrNorm(i) )
            _corrNorm = _threadCorrNorm(i);
        if ( _corrCrit < _threadCorrCrit(i) )
            _corrCrit = _threadCorrCrit(i);
        if ( _residNorm < _threadResidNorm(i) )
            _residNorm = _threadResidNorm(i);
        if ( _residCrit < _threadResidCrit(i) )
            _residCrit = _threadResidCrit(i);
    }

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
RealMatrix LInf::giveConvergenceData()
{
    RealMatrix convData(2, 2);

        convData(0, 0) = _corrNorm;
        convData(0, 1) = _corrCrit;
        convData(1, 0) = _residNorm;
        convData(1, 1) = _residCrit;
    
    return convData;
}
// ----------------------------------------------------------------------------------------
void LInf::initialize( int dofGrpNum )
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
void LInf::processLocalResidualContribution( double contrib, int threadNum )
{
    double val = std::fabs(contrib);
    if ( val > _absTolRes && val > _threadResidCrit(threadNum) )
        _threadResidCrit(threadNum) = val;
}
// ----------------------------------------------------------------------------------------
void LInf::processLocalResidualContribution( const RealVector& contrib, const std::vector<int>& dofGrp, int threadNum )
{
    for ( int i = 0; i < contrib.dim(); i++ )
    {
        double val = std::fabs(contrib(i));
        if ( dofGrp[i] == _dofGrpNum && val > _absTolRes )
        {
            if ( val > _threadResidCrit(threadNum) )
                _threadResidCrit(threadNum) = val;
        }
    }
}
// ----------------------------------------------------------------------------------------
void LInf::readDataFromFile( FILE* fp )
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
void LInf::reportConvergenceStatus()
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
void LInf::resetResidualCriteria()
{
    _threadContribCount.init(_nThreads );
    _threadResidCrit.init(_nThreads );
}