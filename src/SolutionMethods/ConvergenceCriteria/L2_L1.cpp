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
bool L2_L1::checkConvergenceOf( const std::vector<RealVector>& resid
                              , const std::vector<int>& subsysNumbers
                              , const std::vector<Dof*>& dof )
{
    std::chrono::time_point<std::chrono::system_clock> tic, toc;
    std::chrono::duration<double> tictoc;
    tic = std::chrono::high_resolution_clock::now();

    _corrCritPerThread.init(_nThreads, _nDofGroups);
    _corrNormPerThread.init(_nThreads, _nDofGroups);
    _dofGrpCountPerThread.init(_nThreads, _nDofGroups);
    _residNormPerThread.init(_nThreads, _nDofGroups);

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
            int grpIdx = this->giveIndexForDofGroup(grpNum);
            
            // Find subsystem of DOF
            int ssNum = analysisModel().dofManager().giveSubsystemNumberFor(dof[i]);
            int ssIdx = this->giveIndexForSubsystem(ssNum, subsysNumbers);
            
            if ( grpIdx >= 0 && ssIdx >= 0 )
            {
                // Find equation number of DOF
                int eqNo = analysisModel().dofManager().giveEquationNumberAt(dof[i]);

                // Contribution to L2-norm of corrections
                double corrVal = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof[i], correction);
                _corrNormPerThread(threadNum, grpIdx) += corrVal*corrVal;

                // Contribution to L2-norm of incremental solution
                double incVal = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof[i], incremental_value);
                _corrCritPerThread(threadNum, grpIdx) += incVal*incVal;
                
                // Contribution to L2-norm of residual
                double rVal = resid[ssIdx](eqNo);
                _residNormPerThread(threadNum, grpIdx) += rVal*rVal;

                // Increase DOF group count
                _dofGrpCountPerThread(threadNum, grpIdx) += 1.;
            }
        }
    }

    // Accumulate values from different threads
    _corrNorm.init(_nDofGroups);
    _corrCrit.init(_nDofGroups);
    _residNorm.init(_nDofGroups);
    _residCrit.init(_nDofGroups);
    _dofGrpCount.init(_nDofGroups);
    _contribCount.init(_nDofGroups);

    for ( int i = 0; i < _nThreads; i++ )
        for ( int j = 0; j < _nDofGroups; j++ )
        {
            _corrNorm(j) += _corrNormPerThread(i,j);
            _corrCrit(j) += _corrCritPerThread(i,j);
            _residNorm(j) += _residNormPerThread(i,j);
            _residCrit(j) += _residCritPerThread(i,j);
            _dofGrpCount(j) += _dofGrpCountPerThread(i,j);
            _contribCount(j) += _contribCountPerThread(i,j);
        }

    for ( int i = 0; i < _nDofGroups; i++ )
    {
        // Normalize values
        _corrNorm(i) = std::sqrt(_corrNorm(i)/_dofGrpCount(i));
        _corrCrit(i) = std::sqrt(_corrCrit(i)/_dofGrpCount(i));
        _residNorm(i) = std::sqrt(_residNorm(i)/_dofGrpCount(i));
        _residCrit(i) /= _contribCount(i);

        // Apply relative tolerances
        _corrCrit(i) *= _relTolCor(i);
        _residCrit(i) *= _relTolRes(i);

        // Check against absolute tolerances
        if ( _corrCrit(i) < _absTolCor(i) )
            _corrCrit(i) = _absTolCor(i);
        if ( _residCrit(i) < _absTolRes(i) )
            _residCrit(i) = _absTolRes(i);
    }
    
    bool convergenceStatus = true;
    for ( int i = 0; i < _nDofGroups; i++ )
    {
        if ( _corrNorm(i) > _corrCrit(i) )
            convergenceStatus = false;
        if ( _residNorm(i) > _residCrit(i) )
            convergenceStatus = false;
    }
    
    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    diagnostics().addConvergenceCheckTime(tictoc.count());

    return convergenceStatus;
}
// ----------------------------------------------------------------------------------------
RealMatrix L2_L1::giveConvergenceData()
{
    RealMatrix convData(2*_nDofGroups, 2);

    for ( int i = 0; i < _nDofGroups; i++ )
    {
        convData(2*i, 0) = _corrNorm(i);
        convData(2*i, 1) = _corrCrit(i);
        convData(2*i+1, 0) = _residNorm(i);
        convData(2*i+1, 1) = _residCrit(i);
    }

    return convData;
}
// ----------------------------------------------------------------------------------------
void L2_L1::initialize( int nDofGroups )
{
    _nDofGroups = nDofGroups;
    
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
    
    _relTolCor.init(_nDofGroups);
    _relTolRes.init(_nDofGroups);
    _absTolCor.init(_nDofGroups);
    _absTolRes.init(_nDofGroups);

    _corrNorm.init(_nDofGroups);
    _corrCrit.init(_nDofGroups);
    _residNorm.init(_nDofGroups);
    _residCrit.init(_nDofGroups);
}
// ----------------------------------------------------------------------------------------
void L2_L1::processLocalResidualContribution( RealVector& contrib, std::vector<int>& dofGrp, int threadNum )
{
    for ( int i = 0; i < contrib.dim(); i++ )
    {
        int dgIdx = this->giveIndexForDofGroup(dofGrp[i]);
        if ( dgIdx >= 0 && std::fabs(contrib(i)) > _absTolRes(dgIdx) )
        {
            _residCritPerThread( threadNum, dgIdx ) += std::fabs(contrib(i));
            _contribCountPerThread( threadNum, dgIdx) += 1.;
        }
    }
}
// ----------------------------------------------------------------------------------------
void L2_L1::readDataFromFile( FILE* fp )
{
    _dofGrpNum.assign(_nDofGroups, -1);
    for ( int i = 0; i < _nDofGroups; i++ )
    {
        // Read DOF group number
        _dofGrpNum[i] = getIntegerInputFrom(fp, "Failed to read DOF group number from input file!", _name);

        // Read DOF group convergence parameters
        verifyKeyword(fp, "Parameters", _name);

        // Correction tolerance
        _relTolCor(i) = getRealInputFrom(fp, "Failed to read correction tolerance from input file!", _name);

        // Residual tolerance
        _relTolRes(i) = getRealInputFrom(fp, "Failed to read residual tolerance from input file!", _name);

        // Absolute tolerance for corrections
        _absTolCor(i) = getRealInputFrom(fp, "Failed to read absolute correction tolerance from input file!", _name);

        // Absolute tolerance for residuals
        _absTolRes(i) = getRealInputFrom(fp, "Failed to read absolute residual tolerance from input file!", _name);
    }
}
// ----------------------------------------------------------------------------------------
void L2_L1::reportConvergenceStatus()
{
    // Report status
    std::printf("\n\n    DOF Grp   L2-Norm         Criterion");
    std::printf("\n   -------------------------------------------------------");
    
    for ( int i = 0; i < _nDofGroups; i++ )
    {
        // Convergence of corrections
        std::printf("\n     C  %-6d%e   %e ", _dofGrpNum[i], _corrNorm(i), _corrCrit(i));
        if ( _corrNorm(i) <= _corrCrit(i) )
            std::printf(" << CONVERGED");
        
        // Convergence of residuals
        std::printf("\n     R  %-6d%e   %e ", _dofGrpNum[i], _residNorm(i), _residCrit(i));
        if ( _residNorm(i) <= _residCrit(i) )
            std::printf(" << CONVERGED");
    }
}
// ----------------------------------------------------------------------------------------
void L2_L1::resetResidualCriteria()
{
    _contribCountPerThread.init(_nThreads, _nDofGroups);
    _residCritPerThread.init(_nThreads, _nDofGroups);
}