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

#include "LInf_LInf_StoppingCriterion.hpp"
#include <omp.h>
#include <chrono>
#include <cmath>
#include "Core/AnalysisModel.hpp"
#include "Core/Diagnostics.hpp"
#include "Core/DofManager.hpp"

using namespace broomstyx;

LInf_LInf_StoppingCriterion::LInf_LInf_StoppingCriterion()
{
    _name = "L2_L1_StoppingCriterion";
}

LInf_LInf_StoppingCriterion::~LInf_LInf_StoppingCriterion() {}

// Public methods
// ----------------------------------------------------------------------------------------
bool LInf_LInf_StoppingCriterion::checkConvergenceOf( const std::vector<RealVector>& resid
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
                double corrVal = std::fabs(analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof[i], correction));
                if ( corrVal > _corrNormPerThread(threadNum, grpIdx) )
                    _corrNormPerThread(threadNum, grpIdx) = corrVal;

                // Contribution to L2-norm of incremental solution
                double incVal = std::fabs(analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof[i], incremental_value));
                if ( incVal > _corrCritPerThread(threadNum, grpIdx) )
                    _corrCritPerThread(threadNum, grpIdx) = incVal;
                
                // Contribution to L2-norm of residual
                double rVal = std::fabs(resid[ssIdx](eqNo));
                if ( rVal > _residNormPerThread(threadNum, grpIdx) )
                    _residNormPerThread(threadNum, grpIdx) = rVal;
            }
        }
    }

    // Accumulate values from different threads
    _corrNorm.init(_nDofGroups);
    _corrCrit.init(_nDofGroups);
    _residNorm.init(_nDofGroups);
    _residCrit.init(_nDofGroups);

    for ( int i = 0; i < _nThreads; i++ )
        for ( int j = 0; j < _nDofGroups; j++ )
        {
            if ( _corrNorm(j) < _corrNormPerThread(i,j) )
                _corrNorm(j) = _corrNormPerThread(i,j);
            if ( _corrCrit(j) < _corrCritPerThread(i,j) )
                _corrCrit(j) = _corrCritPerThread(i,j);
            if ( _residNorm(j) < _residNormPerThread(i,j) )
                _residNorm(j) = _residNormPerThread(i,j);
            if ( _residCrit(j) < _residCritPerThread(i,j) )
                _residCrit(j) = _residCritPerThread(i,j);
        }

    for ( int i = 0; i < _nDofGroups; i++ )
    {
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
void LInf_LInf_StoppingCriterion::initialize( int nDofGroups )
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
void LInf_LInf_StoppingCriterion::processLocalResidualContribution( RealVector& contrib, std::vector<int>& dofGrp, int threadNum )
{
    for ( int i = 0; i < contrib.dim(); i++ )
    {
        int dgIdx = this->giveIndexForDofGroup(dofGrp[i]);
        double val = std::fabs(contrib(i));
        if ( dgIdx >= 0 && val > _absTolRes(dgIdx) )
        {
            if ( val > _residCritPerThread( threadNum, dgIdx) )
                _residCritPerThread( threadNum, dgIdx ) = std::fabs(contrib(i));
        }
    }
}
// ----------------------------------------------------------------------------------------
void LInf_LInf_StoppingCriterion::reportConvergenceStatus()
{
    // Report status
    std::printf("\n\n    DOF Grp   Inf-Norm         Criterion");
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
void LInf_LInf_StoppingCriterion::resetResidualCriteria()
{
    _contribCountPerThread.init(_nThreads, _nDofGroups);
    _residCritPerThread.init(_nThreads, _nDofGroups);
}