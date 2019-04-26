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

#include "SolutionMethod.hpp"
#include "../Core/AnalysisModel.hpp"
#include "../Core/DomainManager.hpp"
#include "../Core/NumericsManager.hpp"
#include "../Core/SolutionManager.hpp"
#include "../Numerics/Numerics.hpp"

using namespace broomstyx;

SolutionMethod::SolutionMethod() {}

SolutionMethod::~SolutionMethod() {}
// ---------------------------------------------------------------------------
void SolutionMethod::getCurrentLoadStep()
{
    _loadStep = analysisModel().solutionManager().giveCurrentLoadStep();
}
// ---------------------------------------------------------------------------
void SolutionMethod::imposeConstraintsAt( int stage
                                        , const std::vector<BoundaryCondition>& bndCond
                                        , const TimeData& time )
{
    // Loop through all boundary conditions
    for ( int ibc = 0; ibc < (int)bndCond.size(); ibc++ )
    {
        int boundaryId = analysisModel().domainManager().givePhysicalEntityNumberFor(bndCond[ibc].boundaryName());
        Numerics* numerics = analysisModel().numericsManager().giveNumerics(bndCond[ibc].targetNumerics());
        
        int nBCells = analysisModel().domainManager().giveNumberOfBoundaryCells();
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for ( int iCell = 0; iCell < nBCells; iCell++ )
        {
            Cell* curCell = analysisModel().domainManager().giveBoundaryCell(iCell);
            int label = analysisModel().domainManager().giveLabelOf(curCell);

            if ( label == boundaryId )
            {
                // Specifics of constraint imposition are handled by numerics
                numerics->imposeConstraintAt(curCell, stage, bndCond[ibc],time);
            }
        }

        // Some boundary conditions are actually internal conditions, so we need to loop over domain cells too
        int nDCells = analysisModel().domainManager().giveNumberOfDomainCells();

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for ( int iCell = 0; iCell < nDCells; iCell++ )
        {
            Cell* curCell = analysisModel().domainManager().giveDomainCell(iCell);
            int label = analysisModel().domainManager().giveLabelOf(curCell);

            if ( label == boundaryId )
            {
                // Specifics of constraint imposition are handled by numerics
                numerics->imposeConstraintAt(curCell, stage, bndCond[ibc],time);
            }
        }
    }    
}
// ---------------------------------------------------------------------------
bool broomstyx::SolutionMethod::checkConvergenceOfNumericsAt( int stage )
{
    int nCells = analysisModel().domainManager().giveNumberOfDomainCells();
    int nUnconvergedCells = 0;
    
#ifdef _OPENMP
#pragma omp parallel for reduction(+:nUnconvergedCells)    
#endif
    for (int iCell = 0; iCell < nCells; iCell++)
    {
        Cell* curCell = analysisModel().domainManager().giveDomainCell(iCell);
        int label = analysisModel().domainManager().giveLabelOf(curCell);
        Numerics* numerics = analysisModel().domainManager().giveNumericsForDomain(label);

        bool cellConvergence = numerics->performAdditionalConvergenceCheckAt(curCell, stage);
        if ( cellConvergence == false )
            nUnconvergedCells += 1;
    }
    
    if ( nUnconvergedCells > 0 )
        return false;
    else
        return true;
}