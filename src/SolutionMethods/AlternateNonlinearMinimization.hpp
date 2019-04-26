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

#ifndef ALTERNATENONLINEARMINIMIZATION_HPP
#define	ALTERNATENONLINEARMINIMIZATION_HPP

#include "AlternateMinimization.hpp"
#include <map>
#include <vector>

namespace broomstyx
{
    class AlternateNonlinearMinimization : public AlternateMinimization
    {        
    public:
        AlternateNonlinearMinimization();
        virtual ~AlternateNonlinearMinimization();

        int  computeSolutionFor( int stage
                               , const std::vector<BoundaryCondition>& bndCond
                               , const std::vector<FieldCondition>& fldCond
                               , const TimeData& time ) override;
        
        virtual void readDataFromFile( FILE* fp ) override;

    protected:
        std::vector<int> _maxSubsysIter;
        RealMatrix _tolFactor;
        
        std::vector<bool> _enableLineSearch;
        RealVector _slackTolerance;
        RealVector _initMultiplier;
        RealVector _minMultiplier;
        RealVector _maxMultiplier;
        std::vector<int> _maxLineSearchIter;
        
        bool checkSubsystemConvergenceAt( int stage
                                        , int subsys
                                        , const RealVector& resid
                                        , int globalIterCount
                                        , int localIterCount );
    };    
}

#endif	/* ALTERNATENONLINEARMINIMIZATION_HPP */