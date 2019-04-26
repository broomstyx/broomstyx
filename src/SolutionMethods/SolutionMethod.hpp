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

#ifndef SOLUTIONMETHOD_HPP
#define	SOLUTIONMETHOD_HPP

#include <cstdio>
#include <vector>
#include "../Core/BoundaryCondition.hpp"
#include "../Core/FieldCondition.hpp"
#include "../Core/TimeData.hpp"
#include "../Util/RealVector.hpp"

namespace broomstyx
{  
    class LoadStep;

    class SolutionMethod
    {
    public:
        SolutionMethod();
        virtual ~SolutionMethod();

        void getCurrentLoadStep();
        
        void imposeConstraintsAt( int stage
                                , const std::vector<BoundaryCondition>& bndCond
                                , const TimeData& time );
        
        virtual int  computeSolutionFor( int stage
                                       , const std::vector<BoundaryCondition>& bndCond
                                       , const std::vector<FieldCondition>& fldCond
                                       , const TimeData& time ) = 0;
        
        virtual void formSparsityProfileForStage ( int stage ) = 0;
        virtual void initializeSolvers() = 0;
        virtual void readDataFromFile( FILE* fp ) = 0;

    protected:
        LoadStep* _loadStep;
        std::string _name;
        
        bool checkConvergenceOfNumericsAt( int stage );
    };
}

#endif	/* SOLUTIONMETHOD_HPP */