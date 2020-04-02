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

#ifndef LINEARSTATIC_HPP
#define	LINEARSTATIC_HPP

#include "SolutionMethod.hpp"
#include "SparseMatrix/SparseMatrix.hpp"

namespace broomstyx
{
    class LinearSolver;
    
    class LinearStatic : public SolutionMethod
    {
    public:
        LinearStatic();
        virtual ~LinearStatic();

        int  computeSolutionFor( int stage
                               , const std::vector<BoundaryCondition>& bndCond
                               , const std::vector<FieldCondition>& fldCond
                               , const TimeData& time ) override;
        
        void initializeSolvers() override;
        virtual void formSparsityProfileForStage( int stage ) override;
        void readDataFromFile( FILE* fp ) override;

    protected:
        LinearSolver* _solver;
        SparseMatrix* _spMatrix;
        
        virtual void assembleEquations( int stage
                                      , const std::vector<BoundaryCondition>& bndCond
                                      , const std::vector<FieldCondition>& fldCond
                                      , const TimeData& time
                                      , RealVector& rhs );

        void assembleLeftHandSide( int stage, const TimeData& time );
        
        virtual void assembleRightHandSide( int stage
                                          , const std::vector<BoundaryCondition>& bndCond
                                          , const std::vector<FieldCondition>& fldCond
                                          , const TimeData& time
                                          , RealVector& rhs );
        
        void debugPrintSystem( RealVector& rhs );
    };    
}

#endif	/* LINEARSTATIC_HPP */