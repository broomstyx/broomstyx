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

#ifndef NEWTONRAPHSON_HPP
#define	NEWTONRAPHSON_HPP

#include "SolutionMethod.hpp"
#include <map>
#include <tuple>
#include <vector>

namespace broomstyx
{
    class ConvergenceCriterion;
    class Dof;
    class LinearSolver;
    class SparseMatrix;
    
    class NewtonRaphson : public SolutionMethod
    {
    public:
        NewtonRaphson();
        virtual ~NewtonRaphson();
        
        int  computeSolutionFor( int stage
                               , const std::vector<BoundaryCondition>& bndCond
                               , const std::vector<FieldCondition>& fldCond
                               , const TimeData& time ) override;
        
        virtual void formSparsityProfileForStage( int stage ) override;
        void         initializeSolvers() override;
        virtual void readDataFromFile( FILE* fp ) override;
        
    protected:
        int _nDofGroups;
        std::vector<int> _dofGrpNum;
        RealVector       _dofGrpCount;
        
        bool          _symmetry;
        LinearSolver* _solver;
        SparseMatrix* _spMatrix;
        int           _nUnknowns;
        double        _overRelaxation;
        
        std::vector<ConvergenceCriterion*> _convergenceCriterion;
        int _maxIter;
        int _substepCount;
        bool _abortAtMaxIter;
        
        virtual RealVector assembleLeftHandSide( int stage, const TimeData& time );
        virtual void       assembleJacobian( int stage, const TimeData& time );
        
        virtual RealVector assembleRightHandSide( int stage
                                                , const std::vector<BoundaryCondition>& bndCond
                                                , const std::vector<FieldCondition>&    fldCond
                                                , const TimeData& time );
        
        int giveIndexForDofGroup( int dofGroupNum );
    };
}

#endif	/* NEWTONRAPHSON_HPP */