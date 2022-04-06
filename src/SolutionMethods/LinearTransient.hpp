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

#ifndef LINEARTRANSIENTBE_HPP
#define	LINEARTRANSIENTBE_HPP

#include "LinearStatic.hpp"
#include "SparseMatrix/SparseMatrix.hpp"

namespace broomstyx
{
    class LinearTransient final : public LinearStatic
    {
    public:
        LinearTransient();
        virtual ~LinearTransient();

        void formSparsityProfileForStage( int stage ) override;
        
    private:
        void assembleEquations( int stage
                              , const std::vector<BoundaryCondition>& bndCond
                              , const std::vector<FieldCondition>& fldCond
                              , const TimeData& time
                              , RealVector& rhs ) override;
        
        void assembleLeftHandSide( int stage, const TimeData& time );

        void assembleRightHandSide( int stage
                                  , const std::vector<BoundaryCondition>& bndCond
                                  , const std::vector<FieldCondition>& fldCond
                                  , const TimeData& time
                                  , RealVector& rhs ) override;
    };    
}

#endif	/* LINEARTRANSIENTBE_HPP */
