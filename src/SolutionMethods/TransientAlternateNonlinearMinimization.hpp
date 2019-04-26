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

#ifndef TRANSIENTALTERNATENONLINEARMINIMIZATION_HPP
#define TRANSIENTALTERNATENONLINEARMINIMIZATION_HPP

#include "AlternateNonlinearMinimization.hpp"
#include <map>
#include <vector>

namespace broomstyx
{
    class TransientAlternateNonlinearMinimization final : public AlternateNonlinearMinimization
    {
    public:
        TransientAlternateNonlinearMinimization();
        virtual ~TransientAlternateNonlinearMinimization();

        void formSparsityProfileForStage( int stage ) override;
        void readDataFromFile( FILE* fp ) override;
        
    private:
        RealVector assembleLeftHandSide( int stage
                                       , int subsys
                                       , const TimeData& time ) override;
        
        void assembleJacobian( int stage
                             , int subsys
                             , const TimeData& time ) override;
        
        RealVector assembleRightHandSide( int stage
                                        , int subsys
                                        , const std::vector<BoundaryCondition>& bndCond
                                        , const std::vector<FieldCondition>& fldCond
                                        , const TimeData& time ) override;
        
        std::vector<bool> _isTransient;
    };
}

#endif /* TRANSIENTALTERNATENONLINEARMINIMIZATION_HPP */