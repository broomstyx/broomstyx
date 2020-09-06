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

#ifndef CONVERGENCECRITERION_HPP
#define	CONVERGENCECRITERION_HPP

#include <string>
#include <vector>
#include "Util/RealMatrix.hpp"
#include "Util/RealVector.hpp"

namespace broomstyx
{
    class Dof;
    
    class ConvergenceCriterion
    {
    public:
        ConvergenceCriterion() {}
        virtual ~ConvergenceCriterion() {}

        virtual bool checkConvergenceOf( const RealVector& resid, const std::vector<Dof*>& dof ) = 0;
        virtual RealMatrix giveConvergenceData() = 0;
        virtual void initialize( int dofGrpNum ) = 0;
        virtual void processLocalResidualContribution( const RealVector& contrib, const std::vector<int>& dofGrp, int threadNum ) = 0;
        virtual void readDataFromFile( FILE* fp ) = 0;
        virtual void reportConvergenceStatus() = 0;
        virtual void resetResidualCriteria() = 0;
        
    protected:
        std::string _name;
        int _dofGrpNum;
    };
}

#endif	/* CONVERGENCECRITERION_HPP */