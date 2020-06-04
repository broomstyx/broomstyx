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

#ifndef CONVERGENCECHECKER_HPP
#define	CONVERGENCECHECKER_HPP

#include <string>
#include <vector>
#include "Util/RealVector.hpp"

namespace broomstyx
{
    class Dof;
    
    class ConvergenceChecker
    {
    public:
        ConvergenceChecker();
        virtual ~ConvergenceChecker();

        void initialize( int nDofGroups );
        void readDataFromFile( FILE* fp );
        void reportConvergenceStatus();

        virtual bool checkSolutionConvergence() = 0;
        virtual void computeResidalNorms( const std::vector<RealVector>& resid, const std::vector<Dof*>& dof ) = 0;
        virtual void processResidualContribution( std::vector<Dof*> dof, RealVector contrib ) = 0;
        
    protected:
        std::string _name;
        
        int _nDofGroups;
        std::vector<int> _dofGrpNum;
        RealVector _dofGrpCount;
        
        RealVector _relTolCor;
        RealVector _relTolRes;
        RealVector _absTolCor;
        RealVector _absTolRes;

        RealVector _normCor;
        RealVector _normRes;
        RealVector _critCor;
        RealVector _critRes;
    };
}

#endif	/* CONVERGENCECHECKER_HPP */