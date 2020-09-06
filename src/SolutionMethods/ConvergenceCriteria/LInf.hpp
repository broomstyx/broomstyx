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

#ifndef LINF_HPP
#define	LINF_HPP

#include "ConvergenceCriterion.hpp"
#include "Util/RealMatrix.hpp"

namespace broomstyx
{
    class LInf : public ConvergenceCriterion
    {
    public:
        LInf();
        virtual ~LInf();

        bool checkConvergenceOf( const RealVector& resid, const std::vector<Dof*>& dof ) override;
        RealMatrix giveConvergenceData() override;
        void initialize( int dofGrpNum ) override;
        void processLocalResidualContribution( const RealVector& contrib, const std::vector<int>& dofGrp, int threadNum ) override;
        void readDataFromFile( FILE* fp ) override;
        void reportConvergenceStatus() override;
        void resetResidualCriteria() override;

    private:
        std::string _trackingOption;
        int    _nThreads;
        double _contribCount;
        double _dofGrpCount;

        double _relTolCorr;
        double _relTolRes;
        double _absTolCorr;
        double _absTolRes;

        double _corrNorm;
        double _corrCrit;
        double _residNorm;
        double _residCrit;
        
        RealVector _threadContribCount;
        RealVector _threadCorrCrit;
        RealVector _threadCorrNorm;
        RealVector _threadResidCrit;
        RealVector _threadResidNorm;
    };
}

#endif	/* LINF_HPP */