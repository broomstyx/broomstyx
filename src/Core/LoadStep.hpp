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

#ifndef LOADSTEP_HPP
#define	LOADSTEP_HPP

#include <cstdio>
#include <vector>
#include "BoundaryCondition.hpp"
#include "FieldCondition.hpp"
#include "../Util/RealMatrix.hpp"
#include "../Util/RealVector.hpp"

namespace broomstyx
{    
    class SolutionMethod;

    class LoadStep final
    {
    public:
        LoadStep(int lsNum, int nStg);
        virtual ~LoadStep();

        int  giveLoadStepNum();
        void readDataFrom( FILE* fp );
        void solveYourself();
        void writeConvergenceDataForStage( int stg, RealMatrix& convDat );
        void writeIterationDataForStage( int stg, double time, int nIter );

    private:
        struct ProcessData
        {
            std::string domainTag;
            std::string directive;
        };
        
        int    _loadStepNum;
        int    _nStages;
        int    _maxSubsteps;
        double _startTime;
        double _endTime;
        double _dtime;
        
        int _writeInterval;

        std::string _name;
        
        std::vector<ProcessData> _preProcess;
        std::vector<ProcessData> _postProcess;
        
        std::vector<BoundaryCondition> _boundaryCondition;
        std::vector<FieldCondition>    _fieldCondition;
        std::vector<SolutionMethod*>   _solutionMethod;
        
        std::vector<FILE*> _convDatFile;
        std::vector<int>   _convDatCount;
        std::vector<FILE*> _iterDatFile;
        std::vector<int>   _iterDatCount;

        void findConstrainedDofs();
        void performPrefinalCalculationsAtCells();
    };
}

#endif	/* LOADSTEP_HPP */