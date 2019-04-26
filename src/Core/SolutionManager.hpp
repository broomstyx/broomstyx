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

#ifndef SOLUTIONMANAGER_HPP
#define	SOLUTIONMANAGER_HPP

#include <cstdio>
#include <set>
#include <vector>
#include <string>
#include "InitialCondition.hpp"

namespace broomstyx
{
    class LinearSolver;
    class LoadStep;
    class Node;
    class SolutionMethod;
    class UserFunction;

    class SolutionManager final
    {
        friend class AnalysisModel;
        
    public:
        void          commenceSolution();
        LoadStep*     giveCurrentLoadStep();
        int           giveNumberOfSolutionStages();
        std::vector<int>
                      giveRegisteredSolutionStages();
        void          imposeInitialConditions();
        UserFunction* makeNewUserFunction(std::string name);
        void          readInitialConditionsFrom( FILE* fp );
        void          readLoadStepsFrom( FILE* fp );
        void          readNumberOfStagesFrom( FILE* fp );
        void          registerStage( int stage, std::string tag );
        void          reportRegisteredStages();

    private:
        int _nStages;
        std::set<int> _stage;
        std::vector<LoadStep*> _loadStep;
        std::vector<InitialCondition> _initCond;
        LoadStep* _curLoadStep;
        std::vector<UserFunction*> _userFunction;
        
        std::string _name;
        
        SolutionManager();
        virtual ~SolutionManager();
    };
}

#endif	/* SOLUTIONMANAGER_HPP */