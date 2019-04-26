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

#ifndef ANALYSISMODEL_HPP
#define ANALYSISMODEL_HPP

#include <cstdio>
#include <string>

namespace broomstyx
{
    class DofManager;
    class DomainManager;
    class MaterialManager;
    class MeshReader;
    class NumericsManager;
    class OutputManager;
    class SolutionManager;

    class AnalysisModel final
    {
        friend AnalysisModel& analysisModel();
        
    public:
        
        void initializeYourself( std::string filename );
        void solveYourself();
        
        DofManager&      dofManager();
        DomainManager&   domainManager();
        MaterialManager& materialManager();
        NumericsManager& numericsManager();
        OutputManager&   outputManager();
        SolutionManager& solutionManager();
        
        MeshReader& meshReader();

    private:
        std::string _inputFilename;
        std::string _meshFilename;
        
        DofManager*      _dofManager;
        DomainManager*   _domainManager;
        MaterialManager* _materialManager;
        NumericsManager* _numericsManager;
        OutputManager*   _outputManager;
        SolutionManager* _solutionManager;

        MeshReader*      _meshReader;

        AnalysisModel();
        virtual ~AnalysisModel();
        
        void readInputFile();
        void readMeshReaderFrom( FILE* fp );
    };
    
    AnalysisModel& analysisModel();
}

#endif /* ANALYSISMODEL_HPP */