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

#ifndef GMSH_HPP
#define	GMSH_HPP

#include "OutputWriter.hpp"
#include <vector>
#include <string>

namespace broomstyx
{
    class Gmsh final : public OutputWriter
    {
    public:
        Gmsh();
        virtual ~Gmsh();
        
        void initialize() override;
        void readDataFrom( FILE *fp ) override;
        void writeOutput( double time ) override;
        
    private:
        enum DataType { scalar, vector, tensor };
        
        struct OutputData
        {
            DataType         dataType;
            std::vector<int> field;
            std::string      name;
        };
        
        std::string m_outputFilename;
        
        bool m_isBinary;
        
        int m_writeCounter;

        int m_nNodeData;
        int m_nodeDataScalars;
        int m_nodeDataVectors;
        int m_nodeDataTensors;

        int m_nElemData;
        int m_elemDataScalars;
        int m_elemDataVectors;
        int m_elemDataTensors;

        int m_nElemNodeData;
        int m_elemNodeDataScalars;
        int m_elemNodeDataVectors;
        int m_elemNodeDataTensors;

        std::vector<OutputData> m_nodeData;
        std::vector<OutputData> m_elemData;
        std::vector<OutputData> m_elemNodeData;
        
        int giveElementType( int dim, int nNodes );
    };
}

#endif	/* GMSH_HPP */

