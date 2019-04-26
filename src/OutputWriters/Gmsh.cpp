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

#include "Gmsh.hpp"

#include <chrono>
#include <stdexcept>
#include "../Core/AnalysisModel.hpp"
#include "../Core/ObjectFactory.hpp"
#include "../Core/DomainManager.hpp"
#include "../MeshReaders/MeshReader.hpp"
#include "../Numerics/Numerics.hpp"
#include "../Util/readOperations.hpp"

using namespace broomstyx;

registerBroomstyxObject(OutputWriter, Gmsh)

// Constructor
Gmsh::Gmsh()
{
    // Initialize counters
    m_writeCounter = 0;
    
    m_nodeDataScalars = 0;
    m_nodeDataVectors = 0;
    m_nodeDataTensors = 0;
    
    m_elemDataScalars = 0;
    m_elemDataVectors = 0;
    m_elemDataTensors = 0;
    
    m_elemNodeDataScalars = 0;
    m_elemNodeDataVectors = 0;
    m_elemNodeDataTensors = 0;
}

// Destructor
Gmsh::~Gmsh() {}

// Public methods
void Gmsh::initialize()
{
    // Make directory
    const int dir_err = system("mkdir -p Output_Gmsh");
    if ( dir_err == -1 )
        throw std::runtime_error("Error creating directory 'Output_Gmsh'!\n");
}

void Gmsh::readDataFrom( FILE *fp )
{
    std::string str, key, src = "Gmsh (OutputWriter)";
    
    // Output filename
    verifyDeclaration(fp, key = "FILENAME", src);
    m_outputFilename = getStringInputFrom(fp, "Failed to read Gmsh output filename from input file!", src);
    
    // Output style
    verifyDeclaration(fp, key = "OUTPUT_STYLE", src);
    str = getStringInputFrom(fp, "Failed to read Gmsh output style (binary/ascii) from input file!", src);
    if ( str == "binary" )
        m_isBinary = true;
    else if ( str == "ascii" )
        m_isBinary = false;
    else
        throw std::runtime_error("Invalid output style specification '" + str + "' encountered in input file!"
                + std::string("\nValid options: 'binary' or 'ascii'\nSource: ") + src);
    
    // Node data
    verifyDeclaration(fp, key = "NODE_DATA", src);
    m_nNodeData = getIntegerInputFrom(fp, "Failed to read number of node data output from input file!", src);    
    m_nodeData.assign(m_nNodeData, OutputData());
    
    for ( int i = 0; i < m_nNodeData; i++)
    {
        str = getStringInputFrom(fp, "Failed to read node data type from input file!", src);
        if ( str == "SCALAR" )
        {
            m_nodeData[i].dataType = scalar;
            m_nodeData[i].field.assign(1, 0);
        }
        else if ( str == "VECTOR" )
        {
            m_nodeData[i].dataType = vector;
            m_nodeData[i].field.assign(3, 0);
        }
        else if ( str == "TENSOR" )
        {
            m_nodeData[i].dataType = tensor;
            m_nodeData[i].field.assign(9, 0);
        }
        else
            throw std::runtime_error("Invalid point data type '" + str + "' encountered in input file!\nSource: " + src);
        
        // Data name to use in Gmsh
        m_nodeData[i].name = getStringInputFrom(fp, "Failed reading node data label from input file!", src);
        
        // Nodal fields comprising components of output data
        for (int j = 0; j < (int)m_nodeData[i].field.size(); j++)
            m_nodeData[i].field[j] = getIntegerInputFrom(fp, "Failed to read nodal field number for point data from input file.", src);
    }
    
    // Element data
    verifyDeclaration(fp, key = "ELEMENT_DATA", src);
    m_nElemData = getIntegerInputFrom(fp, "Failed reading number of cell data output from input file!", src);
    
    m_elemData.assign(m_nElemData, OutputData());
    
    for ( int i = 0; i < m_nElemData; i++)
    {
        str = getStringInputFrom(fp, "Failed reading element data type from input file!", src);
        if ( str == "SCALAR" )
        {
            m_elemData[i].dataType = scalar;
            m_elemData[i].field.assign(1, 0);
        }
        else if ( str == "VECTOR" )
        {
            m_elemData[i].dataType = vector;
            m_elemData[i].field.assign(3, 0);
        }
        else if ( str == "TENSOR" )
        {
            m_elemData[i].dataType = tensor;
            m_elemData[i].field.assign(9, 0);
        }
        else
            throw std::runtime_error("Invalid element data type '" + str + "' encountered in input file!\nSource: " + src);
        
        // Data name to use in Paraview
        m_elemData[i].name = getStringInputFrom(fp, "Failed reading cell data label from input file!", src);
        
        // Element fields comprising components of output data
        for (int j = 0; j < (int)m_elemData[i].field.size(); j++)
            m_elemData[i].field[j] = getIntegerInputFrom(fp, "Failed to read cell field number for element data from input file!", src);
    }

    // Element node data
    verifyDeclaration(fp, key = "ELEMENT_NODE_DATA", src);
    m_nElemNodeData = getIntegerInputFrom(fp, "Failed reading number of cell data output from input file!", src);
    
    m_elemNodeData.assign(m_nElemNodeData, OutputData());
    
    for ( int i = 0; i < m_nElemNodeData; i++)
    {
        str = getStringInputFrom(fp, "Failed reading element node data type from input file!", src);
        if ( str == "SCALAR" )
        {
            m_elemNodeData[i].dataType = scalar;
            m_elemNodeData[i].field.assign(1, 0);
        }
        else if ( str == "VECTOR" )
        {
            m_elemNodeData[i].dataType = vector;
            m_elemNodeData[i].field.assign(3, 0);
        }
        else if ( str == "TENSOR" )
        {
            m_elemNodeData[i].dataType = tensor;
            m_elemNodeData[i].field.assign(9, 0);
        }
        else
            throw std::runtime_error("Invalid element node data type '" + str + "' encountered in input file!\nSource: " + src);
        
        // Data name to use in Gmsh
        m_elemNodeData[i].name = getStringInputFrom(fp, "Failed reading element node data label from input file!", src);
        
        // Element fields comprising components of output data
        for (int j = 0; j < (int)m_elemNodeData[i].field.size(); j++)
            m_elemNodeData[i].field[j] = getIntegerInputFrom(fp, "Failed to read cell field number for element node data from input file!", src);
    }
}

void Gmsh::writeOutput( double time )
{
    // build complete string for vtu filename
    std::string mshFilename;
    
    mshFilename = "./Output_Gmsh/" + m_outputFilename + "_" + std::to_string(m_writeCounter) + ".msh";
    
    // Write .msh file
    std::printf("\n  %-40s", "Writing results to file ...");
    
    std::chrono::time_point<std::chrono::system_clock> tic, toc;
    std::chrono::duration<double> tictoc;
    tic = std::chrono::high_resolution_clock::now();
    
    FILE *mshFile = fopen(mshFilename.c_str(), "w");

    // Mesh format
    std::fprintf(mshFile, "$MeshFormat");
    std::fprintf(mshFile, "\n2.2 ");
    if ( m_isBinary )
        std::fprintf(mshFile, "1 ");
    else
        std::fprintf(mshFile, "0 ");
    
    std::fprintf(mshFile, "%d", (int)sizeof(double));
    std::fprintf(mshFile, "\n$EndMeshFormat");

    // Physical Names
    std::fprintf(mshFile, "\n$PhysicalNames");
    int nPhysNames = analysisModel().domainManager().giveNumberOfPhysicalNames();
    std::fprintf(mshFile, "\n%d\n", nPhysNames);
    for (  int i = 0; i < nPhysNames; i++ )
    {
        DomainManager::PhysicalEntity physEnt = analysisModel().domainManager().giveDataForPhysicalEntity(i);
        std::fprintf(mshFile, "%d %d %d %s\n", i+1, physEnt.dimension, physEnt.entityNumber, physEnt.name.c_str());
    }
    std::fprintf(mshFile, "$EndPhysicalNames");

    // Nodes
    std::fprintf(mshFile, "\n$Nodes");
    int nNodes = analysisModel().domainManager().giveNumberOfNodes();
    std::fprintf(mshFile, "\n%d\n", nNodes);

    for ( int i = 0; i < nNodes; i++)
    {
        Node* curNode = analysisModel().domainManager().giveNode(i);
        int nodeId = analysisModel().domainManager().giveIdOf(curNode);
        std::fprintf(mshFile, "%d ", nodeId+1);
        
        RealVector coor;
        coor = analysisModel().domainManager().giveCoordinatesOf(curNode);

        if ( coor.dim() == 2 )
            std::fprintf(mshFile, "%25.15e%25.15e%25.15e\n", coor(0), coor(1), 0.0);
        else
            std::fprintf(mshFile, "%25.15e%25.15e%25.15e\n", coor(0), coor(1), coor(2));
    }
    std::fprintf(mshFile, "$EndNodes\n");
    
    // Elements
    std::fprintf(mshFile, "$Elements\n");
    int nElem = analysisModel().domainManager().giveNumberOfDomainCells();
    std::fprintf(mshFile, "%d\n", nElem);
    
    for (  int i = 0; i < nElem; i++)
    {
        Cell* curCell = analysisModel().domainManager().giveDomainCell(i);
        std::fprintf(mshFile, "%d ", i+1);
        
        int elType = analysisModel().domainManager().giveElementTypeOf(curCell);
        int label = analysisModel().domainManager().giveLabelOf(curCell);
        
        std::fprintf(mshFile, "%d 1 %d", elType, label);
        
        std::vector<Node*> cellNode = analysisModel().domainManager().giveNodesOf(curCell);
        
        // Adjust element nodes for 1-based numbering
        for ( int j = 0; j < (int)cellNode.size(); j++)
            std::fprintf(mshFile, " %d", analysisModel().domainManager().giveIdOf(cellNode[j]) + 1);
        std::fprintf(mshFile, "\n");
    }
    std::fprintf(mshFile, "$EndElements\n");
    
    // A. Node Data
    if ( m_nNodeData > 0 )
    {
        for ( int i = 0; i < m_nNodeData; i++)
        {
            std::fprintf(mshFile, "$NodeData\n");
            // String tags
            std::fprintf(mshFile, "1\n\"%s\"\n", m_nodeData[i].name.c_str());
            
            // Real tags
            std::fprintf(mshFile, "1\n%.15e\n", time);
            
            // Integer tags
            std::fprintf(mshFile, "3\n%d\n", m_writeCounter);
            int nComp;
            switch ( m_nodeData[i].dataType )
            {
                case scalar:
                    nComp = 1;
                    break;
                case vector:
                    nComp = 3;
                    break;
                default:
                    nComp = 9;
            }
            std::fprintf(mshFile, "%d\n", nComp);
            
            int nNodes = analysisModel().domainManager().giveNumberOfNodes();
            std::fprintf(mshFile, "%d\n", nNodes);
            for (  int j = 0; j < nNodes; j++ )
            {
                std::fprintf(mshFile, "%d ", j+1);
                Node* curNode = analysisModel().domainManager().giveNode(j);
                
                if ( m_nodeData[i].dataType == scalar )
                    std::fprintf(mshFile, "%25.15e\n", analysisModel().domainManager().giveFieldValueAt(curNode, m_nodeData[i].field[0]));
                else if ( m_nodeData[i].dataType == vector )
                {
                    std::fprintf(mshFile, "%25.15e%25.15e%25.15e\n",
                            analysisModel().domainManager().giveFieldValueAt(curNode, m_nodeData[i].field[0]),
                            analysisModel().domainManager().giveFieldValueAt(curNode, m_nodeData[i].field[1]),
                            analysisModel().domainManager().giveFieldValueAt(curNode, m_nodeData[i].field[2]));
                }
                else
                {
                    std::fprintf(mshFile, "%25.15e%25.15e%25.15e%25.15e%25.15e%25.15e%25.15e%25.15e%25.15e\n",
                            analysisModel().domainManager().giveFieldValueAt(curNode, m_nodeData[i].field[0]),
                            analysisModel().domainManager().giveFieldValueAt(curNode, m_nodeData[i].field[1]),
                            analysisModel().domainManager().giveFieldValueAt(curNode, m_nodeData[i].field[2]),
                            analysisModel().domainManager().giveFieldValueAt(curNode, m_nodeData[i].field[3]),
                            analysisModel().domainManager().giveFieldValueAt(curNode, m_nodeData[i].field[4]),
                            analysisModel().domainManager().giveFieldValueAt(curNode, m_nodeData[i].field[5]),
                            analysisModel().domainManager().giveFieldValueAt(curNode, m_nodeData[i].field[6]),
                            analysisModel().domainManager().giveFieldValueAt(curNode, m_nodeData[i].field[7]),
                            analysisModel().domainManager().giveFieldValueAt(curNode, m_nodeData[i].field[8]));
                }
            }
            std::fprintf(mshFile, "$EndNodeData\n");
        }
    }
    
    // B. Element Data
    if ( m_nElemData > 0 )
    {
        for ( int i = 0; i < m_nElemData; i++)
        {
            std::fprintf(mshFile, "$ElementData\n");
            // String tags
            std::fprintf(mshFile, "1\n\"%s\"\n", m_elemData[i].name.c_str());
            
            // Real tags
            std::fprintf(mshFile, "1\n%.15e\n", time);
            
            // Integer tags
            std::fprintf(mshFile, "3\n%d\n", m_writeCounter);
            int nComp;
            switch ( m_elemData[i].dataType )
            {
                case scalar:
                    nComp = 1;
                    break;
                case vector:
                    nComp = 3;
                    break;
                default:
                    nComp = 9;
            }
            std::fprintf(mshFile, "%d\n", nComp);
            
            int nCells = analysisModel().domainManager().giveNumberOfDomainCells();
            std::fprintf(mshFile, "%d\n", nCells);
            for (  int j = 0; j < nCells; j++ )
            {
                std::fprintf(mshFile, "%d ", j+1);
                Cell* curCell = analysisModel().domainManager().giveDomainCell(j);
                int label = analysisModel().domainManager().giveLabelOf(curCell);
                Numerics* numerics = analysisModel().domainManager().giveNumericsForDomain(label);
                
                if ( m_elemData[i].dataType == scalar )
                    std::fprintf(mshFile, "%25.15e\n", numerics->giveCellFieldValueAt(curCell, m_elemData[i].field[0]));
                else if ( m_elemData[i].dataType == vector )
                {
                    std::fprintf(mshFile, "%25.15e%25.15e%25.15e\n",
                            numerics->giveCellFieldValueAt(curCell, m_elemData[i].field[0]),
                            numerics->giveCellFieldValueAt(curCell, m_elemData[i].field[1]),
                            numerics->giveCellFieldValueAt(curCell, m_elemData[i].field[2]));
                }
                else
                {
                    std::fprintf(mshFile, "%25.15e%25.15e%25.15e%25.15e%25.15e%25.15e%25.15e%25.15e%25.15e\n",
                            numerics->giveCellFieldValueAt(curCell, m_elemData[i].field[0]),
                            numerics->giveCellFieldValueAt(curCell, m_elemData[i].field[1]),
                            numerics->giveCellFieldValueAt(curCell, m_elemData[i].field[2]),
                            numerics->giveCellFieldValueAt(curCell, m_elemData[i].field[3]),
                            numerics->giveCellFieldValueAt(curCell, m_elemData[i].field[4]),
                            numerics->giveCellFieldValueAt(curCell, m_elemData[i].field[5]),
                            numerics->giveCellFieldValueAt(curCell, m_elemData[i].field[6]),
                            numerics->giveCellFieldValueAt(curCell, m_elemData[i].field[7]),
                            numerics->giveCellFieldValueAt(curCell, m_elemData[i].field[8]));
                }
            }
            std::fprintf(mshFile, "$EndElementData\n");
        }
    }
    
    // C. Element Node Data
    if ( m_nElemNodeData > 0 )
    {
        for ( int i = 0; i < m_nElemNodeData; i++)
        {
            std::fprintf(mshFile, "$ElementNodeData\n");
            // String tags
            std::fprintf(mshFile, "1\n\"%s\"\n", m_elemNodeData[i].name.c_str());
            
            // Real tags
            std::fprintf(mshFile, "1\n%.15e\n", time);
            
            // Integer tags
            std::fprintf(mshFile, "3\n%d\n", m_writeCounter);
            int nComp;
            switch ( m_elemNodeData[i].dataType )
            {
                case scalar:
                    nComp = 1;
                    break;
                case vector:
                    nComp = 3;
                    break;
                default:
                    nComp = 9;
            }
            std::fprintf(mshFile, "%d\n", nComp);
            
            int nCells = analysisModel().domainManager().giveNumberOfDomainCells();
            std::fprintf(mshFile, "%d\n", nCells);
            
            for (  int j = 0; j < nCells; j++ )
            {
                Cell* curCell = analysisModel().domainManager().giveDomainCell(j);
                int nCellNodes = analysisModel().domainManager().giveNumberOfNodesOf(curCell);
                std::fprintf(mshFile, "%d %d ", j+1, nCellNodes);
                
                int label = analysisModel().domainManager().giveLabelOf(curCell);
                Numerics* numerics = analysisModel().domainManager().giveNumericsForDomain(label);
                
                if ( m_elemNodeData[i].dataType == scalar )
                {
                    RealVector cellNodeVals = numerics->giveCellNodeFieldValuesAt(curCell, m_elemNodeData[i].field[0]);
                    
                    for (  int k = 0; k < nCellNodes; k++ )
                        std::fprintf(mshFile, "%25.15e", cellNodeVals(k));
                    std::fprintf(mshFile, "\n");
                }
                else if ( m_elemNodeData[i].dataType == vector )
                {
                    std::vector<RealVector> cellNodeVals;
                    cellNodeVals.assign(3, RealVector());
                    for (  int p = 0; p < 3; p++ )
                        cellNodeVals[p] = numerics->giveCellNodeFieldValuesAt(curCell, m_elemNodeData[i].field[p]);
                    
                    for (  int k = 0; k < nCellNodes; k++ )
                        for ( int p = 0; p < 3; p++ )
                            std::fprintf(mshFile, "%25.15e", cellNodeVals[p](k));
                    std::fprintf(mshFile, "\n");
                }
                else
                {
                    std::vector<RealVector> cellNodeVals;
                    cellNodeVals.assign(9, RealVector());
                    for (  int p = 0; p < 9; p++ )
                        cellNodeVals[p] = numerics->giveCellNodeFieldValuesAt(curCell, m_elemNodeData[i].field[p]);
                    
                    for (  int k = 0; k < nCellNodes; k++ )
                        for ( int p = 0; p < 9; p++ )
                            std::fprintf(mshFile, "%25.15e", cellNodeVals[p](k));
                    std::fprintf(mshFile, "\n");
                }
            }
            std::fprintf(mshFile, "$EndElementNodeData\n");
        }
    }
    

    std::fclose(mshFile);
    
    // Increment file count
    m_writeCounter += 1;

    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    std::printf("done (time = %f sec.)\n", tictoc.count());
}