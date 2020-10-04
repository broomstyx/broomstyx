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

#include "GmshReader.hpp"

#include <chrono>
#include <cstdio>
#include <vector>
#include <string>
#include <omp.h>

#include "Core/AnalysisModel.hpp"
#include "Core/DomainManager.hpp"
#include "Core/DomainManager.hpp"
#include "Core/ObjectFactory.hpp"
#include "Util/RealVector.hpp"
#include "Util/readOperations.hpp"

using namespace broomstyx;

registerBroomstyxObject(MeshReader, GmshReader)

// Constructor
GmshReader::GmshReader() {}

// Destructor
GmshReader::~GmshReader() {}

// Public methods
void GmshReader::readMeshFile( std::string filename ) {
    FILE* fp = fopen(filename.c_str(), "r");
    if ( !fp )
        throw std::runtime_error("\nSpecified mesh file '" + filename + "' not found!");
    
    std::printf("  %-40s", "Reading mesh file ...");
    std::fflush(stdout);
    std::chrono::time_point<std::chrono::system_clock> tic, toc;
    std::chrono::duration<double> tictoc;
    tic = std::chrono::high_resolution_clock::now();
    
    std::string str, errmsg, src = "GmshReader (MeshReader)";
    
    // Mesh format
    verifyDeclaration(fp, "$MeshFormat", src);
    double meshFormat = getRealInputFrom(fp, "Failed to read mesh format from mesh file!", src);
    int isBinary = getIntegerInputFrom(fp, "Failed to read ASCII/Binary flag from mesh file!", src);
    int sizeOfDouble = getIntegerInputFrom(fp, "Failed to read size of double from mesh file!", src);
    verifyDeclaration(fp, "$EndMeshFormat", src);
    
    // Physical names
    verifyDeclaration(fp, "$PhysicalNames", src);
    int nPhysEnt = getIntegerInputFrom(fp, "Failed to read number of physical entities from mesh file!", src);
    
    for ( int i = 0; i < nPhysEnt; i++)
    {
        int dimension = getIntegerInputFrom(fp, "Failed to read physical dimension from mesh file!", src);
        int entityNumber = getIntegerInputFrom(fp, "Failed to read physical entity number from mesh file!", src);
        std::string name = getStringInputFrom(fp, "Failed to read physical name from mesh file!", src);
        
        analysisModel().domainManager().createPhysicalEntity(dimension, entityNumber, name);
    }
    verifyDeclaration(fp, "$EndPhysicalNames", src);
    
    // Nodes
    verifyDeclaration(fp, "$Nodes", src);
    int nNodes = getIntegerInputFrom(fp, "Failed to read number of nodes from mesh file!", src);
    
    RealVector coor(3);
    for ( int i = 0; i < nNodes; i++)
    {
        int nodeNum = getIntegerInputFrom(fp, "Failed to read node number from mesh file!", src);
        
        if ( nodeNum != i+1 )
            throw std::runtime_error("Node number" + std::to_string(nodeNum) + "skipped in mesh file!\nSource: " + src);

        // Read nodal coordinates
        coor(0) = getRealInputFrom(fp, "Failed reading nodal x-coor. from mesh file", src);
        coor(1) = getRealInputFrom(fp, "Failed reading nodal y-coor. from mesh file", src);
        coor(2) = getRealInputFrom(fp, "Failed reading nodal z-coor. from mesh file", src);
        
        analysisModel().domainManager().makeNewNodeAt(coor);
        
    }
    verifyDeclaration(fp, "$EndNodes", src);
    analysisModel().domainManager().countNodes();
    
    // Elements
    verifyDeclaration(fp, "$Elements", src);
    int nElements = getIntegerInputFrom(fp, "Failed to read number of elements from mesh file!", src);
    
    int elementNumber, elementType;
    std::vector<int> elementTag, elementNode;
    
    for ( int i = 0; i < nElements; i++) {
        // Read element number
        elementNumber = getIntegerInputFrom(fp,
                errmsg = "Failed reading element number from mesh file!", src);
        
        if ( elementNumber != i+1 )
            throw std::runtime_error("Element number" 
                    + std::to_string(elementNumber)
                    + "skipped in mesh file!\nSource: " + src);
        
        // Read element type and determine number of element nodes
        elementType = getIntegerInputFrom(fp, errmsg = "Failed reading element type from mesh file!", src);
        
        // Read number of tags
        int nTags = getIntegerInputFrom(fp,
                errmsg = "Failed to read number of element tags from mesh file!"
                , src);
        
        elementTag.assign(nTags, 0);
        for ( int j = 0; j < nTags; j++)
            elementTag[j] = getIntegerInputFrom(fp, "Failed reading element tag from mesh file!", src);
        
        // Read number of nodes
        int nElementNodes = numberOfNodesForElementType(elementType);
        elementNode.assign(nElementNodes, 0);
        
        for ( int j = 0; j < nElementNodes; j++)
        {
            elementNode[j] = getIntegerInputFrom(fp, "Failed reading element nodes from mesh file!", src);
            
            // Adjust for 0-based nodal numbering
            --elementNode[j];
        }
        
        // Create cell
        Cell* curCell = analysisModel().domainManager().makeNewCellWithLabel( elementTag[0] );
        
        analysisModel().domainManager().setElementTypeOf(curCell, elementType);
        if ( nTags > 2 )
        {
            int nPartitionTags = elementTag[2];
            int elemPartition = elementTag[3];
            
            if ( nPartitionTags == 1 )
                analysisModel().domainManager().setPartitionOf(curCell, elemPartition);
            else
            {
                analysisModel().domainManager().setPartitionOf(curCell, 0);
            
                std::vector<int> halo;
                halo.assign(nPartitionTags - 1, 0);
                for ( int j = 0; j < nPartitionTags - 1; j++)
                    halo[j] = -elementTag[ 4+j ];
                analysisModel().domainManager().setHaloOf(curCell, halo);
            }
        }
        else
        {
            // Cells are given a partition number of 0 if no partitioning information is given in the mesh file.
            analysisModel().domainManager().setPartitionOf(curCell, 0);
        }
        analysisModel().domainManager().setNodesOf(curCell, elementNode);
    }
    verifyDeclaration(fp, "$EndElements", src);
    
    analysisModel().domainManager().countBoundaryCells();
    analysisModel().domainManager().countDomainCells();
    analysisModel().domainManager().formDomainPartitions();
    
    fclose(fp);
    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    std::printf("done (time = %f sec.)\n\n", tictoc.count());
    
    std::printf("    Meshformat = %.1f, binary flag = %d, size of double = %d\n\n", meshFormat, isBinary, sizeOfDouble);
    
    analysisModel().domainManager().reportStatus();
    
    std::printf("\n");
}

std::vector<int> GmshReader::giveFaceNodeNumbersForElementType( int elType, int face )
{
    std::vector<int> faceNodeNum;
    
    switch ( elType )
    {
        case 2: // 3-node triangle
            switch ( face )
            {
                case 0:
                    faceNodeNum = {0, 1};
                    break;
                case 1:
                    faceNodeNum = {1, 2};
                    break;
                case 2:
                    faceNodeNum = {2, 0};
                    break;
                default:
                    std::printf("\nERROR: Unrecognized face number '%d' for element type '%d'!", face, elType );
                    std::printf("\nSource: GmshReader (MeshReader)\n");
                    throw std::runtime_error("\n");
            }
            break;
        case 4: // 4-node tetrahedron
            switch ( face )
            {
                case 0:
                    faceNodeNum = {0, 3, 2};
                    break;
                case 1:
                    faceNodeNum = {0, 1, 3};
                    break;
                case 2:
                    faceNodeNum = {0, 2, 1};
                    break;
                case 3:
                    faceNodeNum = {1, 2, 3};
                    break;
                default:
                    std::printf("\nERROR: Unrecognized face number '%d' for element type '%d'!", face, elType );
                    std::printf("\nSource: GmshReader (MeshReader)\n");
                    throw std::runtime_error("\n");
            }
            break;
        case 9: // 6-node triangle
            switch ( face )
            {
                case 0:
                    faceNodeNum = {0, 3, 1};
                    break;
                case 1:
                    faceNodeNum = {1, 4, 2};
                    break;
                case 2:
                    faceNodeNum = {2, 5, 0};
                    break;
                default:
                    std::printf("\nERROR: Unrecognized face number '%d' for element type '%d'!", face, elType );
                    std::printf("\nSource: GmshReader (MeshReader)\n");
                    throw std::runtime_error("\n");
            }
            break;
        case 16: // 8-node quadrilateral
            switch ( face )
            {
                case 0:
                    faceNodeNum = {0, 4, 1};
                    break;
                case 1:
                    faceNodeNum = {1, 5, 2};
                    break;
                case 2:
                    faceNodeNum = {2, 6, 3};
                    break;
                case 3:
                    faceNodeNum = {3, 7, 0};
                    break;
                default:
                    std::printf("\nERROR: Unrecognized face number '%d' for element type '%d'!", face, elType );
                    std::printf("\nSource: GmshReader (MeshReader)\n");
                    throw std::runtime_error("\n");
            }
            break;
        default:
            std::printf("\nERROR: Cannot generate faces for unrecognized element type '%d'!", elType);
            std::printf("\nSource: GmshReader (MeshReader)\n");
            throw std::runtime_error("\n");
    }
    
    return faceNodeNum;
}

int GmshReader::giveNumberOfFacesForElementType(int elType)
{
    int nFaces;
    
    switch ( elType)
    {
        case 2: // 3-node triangle
            nFaces = 3;
            break;
        case 4: // 4-node tetrahedron
            nFaces = 4;
            break;
        case 9: // 6-node triangle
            nFaces = 3;
            break;
        case 16: // 8-node quadrilateral
            nFaces = 4;
            break;
        default:
            std::printf("\nERROR: Cannot give number of faces for element type '%d'!", elType);
            std::printf("\nSource: GmshReader\n");
            throw std::runtime_error("\n");
    }
    
    return nFaces;
}

// Private methods
int GmshReader::numberOfNodesForElementType( int elType ) 
{
    int nNodes;

    switch (elType)
    {
        case 1: // 2-node line
            nNodes = 2;
            break;
        case 2: // 3-node triangle
            nNodes = 3;
            break;
        case 3: // 4-node quadrangle
            nNodes = 4;
            break;
        case 4: // 4-node tetrahedron
            nNodes = 4;
            break;
        case 5: // 8-node hexahedron
            nNodes = 8;
            break;
        case 6: // 6-node prism
            nNodes = 6;
            break;
        case 7: // 5-node pyramid
            nNodes = 5;
            break;
        case 8: // 3-node (2nd order) line
            nNodes = 3;
            break;
        case 9: // 6-node (2nd order) triangle
            nNodes = 6;
            break;
        case 10: // 9-node (2nd order) quadrangle
            nNodes = 9;
            break;
        case 11: // 10-node (2nd order) tetrahedron
            nNodes = 10;
            break;
        case 12: // 27-node (2nd order) hexahedron
            nNodes = 27;
            break;
        case 13: // 18-node (2nd order) prism
            nNodes = 18;
            break;
        case 14: // 14-node (2nd order) pyramid
            nNodes = 14;
            break;
        case 15: // 1-node point
            nNodes = 1;
            break;
        case 16: // 8-node (2nd order) quadrangle
            nNodes = 8;
            break;
        case 17: // 20-node (2nd order) hexahedron
            nNodes = 20;
            break;
        case 18: // 15-node (2nd order) prism
            nNodes = 15;
            break;
        case 19: // 13-node (2nd order) pyramid
            nNodes = 13;
            break;
        case 20: // 9-node (3rd order incomplete) triangle
            nNodes = 9;
            break;
        case 21: // 10-node (3rd order) triangle
            nNodes = 10;
            break;
        case 22: // 12-node (4th order incomplete) triangle
            nNodes = 12;
            break;
        case 23: // 15-node (4th order) triangle
            nNodes = 15;
            break;
        case 24: // 15-node (5th order incomplete) triangle
            nNodes = 15;
            break;
        case 25: // 21-node (5th order) triangle
            nNodes = 21;
            break;
        case 26: // 4-node (3rd order) line
            nNodes = 4;
            break;
        case 27: // 5-node (4th order) line
            nNodes = 5;
            break;
        case 28: // 6-node (5th order) line
            nNodes = 6;
            break;
        case 29: // 20-node (3rd order) tetrahedron
            nNodes = 20;
            break;
        case 30: // 35-node (4th order) tetrahedron
            nNodes = 35;
            break;
        case 31: // 56-node (5th order) tetrahedron
            nNodes = 56;
            break;
        case 92: // 64-node (3rd order) hexahedron 
            nNodes = 64;
            break;
        case 93: // 125-node (4th order) hexahedron
            nNodes = 125;
            break;
        default:
            std::printf("\nERROR: Unrecognized element type '%d'!", elType);
            std::printf("\nSource: GmshReader\n");
            throw std::runtime_error("\n");
    }
    
    return nNodes;
}