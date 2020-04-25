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
#include <config.h>
#include <cstdio>

#ifdef USING_DUNE_GRID_BACKEND

#include "DuneGrid_GmshReader.hpp"

#include <Core/AnalysisModel.hpp>
#include <Core/Cell.hpp>
#include <Core/DomainManager.hpp>
#include <Core/ObjectFactory.hpp>
#include <Util/readOperations.hpp>

#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/io/file/gmshreader.hh>

using namespace broomstyx;

registerBroomstyxObject(MeshReader, DuneGrid_GmshReader)

// Constructor
DuneGrid_GmshReader::DuneGrid_GmshReader() {}

// Destructor
DuneGrid_GmshReader::~DuneGrid_GmshReader() {}

// Public methods
// -------------------------------------------------------------------------------
void DuneGrid_GmshReader::readMeshFile( std::string filename )
{
    // First pass: read physical entity names from mesh file

    FILE* fp = std::fopen(filename.c_str(), "r");
    if ( !fp )
        throw std::runtime_error("\nSpecified mesh file '" + filename + "' not found!");
    
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
    fclose(fp);

    // Second pass: Read mesh
    
    std::vector<int> physEntDomain;
    std::vector<int> physEntBoundary;

    Dune::GridFactory<GridType> factory;
    Dune::GmshReader<GridType>::read(factory, filename, physEntBoundary, physEntDomain);
    analysisModel().domainManager()._grid = factory.createGrid();
    *(analysisModel().domainManager()._gridView) = analysisModel().domainManager()._grid->leafGridView();

    // Make broomstyx node objects
    for ( auto& vertex : vertices(*(analysisModel().domainManager()._gridView)) )
    {
        VertexSeedType seed = vertex.seed();
        analysisModel().domainManager().makeNewNodeFrom(seed);
    }
    analysisModel().domainManager().countNodes();

    // Make broomstyx domain cell objects
    for ( auto& element : elements(*(analysisModel().domainManager()._gridView)) )
    {
        int idx = factory.insertionIndex(element);
        if ( idx >= physEntDomain.size() )
            throw std::runtime_error("ERROR: Encountered domain cell without an assigned physical number!");

        DomainSeedType seed = element.seed();
        int cellLabel = physEntDomain[idx];
        Cell* curDomainCell = analysisModel().domainManager().makeNewDomainCellFrom(seed, cellLabel);
        
        // Get vertices of domain cell
        int nDomCellVtx = element.subEntities(GridType::dimension);

        int elemType = this->giveElementTypeFor(nDomCellVtx);
        analysisModel().domainManager().setElementTypeOf(curDomainCell, elemType);
        
        std::vector<VertexSeedType> domCellVertexSeeds;
        for ( int i = 0; i < nDomCellVtx; i++ )
            domCellVertexSeeds.push_back((element.subEntity<GridType::dimension>(i)).seed());

        analysisModel().domainManager().setNodesOf(curDomainCell, domCellVertexSeeds);

        // Make broomstyx boundary cell objects if domain cell lies on the boundary
        for ( auto& intersection : intersections(*(analysisModel().domainManager()._gridView), element) )
        {
            if ( factory.wasInserted(intersection) )
            {
                int idx = factory.insertionIndex(intersection);
                cellLabel = physEntBoundary[idx];
                Cell* curBoundaryCell = analysisModel().domainManager().makeNewBoundaryCellFrom(seed, cellLabel);

                // Get vertices of boundary cell
                int insIdx = intersection.indexInInside();
                // auto facet = element.subEntity<1>(insIdx);
                // int nBndCellVtx = facet.subEntities(GridType::dimension);

                std::vector<int> faceNodeNumbers = this->giveFaceNodeNumbersForElementType(elemType, insIdx);
                std::vector<VertexSeedType> bndCellVertexSeeds;

                for ( int i = 0; i < (int)faceNodeNumbers.size(); i++ )
                    bndCellVertexSeeds.push_back(domCellVertexSeeds[faceNodeNumbers[i]]);

                analysisModel().domainManager().setNodesOf(curBoundaryCell, bndCellVertexSeeds);
            }
        }
    }

    std::printf("\n\nGRID HAS BEEN CREATED!\n\n");
    std::fflush(stdout);
    throw std::runtime_error("BREAKPOINT!\n");
}
// -------------------------------------------------------------------------------
std::vector<int> DuneGrid_GmshReader::giveFaceNodeNumbersForElementType( int elType, int face )
{
    // Uses Dune system of numbering

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
                    faceNodeNum = {2, 0};
                    break;
                case 2:
                    faceNodeNum = {1, 2};
                    break;
                default:
                    std::printf("\nERROR: Unrecognized face number '%d' for element type '%d'!", face, elType );
                    std::printf("\nSource: DuneGrid_GmshReader\n");
                    throw std::runtime_error("\n");
            }
            break;
        case 3: // 4-node quad
            switch ( face )
            {
                case 0:
                    faceNodeNum = {2, 0};
                    break;
                case 1:
                    faceNodeNum = {1, 3};
                    break;
                case 2:
                    faceNodeNum = {0, 1};
                    break;
                case 3:
                    faceNodeNum = {3, 2};
                    break;
                default:
                    std::printf("\nERROR: Unrecognized face number '%d' for element type '%d'!", face, elType );
                    std::printf("\nSource: DuneGrid_GmshReader\n");
                    throw std::runtime_error("\n");
            }
            break;
        case 4: // 4-node tetrahedron
            switch ( face )
            {
                case 0:
                    faceNodeNum = {0, 2, 1};
                    break;
                case 1:
                    faceNodeNum = {0, 1, 3};
                    break;
                case 2:
                    faceNodeNum = {0, 3, 2};
                    break;
                case 3:
                    faceNodeNum = {1, 2, 3};
                    break;
                default:
                    std::printf("\nERROR: Unrecognized face number '%d' for element type '%d'!", face, elType );
                    std::printf("\nSource: DuneGrid_GmshReader\n");
                    throw std::runtime_error("\n");
            }
            break;
        case 5: // 8-node hexahedron
            switch ( face )
            {
                case 0:
                    faceNodeNum = {0, 4, 6, 2};
                    break;
                case 1:
                    faceNodeNum = {1, 3, 7, 5};
                    break;
                case 2:
                    faceNodeNum = {0, 1, 5, 4};
                    break;
                case 3:
                    faceNodeNum = {2, 6, 7, 3};
                    break;
                case 4:
                    faceNodeNum = {0, 2, 3, 1};
                    break;
                case 5:
                    faceNodeNum = {4, 5, 6, 7};
                    break;
                default:
                    std::printf("\nERROR: Unrecognized face number '%d' for element type '%d'!", face, elType );
                    std::printf("\nSource: DuneGrid_GmshReader\n");
                    throw std::runtime_error("\n");
            }
            break;
        default:
            std::printf("\nERROR: Cannot generate faces for unrecognized element type '%d'!", elType);
            std::printf("\nSource: DuneGrid_GmshReader\n");
            throw std::runtime_error("\n");
    }
    
    return faceNodeNum;
}
// -------------------------------------------------------------------------------
int DuneGrid_GmshReader::giveNumberOfFacesForElementType( int elType )
{
    int nFaces;
    
    switch ( elType)
    {
        case 2: // 3-node triangle
        case 9: // 6-node triangle
            nFaces = 3;
            break;
        case 3: // 4-node quad
        case 4: // 4-node tetrahedron
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
// ---------------------------------------------------------------------------------------------
int DuneGrid_GmshReader::giveElementTypeFor( int nVertices )
{
    // Give gmsh element type for cell based on number of vertices
    
    int elType;

    if ( GridType::dimension == 2 )
    {
        switch (nVertices)
        {
            case 3:
                elType = 2;
                break;
            case 4:
                elType = 3;
                break;
            default:
                std::printf("\nERROR: %d-dimensional element with %d vertices is not recognized!", GridType::dimension, elType);
                std::printf("\nSource: DuneGrid_GmshReader\n");
                throw std::runtime_error("\n");
        }
    }
    else if ( GridType::dimension == 3 )
    {
        switch (nVertices)
        {
            case 4:
                elType = 4;
                break;
            case 8:
                elType = 5;
                break;
            default:
                std::printf("\nERROR: %d-dimensional element with %d vertices is not recognized!", GridType::dimension, elType);
                std::printf("\nSource: DuneGrid_GmshReader\n");
                throw std::runtime_error("\n");
        }
    }
    else
    {
        std::printf("\nERROR: Grid with dimension = %d cannot be handled!", GridType::dimension);
        std::printf("\nSource: DuneGrid_GmshReader\n");
        throw std::runtime_error("\n");
    }
}
// ---------------------------------------------------------------------------------------------
int DuneGrid_GmshReader::numberOfNodesForElementType( int elType ) 
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

#endif /* USING_DUNE_GRID_BACKEND */