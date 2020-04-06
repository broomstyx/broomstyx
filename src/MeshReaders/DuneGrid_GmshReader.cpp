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
    // Note -- type alias definition for 'GridType' is given in DomainManager.hpp
    
    Dune::GridFactory<GridType> factory;
    Dune::GmshReader<GridType>::read(
        factory, 
        filename, 
        analysisModel().domainManager()._physicalEntityOfBoundaryCell,
        analysisModel().domainManager()._physicalEntityOfDomainCell);
        
    analysisModel().domainManager()._grid = factory.createGrid();

    // --- CODE TESTING ---
    auto gridView = analysisModel().domainManager()._grid->leafGridView();
    

    std::printf("\nBoundary tags = %d\n", (int)analysisModel().domainManager()._physicalEntityOfBoundaryCell.size());
    for ( int i = 0; i < analysisModel().domainManager()._physicalEntityOfBoundaryCell.size(); i++ )
        std::printf("%d\n", analysisModel().domainManager()._physicalEntityOfBoundaryCell[i]);

    std::printf("\nDomain tags = %d\n", (int)analysisModel().domainManager()._physicalEntityOfDomainCell.size());
    for ( int i = 0; i < analysisModel().domainManager()._physicalEntityOfDomainCell.size(); i++ )
        std::printf("%d\n", analysisModel().domainManager()._physicalEntityOfDomainCell[i]);

    std::printf("\n\nGRID HAS BEEN CREATED!\n\n");
    std::fflush(stdout);
    throw std::runtime_error("BREAKPOINT!\n");
}
// -------------------------------------------------------------------------------
std::vector<int> DuneGrid_GmshReader::giveFaceNodeNumbersForElementType( int elType, int face )
{
    throw std::runtime_error("ERROR: Call to unimplemented method detected!\nSource: DuneGrid_GmshReader\n\n");
}
// -------------------------------------------------------------------------------
int DuneGrid_GmshReader::giveNumberOfFacesForElementType( int elType )
{
    throw std::runtime_error("ERROR: Call to unimplemented method detected!\nSource: DuneGrid_GmshReader\n\n");
}

#endif /* USING_DUNE_GRID_BACKEND */