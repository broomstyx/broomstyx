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

#ifdef USING_DUNE_GRID_BACKEND

#include "DuneGrid_GmshReader.hpp"

#include <Core/AnalysisModel.hpp>
#include <Core/DomainManager.hpp>
#include <Core/ObjectFactory.hpp>
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
    // Note: type alias definition for 'GridType' is given in DomainManager.hpp
    
    Dune::GridFactory<GridType> factory;
    std::vector<int> boundarySegmentToPhysicalEntity;
    std::vector<int> elementToPhysicalEntity;

    Dune::GmshReader<GridType>::read(
        factory, 
        filename, 
        analysisModel().domainManager()._physicalEntityOfBoundaryCell,
        analysisModel().domainManager()._physicalEntityOfDomainCell);
        
    analysisModel().domainManager()._grid = factory.createGrid();
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