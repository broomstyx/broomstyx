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
#ifndef DUNEGRID_GMSHREADER_HPP
#define	DUNEGRID_GMSHREADER_HPP

#include "MeshReader.hpp"

#ifdef USING_DUNE_GRID_BACKEND
#include <cstdio>
#include <string>
#include <vector>

typedef typename Dune::GridSelector::GridType GridType;
typedef typename GridType::template Codim<0>::EntitySeed DomainSeedType;

namespace broomstyx
{    
    class AnalysisModel;
    class Cell;

    class DuneGrid_GmshReader : public MeshReader
    {
    public:
        DuneGrid_GmshReader();
        virtual ~DuneGrid_GmshReader();
        
        // Disable copy constructor and assignment operator
        DuneGrid_GmshReader( const DuneGrid_GmshReader& ) = delete;
        DuneGrid_GmshReader& operator=( const DuneGrid_GmshReader& ) = delete;

        void readMeshFile( std::string filename );
        std::vector<int> giveFaceNodeNumbersForElementType( int elType, int face );
        int  giveNumberOfFacesForElementType( int elType );

    private:
        int giveElementTypeFor( int nVertices );
        int numberOfNodesForElementType( int elType );
    };
}

#endif  /* USING_DUNE_GRID_BACKEND */
#endif	/* DUNEGRID_GMSHREADER_HPP */