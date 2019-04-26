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

#ifndef GMSHREADER_HPP
#define	GMSHREADER_HPP

#include "MeshReader.hpp"

namespace broomstyx 
{
    class AnalysisModel;

    class GmshReader : public MeshReader
    {
    public:
        GmshReader();
        virtual ~GmshReader();

        void readMeshFile( std::string filename );
        std::vector<int> 
             giveFaceNodeNumbersForElementType( int elType, int face );
        int  giveNumberOfFacesForElementType( int elType );

    private:
        int numberOfNodesForElementType( int elType );
    };
}

#endif	/* GMSH_FORMAT_HPP */