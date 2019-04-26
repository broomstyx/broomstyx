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

#include "Cell.hpp"

#include "AnalysisModel.hpp"
#include "Dof.hpp"
#include "Node.hpp"

using namespace broomstyx;

Cell::Cell()
    : numericsStatus(nullptr)
    , _partition(-1)
    , _isPartOfDomain(false)
{}

Cell::~Cell() {}

int Cell::id()
{
    return _id;
}

void Cell::showInfo()
{
    std::printf("\n  Cell ID = %d, label = %d, dim = %d\n", _id, _label, _dim);
    std::printf("  Cell nodes: ");
    
    for ( int i = 0; i < (int)_node.size(); ++i )
        std::printf("%d ", _node[i]->id());
    
    std::printf("\n  Neighbors: ");
    for ( int i = 0; i < (int)_neighbor.size(); ++i )
        if ( _neighbor[i] )
            std::printf("%d ", _neighbor[i]->id());
        else
            std::printf("null ");
    std::printf("\n");
    
    if ( !_isPartOfDomain )
    {
        std::printf("  Number of associated domain cells = %d\n", (int)_assocDomCell.size());
        std::printf("  Associated domain cell = ");
        for ( auto it = _assocDomCell.begin(); it != _assocDomCell.end(); ++it )
            std::printf("%d ", (*it)->id());
        std::printf("\n");
    }
}