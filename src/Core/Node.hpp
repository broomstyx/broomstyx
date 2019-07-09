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

#ifndef NODE_HPP
#define	NODE_HPP

#include <set>
#include <vector>
#include "../Util/RealVector.hpp"

namespace broomstyx
{
    class Dof;
    class Cell;

    class Node
    {
        friend class DomainManager;
        friend class DofManager;

    public:
        Node();
        virtual ~Node();
        
        // Disable copy constructor and assignment operator
        Node( const Node& ) = delete;
        Node& operator=( const Node& ) = delete;
        
        int id();
        
    private:
        int _id;

        RealVector _coordinates;

        std::set<Cell*> _attachedBndCell;
        std::set<Cell*> _attachedDomCell;

        RealVector _fieldVal;
        std::vector<Dof*> _dof;

        bool _isActive;

    };
}

#endif	/* NODE_HPP */