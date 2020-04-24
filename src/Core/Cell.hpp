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

#ifndef CELL_HPP
#define CELL_HPP

#include <config.h>

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <Util/RealVector.hpp>

#ifdef USING_DUNE_GRID_BACKEND
#include "DuneExtensions.hpp"
#endif

namespace broomstyx
{    
    class Dof;
    class Node;
    class NumericsStatus;

    class Cell
    {
        friend class DomainManager;
        friend class DofManager;

    public:
        Cell()
            : numericsStatus(nullptr)
            , _partition(-1)
            , _isPartOfDomain(false)
        {}

        virtual ~Cell() {}
        
        // Disable copy constructor and assignment operator
        Cell( const Cell& ) = delete;
        Cell& operator=( const Cell& ) = delete;

        RealVector cellData;
        NumericsStatus* numericsStatus;
        
    private:
        int _elType;
        int _label;
        int _dim;
        int _id;
        int _partition;

        bool _isPartOfDomain;
        std::vector<Node*> _node;
        std::vector<Dof*> _dof;

#ifdef USING_DUNE_GRID_BACKEND
        DomainSeedType _domainSeed;
        BoundarySeedType _boundarySeed;
#endif

        std::vector<Cell*> _neighbor;
        std::vector<Cell*> _face;
        std::vector<int> _faceOrient;

        // Used by boundary cells
        std::vector<Cell*> _assocDomCell;
        std::vector<int> _halo;

    };
}

#endif /* CELL_HPP */