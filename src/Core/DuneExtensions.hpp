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

/*
  In order to make use of dune-grid as a backend grid manager, broomstyx
  must be compiled using dunecontrol from dune-common, and the following
  variables must be defined as part of CMAKE_FLAGS in the configuration
  file 'config.opts':

    USE_DUNE_GRID = ON
    DUNE_GRID_GRIDTYPE_SELECTOR = ON
    GRIDDIM = <grid dimension to be used>
    DWORLDDIM = <world dimension to be used>
    GRIDTYPE = <specific grid implementation to be used>
*/

#include <config.h>

#ifdef USING_DUNE_GRID_BACKEND
    using GridType = typename Dune::GridSelector::GridType;
    using LeafGridView = typename GridType::LeafGridView;
    using DomainSeedType = typename GridType::template Codim<0>::EntitySeed;
    using BoundarySeedType = typename GridType::template Codim<1>::EntitySeed;
    using VertexSeedType = typename GridType::template Codim<GridType::dimension>::EntitySeed;
#endif