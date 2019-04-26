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

#include "Density.hpp"
#include "../../../Core/ObjectFactory.hpp"
#include "../../../Util/readOperations.hpp"

using namespace broomstyx;

registerBroomstyxObject(Material, Density)

Density::Density()
{
    _name = "Density";
}

Density::~Density() {}

// Public methods
double Density::giveMaterialVariable( const std::string& str, const MaterialStatus* matStatus )
{
    if ( str == "Density" )
        return _rho;
    else
        throw std::runtime_error("Request for unrecognized material variable '" + str + "' made to material class " + _name);
}

double Density::giveParameter( const std::string& str )
{
    if ( str == "Density" )
        return _rho;
    else
        throw std::runtime_error("Request for unrecognized parameter '" + str + "' made to material class " + _name);
}

void Density::readParamatersFrom( FILE* fp )
{
    _rho = getRealInputFrom(fp, "Failed to read material density from input file!", _name);
}