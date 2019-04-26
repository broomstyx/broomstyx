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

#include "Material.hpp"
#include <stdexcept>

using namespace broomstyx;

// Material Status
MaterialStatus::MaterialStatus() {}

MaterialStatus::~MaterialStatus() {}

// Material
Material::Material() {}

Material::~Material() {}

// Public methods
// ----------------------------------------------------------------------------
MaterialStatus* Material::createMaterialStatus() { return nullptr; }
// ----------------------------------------------------------------------------
void Material::destroy( MaterialStatus*& matStatus ) {}
// ----------------------------------------------------------------------------
void Material::initialize() {}
// ----------------------------------------------------------------------------
void Material::readParamatersFrom( FILE* fp ) {}
// ----------------------------------------------------------------------------
void Material::updateStatusFrom( const RealVector& conState, MaterialStatus* matStatus ) {}
// ----------------------------------------------------------------------------
void Material::updateStatusFrom( const RealVector& conState, MaterialStatus* matStatus, const std::string& label ) {}
// ----------------------------------------------------------------------------
double Material::givePotentialFrom(const RealVector& conState, const MaterialStatus* matStatus)
{
    this->error_unimplemented("givePotentialFrom(...)");
    
    // Return statement is only for suppressing compilation warnings
    return 0.;
}
// ----------------------------------------------------------------------------
RealVector Material::giveForceFrom( const RealVector& conState, const MaterialStatus* matStatus )
{
    this->error_unimplemented("giveForceFrom( ... )");
    
    // Return statement is only for suppressing compilation warnings
    RealVector dummy;
    return dummy;
}
// ----------------------------------------------------------------------------
RealVector Material::giveForceFrom( const RealVector& conState, const MaterialStatus* matStatus, const std::string& label )
{
    this->error_unimplemented("giveForceFrom( ..., label )");
    
    // Return statement is only for suppressing compilation warnings
    RealVector dummy;
    return dummy;
}
// ----------------------------------------------------------------------------
RealMatrix Material::giveModulusFrom( const RealVector& conState, const MaterialStatus* matStatus )
{
    this->error_unimplemented("giveModulusFrom( ... )");
    
    // Return statement is only for suppressing compilation warnings
    RealMatrix dummy;
    return dummy;
}
// ----------------------------------------------------------------------------
RealMatrix Material::giveModulusFrom( const RealVector& conState, const MaterialStatus* matStatus, const std::string& label )
{
    this->error_unimplemented("giveModulusFrom( ..., label )");
    
    // Return statement is only for suppressing compilation warnings
    RealMatrix dummy;
    return dummy;
}
// ----------------------------------------------------------------------------
double Material::giveMaterialVariable( const std::string& label
                                     , const MaterialStatus* matStatus )
{
    this->error_unimplemented("giveMaterialVariable( ... )");
    
    // Return statement is only for suppressing compilation warnings
    return 0.;
}
// ----------------------------------------------------------------------------
double Material::giveParameter( const std::string& str )
{
    this->error_unimplemented("giveParameter( ... )");
    
    // Return statement is only for suppressing compilation warnings
    return 0.;
}

// Protected methods
// ----------------------------------------------------------------------------
void Material::error_unimplemented( const std::string& method )
{
    throw std::runtime_error("\n\nError: Call to unimplemented method '" + _name + "::" + method + "' detected!");
}