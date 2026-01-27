/*
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

#include "CubicElasticity.hpp"

#include "Core/ObjectFactory.hpp"
#include "Util/linearAlgebra.hpp"
#include "Util/readOperations.hpp"

using namespace broomstyx;

registerBroomstyxObject( Material, CubicElasticity )

// Constructor
CubicElasticity::CubicElasticity()
{
    _analysisMode = TwoDimensional;
    _C11 = 0.;
    _C12 = 0.;
    _C44 = 0.;
    _name = "CubicElasticity";
}

// Destructor
CubicElasticity::~CubicElasticity() = default;

// Public Methods
// ----------------------------------------------------------------------------
double CubicElasticity::givePotentialFrom( const RealVector& conState, const MaterialStatus* matStatus )
{
    this->checkSizeOf( conState );
    
    RealVector stress = this->giveForceFrom( conState, matStatus );
    
    return 0.5 * stress.dot( conState );
}
// ----------------------------------------------------------------------------
RealVector CubicElasticity::giveForceFrom( const RealVector& conState, const MaterialStatus* matStatus )
{
    this->checkSizeOf( conState );
    
    RealMatrix modulus = this->giveModulusFrom( conState, matStatus );
    RealVector conForce;
    conForce = modulus * conState;
    
    return conForce;
}
// -------------------------------------------------------------------------------------
RealMatrix CubicElasticity::giveModulusFrom( const RealVector& conState, const MaterialStatus* matStatus )
{
    RealMatrix conMod;
    
    if ( _analysisMode == TwoDimensional )
    {
        conMod = { { _C11, _C12,   0. },
                   { _C12, _C11,   0. },
                   {   0.,   0., _C44 } };
    }
    else // _analysisMode == ThreeDimensional
    {
        conMod = { { _C11, _C12, _C12,   0.,   0.,   0. },
                   { _C12, _C11, _C12,   0.,   0.,   0. },
                   { _C12, _C12, _C11,   0.,   0.,   0. },
                   {   0.,   0.,   0., _C44,   0.,   0. },
                   {   0.,   0.,   0.,   0., _C44,   0. },
                   {   0.,   0.,   0.,   0.,   0., _C44 } };
    }

    return conMod;
}
// -------------------------------------------------------------------------------------
void CubicElasticity::readParamatersFrom( FILE* fp )
{
    std::string mode = getStringInputFrom( fp, "Failed to read analysis mode from input file!", _name );
    if ( mode == "2D" )
        _analysisMode = TwoDimensional;
    else if ( mode == "3D" )
        _analysisMode = ThreeDimensional;
    else
        throw std::runtime_error( "ERROR: Invalid analysis mode '" + mode + "' encountered in input file!\nSource: " + _name);

    _C11 = getRealInputFrom( fp, "Failed to read C_11 modulus from input file!", _name );
    _C12 = getRealInputFrom( fp, "Failed to read C_12 modulus from input file!", _name );
    _C44 = getRealInputFrom( fp, "Failed to read C_44 modulus from input file!", _name );
}

// Private methods
// -------------------------------------------------------------------------------------
void CubicElasticity::checkSizeOf( const RealVector& conState )
{
    int reqSize;
    if ( _analysisMode == TwoDimensional )
        reqSize = 3;
    else
        reqSize = 6;

    if ( conState.dim() != reqSize )
        throw std::runtime_error( "Invalid size of vector 'conState' detected!\nSource: " + _name );
}