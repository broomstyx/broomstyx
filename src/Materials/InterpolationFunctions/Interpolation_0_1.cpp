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

#include "Interpolation_0_1.hpp"

#include "Core/ObjectFactory.hpp"
#include <cmath>

using namespace broomstyx;

registerBroomstyxObject( Material, Interpolation_0_1 )

// Constructor
Interpolation_0_1::Interpolation_0_1()
{
    _name = "Interpolation_0_1";
}

// Destructor
Interpolation_0_1::~Interpolation_0_1() = default;

// Public Methods
// ----------------------------------------------------------------------------
double Interpolation_0_1::givePotentialFrom( const RealVector& conState, const MaterialStatus* matStatus )
{
    this->checkSizeOf( conState );
    double c = conState( 0 );

    return 6. * std::pow( c, 5. ) - 15. * std::pow( c, 4. ) + 10. * std::pow( c, 3. );
}
// ----------------------------------------------------------------------------
RealVector Interpolation_0_1::giveForceFrom( const RealVector& conState, const MaterialStatus* matStatus )
{
    this->checkSizeOf( conState );
    double c = conState( 0 );

    RealVector conForce( 1 );
    conForce( 0 ) = 30. * std::pow( c, 4. ) - 60. * std::pow( c, 3. ) + 30. * std::pow( c, 2. );

    return conForce;
}
// -------------------------------------------------------------------------------------
RealMatrix Interpolation_0_1::giveModulusFrom( const RealVector& conState, const MaterialStatus* matStatus )
{
    this->checkSizeOf( conState );
    double c = conState( 0 );

    RealMatrix conMod( 1,1 );
    conMod( 0,0 ) = 120. * std::pow( c, 3. ) - 180. * c * c + 60. * c;
    
    return conMod;
}

// Private methods
// -------------------------------------------------------------------------------------
void Interpolation_0_1::checkSizeOf( const RealVector& conState )
{
    if ( conState.dim() != 1 )
        throw std::runtime_error( "Invalid size of vector 'conState' detected!\nSource: " + _name );
}