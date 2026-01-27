/*
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

// ---------------------------------------------------------------------------
//   Declaration in input file:
//
//   <n>  Interpolation_0_1
//
//   where n  = material ID
// ---------------------------------------------------------------------------
//   Requirements on function arguments
//
//   conState = { c }
//
//   where c is a real number between 0 and 1.
//
// ---------------------------------------------------------------------------

#ifndef INTERPOLATION_0_1
#define	INTERPOLATION_0_1

#include "Materials/Material.hpp"

namespace broomstyx
{
    class Interpolation_0_1 final : public Material
    {
    public:
        Interpolation_0_1();
        ~Interpolation_0_1() override;

        double     givePotentialFrom( const RealVector& conState, const MaterialStatus* matStatus ) override;
        RealVector giveForceFrom( const RealVector& conState, const MaterialStatus* matStatus ) override;
        RealMatrix giveModulusFrom( const RealVector& conState, const MaterialStatus* matStatus ) override;

    private:
        void checkSizeOf( const RealVector& conState );
    };
}

#endif	/* INTERPOLATION_0_1 */