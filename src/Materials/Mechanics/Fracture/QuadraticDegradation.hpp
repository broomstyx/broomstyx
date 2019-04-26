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

//
//   Material class 'QuadraticDegradation' implements the standard quadratic
//   degradation function that is used in variational (phase-field) models for 
//   fracture. The degradation function has the form
//
//                          2
//            g(x) = (1 - x)
//
//   Vector 'conState' should have exactly one component, which is the
//   value of the phase-field that will be used to calculate the degradation
//   function and its derivatives.
//
// ---------------------------------------------------------------------------
//   Declaration in input file:
//
//   <id>  QuadraticDegradation
//
//   where  id = material ID
//   
// ---------------------------------------------------------------------------

#ifndef QUADRATICDEGRADATION_HPP
#define	QUADRATICDEGRADATION_HPP

#include "../../Material.hpp"

namespace broomstyx
{
    class QuadraticDegradation : public Material
    {
    public:
        QuadraticDegradation();
        virtual ~QuadraticDegradation();

        double     givePotentialFrom( const RealVector& conState, const MaterialStatus* matStatus ) override;
        RealVector giveForceFrom( const RealVector& conState, const MaterialStatus* matStatus ) override;
        RealMatrix giveModulusFrom( const RealVector& conState, const MaterialStatus* matStatus ) override;
        
    private:
    };
}

#endif	/* QUADRATICDEGRADATION_HPP */