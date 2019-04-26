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
//   Material class 'CorrectedExponentialDegradation' implements a 1-parameter energy 
//   degradation function that is used in variational (phase-field) models for 
//   fracture. The degradation function has the form
//
//                                      n
//                             -k(1 - x)
//                        1 - e             
//            g(x) =   __________________
//                                 -k
//                            1 - e
//
//   where k = k(n) defined as
//
//                             n-2            n-1
//              (n - 1)y(1 - y)     +  (1 - y)
//     k(n) = _____________________________________
//                               2n-2 
//                    n y (1 - y)
//
//                               2 
//             -(n + 1) + sqrt(5n  + 2n - 7)
//        y = _______________________________
//                      2
//                   2(n  - 2)
//
//
//   Material function flags:
//      
//      0 : energy degradation function
//      1 : 1st derivative of degradation function
//      2 : 2nd derivative of degradation function
//
//   Vector 'conState' should have exactly one component, which is the
//   value of the phase-field that will be used to calculate the degradation
//   function and its derivatives.
//
// ---------------------------------------------------------------------------
//   Declaration in input file:
//
//   <id>  CorrectedExponentialDegradation <n>
//
//   where  id = material ID
//           n = exponent (real number greater or equal to 2)
//   
// ---------------------------------------------------------------------------
//   Requirements on function arguments
//
//   conState = {phaseField}
//
//   Argument 'matData' is not used.
//
// ---------------------------------------------------------------------------

#ifndef CORRECTEDEXPONENTIALDEGRADATION_HPP
#define	CORRECTEDEXPONENTIALDEGRADATION_HPP

#include "../../Material.hpp"

namespace broomstyx
{
    class CorrectedExponentialDegradation : public Material
    {
    public:
        CorrectedExponentialDegradation();
        virtual ~CorrectedExponentialDegradation();
        
        double     givePotentialFrom( const RealVector& conState, const MaterialStatus* matStatus ) override;
        RealVector giveForceFrom( const RealVector& conState, const MaterialStatus* matStatus ) override;
        RealMatrix giveModulusFrom( const RealVector& conState, const MaterialStatus* matStatus ) override;
        void       readParamatersFrom( FILE* fp ) override;
        
    private:
        double _n;
        double _k;
        double _a2;
        double _a3;
        double _w;
    };
}

#endif	/* CORRECTEDEXPONENTIALDEGRADATION_HPP */