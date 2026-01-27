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
//   <n>  CubicElasticity <AM> <C11> <C12> <C44>
//
//   where n  = material ID
//         AM = Analysis mode: '2D'/'3D'
//         C11, C12, C44 = Elasticity tensor components
//
// ---------------------------------------------------------------------------
//   Requirements on function arguments
//
//   For 2D:
//      conState = {eps_11, eps_22, ga_12}^T
//
//   For 3D:
//      conState = {eps_11, eps_22, eps_33, ga_23, ga_13, ga_12}^T
//
// ---------------------------------------------------------------------------

#ifndef CUBICELASTICITY_HPP
#define	CUBICELASTICITY_HPP

#include "Materials/Material.hpp"

namespace broomstyx
{
    class CubicElasticity final : public Material
    {
    public:
        CubicElasticity();
        ~CubicElasticity() override;

        double     givePotentialFrom( const RealVector& conState, const MaterialStatus* matStatus ) override;
        RealVector giveForceFrom( const RealVector& conState, const MaterialStatus* matStatus ) override;
        RealMatrix giveModulusFrom( const RealVector& conState, const MaterialStatus* matStatus ) override;
        void       readParamatersFrom( FILE* fp ) override;

    private:
        enum AnalysisMode
        {
            TwoDimensional,
            ThreeDimensional
        };

        AnalysisMode _analysisMode;
        
        double _C11;
        double _C12;
        double _C44;

        void checkSizeOf( const RealVector& conState );
    };    
}

#endif	/* CUBICELASTICITY_HPP */