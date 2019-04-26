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

// ---------------------------------------------------------------------------
//   Declaration in input file:
//
//   <n>  LinearIsotropicElasticity <AM> <E> <nu>
//
//   where n  = material ID
//         AM = Analysis mode: 'PlaneStress'/'PlaneStrain'/'Axisymmetric'/'3D'
//          E = Young's modulus
//         nu = Poisson ratio
//
// ---------------------------------------------------------------------------
//   Requirements on function arguments
//
//   For plane stress:
//      conState = {eps_11, eps_22, ga_12}^T
//
//   For plane strain, axisymmetry:
//      conState = {eps_11, eps_22, eps_33, ga_12}^T
//
//   For 3D:
//      conState = {eps_11, eps_22, eps_33, ga_23, ga_13, ga_12}^T
//
//   For Torsion:
//      conState = {ga_13, ga_23}^T 
//
// ---------------------------------------------------------------------------

#ifndef LINEARISOTROPICELASTICITY_HPP
#define	LINEARISOTROPICELASTICITY_HPP

#include "../../Material.hpp"

namespace broomstyx
{
    class LinearIsotropicElasticity final : public Material
    {
    public:
        LinearIsotropicElasticity();
        virtual ~LinearIsotropicElasticity();

        double     givePotentialFrom( const RealVector& conState, const MaterialStatus* matStatus ) override;
        RealVector giveForceFrom( const RealVector& conState, const MaterialStatus* matStatus ) override;
        RealMatrix giveModulusFrom( const RealVector& conState, const MaterialStatus* matStatus ) override;
        double     giveParameter( const std::string& str ) override;
        void       readParamatersFrom( FILE* fp ) override;

    private:
        int _analysisMode;
        
        double _E;
        double _nu;
        double _G;
    };    
}

#endif	/* LINEARISOTROPICELASTICITY_HPP */