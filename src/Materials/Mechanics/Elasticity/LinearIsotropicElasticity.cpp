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

#include "LinearIsotropicElasticity.hpp"

#include <cstring>
#include "../../../Core/ObjectFactory.hpp"
#include "../../../Util/linearAlgebra.hpp"
#include "../../../Util/readOperations.hpp"

using namespace broomstyx;

registerBroomstyxObject(Material, LinearIsotropicElasticity)

// Constructor
LinearIsotropicElasticity::LinearIsotropicElasticity()
{
    _E = 0.;
    _nu = 0.;
    _G = 0.;
    _name = "LinearIsotropicElasticity";
}

// Destructor
LinearIsotropicElasticity::~LinearIsotropicElasticity() {}

// Public Methods
// ----------------------------------------------------------------------------
double LinearIsotropicElasticity::givePotentialFrom( const RealVector& conState, const MaterialStatus* matStatus )
{
    int reqSize;
    switch ( _analysisMode )
    {
        case 1:
            reqSize = 3;
            break;
        case 2:
        case 3:
            reqSize = 4;
            break;
        case 4:
            reqSize = 6;
            break;
        default:
            reqSize = 2;
    }
    if ( conState.dim() != reqSize )
        throw std::runtime_error("Invalid size of vector 'conState' detected!\nSource: " + _name);
    
    RealVector stress = this->giveForceFrom(conState, matStatus);
    
    return 0.5*stress.dot(conState);
}
// ----------------------------------------------------------------------------
RealVector LinearIsotropicElasticity::giveForceFrom( const RealVector& conState, const MaterialStatus* matStatus )
{
    int reqSize;
    switch ( _analysisMode )
    {
        case 1:
            reqSize = 3;
            break;
        case 2:
        case 3:
            reqSize = 4;
            break;
        case 4:
            reqSize = 6;
            break;
        default:
            reqSize = 2;
    }
    if ( conState.dim() != reqSize )
        throw std::runtime_error("Invalid size of vector 'conState' detected!\nSource: " + _name);
    
    RealMatrix modulus = this->giveModulusFrom(conState, matStatus);
    RealVector conForce;
    conForce = modulus*conState;
    
    return conForce;
}

RealMatrix LinearIsotropicElasticity::giveModulusFrom( const RealVector& conState, const MaterialStatus* matStatus )
{
    double c1, c2, lambda;
    RealMatrix conMod;
    
    lambda = _nu*_E/((1. + _nu)*(1. - 2.*_nu));
    
    switch ( _analysisMode )
    {
        case 1:  // Plane stress
            c1 = _E/(1 - _nu*_nu);
            c2 = _nu*c1;
            
            conMod = {{c1, c2, 0},
                      {c2, c1, 0},
                      {0,  0, _G}};
            break;
            
        case 2:  // Plane strain
        case 3:  // Axisymmetry
            c1 = lambda + 2.*_G;
            c2 = lambda;
            
            conMod = {{c1, c2, c2, 0},
                      {c2, c1, c2, 0},
                      {c2, c2, c1, 0},
                      {0,  0,  0, _G}};
            break;
            
        case 4:  // Full 3D
            c1 = lambda + 2*_G;
            c2 = lambda;
            
            conMod = {{c1, c2, c2, 0,  0,  0},
                      {c2, c1, c2, 0,  0,  0},
                      {c2, c2, c1, 0,  0,  0},
                      {0,  0,  0, _G,  0,  0},
                      {0,  0,  0,  0, _G,  0},
                      {0,  0,  0,  0,  0, _G}};
            break;
            
        default: // Torsion
            conMod = {{_G, 0},
                      {0, _G}};
    }
    
    return conMod;
}

double LinearIsotropicElasticity::giveParameter( const std::string& str )
{
    if ( str == "YoungModulus" )
        return _E;
    else if ( str == "PoissonRatio" )
        return _nu;
    else if ( str == "ShearModulus" )
        return _G;
    else
        throw std::runtime_error("Request for unrecognized parameter '" + str + "' made to material class " + _name);
}

void LinearIsotropicElasticity::readParamatersFrom( FILE* fp )
{
    std::string mode = getStringInputFrom(fp, "Failed to read analysis mode from input file!", _name);
    if ( mode == "PlaneStress" )
        _analysisMode = 1;
    else if ( mode == "PlaneStrain" )
        _analysisMode = 2;
    else if ( mode == "Axisymmetric" )
        _analysisMode = 3;
    else if ( mode == "3D" )
        _analysisMode = 4;
    else if ( mode == "Torsion" )
        _analysisMode = 5;
    else
        throw std::runtime_error("ERROR: Invalid value '" + mode + "' specified for analysis mode in input file!\nSource: " + _name);
    
    switch (_analysisMode)
    {
        case 5:
            _G = getRealInputFrom(fp, "Failed to read shear modulus from input file!", _name);
            break;
        default:
            _E = getRealInputFrom(fp, "Failed to read Young's modulus from input file!", _name);
            _nu = getRealInputFrom(fp, "Failed to read Poisson's ratio from input file!", _name);
            _G = _E/(2.*(1. + _nu));
    }
}