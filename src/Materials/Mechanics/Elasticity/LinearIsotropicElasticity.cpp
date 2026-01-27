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

#include "LinearIsotropicElasticity.hpp"

#include "Core/ObjectFactory.hpp"
#include "Util/linearAlgebra.hpp"
#include "Util/readOperations.hpp"

using namespace broomstyx;

registerBroomstyxObject(Material, LinearIsotropicElasticity)

// Constructor
LinearIsotropicElasticity::LinearIsotropicElasticity()
{
    _analysisMode = TwoDimensional;
    _E = 0.;
    _nu = 0.;
    _G = 0.;
    _K = 0.;
    _name = "LinearIsotropicElasticity";
}

// Destructor
LinearIsotropicElasticity::~LinearIsotropicElasticity() = default;

// Public Methods
// ----------------------------------------------------------------------------
double LinearIsotropicElasticity::givePotentialFrom( const RealVector& conState, const MaterialStatus* matStatus )
{
    this->checkSizeOf( conState );
    
    RealVector stress = this->giveForceFrom( conState, matStatus );
    
    return 0.5 * stress.dot( conState );
}
double LinearIsotropicElasticity::givePotentialFrom( const RealVector& conState, const MaterialStatus* matStatus, const std::string& label )
{
    double potential;

    this->checkSizeOf( conState );

    if ( label == "Volumetric" )
    {
        RealVector volStrainVec = _Pvol * conState;
        RealVector volStressVec = this->giveForceFrom( volStrainVec, matStatus );

        potential = 0.5 * volStressVec.dot( volStrainVec );
    }
    else if ( label == "Deviatoric" )
    {
        RealVector devStrainVec = conState - _Pvol * conState;
        RealVector devStressVec = this->giveForceFrom( devStrainVec, matStatus );

        potential = 0.5 * devStressVec.dot( devStrainVec );
    }
    else
        throw std::runtime_error( "Invalid label '" + label + "' encountered for calculation of potential!\nSource: " + _name );

    return potential;
}
// ----------------------------------------------------------------------------
RealVector LinearIsotropicElasticity::giveForceFrom( const RealVector& conState, const MaterialStatus* matStatus )
{
    this->checkSizeOf( conState );
    
    RealMatrix modulus = this->giveModulusFrom( conState, matStatus );
    RealVector conForce;
    conForce = modulus * conState;
    
    return conForce;
}
// ----------------------------------------------------------------------------
RealVector LinearIsotropicElasticity::giveForceFrom( const RealVector& conState, const MaterialStatus* matStatus, const std::string& label )
{
    this->checkSizeOf( conState );

    RealMatrix modulus = this->giveModulusFrom( conState, matStatus );
    RealVector conForce;

    if ( label == "Volumetric" )
        conForce = modulus * ( _Pvol * conState );
    else if ( label == "Deviatoric" )
        conForce = modulus * ( ( _I - _Pvol ) * conState );
    else
        throw std::runtime_error( "Invalid label '" + label + "' encountered for calculation of force!\nSource: " + _name );

    return conForce;
}
// -------------------------------------------------------------------------------------
RealMatrix LinearIsotropicElasticity::giveModulusFrom( const RealVector& conState, const MaterialStatus* matStatus )
{
    double c1, c2, lambda;
    RealMatrix conMod;
    
    lambda = _nu * _E / ( ( 1. + _nu ) * ( 1. - 2. * _nu ) );
    
    if ( _analysisMode == TwoDimensional || _analysisMode == PlaneStress )
    {
        c1 = _E / ( 1. - _nu * _nu );
        c2 = _nu * c1;

        conMod = { { c1, c2, 0 },
                   { c2, c1, 0 },
                   { 0,  0,  _G } };
    }
    else if ( _analysisMode == PlaneStrain || _analysisMode == Axisymmetric )
    {
        c1 = lambda + 2. * _G;
        c2 = lambda;

        conMod = { { c1, c2, c2, 0 },
                   { c2, c1, c2, 0 },
                   { c2, c2, c1, 0 },
                   { 0,  0,  0, _G } };
    }
    else if ( _analysisMode == ThreeDimensional )
    {
        c1 = lambda + 2 * _G;
        c2 = lambda;

        conMod = { { c1, c2, c2, 0,  0,  0 },
                   { c2, c1, c2, 0,  0,  0 },
                   { c2, c2, c1, 0,  0,  0 },
                   { 0,  0,  0, _G,  0,  0 },
                   { 0,  0,  0,  0, _G,  0 },
                   { 0,  0,  0,  0,  0, _G } };
    }
    else // _analysisMode == Torsion
    {
        conMod = { { _G, 0 },
                   { 0, _G } };
    }

    return conMod;
}
// -------------------------------------------------------------------------------------
RealMatrix LinearIsotropicElasticity::giveModulusFrom( const RealVector& conState, const MaterialStatus* matStatus, const std::string& label )
{
    RealMatrix modulus = this->giveModulusFrom( conState, matStatus );
    RealMatrix conMod;

    if ( label == "Volumetric" )
        conMod = modulus * _Pvol;
    else if ( label == "Deviatoric" )
        conMod = modulus * ( _I - _Pvol );
    else
        throw std::runtime_error( "Invalid label '" + label + "' encountered for calculation of modulus!\nSource: " + _name );

    return conMod;
}
// -------------------------------------------------------------------------------------
double LinearIsotropicElasticity::giveParameter( const std::string& str )
{
    if ( str == "YoungModulus" )
        return _E;
    else if ( str == "PoissonRatio" )
        return _nu;
    else if ( str == "ShearModulus" )
        return _G;
    else if ( str == "BulkModulus" )
        return _K;
    else
        throw std::runtime_error( "Request for unrecognized parameter '" + str + "' made to material class " + _name );
}
// -------------------------------------------------------------------------------------
void LinearIsotropicElasticity::readParamatersFrom( FILE* fp )
{
    std::string mode = getStringInputFrom( fp, "Failed to read analysis mode from input file!", _name );
    if ( mode == "2D" )
        _analysisMode = TwoDimensional;
    else if ( mode == "PlaneStress" )
        _analysisMode = PlaneStress;
    else if ( mode == "PlaneStrain" )
        _analysisMode = PlaneStrain;
    else if ( mode == "Axisymmetric" )
        _analysisMode = Axisymmetric;
    else if ( mode == "3D" )
        _analysisMode = ThreeDimensional;
    else if ( mode == "Torsion" )
        _analysisMode = Torsion;
    else
        throw std::runtime_error( "ERROR: Invalid analysis mode '" + mode + "' encountered in input file!\nSource: " + _name);

    if ( _analysisMode == Torsion )
        _G = getRealInputFrom( fp, "Failed to read shear modulus from input file!", _name );
    else
    {
        _E = getRealInputFrom( fp, "Failed to read Young's modulus from input file!", _name );
        _nu = getRealInputFrom( fp, "Failed to read Poisson's ratio from input file!", _name );
        _G = _E / ( 2.* ( 1. + _nu ) );

        if ( _analysisMode == TwoDimensional )
            _K = _E / ( 2. * ( 1. - _nu ) );
        else
            _K = _E / ( 3. * ( 1. - 2. * _nu ) );
    }

    // Declare values of _Pvol and _I
    if ( _analysisMode == TwoDimensional || _analysisMode == PlaneStress )
    {
        _I = { { 1., 0., 0. },
               { 0., 1., 0. },
               { 0., 0., 1. } };

        double f = 0.5;
        _Pvol = { { f,  f,  0. },
                  { f,  f,  0. },
                  { 0., 0., 0. } };
    }
    else if ( _analysisMode == PlaneStrain || _analysisMode == Axisymmetric )
    {
        _I = { { 1., 0., 0., 0. },
               { 0., 1., 0., 0. },
               { 0., 0., 1., 0. },
               { 0., 0., 0., 1. } };

        double f = 1. / 3.;
        _Pvol = { { f,  f,  f,  0. },
                  { f,  f,  f,  0. },
                  { f,  f,  f,  0. },
                  { 0., 0., 0., 0. } };
    }
    else if ( _analysisMode == ThreeDimensional )
    {
        _I = { { 1., 0., 0., 0., 0., 0. },
               { 0., 1., 0., 0., 0., 0. },
               { 0., 0., 1., 0., 0., 0. },
               { 0., 0., 0., 1., 0., 0. },
               { 0., 0., 0., 0., 1., 0. },
               { 0., 0., 0., 0., 0., 1. } };

        double f = 1. / 3.;
        _Pvol = { { f,  f,  f,  0., 0., 0. },
                  { f,  f,  f,  0., 0., 0. },
                  { f,  f,  f,  0., 0., 0. },
                  { 0., 0., 0., 0., 0., 0. },
                  { 0., 0., 0., 0., 0., 0. },
                  { 0., 0., 0., 0., 0., 0. } };
    }
}

// Private methods
// -------------------------------------------------------------------------------------
void LinearIsotropicElasticity::checkSizeOf( const RealVector& conState )
{
    int reqSize;
    switch ( _analysisMode )
    {
        case TwoDimensional:
        case PlaneStress:
            reqSize = 3;
            break;
        case PlaneStrain:
        case Axisymmetric:
            reqSize = 4;
            break;
        case ThreeDimensional:
            reqSize = 6;
            break;
        default: // Torsion
            reqSize = 2;
    }

    if ( conState.dim() != reqSize )
        throw std::runtime_error( "Invalid size of vector 'conState' detected!\nSource: " + _name );
}