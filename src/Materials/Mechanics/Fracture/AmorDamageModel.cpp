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

#include "AmorDamageModel.hpp"
#include "Core/ObjectFactory.hpp"
#include "Util/readOperations.hpp"
#include "Util/linearAlgebra.hpp"

using namespace broomstyx;

registerBroomstyxObject( Material, AmorDamageModel )

// Material status
MaterialStatus_AmorDamageModel::MaterialStatus_AmorDamageModel()
    : _materialStatus { nullptr, nullptr }
    , _volElasticEnergy( 0. )
    , _devElasticEnergy( 0. )
{}

// Constructor
AmorDamageModel::AmorDamageModel()
    : _analysisMode( TwoDimensional )
    , _elasticityModel( nullptr )
    , _degradationFunction( nullptr )
    , _kT( 0. )
    , _kC( 0. )
{
    _name = "AmorDamageModel";
}

// Public methods
// ----------------------------------------------------------------------------
MaterialStatus* AmorDamageModel::createMaterialStatus()
{
    MaterialStatus* matStatus = new MaterialStatus_AmorDamageModel();
    auto mst = this->accessMaterialStatus( matStatus );
    
    mst->_materialStatus[ 0 ] = _elasticityModel->createMaterialStatus();
    mst->_materialStatus[ 1 ] = _degradationFunction->createMaterialStatus();
    
    return matStatus;
}
// ----------------------------------------------------------------------------
void AmorDamageModel::destroy( MaterialStatus*& matStatus )
{
    auto mst = this->accessMaterialStatus( matStatus );
    delete mst->_materialStatus[ 0 ];
    delete mst->_materialStatus[ 1 ];
}
// ----------------------------------------------------------------------------
double AmorDamageModel::givePotentialFrom( const RealVector& conState, const MaterialStatus* matStatus )
{
    auto mst = this->accessConstMaterialStatus( matStatus );

    RealVector strain, phi;
    std::tie( strain, phi ) = this->retrieveStrainAndPhasefieldFrom( conState );

    // Strain decomposition
    RealVector volStrain;
    volStrain = _Pvol * strain;

    double degFcn = _degradationFunction->givePotentialFrom( phi, mst->_materialStatus[ 1 ] );
    double volEnergy = _elasticityModel->givePotentialFrom( strain, mst->_materialStatus[ 0 ], "Volumetric" );
    double devEnergy = _elasticityModel->givePotentialFrom( strain, mst->_materialStatus[ 0 ], "Deviatoric" );
    
    double potential;
    if ( volStrain( 0 ) > 0. )
        potential = ( _kT + ( 1. - _kT ) * degFcn ) * ( volEnergy + devEnergy );
    else
        potential = ( _kC + ( 1. - _kC ) * degFcn ) * devEnergy + volEnergy;
    
    return potential;
}
// ----------------------------------------------------------------------------
double AmorDamageModel::givePotentialFrom( const RealVector& conState, const MaterialStatus* matStatus, const std::string& label )
{
    auto mst = this->accessConstMaterialStatus( matStatus );

    RealVector strain, phi;
    std::tie( strain, phi ) = this->retrieveStrainAndPhasefieldFrom( conState );

    // Strain decomposition
    RealVector volStrain;
    volStrain = _Pvol * strain;

    double potential;
    if ( label == "Volumetric" )
    {
        double volEnergy = _elasticityModel->givePotentialFrom( strain, mst->_materialStatus[ 0 ], "Volumetric" );
        if ( volStrain( 0 ) > 0. )
        {
            double degFcn = _degradationFunction->givePotentialFrom( phi, mst->_materialStatus[ 1 ] );
            potential = ( _kT + ( 1. - _kT ) * degFcn ) * volEnergy;
        }
        else
            potential = volEnergy;
    }
    else if ( label == "Deviatoric" )
    {
        double degFcn = _degradationFunction->givePotentialFrom( phi, mst->_materialStatus[ 1 ] );
        double devEnergy = _elasticityModel->givePotentialFrom( strain, mst->_materialStatus[ 0 ], "Deviatoric" );
        if ( volStrain( 0 ) > 0. )
            potential = ( _kT + ( 1. - _kT ) * degFcn ) * devEnergy;
        else
            potential = ( _kC + ( 1. - _kC ) * degFcn ) * devEnergy;
    }

    return potential;
}
// ----------------------------------------------------------------------------
RealVector AmorDamageModel::giveForceFrom( const RealVector& conState, const MaterialStatus* matStatus, const std::string& label )
{
    RealVector conForce;
    
    auto mst = this->accessConstMaterialStatus(matStatus);

    RealVector strain, phi;
    std::tie( strain, phi ) = this->retrieveStrainAndPhasefieldFrom( conState );
    // Strain decomposition
    RealVector volStrain;
    volStrain = _Pvol * strain;

    if ( label == "Mechanics" )
    {
        RealVector volStress = _elasticityModel->giveForceFrom( strain, mst->_materialStatus[ 0 ], "Volumetric" );
        RealVector devStress = _elasticityModel->giveForceFrom( strain, mst->_materialStatus[ 0 ], "Deviatoric" );
        double degFcn = _degradationFunction->givePotentialFrom( phi, mst->_materialStatus[ 1 ] );
        
        if ( volStrain( 0 ) > 0. )
            conForce = ( _kT + ( 1. - _kT ) * degFcn ) * ( volStress + devStress );
        else
            conForce = ( _kC + ( 1. - _kC ) * degFcn ) * devStress + volStress;
    }
    else if ( label == "Mechanics_Volumetric" )
    {
        RealVector volStress = _elasticityModel->giveForceFrom( strain, mst->_materialStatus[ 0 ], "Volumetric" );
        if ( volStrain( 0 ) > 0. )
        {
            double degFcn = _degradationFunction->givePotentialFrom( phi, mst->_materialStatus[ 1 ] );
            conForce = ( _kT + ( 1. - _kT ) * degFcn ) * volStress;
        }
        else
            conForce = volStress;
    }
    else if ( label == "Mechanics_Deviatoric" )
    {
        RealVector devStress = _elasticityModel->giveForceFrom( strain, mst->_materialStatus[ 0 ], "Deviatoric" );
        double degFcn = _degradationFunction->givePotentialFrom( phi, mst->_materialStatus[ 1 ] );
        if ( volStrain( 0 ) > 0. )
            conForce = ( _kT + ( 1. - _kT ) * degFcn ) * devStress;
        else
            conForce = ( _kC + ( 1. - _kC ) * degFcn ) * devStress;
    }
    else if ( label == "PhaseField" )
    {
        RealVector DdegFcn = _degradationFunction->giveForceFrom( phi, mst->_materialStatus[ 1 ] );
        if ( volStrain( 0 ) > 0. )
            conForce = { ( 1. - _kT ) * DdegFcn( 0 ), mst->_volElasticEnergy + mst->_devElasticEnergy };
        else
            conForce = { ( 1. - _kC ) * DdegFcn( 0 ), mst->_devElasticEnergy };
    }
    else if ( label == "PhaseField_Volumetric" )
    {
        RealVector DdegFcn = _degradationFunction->giveForceFrom( phi, mst->_materialStatus[ 1 ] );
        if ( volStrain( 0 ) > 0. )
            conForce = { ( 1. - _kT ) * DdegFcn( 0 ), mst->_volElasticEnergy };
        else
            conForce = { ( 1. - _kC ) * DdegFcn( 0 ), 0. };
    }
    else if ( label == "PhaseField_Deviatoric" )
    {
        RealVector DdegFcn = _degradationFunction->giveForceFrom( phi, mst->_materialStatus[ 1 ] );
        if ( volStrain( 0 ) > 0. )
            conForce = { ( 1. - _kT ) * DdegFcn( 0 ), mst->_devElasticEnergy };
        else
            conForce = { ( 1. - _kC ) * DdegFcn( 0 ), mst->_devElasticEnergy };
    }
    else
        throw std::runtime_error( "Error: Invalid subsytem label '" + label +
                "' encountered in constitutive force calculation!\nSource: " + _name );
    
    return conForce;
}
// ----------------------------------------------------------------------------
RealMatrix AmorDamageModel::giveModulusFrom( const RealVector& conState, const MaterialStatus* matStatus, const std::string& label )
{
    RealMatrix conMod;
    
    auto mst = this->accessConstMaterialStatus( matStatus );

    RealVector strain, phi;
    std::tie( strain, phi ) = this->retrieveStrainAndPhasefieldFrom( conState );
    
    // Strain decomposition
    RealVector volStrain;
    volStrain = _Pvol * strain;

    if ( label == "Mechanics" )
    {
        RealMatrix volModulus = _elasticityModel->giveModulusFrom( strain, mst->_materialStatus[ 0 ], "Volumetric" );
        RealMatrix devModulus = _elasticityModel->giveModulusFrom( strain, mst->_materialStatus[ 0 ], "Deviatoric" );
        double degFcn = _degradationFunction->givePotentialFrom( phi, mst->_materialStatus[ 1 ] );
        
        if ( volStrain( 0 ) > 0. )
            conMod = ( _kT + ( 1. - _kT ) * degFcn ) * ( volModulus + devModulus );
        else
            conMod = ( _kC + ( 1. - _kC ) * degFcn ) * devModulus + volModulus;
    }
    else if ( label == "Mechanics_Volumetric" )
    {
        RealMatrix volModulus = _elasticityModel->giveModulusFrom( strain, mst->_materialStatus[ 0 ], "Volumetric" );

        if ( volStrain( 0 ) > 0. )
        {
            double degFcn = _degradationFunction->givePotentialFrom( phi, mst->_materialStatus[ 1 ] );
            conMod = ( _kT + ( 1. - _kT ) * degFcn ) * volModulus;
        }
        else
            conMod = volModulus;
    }
    else if ( label == "Mechanics_Deviatoric" )
    {
        RealMatrix devModulus = _elasticityModel->giveModulusFrom( strain, mst->_materialStatus[ 0 ], "Deviatoric" );
        double degFcn = _degradationFunction->givePotentialFrom( phi, mst->_materialStatus[ 1 ] );
        if ( volStrain( 0 ) > 0. )
            conMod = ( _kT + ( 1. - _kT ) * degFcn ) * devModulus;
        else
            conMod = ( _kC + ( 1. - _kC ) * degFcn ) * devModulus;
    }
    else if ( label == "PhaseField" )
    {
        RealMatrix DDdegFcn = _degradationFunction->giveModulusFrom( phi, mst->_materialStatus[ 1 ] );
        if ( volStrain( 0 ) > 0. )
            conMod = { { ( 1. - _kT ) * DDdegFcn( 0,0 ), 0. },
                       { 0., mst->_volElasticEnergy + mst->_devElasticEnergy } };
        else
            conMod = { { ( 1. - _kC ) * DDdegFcn( 0,0 ), 0. },
                       { 0., mst->_devElasticEnergy } };
    }
    else if ( label == "PhaseField_Volumetric" )
    {
        RealMatrix DDdegFcn = _degradationFunction->giveModulusFrom( phi, mst->_materialStatus[ 1 ] );
        if ( volStrain( 0 ) > 0. )
            conMod = { { ( 1. - _kT ) * DDdegFcn( 0,0 ), 0. },
                       { 0., mst->_volElasticEnergy } };
        else
            conMod = { { ( 1. - _kC ) * DDdegFcn( 0,0 ), 0. },
                       { 0., 0. } };
    }
    else if ( label == "PhaseField_Deviatoric" )
    {
        RealMatrix DDdegFcn = _degradationFunction->giveModulusFrom( phi, mst->_materialStatus[ 1 ] );
        if ( volStrain( 0 ) > 0. )
            conMod = { { ( 1. - _kT ) * DDdegFcn( 0,0 ), 0. },
                       { 0., mst->_devElasticEnergy } };
        else
            conMod = { { ( 1. - _kC ) * DDdegFcn( 0,0 ), 0. },
                       { 0., mst->_devElasticEnergy } };
    }
    else
        throw std::runtime_error( "Error: Invalid subsytem label '" + label +
                "' encountered in constitutive force calculation!\nSource: " + _name );
    
    return conMod;
}
// ----------------------------------------------------------------------------
void AmorDamageModel::readParamatersFrom( FILE* fp )
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
    else
        throw std::runtime_error( "Invalid analysis mode '" + mode + "' encountered in input file!\nSource: " + _name );

    std::string elasticityModel = getStringInputFrom( fp, "Failed to read elasticity model from input file!", _name );
    _elasticityModel = objectFactory().instantiateMaterial( elasticityModel );
    _elasticityModel->readParamatersFrom( fp );
    
    std::string degFcn = getStringInputFrom( fp, "Failed to read degradation function model from input file!", _name );
    _degradationFunction = objectFactory().instantiateMaterial( degFcn );
    _degradationFunction->readParamatersFrom( fp );

    verifyKeyword( fp, "Stabilization", _name );
    _kT = getRealInputFrom( fp, "Failed to read tension stabilization constant from input file!", _name );
    _kC = getRealInputFrom( fp, "Failed to read compression stabilization constant from input file!", _name );

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
    else
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
// ----------------------------------------------------------------------------
void AmorDamageModel::updateStatusFrom( const RealVector& conState, MaterialStatus* matStatus )
{
    RealVector strain, phi;

    std::tie( strain,phi ) = retrieveStrainAndPhasefieldFrom( conState );

    RealVector volStrain = _Pvol * strain;

    auto mst = this->accessMaterialStatus( matStatus );
    _elasticityModel->updateStatusFrom( strain, mst->_materialStatus[ 0 ] );
    _degradationFunction->updateStatusFrom( phi, mst->_materialStatus[ 1 ] );

    mst->_devElasticEnergy = _elasticityModel->givePotentialFrom( strain, mst->_materialStatus[ 0 ], "Deviatoric" );
    if ( volStrain( 0 ) > 0. )
        mst->_volElasticEnergy = _elasticityModel->givePotentialFrom( strain, mst->_materialStatus[ 0 ], "Volumetric" );
    else
        mst->_volElasticEnergy = 0.;
}
// ----------------------------------------------------------------------------
void AmorDamageModel::updateStatusFrom( const RealVector& conState, MaterialStatus* matStatus, const std::string& label )
{
    auto mst = this->accessMaterialStatus( matStatus );
    RealVector strain, phi;
    std::tie( strain, phi ) = retrieveStrainAndPhasefieldFrom( conState );

    // Strain decomposition
    RealVector volStrain, devStrain;
    volStrain = _Pvol * strain;

    if ( label == "Mechanics" )
        _elasticityModel->updateStatusFrom( strain, mst->_materialStatus[ 0 ] );
    else if ( label == "PhaseField" )
    {
        _elasticityModel->updateStatusFrom( strain, mst->_materialStatus[ 0 ] );
        _degradationFunction->updateStatusFrom( phi, mst->_materialStatus[ 1 ] );

        mst->_devElasticEnergy = _elasticityModel->givePotentialFrom( strain, mst->_materialStatus[ 0 ], "Deviatoric" );
        if ( volStrain( 0 ) > 0. )
            mst->_volElasticEnergy = _elasticityModel->givePotentialFrom( strain, mst->_materialStatus[ 0 ], "Volumetric" );
        else
            mst->_volElasticEnergy = 0.;
    }
    else
        throw std::runtime_error( "Error: Invalid size for constitutive state vector used for material status update!\nSource: " + _name );
}

// Private methods
// ----------------------------------------------------------------------------
MaterialStatus_AmorDamageModel*
AmorDamageModel::accessMaterialStatus( MaterialStatus* matStatus )
{
    auto mst = dynamic_cast<MaterialStatus_AmorDamageModel*>( matStatus );
    if ( !mst )
        throw std::runtime_error( "Error: Unable to access material status!\nSource: " + _name );
    
    return mst;
}
// ----------------------------------------------------------------------------
const MaterialStatus_AmorDamageModel*
AmorDamageModel::accessConstMaterialStatus( const MaterialStatus* matStatus )
{
    auto mst = dynamic_cast<const MaterialStatus_AmorDamageModel*>( matStatus );
    if ( !mst )
        throw std::runtime_error( "Error: Unable to access material status!\nSource: " + _name );
    
    return mst;
}
// ----------------------------------------------------------------------------
std::tuple< RealVector, RealVector >
AmorDamageModel::retrieveStrainAndPhasefieldFrom( const broomstyx::RealVector& conState )
{
    RealVector strain, phi;
    if ( _analysisMode == TwoDimensional || _analysisMode == PlaneStress )
    {
        strain = { conState( 0 ), conState( 1 ), conState( 2 ) };
        phi = { conState( 3 ) };
    }
    else if ( _analysisMode == PlaneStrain || _analysisMode == Axisymmetric )
    {
        strain = { conState( 0 ), conState( 1 ), conState( 2 ), conState( 3 ) };
        phi = { conState( 4 ) };
    }
    else
    {
        strain = { conState( 0 ), conState( 1 ), conState( 2 ), conState( 3 ), conState( 4 ), conState( 5 ) };
        phi = { conState( 6 ) };
    }

    return std::make_tuple( std::move( strain), phi );
}