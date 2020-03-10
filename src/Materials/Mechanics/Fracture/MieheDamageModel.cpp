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

#include "MieheDamageModel.hpp"
#include "Core/ObjectFactory.hpp"
#include "Util/readOperations.hpp"
#include "Util/linearAlgebra.hpp"
#include <cmath>

#define PI 3.14159265358979323846
#define ZERO_THRESHOLD 1.0e-14

using namespace broomstyx;

registerBroomstyxObject(Material, MieheDamageModel)

// Material status
MaterialStatus_MieheDamageModel::MaterialStatus_MieheDamageModel()
    : _materialStatus(nullptr)
    , _prinStrain({ZERO_THRESHOLD, -ZERO_THRESHOLD, 0})
    , _M({{1, 0, 0}, {0, 1, 0}, {0, 0, 1}})
    , _historyField(0.)
    , _drivingForce(0.)
{}

MaterialStatus_MieheDamageModel::~MaterialStatus_MieheDamageModel() {}

// Constructor
MieheDamageModel::MieheDamageModel()
{}

// Destructor
MieheDamageModel::~MieheDamageModel() {}

// Public methods
// ----------------------------------------------------------------------------
MaterialStatus* MieheDamageModel::createMaterialStatus()
{
    MaterialStatus* matStatus = new MaterialStatus_MieheDamageModel();
    auto mst = this->accessMaterialStatus(matStatus);    
    mst->_materialStatus = _degradationFunction->createMaterialStatus();
    return matStatus;
}
// ----------------------------------------------------------------------------
void MieheDamageModel::destroy( MaterialStatus*& matStatus )
{
    auto mst = this->accessMaterialStatus(matStatus);
    delete mst->_materialStatus;
}
// ----------------------------------------------------------------------------
double MieheDamageModel::givePotentialFrom( const RealVector& conState, const MaterialStatus* matStatus )
{
    auto mst = this->accessConstMaterialStatus(matStatus);
    RealVector strain, phi;
    std::tie(strain, phi) = this->getStrainAndPhaseFieldFrom(conState);
        
    double degFcn = _degradationFunction->givePotentialFrom(phi, mst->_materialStatus);
    double potential = 0.;
    double volStrain = mst->_prinStrain(0) + mst->_prinStrain(1) + mst->_prinStrain(2);
    
    if ( volStrain > ZERO_THRESHOLD )
        potential += 0.5*degFcn*_lambda*volStrain*volStrain;
    else
        potential += 0.5*_lambda*volStrain*volStrain;
    
    for ( int i = 0; i < 3; i++ )
        if ( mst->_prinStrain(i) > ZERO_THRESHOLD )
            potential += degFcn*_G*mst->_prinStrain(i)*mst->_prinStrain(i);
        else
            potential += _G*mst->_prinStrain(i)*mst->_prinStrain(i);
    
    return potential;
}
// ----------------------------------------------------------------------------
RealVector MieheDamageModel::giveForceFrom( const RealVector&     conState
                                          , const MaterialStatus* matStatus
                                          , const std::string&    label )
{
    RealVector conForce;

    RealVector strain, phi;
    std::tie(strain, phi) = this->getStrainAndPhaseFieldFrom(conState);
    
    auto mst = this->accessConstMaterialStatus(matStatus);
    if ( label == "Mechanics" )
    {
        // Calculate stress components in principal orientation
        double degFcn = _degradationFunction->givePotentialFrom(phi, mst->_materialStatus);
        double epsV = mst->_prinStrain(0) + mst->_prinStrain(1) + mst->_prinStrain(2);

        double volComp;
        if ( epsV > ZERO_THRESHOLD )
            volComp = degFcn*_lambda*epsV;
        else
            volComp = _lambda*epsV;
        
        RealMatrix sigmaP({{volComp, 0., 0.},
                           {0., volComp, 0.},
                           {0., 0., volComp}});
        
        for ( int i = 0; i < 3; i++ )
            if ( mst->_prinStrain(i) > ZERO_THRESHOLD )
                sigmaP(i,i) += 2.*_G*degFcn*mst->_prinStrain(i);
            else
                sigmaP(i,i) += 2.*_G*mst->_prinStrain(i);
        
        // Rotate to x-y axes
        RealMatrix Rmat({{mst->_M(0,0), mst->_M(0,1), 0.},
                         {-mst->_M(0,1), mst->_M(0,0), 0.},
                         {0., 0., 1.}});

        RealMatrix sigma = Rmat*sigmaP*trp(Rmat);

        conForce = {sigma(0,0), sigma(1,1), sigma(2,2), sigma(0,1)};
    }
    else if ( label == "PhaseField" )
    {
        RealVector DdegFcn = _degradationFunction->giveForceFrom(phi, mst->_materialStatus);
        conForce = DdegFcn*mst->_drivingForce;
    }
    else
        throw std::runtime_error("Error: Invalid label '" + label + "' encountered in constitutive force calculation!\nSource: " + _name);
    
    return conForce;
}
// ----------------------------------------------------------------------------
double MieheDamageModel::giveMaterialVariable( const std::string& label, const MaterialStatus* matStatus )
{
    auto mst = this->accessConstMaterialStatus(matStatus);

    if ( label == "eigVec1_1" )
        return mst->_M(0,0);
    else if ( label == "eigVec1_2" )
        return mst->_M(1,0);
    else if ( label == "eigVec2_1" )
        return mst->_M(0,1);
    else if ( label == "eigVec2_2" )
        return mst->_M(1,1);
    else
        throw std::runtime_error("Invalid tag '" + label + "' supplied in field output request made to material '" + _name + "'!");
}
// ----------------------------------------------------------------------------
RealMatrix MieheDamageModel::giveModulusFrom( const RealVector&     conState
                                            , const MaterialStatus* matStatus
                                            , const std::string&    label )
{
    RealMatrix conMod;
    
    RealVector strain, phi;
    std::tie(strain, phi) = this->getStrainAndPhaseFieldFrom(conState);

    auto mst = this->accessConstMaterialStatus(matStatus);

    if ( label == "Mechanics" )
    {
        if ( _analysisMode == 2 /* Plane strain */)
        {
            conMod.init(4,4);
            
            RealMatrix M = mst->_M;
            double degFcn = _degradationFunction->givePotentialFrom(phi, mst->_materialStatus);
            double epsV = mst->_prinStrain(0) + mst->_prinStrain(1) + mst->_prinStrain(2);

            double coefA;
            if ( epsV > ZERO_THRESHOLD )
                coefA = degFcn*_lambda;
            else
                coefA = _lambda;

            for ( int a = 0; a < 3; a++ )
            {
                double coefB;
                if ( mst->_prinStrain(a) > ZERO_THRESHOLD )
                    coefB = 2.*degFcn*_G;
                else
                    coefB = 2.*_G;
                
                for ( int b = 0; b < 3; b++ )
                {
                    RealMatrix nprod({{M(0,a)*M(0,a)*M(0,b)*M(0,b), M(0,a)*M(0,a)*M(1,b)*M(1,b), M(0,a)*M(0,a)*M(2,b)*M(2,b), M(0,a)*M(0,a)*M(0,b)*M(1,b)},
                                      {M(1,a)*M(1,a)*M(0,b)*M(0,b), M(1,a)*M(1,a)*M(1,b)*M(1,b), M(1,a)*M(1,a)*M(2,b)*M(2,b), M(1,a)*M(1,a)*M(0,b)*M(1,b)},
                                      {M(2,a)*M(2,a)*M(0,b)*M(0,b), M(2,a)*M(2,a)*M(1,b)*M(1,b), M(2,a)*M(2,a)*M(2,b)*M(2,b), M(2,a)*M(2,a)*M(0,b)*M(1,b)},
                                      {M(0,a)*M(1,a)*M(0,b)*M(0,b), M(0,a)*M(1,a)*M(1,b)*M(1,b), M(0,a)*M(1,a)*M(2,b)*M(2,b), M(0,a)*M(1,a)*M(0,b)*M(1,b)}});
                    
                    if ( a == b )
                        conMod += (coefA + coefB)*nprod;
                    else
                        conMod += coefA*nprod;
                }
            }

            /* Loop index for a and b run from 0--1 instead of 0--2 since for plane strain the eigenvector associated
               with the 3rd principal strain (i.e. {0, 0, 1}) results in zero contributions when a = 2 or b = 2.
            */
            for ( int a = 0; a < 2; a++ ) 
            {
                for ( int b = 0; b < 2; b++ )
                    if ( a != b )
                    {
                        double coefC = 0.;
                        if ( mst->_prinStrain(a) > ZERO_THRESHOLD )
                            coefC += 2.*degFcn*_G*mst->_prinStrain(a)/(mst->_prinStrain(a) - mst->_prinStrain(b));
                        else
                            coefC += 2.*_G*mst->_prinStrain(a)/(mst->_prinStrain(a) - mst->_prinStrain(b));

                        if ( mst->_prinStrain(b) > ZERO_THRESHOLD )
                            coefC -= 2.*degFcn*_G*mst->_prinStrain(b)/(mst->_prinStrain(a) - mst->_prinStrain(b));
                        else
                            coefC -= 2.*_G*mst->_prinStrain(b)/(mst->_prinStrain(a) - mst->_prinStrain(b));

                        RealMatrix nprod1({{M(0,a)*M(0,b)*M(0,a)*M(0,b), M(0,a)*M(0,b)*M(1,a)*M(1,b), M(0,a)*M(0,b)*M(2,a)*M(2,b), M(0,a)*M(0,b)*M(0,a)*M(1,b)},
                                           {M(1,a)*M(1,b)*M(0,a)*M(0,b), M(1,a)*M(1,b)*M(1,a)*M(1,b), M(1,a)*M(1,b)*M(2,a)*M(2,b), M(1,a)*M(1,b)*M(0,a)*M(1,b)},
                                           {M(2,a)*M(2,b)*M(0,a)*M(0,b), M(2,a)*M(2,b)*M(1,a)*M(1,b), M(2,a)*M(2,b)*M(2,a)*M(2,b), M(2,a)*M(2,b)*M(0,a)*M(1,b)},
                                           {M(0,a)*M(1,b)*M(0,a)*M(0,b), M(0,a)*M(1,b)*M(1,a)*M(1,b), M(0,a)*M(1,b)*M(2,a)*M(2,b), M(0,a)*M(1,b)*M(0,a)*M(1,b)}});
                        
                        RealMatrix nprod2({{M(0,a)*M(0,b)*M(0,b)*M(0,a), M(0,a)*M(0,b)*M(1,b)*M(1,a), M(0,a)*M(0,b)*M(2,b)*M(2,a), M(0,a)*M(0,b)*M(0,b)*M(1,a)},
                                           {M(1,a)*M(1,b)*M(0,b)*M(0,a), M(1,a)*M(1,b)*M(1,b)*M(1,a), M(1,a)*M(1,b)*M(2,b)*M(2,a), M(1,a)*M(1,b)*M(0,b)*M(1,a)},
                                           {M(2,a)*M(2,b)*M(0,b)*M(0,a), M(2,a)*M(2,b)*M(1,b)*M(1,a), M(2,a)*M(2,b)*M(2,b)*M(2,a), M(2,a)*M(2,b)*M(0,b)*M(1,a)},
                                           {M(0,a)*M(1,b)*M(0,b)*M(0,a), M(0,a)*M(1,b)*M(1,b)*M(1,a), M(0,a)*M(1,b)*M(2,b)*M(2,a), M(0,a)*M(1,b)*M(0,b)*M(1,a)}});

                        conMod += 0.5*coefC*(nprod1 + nprod2);
                    }
            }
        }
        else
            throw std::runtime_error("Error: Cannot calculate constitutive modulus for current analysis mode!\nSource: " + _name);
    }
    else if ( label == "PhaseField" )
    {
        RealMatrix DDdegFcn = _degradationFunction->giveModulusFrom(phi, mst->_materialStatus);
        conMod = DDdegFcn*mst->_drivingForce;
    }
    else
        throw std::runtime_error("Error: Invalid label '" + label + "' encountered in constitutive modulus calculation!\nSource: " + _name);
    
    return conMod;
}
// ----------------------------------------------------------------------------
void MieheDamageModel::readParamatersFrom( FILE* fp )
{
    std::string mode = getStringInputFrom(fp, "Failed to read analysis mode from input file!", _name);
    if ( mode == "PlaneStrain" )
        _analysisMode = 2;
    else
        throw std::runtime_error("ERROR: Invalid value '" + mode + "' specified for analysis mode in input file!\nSource: " + _name);
    
    _E = getRealInputFrom(fp, "Failed to read Young's modulus from input file!", _name);
    _nu = getRealInputFrom(fp, "Failed to read Poisson's ratio from input file!", _name);
    _lambda = _nu*_E/((1. + _nu)*(1. - 2*_nu));
    _G = _E/(2.*(1. + _nu));
    
    std::string degFcn = getStringInputFrom(fp, "Failed to read degradation function model from input file!", _name);
    _degradationFunction = objectFactory().instantiateMaterial(degFcn);
    _degradationFunction->readParamatersFrom(fp);
    
    verifyKeyword(fp, "IrreversibilityThreshold", _name);
    _phiIrrev = getRealInputFrom(fp, "Failed to read irreversibility threshold from input file!", _name);
}
// ----------------------------------------------------------------------------
void MieheDamageModel::updateStatusFrom( const RealVector& conState, MaterialStatus* matStatus )
{
    RealVector strain, phi;
    std::tie(strain, phi) = this->getStrainAndPhaseFieldFrom(conState);

    auto mst = this->accessMaterialStatus(matStatus);

    // Spectral decomposition
    if ( _analysisMode == 2 /* Plane strain */ )
    {
        double epsC = 0.5*(strain(0) + strain(1));
        double epsD = strain(0) - strain(1);

        // Perturb in case of equal principal strains
        if ( std::fabs(epsD) < ZERO_THRESHOLD )
            epsD = ZERO_THRESHOLD;

        // Eigenvectors
        double theta1 = 0.5*atan(strain(3)/epsD);
        double theta2 = theta1 + PI/2;
        mst->_M = {{ cos(theta1), cos(theta2), 0. },
                   { sin(theta1), sin(theta2), 0. },
                   { 0.,          0.,          1. }};
        
        // Eigenvalues
        mst->_prinStrain = { epsC + 0.5*epsD*cos(2.*theta1) + 0.5*strain(3)*sin(2.*theta1),
                             epsC + 0.5*epsD*cos(2.*theta2) + 0.5*strain(3)*sin(2.*theta2),
                             0. };
    }
    else
        throw std::runtime_error("Error: Cannot update material status for current analysis mode!\nSource: " + _name);

    _degradationFunction->updateStatusFrom(phi, mst->_materialStatus);
    
    // Phase-field driving force
    double volStrain = mst->_prinStrain(0) + mst->_prinStrain(1) + mst->_prinStrain(2);
    mst->_drivingForce = 0.;

    if ( volStrain > ZERO_THRESHOLD )
        mst->_drivingForce += 0.5*_lambda*volStrain*volStrain;
    
    for ( int i = 0; i < 3; i++ )
        if ( mst->_prinStrain(i) > ZERO_THRESHOLD )
            mst->_drivingForce += _G*mst->_prinStrain(i)*mst->_prinStrain(i);
    
    // Enforce irreversibility
    if ( mst->_drivingForce < mst->_historyField && phi(0) > _phiIrrev )
        mst->_drivingForce = mst->_historyField;
}
// ----------------------------------------------------------------------------
void MieheDamageModel::updateStatusFrom( const RealVector& conState, MaterialStatus* matStatus, const std::string& label )
{
    auto mst = this->accessMaterialStatus(matStatus);
    
    if ( label == "InitializeHistoryField" )
    {
        /***********************************
         *
         *  conState(0) = phi
         *  conState(1) = g'(phi)*Psi0
         * 
         ***********************************/
        
        RealVector phi({conState(0)});
        RealVector Dgphi = _degradationFunction->giveForceFrom(phi, mst->_materialStatus);
        
        // Set history field
        mst->_historyField = conState(1)/Dgphi(0);
    }
    else if ( label == "FinalizeHistoryField" )
    {
        RealVector strain, phi;
        std::tie(strain, phi) = this->getStrainAndPhaseFieldFrom(conState);

        _degradationFunction->updateStatusFrom(phi, mst->_materialStatus);
            
        // Phase-field driving force
        double volStrain = mst->_prinStrain(0) + mst->_prinStrain(1) + mst->_prinStrain(2);
        mst->_drivingForce = 0.;

        if ( volStrain > ZERO_THRESHOLD )
            mst->_drivingForce += 0.5*_lambda*volStrain*volStrain;
        
        for ( int i = 0; i < 3; i++ )
            if ( mst->_prinStrain(i) > ZERO_THRESHOLD )
                mst->_drivingForce += _G*mst->_prinStrain(i)*mst->_prinStrain(i);
        
        // Enforce irreversibility
        if ( mst->_drivingForce > mst->_historyField && phi(0) > _phiIrrev )
            mst->_historyField = mst->_drivingForce;
    }
    else
        throw std::runtime_error("Error: Invalid label '" + label + "' encountered during material status update!\nSource: " + _name);
}