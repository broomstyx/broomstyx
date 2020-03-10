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

registerBroomstyxObject(Material, AmorDamageModel)

// Material status
MaterialStatus_AmorDamageModel::MaterialStatus_AmorDamageModel()
    : _materialStatus {nullptr, nullptr}
    , _historyField(0.)
    , _elasticEnergy(0.)
{}

MaterialStatus_AmorDamageModel::~MaterialStatus_AmorDamageModel() {}

// Constructor
AmorDamageModel::AmorDamageModel()
    : _elasticityModel(nullptr)
    , _degradationFunction(nullptr)
    , _phiIrrev(0.)
    , _I({{1., 0., 0., 0.},
          {0., 1., 0., 0.},
          {0., 0., 1., 0.},
          {0., 0., 0., 1.}})
    , _P({{0.5, 0.5, 0., 0.},
          {0.5, 0.5, 0., 0.},
          {0., 0., 0., 0.},
          {0., 0., 0., 0.}})
{}

// Destructor
AmorDamageModel::~AmorDamageModel() {}

// Public methods
// ----------------------------------------------------------------------------
MaterialStatus* AmorDamageModel::createMaterialStatus()
{
    MaterialStatus* matStatus = new MaterialStatus_AmorDamageModel();
    auto mst = this->accessMaterialStatus(matStatus);
    
    mst->_materialStatus[0] = _elasticityModel->createMaterialStatus();
    mst->_materialStatus[1] = _degradationFunction->createMaterialStatus();
    
    return matStatus;
}
// ----------------------------------------------------------------------------
void AmorDamageModel::destroy( MaterialStatus*& matStatus )
{
    auto mst = this->accessMaterialStatus(matStatus);
    delete mst->_materialStatus[0];
    delete mst->_materialStatus[1];
}
// ----------------------------------------------------------------------------
double AmorDamageModel::givePotentialFrom( const RealVector& conState, const MaterialStatus* matStatus )
{
    auto mst = this->accessConstMaterialStatus(matStatus);
    
    // Split constitutive state
    RealVector strain({conState(0), conState(1), conState(2), conState(3)});
    RealVector phi({conState(4)});
    
    // Strain decomposition
    RealVector volStrain, devStrain;
    volStrain = _P*strain;
    devStrain = strain - volStrain;
    
    double degFcn = _degradationFunction->givePotentialFrom(phi, mst->_materialStatus[1]);
    double volEnergy = _elasticityModel->givePotentialFrom(volStrain, mst->_materialStatus[0]);
    double devEnergy = _elasticityModel->givePotentialFrom(devStrain, mst->_materialStatus[0]);
    
    double potential;
    if ( strain(0) + strain(1) > 0. )
        potential = degFcn*(volEnergy + devEnergy);
    else
        potential = degFcn*devEnergy + volEnergy;
    
    return potential;
}
// ----------------------------------------------------------------------------
RealVector AmorDamageModel::giveForceFrom( const RealVector&     conState
                                         , const MaterialStatus* matStatus
                                         , const std::string&    label )
{
    RealVector conForce;
    
    auto mst = this->accessConstMaterialStatus(matStatus);
    
    // Split constitutive state
    RealVector strain({conState(0), conState(1), conState(2), conState(3)});
    RealVector phi({conState(4)});
    
    // Strain decomposition
    RealVector volStrain, devStrain;
    volStrain = _P*strain;
    devStrain = strain - volStrain;
    
    if ( label == "Mechanics" )
    {
        RealVector volStress = _elasticityModel->giveForceFrom(volStrain, mst->_materialStatus[0]);
        RealVector devStress = _elasticityModel->giveForceFrom(devStrain, mst->_materialStatus[0]);
        double degFcn = _degradationFunction->givePotentialFrom(phi, mst->_materialStatus[1]);
        
        if ( strain(0) + strain(1) > 0. )
            conForce = degFcn*(volStress + devStress);
        else
            conForce = degFcn*devStress + volStress;
    }
    else if ( label == "PhaseField" )
    {
        RealVector DdegFcn = _degradationFunction->giveForceFrom(phi, mst->_materialStatus[1]);
        conForce = DdegFcn*mst->_elasticEnergy;
    }
    else
        throw std::runtime_error("Error: Invalid subsytem label '" + label +
                "' encountered in constitutive force calculation!\nSource: " + _name);
    
    return conForce;
}
// ----------------------------------------------------------------------------
RealMatrix AmorDamageModel::giveModulusFrom( const RealVector&     conState
                                           , const MaterialStatus* matStatus
                                           , const std::string&    label )
{
    RealMatrix conMod;
    
    auto mst = this->accessConstMaterialStatus(matStatus);
    
    // Split constitutive state
    RealVector strain({conState(0), conState(1), conState(2), conState(3)});
    RealVector phi({conState(4)});
    
    // Strain decomposition
    RealVector volStrain, devStrain;
    volStrain = _P*strain;
    devStrain = strain - volStrain;
    
    if ( label == "Mechanics" )
    {
        RealMatrix volModulus = _elasticityModel->giveModulusFrom(volStrain, mst->_materialStatus[0]);
        RealMatrix devModulus = _elasticityModel->giveModulusFrom(devStrain, mst->_materialStatus[0]);
        double degFcn = _degradationFunction->givePotentialFrom(phi, mst->_materialStatus[1]);
        
        if ( strain(0) + strain(1) > 0. )
            conMod = degFcn*(volModulus*_P + devModulus*(_I - _P));
        else
            conMod = degFcn*devModulus*(_I - _P) + volModulus*_P;
    }
    else if ( label == "PhaseField" )
    {
        RealMatrix DDdegFcn = _degradationFunction->giveModulusFrom(phi, mst->_materialStatus[1]);
        conMod = DDdegFcn*mst->_elasticEnergy;
    }
    else
        throw std::runtime_error("Error: Invalid subsytem label '" + label +
                "' encountered in constitutive force calculation!\nSource: " + _name);
    
    return conMod;
}
// ----------------------------------------------------------------------------
void AmorDamageModel::readParamatersFrom( FILE* fp )
{
    std::string elasticityModel = getStringInputFrom(fp, "Failed to read elasticity model from input file!", _name);
    _elasticityModel = objectFactory().instantiateMaterial(elasticityModel);
    _elasticityModel->readParamatersFrom(fp);
    
    std::string degFcn = getStringInputFrom(fp, "Failed to read degradation function model from input file!", _name);
    _degradationFunction = objectFactory().instantiateMaterial(degFcn);
    _degradationFunction->readParamatersFrom(fp);
    
    verifyKeyword(fp, "IrreversibilityThreshold", _name);
    _phiIrrev = getRealInputFrom(fp, "Failed to read irreversibility threshold from input file!", _name);
}
// ----------------------------------------------------------------------------
void AmorDamageModel::updateStatusFrom( const RealVector& conState, MaterialStatus* matStatus )
{
    if ( conState.dim() == 5 )
    {
        RealVector strain({conState(0), conState(1), conState(2), conState(3)});
        RealVector phi({conState(4)});
        
        // Strain decomposition
        RealVector volStrain, devStrain;
        volStrain = _P*strain;
        devStrain = strain - volStrain;
        
        auto mst = this->accessMaterialStatus(matStatus);
        _elasticityModel->updateStatusFrom(strain, mst->_materialStatus[0]);
        _degradationFunction->updateStatusFrom(phi, mst->_materialStatus[1]);
        
        double volEnergy, devEnergy;
        volEnergy = _elasticityModel->givePotentialFrom(volStrain, mst->_materialStatus[0]);
        devEnergy = _elasticityModel->givePotentialFrom(devStrain, mst->_materialStatus[0]);
        
        if ( strain(0) + strain(1) > 0. )
            mst->_elasticEnergy = volEnergy + devEnergy;
        else
            mst->_elasticEnergy = devEnergy;
        
        // Enforce irreversibility
        if ( mst->_elasticEnergy < mst->_historyField && phi(0) > _phiIrrev )
            mst->_elasticEnergy = mst->_historyField;
    }
    else
        throw std::runtime_error("Error: Invalid size for constitutive state vector used for material status update!\nSource: " + _name);
}
// ----------------------------------------------------------------------------
void AmorDamageModel::updateStatusFrom( const RealVector& conState, MaterialStatus* matStatus, const std::string& label )
{
    auto mst = this->accessMaterialStatus(matStatus);
    
    if ( label == "Mechanics" )
    {
        RealVector strain({conState(0), conState(1), conState(2), conState(3)});
        _elasticityModel->updateStatusFrom(strain, mst->_materialStatus[0]);
    }
    else if ( label == "PhaseField" )
    {
        RealVector strain({conState(0), conState(1), conState(2), conState(3)});
        RealVector phi({conState(4)});
        _degradationFunction->updateStatusFrom(phi, mst->_materialStatus[1]);
        
        // Strain decomposition
        RealVector volStrain, devStrain;
        volStrain = _P*strain;
        devStrain = strain - volStrain;
        
        double volEnergy, devEnergy;
        volEnergy = _elasticityModel->givePotentialFrom(volStrain, mst->_materialStatus[0]);
        devEnergy = _elasticityModel->givePotentialFrom(devStrain, mst->_materialStatus[0]);
        
        if ( strain(0) + strain(1) > 0. )
            mst->_elasticEnergy = volEnergy + devEnergy;
        else
            mst->_elasticEnergy = devEnergy;
        
        // Enforce irreversibility
        if ( mst->_elasticEnergy < mst->_historyField && phi(0) > _phiIrrev )
            mst->_elasticEnergy = mst->_historyField;
    }
    else if ( label == "InitializeHistoryField" )
    {
        /***********************************
         *
         *  conState(0) = phi
         *  conState(1) = g'(phi)*Psi0
         * 
         ***********************************/
        
        auto mst = this->accessMaterialStatus(matStatus);
        RealVector phi({conState(0)});
        RealVector Dgphi = _degradationFunction->giveForceFrom(phi, mst->_materialStatus[1]);
        
        // Set history field
        mst->_historyField = conState(1)/Dgphi(0);
    }
    else if ( label == "FinalizeHistoryField" )
    {
        RealVector strain({conState(0), conState(1), conState(2), conState(3)});
        RealVector phi({conState(4)});
        
        auto mst = this->accessMaterialStatus(matStatus);
        
        // Strain decomposition
        RealVector volStrain, devStrain;
        volStrain = _P*strain;
        devStrain = strain - volStrain;
        
        double volEnergy, devEnergy;
        volEnergy = _elasticityModel->givePotentialFrom(volStrain, mst->_materialStatus[0]);
        devEnergy = _elasticityModel->givePotentialFrom(devStrain, mst->_materialStatus[0]);
        
        if ( strain(0) + strain(1) > 0. )
            mst->_elasticEnergy = volEnergy + devEnergy;
        else
            mst->_elasticEnergy = devEnergy;
        
        if ( mst->_elasticEnergy > mst->_historyField && phi(0) > _phiIrrev )
            mst->_historyField = mst->_elasticEnergy;
    }
    else
        throw std::runtime_error("Error: Invalid size for constitutive state vector used for material status update!\nSource: " + _name);
}

// Private methods
// ----------------------------------------------------------------------------
MaterialStatus_AmorDamageModel*
AmorDamageModel::accessMaterialStatus( MaterialStatus* matStatus )
{
    auto mst = dynamic_cast<MaterialStatus_AmorDamageModel*>(matStatus);
    if ( !mst )
        throw std::runtime_error("Error: Unable to access material status!\nSource: " + _name);
    
    return mst;
}
// ----------------------------------------------------------------------------
const MaterialStatus_AmorDamageModel*
AmorDamageModel::accessConstMaterialStatus( const MaterialStatus* matStatus )
{
    auto mst = dynamic_cast<const MaterialStatus_AmorDamageModel*>(matStatus);
    if ( !mst )
        throw std::runtime_error("Error: Unable to access material status!\nSource: " + _name);
    
    return mst;
}