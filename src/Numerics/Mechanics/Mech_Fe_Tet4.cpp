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

#include "Mech_Fe_Tet4.hpp"
#include <cmath>
#include <string>
#include <vector>
#include "Core/AnalysisModel.hpp"
#include "Core/Cell.hpp"
#include "Core/ObjectFactory.hpp"
#include "Core/DomainManager.hpp"
#include "Materials/Material.hpp"
#include "Util/linearAlgebra.hpp"
#include "BasisFunctions/Triangle_P1.hpp"

using namespace broomstyx;

registerBroomstyxObject(Numerics, Mech_Fe_Tet4)

// Numerics status
NumericsStatus_Mech_Fe_Tet4::NumericsStatus_Mech_Fe_Tet4()
    : _strain(RealVector(6))
    , _stress(RealVector(6))
    , _gradU(RealMatrix(3,3))
    , _materialStatus {nullptr, nullptr}
{}

NumericsStatus_Mech_Fe_Tet4::~NumericsStatus_Mech_Fe_Tet4()
{}

// Constructor
Mech_Fe_Tet4::Mech_Fe_Tet4()
{   
    _dim = 3;
    _dofPerNode = 3;
    _nNodes = 4;
    _nMaterials = 2;
    _nStages = 1;
    _nSubsystems = 1;
    _name = "Mech_Fe_Tet4";
    
    // Lone Gauss point coordinate is at natural coordinate (1/4, 1/4, 1/4), weight = 1/6
    _gpNatCoor = {0.25, 0.25, 0.25};
    _wt = 1./6.;
    
    // Pre-calculate shape functions and derivatives at Gauss point
    _basisFunctionValues = _basisFunction.giveBasisFunctionsAt(_gpNatCoor);
    std::vector<RealVector> dpsiNat = _basisFunction.giveBasisFunctionDerivativesAt(_gpNatCoor);
    _basisFunctionDerivatives.init(3,4);
    for ( int i = 0; i < 4; i++)
    {
        _basisFunctionDerivatives(0,i) = dpsiNat[0](i);
        _basisFunctionDerivatives(1,i) = dpsiNat[1](i);
        _basisFunctionDerivatives(2,i) = dpsiNat[2](i);
    }
}

// Destructor
Mech_Fe_Tet4::~Mech_Fe_Tet4()
{}

// Private methods
// ----------------------------------------------------------------------------
void Mech_Fe_Tet4::deleteNumericsAt(Cell* targetCell)
{
    auto cns = this->getNumericsStatusAt(targetCell);
    auto material = this->giveMaterialSetFor(targetCell);
    
    material[0]->destroy(cns->_materialStatus[0]);
    material[1]->destroy(cns->_materialStatus[1]);
    
    delete cns;
}
// ----------------------------------------------------------------------------
void Mech_Fe_Tet4::finalizeDataAt( Cell* targetCell )
{
    // Retrieve numerics status at cell
    auto cns = this->getNumericsStatusAt(targetCell);
    
    // Retrieve nodal DOFs local to element
    std::vector<Dof*> dof = this->giveNodalDofsAt(targetCell);
    
    // Retrieve DOF values
    RealVector uVec = this->giveLocalDisplacementsAt(dof, converged_value);
    
    // Construct displacement gradient
    RealMatrix bmatGradU, uMat(4,3);
    
    uMat(0,0) = uVec(0); uMat(0,1) = uVec(1);  uMat(0,2) = uVec(2);
    uMat(1,0) = uVec(3); uMat(1,1) = uVec(4);  uMat(1,2) = uVec(5);
    uMat(2,0) = uVec(6); uMat(2,1) = uVec(7);  uMat(2,2) = uVec(8);
    uMat(3,0) = uVec(9); uMat(3,1) = uVec(10); uMat(3,2) = uVec(11);
    bmatGradU = this->giveGradBmatAt(targetCell);
    cns->_gradU = bmatGradU*uMat;
    
    // Construct strain
    RealMatrix bmat = this->giveBmatAt(targetCell);
    cns->_strain = bmat*uVec;
    
    // Retrieve material set for element
    std::vector<Material*> material = this->giveMaterialSetFor(targetCell);

    // Get stress
    cns->_stress = material[1]->giveForceFrom(cns->_strain, cns->_materialStatus[1]);
}
// ----------------------------------------------------------------------------
double Mech_Fe_Tet4::giveCellFieldValueAt( Cell* targetCell, int fieldNum )
{
    RealVector val;
    std::string fieldTag;
    
    try
    {
        fieldTag = _cellFieldOutput.at(fieldNum);
    }
    catch (std::exception& e)
    {
        fieldTag = "unassigned";
    }
    
    std::tie(val, std::ignore) = this->giveFieldOutputAt(targetCell, fieldTag);
    
    return val(0);
}
// ----------------------------------------------------------------------------
RealVector Mech_Fe_Tet4::giveCellNodeFieldValuesAt( Cell* targetCell, int fieldNum )
{
    double val = this->giveCellFieldValueAt(targetCell, fieldNum);
    return RealVector({val, val, val, val});
}
// ----------------------------------------------------------------------------
std::vector<RealVector> Mech_Fe_Tet4::giveEvaluationPointsFor( Cell *targetCell)
{
    std::vector<Node*> node = analysisModel().domainManager().giveNodesOf(targetCell);
    
    std::vector<RealVector> coor(1, RealVector());
    coor[0].init(3);
    
    for ( int i = 0; i < 4; i++)
    {
        RealVector nodeCoor = analysisModel().domainManager().giveCoordinatesOf(node[i]);
        coor[0](0) += nodeCoor(0)/4.0;
        coor[0](1) += nodeCoor(1)/4.0;
        coor[0](2) += nodeCoor(2)/4.0;
    }
    
    return coor;
}
// ----------------------------------------------------------------------------
std::tuple< RealVector, RealVector >
Mech_Fe_Tet4::giveFieldOutputAt( Cell* targetCell, const std::string& fieldTag )
{
    RealVector fieldVal(1), weight(1);
    
    auto cns = this->getNumericsStatusAt(targetCell);
    weight(0) = _wt*cns->_Jdet;
    
    std::vector<Material*> material = this->giveMaterialSetFor(targetCell);
    
    if ( fieldTag == "unassigned" )
        fieldVal(0) = 0.;
    else if ( fieldTag == "s_xx" )
        fieldVal(0) = cns->_stress(0);
    else if ( fieldTag == "s_yy" )
        fieldVal(0) = cns->_stress(1);
    else if ( fieldTag == "s_zz" )
        fieldVal(0) = cns->_stress(2);
    else if ( fieldTag == "s_yz" )
        fieldVal(0) = cns->_stress(3);
    else if ( fieldTag == "s_xz" )
        fieldVal(0) = cns->_stress(4);
    else if ( fieldTag == "s_xy" )
        fieldVal(0) = cns->_stress(5);
    else if ( fieldTag == "ux_x" )
        fieldVal(0) = cns->_gradU(0,0);
    else if ( fieldTag == "uy_x" )
        fieldVal(0) = cns->_gradU(0,1);
    else if ( fieldTag == "uz_x" )
        fieldVal(0) = cns->_gradU(0,2);
    else if ( fieldTag == "ux_y" )
        fieldVal(0) = cns->_gradU(1,0);
    else if ( fieldTag == "uy_y" )
        fieldVal(0) = cns->_gradU(1,1);
    else if ( fieldTag == "uz_y" )
        fieldVal(0) = cns->_gradU(1,2);
    else if ( fieldTag == "ux_z" )
        fieldVal(0) = cns->_gradU(2,0);
    else if ( fieldTag == "uy_z" )
        fieldVal(0) = cns->_gradU(2,1);
    else if ( fieldTag == "uz_z" )
        fieldVal(0) = cns->_gradU(2,2);
    else if ( fieldTag == "g_yz" )
        fieldVal(0) = cns->_gradU(1,2) + cns->_gradU(2,1);
    else if ( fieldTag == "g_xz" )
        fieldVal(0) = cns->_gradU(2,0) + cns->_gradU(0,2);
    else if ( fieldTag == "g_xy" )
        fieldVal(0) = cns->_gradU(1,0) + cns->_gradU(0,1);
    else if ( fieldTag == "ep_xx" )
    {
        try
        {
            fieldVal(0) = material[1]->giveMaterialVariable("plasticStrain_xx", cns->_materialStatus[1]);
        }
        catch (...)
        {
            fieldVal(0) = 0.;
        }
    }
    else if ( fieldTag == "ep_yy" )
    {
        try
        {
            fieldVal(0) = material[1]->giveMaterialVariable("plasticStrain_yy", cns->_materialStatus[1]);
        }
        catch (...)
        {
            fieldVal(0) = 0.;
        }
    }
    else if ( fieldTag == "ep_zz" )
    {
        try
        {
            fieldVal(0) = material[1]->giveMaterialVariable("plasticStrain_zz", cns->_materialStatus[1]);
        }
        catch (...)
        {
            fieldVal(0) = 0.;
        }
    }
    else if ( fieldTag == "gp_yz" )
    {
        try
        {
            fieldVal(0) = material[1]->giveMaterialVariable("plasticStrain_yz", cns->_materialStatus[1]);
        }
        catch (...)
        {
            fieldVal(0) = 0.;
        }
    }
    else if ( fieldTag == "ep_xz" )
    {
        try
        {
            fieldVal(0) = material[1]->giveMaterialVariable("plasticStrain_xz", cns->_materialStatus[1]);
        }
        catch (...)
        {
            fieldVal(0) = 0.;
        }
    }
    else if ( fieldTag == "ep_xy" )
    {
        try
        {
            fieldVal(0) = material[1]->giveMaterialVariable("plasticStrain_xy", cns->_materialStatus[1]);
        }
        catch (...)
        {
            fieldVal(0) = 0.;
        }
    }
    else if ( fieldTag == "ene" )
        fieldVal(0) = 0.5*(cns->_stress(0)*cns->_gradU(0,0) +
                           cns->_stress(1)*cns->_gradU(1,1) +
                           cns->_stress(2)*cns->_gradU(2,2) +
                           cns->_stress(3)*(cns->_gradU(1,2) + cns->_gradU(2,1)) +
                           cns->_stress(3)*(cns->_gradU(0,2) + cns->_gradU(2,0)) +
                           cns->_stress(3)*(cns->_gradU(0,1) + cns->_gradU(1,0)));
    else
        throw std::runtime_error("Invalid tag '" + fieldTag + "' supplied in field output request made to numerics '" + _name + "'!");
    
    return std::make_tuple(std::move(fieldVal), std::move(weight));
}
// ----------------------------------------------------------------------------
std::tuple< std::vector<Dof*>
          , std::vector<Dof*>
          , RealVector >
Mech_Fe_Tet4::giveStaticCoefficientMatrixAt( Cell*           targetCell
                                           , int             stage
                                           , int             subsys
                                           , const TimeData& time )
{
    std::vector<Dof*> rowDof, colDof;
    RealVector coefVal;
    
    if ( stage == _stage[0] && (subsys == _subsystem[0] || subsys == UNASSIGNED) )
    {
        // Retrieve numerics status
        auto cns = this->getNumericsStatusAt(targetCell);
        
        // Retrieve nodal DOFs local to element
        std::vector<Dof*> dof = this->giveNodalDofsAt(targetCell);
        
        // Retrieve material set for element
        std::vector<Material*> material = this->giveMaterialSetFor(targetCell);
        
        // Get tangent modulus
        RealMatrix cmat = material[1]->giveModulusFrom(cns->_strain, cns->_materialStatus[1]);

        // Calculate stiffness matrix
        RealMatrix bmat = this->giveBmatAt(targetCell);
        RealMatrix kmat = _wt*cns->_Jdet*trp(bmat)*cmat*bmat;
        
        rowDof.assign(144, nullptr);
        colDof.assign(144, nullptr);
        coefVal.init(144);
        
        int counter = 0;
        for ( int i = 0; i < 12; i++)
            for ( int j = 0; j < 12; j++)
            {
                rowDof[counter] = dof[i];
                colDof[counter] = dof[j];
                coefVal(counter) = kmat(i,j);
                ++counter;
            }
    }
    
    return std::make_tuple(std::move(rowDof), std::move(colDof), std::move(coefVal));
}
// ----------------------------------------------------------------------------
std::tuple< std::vector<Dof*>
          , RealVector >
Mech_Fe_Tet4::giveStaticLeftHandSideAt(Cell*            targetCell
                                      , int             stage
                                      , int             subsys
                                      , const TimeData& time )
{
    std::vector<Dof*> rowDof;
    RealVector lhs;
    
    if ( stage == _stage[0] && (subsys == _subsystem[0] || subsys == UNASSIGNED) )
    {
        // Retrieve numerics status
        auto cns = this->getNumericsStatusAt(targetCell);
        
        // Local displacements
        rowDof = this->giveNodalDofsAt(targetCell);
        RealVector u = this->giveLocalDisplacementsAt(rowDof, current_value);
            
        // Compute local strains
        RealMatrix bmat = giveBmatAt(targetCell);
        cns->_strain = bmat*u;
        
        // Update material state and compute stress
        std::vector<Material*> material = this->giveMaterialSetFor(targetCell);
        material[1]->updateStatusFrom(cns->_strain, cns->_materialStatus[1]);
        cns->_stress = material[1]->giveForceFrom(cns->_strain, cns->_materialStatus[1]);

        // Calculate lhs
        lhs = _wt*cns->_Jdet*trp(bmat)*cns->_stress;
    }
    
    return std::make_tuple(std::move(rowDof), std::move(lhs));
}
// ----------------------------------------------------------------------------
std::tuple< std::vector<Dof*>
          , RealVector >
Mech_Fe_Tet4::giveStaticRightHandSideAt( Cell*                    targetCell
                                              , int                      stage
                                              , int                      subsys
                                              , const BoundaryCondition& bndCond
                                              , const TimeData&          time )
{   
    std::vector<Dof*> rowDof;
    RealVector rhs;
    
    if ( stage == _stage[0] && (subsys == _subsystem[0] || subsys == UNASSIGNED) )
    {
        if ( bndCond.conditionType() == "Traction" )
        {
            // Retrieve nodes of boundary element
            std::vector<Node*> node = analysisModel().domainManager().giveNodesOf(targetCell);
            
            // ----------------------------------------------------------------
            // Note: Boundary element must be a 2-node line
            // ----------------------------------------------------------------
            
            if ( (int)node.size() != 3 )
                throw std::runtime_error("Error: Traction boundary condition for '" + _name + "' requires 3-node boundary elements!");
            
            RealVector coor0, coor1, coor2, dr0, dr1, dr2;
            coor0 = analysisModel().domainManager().giveCoordinatesOf(node[0]);
            coor1 = analysisModel().domainManager().giveCoordinatesOf(node[1]);
            coor2 = analysisModel().domainManager().giveCoordinatesOf(node[2]);

            // Calculate area using Heron's formula
            dr0 = coor1 - coor0;
            dr1 = coor2 - coor1;
            dr2 = coor0 - coor2;
            double a = std::sqrt(dr0.dot(dr0));
            double b = std::sqrt(dr1.dot(dr1));
            double c = std::sqrt(dr2.dot(dr2));
            double s = 0.5*(a + b + c);
            double area = std::sqrt(s*(s - a)*(s - b)*(s - c));

            // Determine proper value of boundary condition
            RealVector center;
            center = (coor0 + coor1 + coor2)/3.;
            double bcVal = bndCond.valueAt(center, time);

            // Construct local RHS vector (only for relevant DOFs)
            rhs.init(3);
            rhs(0) = area*bcVal/3.;
            rhs(1) = area*bcVal/3.;
            rhs(2) = area*bcVal/3.;

            // Construct global address vector
            rowDof.assign(3, nullptr);

            int dofNum = analysisModel().dofManager().giveIndexForNodalDof(bndCond.targetDof());
            rowDof[0] = analysisModel().domainManager().giveNodalDof(dofNum, node[0]);
            rowDof[1] = analysisModel().domainManager().giveNodalDof(dofNum, node[1]);
            rowDof[2] = analysisModel().domainManager().giveNodalDof(dofNum, node[2]);
        }
        else if ( bndCond.conditionType() == "ConcentratedForce" ) // Boundary conditions for 1-node point
        {
            // Retrieve nodes of boundary element
            std::vector<Node*> node = analysisModel().domainManager().giveNodesOf(targetCell);
            
            // ----------------------------------------------------------------
            // Note: Boundary element must be a 1-node point
            // ----------------------------------------------------------------
            
            if ( (int)node.size() != 1 )
                throw std::runtime_error("Error: Concentrated force boundary condition for '" + _name + "' requires 1-node boundary elements!");
            
            // Determine proper value of boundary condition
            RealVector coor = analysisModel().domainManager().giveCoordinatesOf(node[0]);
            double bcVal = bndCond.valueAt(coor, time);

            // Construct local RHS vector (only for relevant DOFs)
            rhs.init(1);
            rhs(0) = bcVal;

            // Construct global address vector
            rowDof.assign(1, nullptr);

            int dofNum = analysisModel().dofManager().giveIndexForNodalDof(bndCond.targetDof());
            rowDof[0] = analysisModel().domainManager().giveNodalDof(dofNum, node[0]);
        }
    }
    
    return std::make_tuple(std::move(rowDof), std::move(rhs));
}
// ----------------------------------------------------------------------------
std::tuple< std::vector<Dof*>
          , RealVector >
Mech_Fe_Tet4::giveStaticRightHandSideAt( Cell*                 targetCell
                                              , int                   stage
                                              , int                   subsys
                                              , const FieldCondition& fldCond
                                              , const TimeData&       time )
{
    std::vector<Dof*> rowDof;
    RealVector rhs;
    
    if ( stage == _stage[0] && (subsys == _subsystem[0] || subsys == UNASSIGNED) )
    {
        // Evaluate body forces arising from field conditions
        std::string cndType = fldCond.conditionType();
        if ( cndType == "Acceleration_X" || cndType == "Acceleration_Y" || cndType == "Acceleration_Z" )
        {
            auto cns = this->getNumericsStatusAt(targetCell);
            
            // Get element density
            std::vector<Material*> material = this->giveMaterialSetFor(targetCell);
            double rho = material[0]->giveMaterialVariable("Density", cns->_materialStatus[0]);

            // Retrieve nodal DOFs for element
            std::vector<Dof*> dof = this->giveNodalDofsAt(targetCell);
            
            // Note that field condition is assumed to be piecewise linear over each cell
            std::vector<RealVector> coor = this->giveEvaluationPointsFor(targetCell);
            double accel = fldCond.valueAt(coor[0], time);
            double force = accel*rho*_wt*cns->_Jdet/4.0;
            
            rhs = {force, force, force, force};
            
            rowDof.assign(4, nullptr);
            if ( cndType == "Acceleration_X" )
            {
                rowDof[0] = dof[0];
                rowDof[1] = dof[3];
                rowDof[2] = dof[6];
                rowDof[3] = dof[9];
            }
            else if ( cndType == "Acceleration_Y" )
            {
                rowDof[0] = dof[1];
                rowDof[1] = dof[4];
                rowDof[2] = dof[7];
                rowDof[3] = dof[10];
            }
            else // Acceleration_Z
            {
                rowDof[0] = dof[2];
                rowDof[1] = dof[5];
                rowDof[2] = dof[8];
                rowDof[3] = dof[11];
            }
        }        
    }
    
    return std::make_tuple(std::move(rowDof), std::move(rhs));
}
// ----------------------------------------------------------------------------
void Mech_Fe_Tet4::imposeConstraintAt( Cell*                    targetCell
                                     , int                      stage
                                     , const BoundaryCondition& bndCond
                                     , const TimeData&          time )
{                    
    // Only for essential BCs on nodal DOFs
    if ( bndCond.conditionType() == "NodalConstraint" )
    {
        // Retrieve nodes of boundary element
        std::vector<Node*> node = analysisModel().domainManager().giveNodesOf(targetCell);
        int dofNum = analysisModel().dofManager().giveIndexForNodalDof(bndCond.targetDof());
        
        for ( int j = 0; j < (int)node.size(); j++)
        {
            Dof* targetDof = analysisModel().domainManager().giveNodalDof(dofNum, node[j]);
            RealVector coor = analysisModel().domainManager().giveCoordinatesOf(node[j]);
            double bcVal = bndCond.valueAt(coor, time);
            analysisModel().dofManager().setConstraintValueAt(targetDof, bcVal);
        }
    }
}
// ----------------------------------------------------------------------------
void Mech_Fe_Tet4::initializeMaterialsAt( Cell* targetCell )
{
    auto cns = this->getNumericsStatusAt(targetCell);
    std::vector<Material*> material = this->giveMaterialSetFor(targetCell);
    
    cns->_materialStatus[0] = material[0]->createMaterialStatus();
    cns->_materialStatus[1] = material[1]->createMaterialStatus();
}
// ----------------------------------------------------------------------------
void Mech_Fe_Tet4::initializeNumericsAt( Cell* targetCell )
{
    // Check that cell is a 4-node tetrahedron
    int nNodes = analysisModel().domainManager().giveNumberOfNodesOf(targetCell);
    if ( nNodes != 4 )
        throw std::runtime_error("ERROR: Domain cell with " + std::to_string(nNodes) + " nodes detected! 'Mech_Fe_Tet4' requires 4-node tetrahedra.\n");

    targetCell->numericsStatus = new NumericsStatus_Mech_Fe_Tet4();
    auto cns = this->getNumericsStatusAt(targetCell);
    
    // Pre-calculate det(J) and inv(J);
    RealMatrix Jmat = this->giveJacobianMatrixAt(targetCell, _gpNatCoor);
    cns->_Jdet = Jmat(0,0)*Jmat(1,1)*Jmat(2,2) + Jmat(0,1)*Jmat(1,2)*Jmat(2,0) + Jmat(0,2)*Jmat(1,0)*Jmat(2,1)
               - Jmat(2,0)*Jmat(1,1)*Jmat(0,2) - Jmat(2,1)*Jmat(1,2)*Jmat(0,0) - Jmat(2,2)*Jmat(1,0)*Jmat(0,1);
    
    // Sanity check
    if ( cns->_Jdet <= 0 )
        throw std::runtime_error("Calculation of negative volume detected!\nSource: " + _name);
    
    cns->_JmatInv = inv(Jmat);
}
// ----------------------------------------------------------------------------
void Mech_Fe_Tet4::setDofStagesAt( Cell* targetCell )
{
    // Element uses only one stage i.e. stage[0], and all nodal degrees of 
    // freedom get assigned to this stage
    std::vector<Node*> node = analysisModel().domainManager().giveNodesOf(targetCell);
    
    for ( int i = 0; i < (int)node.size(); i++)
    {
        Dof *dof_x, *dof_y, *dof_z;
        
        dof_x = analysisModel().domainManager().giveNodalDof(_nodalDof[0],node[i]);
        dof_y = analysisModel().domainManager().giveNodalDof(_nodalDof[1],node[i]);
        dof_z = analysisModel().domainManager().giveNodalDof(_nodalDof[2],node[i]);
        
        analysisModel().dofManager().setStageFor(dof_x, _stage[0]);
        analysisModel().dofManager().setStageFor(dof_y, _stage[0]);
        analysisModel().dofManager().setStageFor(dof_z, _stage[0]);
    }
}

// Private methods
// ---------------------------------------------------------------------------
NumericsStatus_Mech_Fe_Tet4*
Mech_Fe_Tet4::getNumericsStatusAt( Cell* targetCell )
{
    auto cns = dynamic_cast<NumericsStatus_Mech_Fe_Tet4*>(targetCell->numericsStatus);
    if ( !cns )
        throw std::runtime_error("Error: Unable to retrieve numerics status at cell!\nSource: " + _name);
    
    return cns;
}
// ----------------------------------------------------------------------------
RealMatrix Mech_Fe_Tet4::giveBmatAt( Cell* targetCell )
{
    auto cns = this->getNumericsStatusAt(targetCell);
    RealMatrix dpsi = cns->_JmatInv*_basisFunctionDerivatives;
    
    RealMatrix bmat(6,12);
    for ( int i = 0; i < 4; i++ )
    {
        bmat(0,3*i)   = dpsi(0,i);
        bmat(1,3*i+1) = dpsi(1,i);
        bmat(2,3*i+2) = dpsi(2,i);
        bmat(3,3*i+1) = dpsi(2,i); bmat(3,3*i+2) = dpsi(1,i);
        bmat(4,3*i)   = dpsi(2,i); bmat(4,3*i+2) = dpsi(0,i);
        bmat(5,3*i)   = dpsi(1,i); bmat(5,3*i+1) = dpsi(0,i);
    }
    
    return bmat;
}
// ----------------------------------------------------------------------------
RealMatrix Mech_Fe_Tet4::giveGradBmatAt( Cell* targetCell )
{
    auto cns = this->getNumericsStatusAt(targetCell);
    return cns->_JmatInv*_basisFunctionDerivatives;
}
// ----------------------------------------------------------------------------
RealMatrix Mech_Fe_Tet4::giveJacobianMatrixAt( Cell* targetCell, const RealVector& natCoor )
{
    std::vector<Node*> node = analysisModel().domainManager().giveNodesOf(targetCell);
#ifndef NDEBUG
    if ( (int)node.size() != 4 )
        throw std::runtime_error("\nError: Basis function requires 4 nodes be specified for calculation of Jacobian matrix!\nSource: " + _name);
#endif
    
    RealMatrix coorMat(4,3);
    
    for ( int i = 0; i < 4; i++ )
    {
        RealVector nodeCoor = analysisModel().domainManager().giveCoordinatesOf(node[i]);
        coorMat(i,0) = nodeCoor(0);
        coorMat(i,1) = nodeCoor(1);
        coorMat(i,2) = nodeCoor(2);
    }
    
   return _basisFunctionDerivatives*coorMat;
}
// ----------------------------------------------------------------------------
RealVector Mech_Fe_Tet4::giveLocalDisplacementsAt( std::vector<Dof*>& dof, ValueType valType)
{    
    RealVector u(12);
    for ( int i = 0; i < 12; i++ )
        u(i) = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof[i], valType);
    return u;
}
// ----------------------------------------------------------------------------
std::vector<Dof*> Mech_Fe_Tet4::giveNodalDofsAt( Cell* targetCell )
{
    std::vector<Node*> node = analysisModel().domainManager().giveNodesOf(targetCell);
    std::vector<Dof*> dof(12, nullptr);
    for ( int i = 0; i < 4; i++ )
    {
        dof[3*i]   = analysisModel().domainManager().giveNodalDof(_nodalDof[0], node[i]);
        dof[3*i+1] = analysisModel().domainManager().giveNodalDof(_nodalDof[1], node[i]);
        dof[3*i+2] = analysisModel().domainManager().giveNodalDof(_nodalDof[2], node[i]);
    }
    return dof;
}