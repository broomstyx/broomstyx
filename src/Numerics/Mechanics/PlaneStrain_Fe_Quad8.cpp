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

#include "PlaneStrain_Fe_Quad8.hpp"
#include <cmath>
#include <stdexcept>
#include <string>
#include "Core/AnalysisModel.hpp"
#include "Core/EvalPoint.hpp"
#include "Core/DomainManager.hpp"
#include "Core/ObjectFactory.hpp"
#include "Materials/Material.hpp"
#include "User/UserFunction.hpp"
#include "Util/linearAlgebra.hpp"

#include "IntegrationRules/Legendre_1D.hpp"
#include "IntegrationRules/Legendre_2D_Quad.hpp"
#include "BasisFunctions/Line_P2.hpp"
#include "BasisFunctions/Quadrilateral_S2.hpp"

using namespace broomstyx;

registerBroomstyxObject(Numerics, PlaneStrain_Fe_Quad8)

// Cell numerics status
CellNumericsStatus_PlaneStrain_Fe_Quad8::CellNumericsStatus_PlaneStrain_Fe_Quad8( int nGaussPts )
{
    _nGaussPts = nGaussPts;
    _gp.assign(_nGaussPts, EvalPoint());
    
    for ( int i = 0; i < _nGaussPts; i++ )
        _gp[i].numericsStatus = new EvalPtNumericsStatus_PlaneStrain_Fe_Quad8();
}

CellNumericsStatus_PlaneStrain_Fe_Quad8::~CellNumericsStatus_PlaneStrain_Fe_Quad8()
{
    for ( int i = 0; i < _nGaussPts; i++ )
        delete _gp[i].numericsStatus;
}

// Integration point numerics status
EvalPtNumericsStatus_PlaneStrain_Fe_Quad8::EvalPtNumericsStatus_PlaneStrain_Fe_Quad8()
    : _strain(RealVector(4))
    , _stress(RealVector(4))
    , _gradU(RealMatrix(2,2))
    , _materialStatus {nullptr, nullptr}
{}
    
EvalPtNumericsStatus_PlaneStrain_Fe_Quad8::~EvalPtNumericsStatus_PlaneStrain_Fe_Quad8()
{}

// Constructor
PlaneStrain_Fe_Quad8::PlaneStrain_Fe_Quad8()
{    
    _dim = 2;
    _dofPerCell = 0;
    _dofPerNode = 2;
    
    _nNodes = 8;
    _nMaterials = 2;
    _nStages = 1;
    _nSubsystems = 1;
    
    _name = "PlaneStrain_Fe_Quad8";
    
    _basisFunction = new Quadrilateral_S2();
    _edgeBasisFunction = new Line_P2();
    _integrationRule = new Legendre_2D_Quad(4);
    _edgeIntegrationRule = new Legendre_1D(2);
}

// Destructor
PlaneStrain_Fe_Quad8::~PlaneStrain_Fe_Quad8()
{
    delete _basisFunction;
    delete _edgeBasisFunction;
    delete _integrationRule;
    delete _edgeIntegrationRule;
}

// Public methods
// ----------------------------------------------------------------------------
void PlaneStrain_Fe_Quad8::deleteNumericsAt( Cell* targetCell )
{
    auto cns = this->getNumericsStatusAt(targetCell);
    std::vector<Material*> material = this->giveMaterialSetFor(targetCell);
    
    for ( int i = 0; i < cns->_nGaussPts; i++ )
    {
        auto gpns = this->getNumericsStatusAt(cns->_gp[i]);
        material[0]->destroy(gpns->_materialStatus[0]);
        material[1]->destroy(gpns->_materialStatus[1]);
    }
    
    delete cns;
}
// ----------------------------------------------------------------------------
void PlaneStrain_Fe_Quad8::finalizeDataAt( Cell* targetCell )
{
    // Retrieve nodal DOFs local to element
    std::vector<Dof*> dof = this->giveNodalDofsAt(targetCell);
    
    // Retrieve DOF values
    RealVector uVec = this->giveLocalDisplacementsAt(targetCell, converged_value);
    
    // Displacement values in matrix form
    RealMatrix u(8,2);
    u = {{uVec(0),  uVec(1)},
         {uVec(2),  uVec(3)},
         {uVec(4),  uVec(5)},
         {uVec(6),  uVec(7)},
         {uVec(8),  uVec(9)},
         {uVec(10), uVec(11)},
         {uVec(12), uVec(13)},
         {uVec(14), uVec(15)}};
    
    // Gauss point locations and weights
    std::vector<RealVector> gpLoc;
    RealVector gpWt;
    std::tie(gpLoc,gpWt) = _integrationRule->giveIntegrationPointsAndWeights();
    
    // Retrieve material set for element
    std::vector<Material*> material = this->giveMaterialSetFor(targetCell);
    
    // Cell numerics status
    auto cns = this->getNumericsStatusAt(targetCell);
    
    for ( int i = 0; i < cns->_nGaussPts; i++)
    {
        // Numerics status at Gauss point
        auto gpns = this->getNumericsStatusAt(cns->_gp[i]);
        
        // Shape functions derivatives
        RealMatrix Jmat = this->giveJacobianMatrixAt(targetCell, cns->_gp[i].coordinates);
        
        std::vector<RealVector> dpsiNat = _basisFunction->giveBasisFunctionDerivativesAt(cns->_gp[i].coordinates);
        RealMatrix dpsiNatMat(2,8);
        
        for ( int j = 0; j < 8; j++)
        {
            dpsiNatMat(0,j) = dpsiNat[0](j);
            dpsiNatMat(1,j) = dpsiNat[1](j);
        }
        
        RealMatrix dpsi = inv(Jmat)*dpsiNatMat;
        RealMatrix bmat = this->giveBmatAt(targetCell, cns->_gp[i].coordinates);
        
        // Strain vector
        gpns->_gradU = dpsi*u;
        gpns->_strain = bmat*uVec;
        
        // Stress
        gpns->_stress = material[1]->giveForceFrom(gpns->_strain, gpns->_materialStatus[1]);
    }
}
// ---------------------------------------------------------------------------
RealVector PlaneStrain_Fe_Quad8::giveCellNodeFieldValuesAt( Cell* targetCell, int fieldNum )
{
    RealVector gpVals, wt, nodeVals(8);
    std::string fieldTag;
    
    // Linear extrapolation from Gauss points to nodes (extrapolation matrix is specific to 2x2 integration rule
    // using the arrangement implied in the tensor-product implementation found in class 'Legendre_2d_Quad')

    const double c1 = std::sqrt(3.)/2. + 1.;
    const double c2 = -0.5;
    const double c3 = -std::sqrt(3.)/2. + 1.;
    
    RealMatrix coefMat({{c1, c2, c2, c3},
                        {c2, c1, c3, c2},
                        {c3, c2, c2, c1},
                        {c2, c3, c1, c2}});
    
    try
    {
        fieldTag = _cellFieldOutput.at(fieldNum);
    }
    catch (std::exception& e)
    {
        fieldTag = "unassigned";
    }
    
    // We need a special method for calculating the elastic strain energy
    // which has quadratic behavior within the element
    if ( fieldTag == "ene" )
    {
        RealVector s_xx, s_yy, s_xy, e_xx, e_yy, g_xy;
        std::string subTag;
        
        std::tie(gpVals,wt) = this->giveFieldOutputAt(targetCell, subTag = "s_xx");
        s_xx = coefMat*gpVals;
        std::tie(gpVals,wt) = this->giveFieldOutputAt(targetCell, subTag = "s_yy");
        s_yy = coefMat*gpVals;
        std::tie(gpVals,wt) = this->giveFieldOutputAt(targetCell, subTag = "s_xy");
        s_xy = coefMat*gpVals;
        std::tie(gpVals,wt) = this->giveFieldOutputAt(targetCell, subTag = "ux_x");
        e_xx = coefMat*gpVals;
        std::tie(gpVals,wt) = this->giveFieldOutputAt(targetCell, subTag = "uy_y");
        e_yy = coefMat*gpVals;
        std::tie(gpVals,wt) = this->giveFieldOutputAt(targetCell, subTag = "g_xy");
        g_xy = coefMat*gpVals;
        
        nodeVals(0) = 0.5*(s_xx(0)*e_xx(0) + s_yy(0)*e_yy(0) + s_xy(0)*g_xy(0));
        nodeVals(1) = 0.5*(s_xx(1)*e_xx(1) + s_yy(1)*e_yy(1) + s_xy(1)*g_xy(1));
        nodeVals(2) = 0.5*(s_xx(2)*e_xx(2) + s_yy(2)*e_yy(2) + s_xy(2)*g_xy(2));
        nodeVals(3) = 0.5*(s_xx(3)*e_xx(3) + s_yy(3)*e_yy(3) + s_xy(3)*g_xy(3));
        nodeVals(4) = 0.5*(0.25*(s_xx(0) + s_xx(1))*(e_xx(0) + e_xx(1))
                         + 0.25*(s_yy(0) + s_yy(1))*(e_yy(0) + e_yy(1))
                         + 0.25*(s_xy(0) + s_xy(1))*(g_xy(0) + g_xy(1)));
        nodeVals(5) = 0.5*(0.25*(s_xx(1) + s_xx(2))*(e_xx(1) + e_xx(2))
                         + 0.25*(s_yy(1) + s_yy(2))*(e_yy(1) + e_yy(2))
                         + 0.25*(s_xy(1) + s_xy(2))*(g_xy(1) + g_xy(2)));
        nodeVals(6) = 0.5*(0.25*(s_xx(2) + s_xx(3))*(e_xx(2) + e_xx(3))
                         + 0.25*(s_yy(2) + s_yy(3))*(e_yy(2) + e_yy(3))
                         + 0.25*(s_xy(2) + s_xy(3))*(g_xy(2) + g_xy(3)));
        nodeVals(7) = 0.5*(0.25*(s_xx(3) + s_xx(0))*(e_xx(3) + e_xx(0))
                         + 0.25*(s_yy(3) + s_yy(0))*(e_yy(3) + e_yy(0))
                         + 0.25*(s_xy(3) + s_xy(0))*(g_xy(3) + g_xy(0)));
    }
    else
    {
        std::tie(gpVals,wt) = this->giveFieldOutputAt(targetCell, fieldTag);
        RealVector linVals = coefMat*gpVals;
        
        nodeVals(0) = linVals(0);
        nodeVals(1) = linVals(1);
        nodeVals(2) = linVals(2);
        nodeVals(3) = linVals(3);
        nodeVals(4) = 0.5*(linVals(0) + linVals(1));
        nodeVals(5) = 0.5*(linVals(1) + linVals(2));
        nodeVals(6) = 0.5*(linVals(2) + linVals(3));
        nodeVals(7) = 0.5*(linVals(3) + linVals(0));
    }
    
    return nodeVals;
}
// ---------------------------------------------------------------------------
std::vector<RealVector> PlaneStrain_Fe_Quad8::giveEvaluationPointsFor( Cell *targetCell)
{
    std::vector<Node*> node = analysisModel().domainManager().giveNodesOf(targetCell);
    int nNodes = node.size();
    
    if ( nNodes != 8 )
        throw std::runtime_error("Error: Cell with " + std::to_string(nNodes) + " nodes detected!\n" + _name + " requires 8 nodes per cell.");
    
    int nGaussPts = _integrationRule->giveNumberOfIntegrationPoints();
    std::vector<RealVector> coor(nGaussPts, RealVector());
    std::vector<RealVector> gpLoc;
    RealVector gpWt;
    std::tie(gpLoc, gpWt) = _integrationRule->giveIntegrationPointsAndWeights();
    
    RealVector x(8), y(8), z(8);
    for ( int i = 0; i < 8; i++)
    {
        RealVector nodeCoor = analysisModel().domainManager().giveCoordinatesOf(node[i]);
        x(i) = nodeCoor(0);
        y(i) = nodeCoor(1);
        z(i) = nodeCoor(2);
    }
    
    for ( int i = 0; i < nGaussPts; i++)
    {
        coor[i].init(3);
        RealVector psi = _basisFunction->giveBasisFunctionsAt(gpLoc[i]);
        
        coor[i](0) = x.dot(psi);
        coor[i](1) = y.dot(psi);
        coor[i](2) = z.dot(psi);
    }
    
    return coor;
}
// ---------------------------------------------------------------------------
std::tuple< RealVector
          , RealVector >
PlaneStrain_Fe_Quad8::giveFieldOutputAt( Cell* targetCell, const std::string& fieldTag )
{
    auto cns = this->getNumericsStatusAt(targetCell);
    std::vector<Material*> material = this->giveMaterialSetFor(targetCell);
    
    RealVector fieldVal(cns->_nGaussPts), weight(cns->_nGaussPts);
    for ( int i = 0; i < cns->_nGaussPts; i++)
    {
        auto gpns = this->getNumericsStatusAt(cns->_gp[i]);
        
        if ( fieldTag == "unassigned" )
            fieldVal(i) = 0.;
        else if ( fieldTag == "s_xx" )
            fieldVal(i) = gpns->_stress(0);
        else if ( fieldTag == "s_yy" )
            fieldVal(i) = gpns->_stress(1);
        else if ( fieldTag == "s_zz" )
            fieldVal(i) = gpns->_stress(2);
        else if ( fieldTag == "s_xy" )
            fieldVal(i) = gpns->_stress(3);
        else if ( fieldTag == "ux_x" )
            fieldVal(i) = gpns->_gradU(0,0);
        else if ( fieldTag == "uy_x" )
            fieldVal(i) = gpns->_gradU(0,1);
        else if ( fieldTag == "ux_y" )
            fieldVal(i) = gpns->_gradU(1,0);
        else if ( fieldTag == "uy_y" )
            fieldVal(i) = gpns->_gradU(1,1);
        else if ( fieldTag == "g_xy" )
            fieldVal(i) = gpns->_strain(3);
        else if ( fieldTag == "ep_xx" || fieldTag == "ep_yy" || fieldTag == "ep_zz" || fieldTag == "ep_xy" )
        {
            try
            {
                fieldVal(i) = material[1]->giveMaterialVariable(fieldTag, gpns->_materialStatus[1]);
            }
            catch (...)
            {
                fieldVal(i) = 0.;
            }
        }
        else if ( fieldTag == "ene" )
            fieldVal(i) = 0.5*gpns->_stress.dot(gpns->_strain);
        else
            throw std::runtime_error("Invalid tag '" + fieldTag + "' supplied in field output request made to " + _name + ".");
        
        // Jacobian matrix and determinant
        RealMatrix Jmat = this->giveJacobianMatrixAt(targetCell, cns->_gp[i].coordinates);
        double J = Jmat(0,0)*Jmat(1,1) - Jmat(1,0)*Jmat(0,1);
        
        weight(i) = J*cns->_gp[i].weight;
    }
    
    return std::make_tuple(std::move(fieldVal), std::move(weight));
}
// ---------------------------------------------------------------------------
std::tuple< std::vector<Dof*>
          , std::vector<Dof*>
          , RealVector >
PlaneStrain_Fe_Quad8::giveStaticCoefficientMatrixAt( Cell*           targetCell
                                                  , int             stage
                                                  , int             subsys
                                                  , const TimeData& time )
{
    std::vector<Dof*> rowDof, colDof;
    RealVector coefVal;
    
    if ( stage == _stage[0] && (subsys == _subsystem[0] || subsys == UNASSIGNED) )
    {
        RealMatrix kmat(16,16);
        
        // Retrieve nodal DOFs local to element
        std::vector<Dof*> dof = this->giveNodalDofsAt(targetCell);
        
        // Retrieve material set for element
        std::vector<Material*> material = this->giveMaterialSetFor(targetCell);
        
        // Cell numerics status
        auto cns = this->getNumericsStatusAt(targetCell);
        
        // Loop through Gauss points
        for ( int i = 0; i < cns->_nGaussPts; i++)
        {
            // Jacobian matrix and determinant
            RealMatrix Jmat = this->giveJacobianMatrixAt(targetCell, cns->_gp[i].coordinates);
            double J = Jmat(0,0)*Jmat(1,1) - Jmat(1,0)*Jmat(0,1);
            
            // Bmat
            RealMatrix bmat = this->giveBmatAt(targetCell, cns->_gp[i].coordinates);
            
            // Get tangent modulus
            auto gpns = this->getNumericsStatusAt(cns->_gp[i]);
            RealMatrix cmat = material[1]->giveModulusFrom(gpns->_strain, gpns->_materialStatus[1]);
            
            // Add kmat contribution from Gauss point
            kmat = kmat + trp(bmat)*cmat*(bmat*(J*cns->_gp[i].weight));
        }
        
        rowDof.assign(256, nullptr);
        colDof.assign(256, nullptr);
        coefVal.init(256);
        
        for ( int i = 0; i < 16; i++)
            for ( int j = 0; j < 16; j++)
            {
                rowDof[16*i + j] = dof[i];
                colDof[16*i + j] = dof[j];
                coefVal(16*i + j) = kmat(i,j);
            }
    }
    
    return std::make_tuple(std::move(rowDof), std::move(colDof), std::move(coefVal));
}
// ---------------------------------------------------------------------------
std::tuple< std::vector<Dof*>
          , RealVector >
PlaneStrain_Fe_Quad8::giveStaticLeftHandSideAt( Cell*           targetCell
                                             , int             stage
                                             , int             subsys
                                             , const TimeData& time )
{
    std::vector<Dof*> rowDof;
    RealVector lhs;
    
    if ( stage == _stage[0] && (subsys == _subsystem[0] || subsys == UNASSIGNED) )
    {
        // Retrieve nodal DOFs local to element
        rowDof = this->giveNodalDofsAt(targetCell);
        
        // Retrieve current value of local displacements
        RealVector u = this->giveLocalDisplacementsAt(targetCell, current_value);
            
        // Retrieve material set for element
        std::vector<Material*> material = this->giveMaterialSetFor(targetCell);
        
        // Get numerics status at cell
        auto cns = this->getNumericsStatusAt(targetCell);
        
        // Initialize local LHS
        lhs.init(16);
        
        // Loop through Gauss points
        for ( int i = 0; i < cns->_nGaussPts; i++ )
        {
            // Jacobian matrix and determinant
            RealMatrix Jmat = this->giveJacobianMatrixAt(targetCell, cns->_gp[i].coordinates);
            double J = Jmat(0,0)*Jmat(1,1) - Jmat(1,0)*Jmat(0,1);
            
            // Compute local strains
            auto gpns = this->getNumericsStatusAt(cns->_gp[i]);
            RealMatrix bmat = this->giveBmatAt(targetCell, cns->_gp[i].coordinates);
            gpns->_strain = bmat*u;
            
            // Update material state
            material[1]->updateStatusFrom(gpns->_strain, gpns->_materialStatus[1]);
            
            // Compute stress
            gpns->_stress = material[1]->giveForceFrom(gpns->_strain, gpns->_materialStatus[1]);
            
            // Add lhs contribution from Gauss point
            lhs = lhs + trp(bmat)*(gpns->_stress*(J*cns->_gp[i].weight));
        }
    }
    
    return std::make_tuple(std::move(rowDof), std::move(lhs));
}
// ---------------------------------------------------------------------------
std::tuple< std::vector<Dof*>
          , RealVector >
PlaneStrain_Fe_Quad8::giveStaticRightHandSideAt( Cell* targetCell
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
            // ---------------------------------------------
            // Note: Boundary element must be a 3-node line!
            // ---------------------------------------------

            // Retrieve nodes of boundary element
            std::vector<Node*> node = analysisModel().domainManager().giveNodesOf(targetCell);

            if ( (int)node.size() != 3 )
                throw std::runtime_error("Error: Traction boundary condition for '" + _name + "' requires 3-node boundary elements!");

            // Retrieve x and y values of nodes
            RealVector x(3), y(3), coor;

            for ( int i = 0; i < 3; i++)
            {
                coor = analysisModel().domainManager().giveCoordinatesOf(node[i]);
                x(i) = coor(0);
                y(i) = coor(1);
            }

            // Construct local RHS vector (only for relevant DOFs)
            rhs.init(3);

            // Edge gauss point locations and weights
            int nGaussPts = _edgeIntegrationRule->giveNumberOfIntegrationPoints();
            std::vector<RealVector> gpLoc;
            RealVector gpWt;

            // Cycle through Gauss points
            std::tie(gpLoc,gpWt) = _edgeIntegrationRule->giveIntegrationPointsAndWeights();

            for ( int i = 0; i < nGaussPts; i++)
            {
                // Retrieve shape functions and their derivatives
                RealVector psi = _edgeBasisFunction->giveBasisFunctionsAt(gpLoc[i]);
                std::vector<RealVector> dpsi = _edgeBasisFunction->giveBasisFunctionDerivativesAt(gpLoc[i]);

                // Physical coordinate of Gauss point needed for evaluating BC
                RealVector gpCoor(3);
                gpCoor(0) = psi.dot(x);
                gpCoor(1) = psi.dot(y);
                double bcVal = bndCond.valueAt(gpCoor, time);

                double dxde = dpsi[0].dot(x);
                double dyde = dpsi[0].dot(y);
                double Je = std::sqrt(dxde*dxde + dyde*dyde);

                // Gauss point contribution to RHS
                rhs = rhs + psi*(bcVal*Je*gpWt(i));
            }

            int dofNum = analysisModel().dofManager().giveIndexForNodalDof(bndCond.targetDof());

            rowDof.assign(3, nullptr);
            rowDof[0] = analysisModel().domainManager().giveNodalDof(dofNum, node[0]);
            rowDof[1] = analysisModel().domainManager().giveNodalDof(dofNum, node[1]);
            rowDof[3] = analysisModel().domainManager().giveNodalDof(dofNum, node[2]);
        }
        else if ( bndCond.conditionType() == "ConcentratedForce" )
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

            // Note that bndCond.targetDof() assumes 1-based notation but
            // first element of nodalDof is stored at index 0!
            // So you need to make adjustments
            int dofNum = analysisModel().dofManager().giveIndexForNodalDof(bndCond.targetDof());

            rowDof[0] = analysisModel().domainManager().giveNodalDof(dofNum, node[0]);
        }
    }
    
    return std::make_tuple(std::move(rowDof), std::move(rhs));
}
// ---------------------------------------------------------------------------
std::tuple< std::vector<Dof*>
          , RealVector >
PlaneStrain_Fe_Quad8::giveStaticRightHandSideAt( Cell*                 targetCell
                                              , int                   stage
                                              , int                   subsys
                                              , const FieldCondition& fldCond
                                              , const TimeData&       time )
{
    std::vector<Dof*> rowDof;
    RealVector rhs;
    
    if ( stage == _stage[0] && (subsys == _subsystem[0] || subsys == UNASSIGNED) )
    {
        // Retrieve nodal DOFs for element
        std::vector<Dof*> dof = this->giveNodalDofsAt(targetCell);
        
        // Evaluate body forces arising from field conditions
        std::string cndType = fldCond.conditionType();
        if ( cndType == "Acceleration_X" || cndType == "Acceleration_Y" )
        {
            auto cns = this->getNumericsStatusAt(targetCell);
            
            // Initialize local right hand side
            rhs.init(8);
            
            for ( int i = 0; i < cns->_nGaussPts; i++ )
            {
                double accel = fldCond.valueAt(cns->_gp[i].coordinates, time);
                
                // Jacobian matrix and determinant
                RealMatrix Jmat = this->giveJacobianMatrixAt(targetCell, cns->_gp[i].coordinates);
                double J = Jmat(0,0)*Jmat(1,1) - Jmat(1,0)*Jmat(0,1);
                
                // Get element density
                std::vector<Material*> material = this->giveMaterialSetFor(targetCell);

                // Get numerics status at Gauss point
                auto gpns = this->getNumericsStatusAt(cns->_gp[i]);
                double rho = material[0]->giveMaterialVariable("Density", gpns->_materialStatus[0]);

                // Get shape functions
                RealVector psi = _basisFunction->giveBasisFunctionsAt(cns->_gp[i].coordinates);

                // Add Gauss point contribution to rhs
#pragma GCC ivdep
                for ( int j = 0; j < 8; j++)
                    rhs(j) = rhs(j) + (accel*rho*J*cns->_gp[i].weight)*psi(j);
            }
            
            rowDof.assign(8, nullptr);
            if ( cndType == "Acceleration_X" )
            {
#pragma GCC ivdep
                for ( int i = 0; i < 8; i++ )
                    rowDof[i] = dof[2*i];
            }
            else
            {
#pragma GCC ivdep
                for ( int i = 0; i < 8; i++ )
                    rowDof[i] = dof[2*i+1];
            }
        }        
    }
    
    return std::make_tuple(std::move(rowDof), std::move(rhs));
}
// ---------------------------------------------------------------------------
void PlaneStrain_Fe_Quad8::imposeConstraintAt( Cell* targetCell
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
// ---------------------------------------------------------------------------
void PlaneStrain_Fe_Quad8::initializeMaterialsAt( Cell* targetCell )
{
    auto material = this->giveMaterialSetFor(targetCell);
    auto cns = this->getNumericsStatusAt(targetCell);
    
    for ( int i = 0; i < cns->_nGaussPts; i++ )
    {
        auto gpns = this->getNumericsStatusAt(cns->_gp[i]);
        gpns->_materialStatus[0] = material[0]->createMaterialStatus();
        gpns->_materialStatus[1] = material[1]->createMaterialStatus();
    }
}
// ---------------------------------------------------------------------------
void PlaneStrain_Fe_Quad8::initializeNumericsAt( Cell* targetCell ) 
{
    std::vector<RealVector> gpCoor;
    RealVector gpWt;
    
    int nGaussPts = _integrationRule->giveNumberOfIntegrationPoints();
    std::tie(gpCoor, gpWt) = _integrationRule->giveIntegrationPointsAndWeights();
    
    targetCell->numericsStatus = new CellNumericsStatus_PlaneStrain_Fe_Quad8(nGaussPts);
    auto cns = this->getNumericsStatusAt(targetCell);
    for ( int i = 0; i < nGaussPts; i++ )
    {
        cns->_gp[i].coordinates = gpCoor[i];
        cns->_gp[i].weight = gpWt(i);
    }
}
// ---------------------------------------------------------------------------
void PlaneStrain_Fe_Quad8::setDofStagesAt( Cell* targetCell )
{
    // A. Nodal DOFs
    // Element uses only one stage i.e. stage[0], and
    // both nodal degrees of freedom get assigned to this stage
    std::vector<Node*> node = analysisModel().domainManager().giveNodesOf(targetCell);
    int nNodes = node.size();
    
    if ( nNodes != 8 )
        throw std::runtime_error("Error: Cell with " + std::to_string(nNodes) + " nodes detected!\n" + _name + " requires 8 nodes per cell.");
    
#pragma GCC ivdep
    for ( int i = 0; i < 8; i++)
    {
        Dof *dof_x, *dof_y;
        
        dof_x = analysisModel().domainManager().giveNodalDof(_nodalDof[0], node[i]);
        dof_y = analysisModel().domainManager().giveNodalDof(_nodalDof[1], node[i]);
        
        analysisModel().dofManager().setStageFor(dof_x, _stage[0]);
        analysisModel().dofManager().setStageFor(dof_y, _stage[0]);
    }
        
    // B. Element DOFs
    // This element doesn't use elemental degrees of freedom    
}

// Private methods
// ---------------------------------------------------------------------------
CellNumericsStatus_PlaneStrain_Fe_Quad8*
PlaneStrain_Fe_Quad8::getNumericsStatusAt( Cell* targetCell )
{
    auto cns = dynamic_cast<CellNumericsStatus_PlaneStrain_Fe_Quad8*>(targetCell->numericsStatus);
    if ( !cns )
        throw std::runtime_error("Error: Unable to retrieve numerics status at cell!\nSource: " + _name);
    
    return cns;
}
// ---------------------------------------------------------------------------
EvalPtNumericsStatus_PlaneStrain_Fe_Quad8*
PlaneStrain_Fe_Quad8::getNumericsStatusAt( EvalPoint& gp )
{
    auto gpns = dynamic_cast<EvalPtNumericsStatus_PlaneStrain_Fe_Quad8*>(gp.numericsStatus);
    if ( !gpns )
        throw std::runtime_error("Error: Unable to retrieve numerics status at integration point!\nSource: " + _name);
    
    return gpns;
}
// ---------------------------------------------------------------------------
RealMatrix PlaneStrain_Fe_Quad8::giveJacobianMatrixAt( Cell* targetCell, const RealVector& natCoor )
{
    std::vector<Node*> node = analysisModel().domainManager().giveNodesOf(targetCell);
#ifndef NDEBUG
    if ( (int)node.size() != 8 )
        throw std::runtime_error("\nError: Basis function requires 8 nodes be specified for calculation of Jacobian matrix!\nSource: " + _name);
#endif
    
    std::vector<RealVector> dpsi = _basisFunction->giveBasisFunctionDerivativesAt(natCoor);
    RealVector x(8), y(8);
    
#pragma GCC ivdep
    for ( int i = 0; i < 8; i++ )
    {
        RealVector nodeCoor = analysisModel().domainManager().giveCoordinatesOf(node[i]);
        x(i) = nodeCoor(0);
        y(i) = nodeCoor(1);
    }
    
    RealMatrix jmat(2,2);
    jmat(0,0) = dpsi[0].dot(x);
    jmat(0,1) = dpsi[0].dot(y);
    jmat(1,0) = dpsi[1].dot(x);
    jmat(1,1) = dpsi[1].dot(y);
    
    return jmat;
}
// ----------------------------------------------------------------------------
RealMatrix PlaneStrain_Fe_Quad8::giveBmatAt( Cell* targetCell, const RealVector& natCoor )
{
    RealMatrix dpsi, bmat(4,16);
    
    RealMatrix jmat = this->giveJacobianMatrixAt(targetCell, natCoor);
    std::vector<RealVector> dpsiNat = _basisFunction->giveBasisFunctionDerivativesAt(natCoor);
    
    RealMatrix dpsiNatMat(2,8);

    for ( int i = 0; i < 8; i++)
    {
        dpsiNatMat(0,i) = dpsiNat[0](i);
        dpsiNatMat(1,i) = dpsiNat[1](i);
    }
    
    dpsi = inv(jmat)*dpsiNatMat;
    
    bmat = {
        {dpsi(0,0), 0,         dpsi(0,1), 0,         dpsi(0,2), 0,         dpsi(0,3), 0,
         dpsi(0,4), 0,         dpsi(0,5), 0,         dpsi(0,6), 0,         dpsi(0,7), 0},
        {0,         dpsi(1,0), 0,         dpsi(1,1), 0,         dpsi(1,2), 0,         dpsi(1,3),
         0,         dpsi(1,4), 0,         dpsi(1,5), 0,         dpsi(1,6), 0,         dpsi(1,7)},
        {0,         0,         0,         0,         0,         0,         0,         0,
         0,         0,         0,         0,         0,         0,         0,         0},
        {dpsi(1,0), dpsi(0,0), dpsi(1,1), dpsi(0,1), dpsi(1,2), dpsi(0,2), dpsi(1,3), dpsi(0,3),
         dpsi(1,4), dpsi(0,4), dpsi(1,5), dpsi(0,5), dpsi(1,6), dpsi(0,6), dpsi(1,7), dpsi(0,7)}
    };

    return bmat;
}
// ----------------------------------------------------------------------------
RealVector PlaneStrain_Fe_Quad8::giveLocalDisplacementsAt( Cell* targetCell, ValueType valType )
{    
    std::vector<Dof*> dof = this->giveNodalDofsAt(targetCell);
    RealVector u(16);
    
#pragma GCC ivdep
    for ( int i = 0; i < 16; i++ )
        u(i) = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof[i], valType);

    return u;
}
// ----------------------------------------------------------------------------
std::vector<Dof*> PlaneStrain_Fe_Quad8::giveNodalDofsAt( Cell* targetCell )
{
    std::vector<Node*> node = analysisModel().domainManager().giveNodesOf(targetCell);
    std::vector<Dof*> dof(16, nullptr);

#pragma GCC ivdep
    for ( int i = 0; i < 8; i++ )
    {
        dof[2*i] = analysisModel().domainManager().giveNodalDof(_nodalDof[0], node[i]);
        dof[2*i+1] = analysisModel().domainManager().giveNodalDof(_nodalDof[1], node[i]);
    }
    
    return dof;
}