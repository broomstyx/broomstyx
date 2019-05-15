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

#include "StVenantTorsion_Fe_Tri6.hpp"
#include <cmath>
#include <stdexcept>
#include <string>
#include "../../Core/AnalysisModel.hpp"
#include "../../Core/EvalPoint.hpp"
#include "../../Core/DomainManager.hpp"
#include "../../Core/ObjectFactory.hpp"
#include "../../Materials/Material.hpp"
#include "../../User/UserFunction.hpp"
#include "../../Util/linearAlgebra.hpp"

#include "../../IntegrationRules/Legendre_1D.hpp"
#include "../../IntegrationRules/Legendre_2D_Tri.hpp"
#include "../../BasisFunctions/Line_P2.hpp"
#include "../../BasisFunctions/Triangle_P2.hpp"
#include "../../Util/readOperations.hpp"

#include "../../Core/Node.hpp"

using namespace broomstyx;

registerBroomstyxObject(Numerics, StVenantTorsion_Fe_Tri6)

// Cell numerics status
CellNumericsStatus_StVenantTorsion_Fe_Tri6::CellNumericsStatus_StVenantTorsion_Fe_Tri6( int nGaussPts )
{
    _nodalStrain_zx.init(6);
    _nodalStrain_zy.init(6);
    _nodalStrain_mag.init(6);
    _nodalStress_zx.init(6);
    _nodalStress_zy.init(6);
    _nodalStress_mag.init(6);
    _zDisp.init(6);
    _torqueContrib = 0.;
}

CellNumericsStatus_StVenantTorsion_Fe_Tri6::~CellNumericsStatus_StVenantTorsion_Fe_Tri6() {}

// Constructor
StVenantTorsion_Fe_Tri6::StVenantTorsion_Fe_Tri6()
{
    _dim = 2;
    _dofPerCell = 0;
    _dofPerNode = 1;
    
    _nNodes = 6;
    _nMaterials = 1;
    _nStages = 1;
    _nSubsystems = 1;
    
    _name = "StVenantTorsion_Fe_Tri6";
    
    _basisFunction = new Triangle_P2();
    _integrationRule = new Legendre_2D_Tri(6);
    
    _edgeBasisFunction = new Line_P2();
    _edgeIntegrationRule = new Legendre_1D(3);
    
    _torqueIntegrationRule = new Legendre_2D_Tri(6);
    
    _az = nullptr;
    _kx = nullptr;
    _ky = nullptr;
    _dofGrp = 0;
    _beta = 0.;
}

// Destructor
StVenantTorsion_Fe_Tri6::~StVenantTorsion_Fe_Tri6()
{
    // Print info on calculated center of twist
    double xc = analysisModel().dofManager().giveValueOfPrimaryVariableAt(_kx, converged_value);
    double yc = analysisModel().dofManager().giveValueOfPrimaryVariableAt(_ky, converged_value);
    
    std::printf("  Calculated center of twist\n    x = %.6e\n    y = %.6e\n\n", xc, yc);
    std::fflush(stdout);
    
    delete _basisFunction;
    delete _edgeBasisFunction;
    delete _integrationRule;
    delete _edgeIntegrationRule;
    delete _torqueIntegrationRule;
    
    // Delete lagrange multiplier DOF
    analysisModel().dofManager().destroyNumericsDof(_az);
    analysisModel().dofManager().destroyNumericsDof(_kx);
    analysisModel().dofManager().destroyNumericsDof(_ky);
}

// Public methods
// ----------------------------------------------------------------------------
void StVenantTorsion_Fe_Tri6::deleteNumericsAt( Cell* targetCell )
{
    delete targetCell->numericsStatus;
}
// ----------------------------------------------------------------------------
void StVenantTorsion_Fe_Tri6::finalizeDataAt( Cell* targetCell )
{
    double az = analysisModel().dofManager().giveValueOfPrimaryVariableAt(_az, converged_value);
    double xc = analysisModel().dofManager().giveValueOfPrimaryVariableAt(_kx, converged_value);
    double yc = analysisModel().dofManager().giveValueOfPrimaryVariableAt(_ky, converged_value);
    
    // Retrieve local nodal values
    std::vector<Node*> node = analysisModel().domainManager().giveNodesOf(targetCell);
    std::vector<Dof*> dof = this->giveNodalDofsAt(targetCell);
    RealVector uVec(6), xNode(6), yNode(6);
    
#pragma GCC ivdep
    for ( int i = 0; i < 6; i++ )
    {
        uVec(i) = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof[i], converged_value);
        RealVector nodeCoor = analysisModel().domainManager().giveCoordinatesOf(node[i]);
        xNode(i) = nodeCoor(0);
        yNode(i) = nodeCoor(1);
    }
    
    std::vector<RealVector> nodeNatCoor;
    nodeNatCoor.assign(6, RealVector());
    nodeNatCoor[0] = {0, 0};
    nodeNatCoor[1] = {1, 0};
    nodeNatCoor[2] = {0, 1};
    nodeNatCoor[3] = {0.5, 0};
    nodeNatCoor[4] = {0.5, 0.5};
    nodeNatCoor[5] = {0, 0.5};
    
    auto cns = this->getNumericsStatusAt(targetCell);
    
    // Calculate stress at nodes
    for ( int i = 0; i < 6; i++ )
    {
        RealMatrix bmat = this->giveBmatAt(targetCell, nodeNatCoor[i]);
        RealVector gradU = bmat*uVec;
        cns->_nodalStrain_zx(i) = _beta*(gradU(0) - yNode(i));
        cns->_nodalStrain_zy(i) = _beta*(gradU(1) + xNode(i));
    }
    
    auto material = this->giveMaterialSetFor(targetCell);
    double G = material[0]->giveParameter("ShearModulus");
    
    cns->_nodalStress_zx = G*cns->_nodalStrain_zx;
    cns->_nodalStress_zy = G*cns->_nodalStrain_zy;
    
#pragma GCC ivdep
    for ( int i = 0; i < 6; i++ )
    {
        cns->_nodalStrain_mag(i) = std::sqrt(cns->_nodalStrain_zx(i)*cns->_nodalStrain_zx(i) 
                + cns->_nodalStrain_zy(i)*cns->_nodalStrain_zy(i));
        cns->_nodalStress_mag(i) = std::sqrt(cns->_nodalStress_zx(i)*cns->_nodalStress_zx(i) 
                + cns->_nodalStress_zy(i)*cns->_nodalStress_zy(i));
    }
    
    // Calculate element contribution to torque
    std::vector<RealVector> gpLoc;
    RealVector gpWt;
    std::tie(gpLoc,gpWt) = _torqueIntegrationRule->giveIntegrationPointsAndWeights();
    int nTgp = _torqueIntegrationRule->giveNumberOfIntegrationPoints();
    
    cns->_torqueContrib = 0.;
    for ( int i = 0; i < nTgp; i++ )
    {
        // Shape functions and derivatives
        RealMatrix bmat = this->giveBmatAt(targetCell, gpLoc[i]);
        RealVector nmat = this->giveNmatAt(gpLoc[i]);
        
        // Calculate x and y coordinates of Gauss point
        double x = nmat.dot(xNode);
        double y = nmat.dot(yNode);
        
        // stresses
        RealVector gradU = bmat*uVec;
        double tau_zx = _beta*(gradU(0) - y);
        double tau_zy = _beta*(gradU(1) + x);
        cns->_torqueContrib += (x*tau_zy - y*tau_zx)*gpWt(i);
    }
    
    for ( int i = 0; i < 6; i++ )
        cns->_zDisp(i) = _beta*(uVec(i) + az + xc*yNode(i) - yc*xNode(i));
}
// ---------------------------------------------------------------------------
RealVector StVenantTorsion_Fe_Tri6::giveCellNodeFieldValuesAt( Cell* targetCell, int fieldNum )
{
    RealVector nodeVals(6);
    
    auto cns = this->getNumericsStatusAt(targetCell);
    
    std::string fieldTag;
    try
    {
        fieldTag = _cellFieldOutput.at(fieldNum);
    }
    catch (std::exception& e)
    {
        fieldTag = "unassigned";
    }
    
    if ( fieldTag == "s_zx" )
        nodeVals = cns->_nodalStress_zx;
    else if ( fieldTag == "s_zy" )
        nodeVals = cns->_nodalStress_zy;
    else if ( fieldTag == "s_mag" )
        nodeVals = cns->_nodalStress_mag;
    else if ( fieldTag == "g_zx" )
        nodeVals = cns->_nodalStrain_zx;
    else if ( fieldTag == "g_zy" )
        nodeVals = cns->_nodalStrain_zy;
    else if ( fieldTag == "g_mag" )
        nodeVals = cns->_nodalStrain_mag;
    else if ( fieldTag == "u_z" )
        nodeVals = cns->_zDisp;
    
    return nodeVals;
}
// ---------------------------------------------------------------------------
std::tuple< std::vector<Dof*>
          , std::vector<Dof*>
          , RealVector >
StVenantTorsion_Fe_Tri6::giveStaticCoefficientMatrixAt( Cell*           targetCell
                                                      , int             stage
                                                      , int             subsys
                                                      , const TimeData& time )
{
    std::vector<Dof*> rowDof, colDof;
    RealVector coefVal;
    
    if ( stage == _stage[0] && (subsys == _subsystem[0] || subsys == UNASSIGNED) )
    {
        std::vector<Node*> node = analysisModel().domainManager().giveNodesOf(targetCell);
        std::vector<Dof*> dof(6, nullptr);
        RealVector xNode(6), yNode(6);
        
#pragma GCC ivdep
        for ( int i = 0; i < 6; i++ )
        {
            RealVector coor = analysisModel().domainManager().giveCoordinatesOf(node[i]);
            dof[i] = analysisModel().domainManager().giveNodalDof(_nodalDof[0], node[i]);
            xNode(i) = coor(0);
            yNode(i) = coor(1);
        }
        
        std::vector<Material*> material = this->giveMaterialSetFor(targetCell);
        double G = material[0]->giveParameter("ShearModulus");
        
        std::vector<RealVector> gpLoc;
        RealVector gpWt;
        int nGaussPts = _integrationRule->giveNumberOfIntegrationPoints();
        std::tie(gpLoc,gpWt) = _integrationRule->giveIntegrationPointsAndWeights();
        
        RealMatrix kmat(6,6);
        RealVector pmat_x(6), pmat_y(6), qmat(6), qmat_x(6), qmat_y(6);
        double A = 0., S_x = 0., S_y = 0., I_xx = 0., I_xy = 0., I_yy = 0.;
        
        // Loop through Gauss points
        for ( int i = 0; i < nGaussPts; i++)
        {
            // Jacobian matrix and determinant
            RealMatrix Jmat = this->giveJacobianMatrixAt(targetCell, gpLoc[i]);
            double jac = Jmat(0,0)*Jmat(1,1) - Jmat(1,0)*Jmat(0,1);
            
            // Shape functions and derivatives
            RealVector nmat = this->giveNmatAt(gpLoc[i]);
            RealMatrix bmat = this->giveBmatAt(targetCell, gpLoc[i]);
            
            // Cartesian coordinates of Gauss point
            double x = nmat.dot(xNode);
            double y = nmat.dot(yNode);
            
            // Add kmat contribution from Gauss point
            kmat += G*trp(bmat)*bmat*jac*gpWt(i);
            qmat += nmat*jac*gpWt(i);
            qmat_x += y*nmat*jac*gpWt(i);
            qmat_y -= x*nmat*jac*gpWt(i);
            A += jac*gpWt(i);
            S_x += y*jac*gpWt(i);
            S_y += x*jac*gpWt(i);
            I_xx += y*y*jac*gpWt(i);
            I_xy += x*y*jac*gpWt(i);
            I_yy += x*x*jac*gpWt(i);
        }
        
        rowDof.assign(75, nullptr);
        colDof.assign(75, nullptr);
        coefVal.init(75);
        
        int counter = 0;
        for ( int i = 0; i < 6; i++ )
            for ( int j = 0; j < 6; j++ )
            {
                rowDof[counter] = dof[i];
                colDof[counter] = dof[j];
                coefVal(counter++) = kmat(i,j);
            }
        
        for ( int i = 0; i < 6; i++ )
        {
            rowDof[counter] = _az;
            colDof[counter] = dof[i];
            coefVal(counter++) = qmat(i);
            
            rowDof[counter] = _kx;
            colDof[counter] = dof[i];
            coefVal(counter++) = qmat_x(i);
            
            rowDof[counter] = _ky;
            colDof[counter] = dof[i];
            coefVal(counter++) = qmat_y(i);
        }
        
        rowDof[counter] = _az;
        colDof[counter] = _az;
        coefVal(counter++) = A;
        
        rowDof[counter] = _az;
        colDof[counter] = _kx;
        coefVal(counter++) = S_x;
        
        rowDof[counter] = _kx;
        colDof[counter] = _az;
        coefVal(counter++) = S_x;
        
        rowDof[counter] = _az;
        colDof[counter] = _ky;
        coefVal(counter++) = -S_y;
        
        rowDof[counter] = _ky;
        colDof[counter] = _az;
        coefVal(counter++) = -S_y;
        
        rowDof[counter] = _kx;
        colDof[counter] = _kx;
        coefVal(counter++) = I_xx;
        
        rowDof[counter] = _kx;
        colDof[counter] = _ky;
        coefVal(counter++) = -I_xy;
        
        rowDof[counter] = _ky;
        colDof[counter] = _kx;
        coefVal(counter++) = -I_xy;
        
        rowDof[counter] = _ky;
        colDof[counter] = _ky;
        coefVal(counter++) = I_yy;
    }
    
    return std::make_tuple(std::move(rowDof), std::move(colDof), std::move(coefVal));
}
// ----------------------------------------------------------------------------
std::tuple< std::vector<Dof*>
          , RealVector > 
StVenantTorsion_Fe_Tri6::giveStaticRightHandSideAt( Cell*                    targetCell
                                                  , int                      stage
                                                  , int                      subsys
                                                  , const BoundaryCondition& bndCond
                                                  , const TimeData&          time )
{
    std::vector<Dof*> rowDof;
    RealVector rhs;
    
    if ( stage == _stage[0] && (subsys == _subsystem[0] || subsys == UNASSIGNED) 
                            && bndCond.conditionType() == "TorsionBoundaryTerm" )
    {
        // ---------------------------------------------
        // Note: Boundary element must be a 3-node line!
        // ---------------------------------------------

        // Retrieve nodes of boundary element
        std::vector<Node*> node = analysisModel().domainManager().giveNodesOf(targetCell);

        if ( (int)node.size() != 3 )
            throw std::runtime_error("Incompatible boundary cell encountered for numerics '" + _name + ".");
        
        // Retrieve x and y values of nodes
        RealVector xNode(3), yNode(3);
        
        for ( int i = 0; i < 3; i++ )
        {
            RealVector coor = analysisModel().domainManager().giveCoordinatesOf(node[i]);
            xNode(i) = coor(0);
            yNode(i) = coor(1);
        }
        
        // Construct local RHS vector
        rhs.init(3);
        
        // Edge gauss point locations and weights
        int nGaussPts = _edgeIntegrationRule->giveNumberOfIntegrationPoints();
        std::vector<RealVector> gpLoc;
        RealVector gpWt;
        
        // Cycle through Gauss points
        std::tie(gpLoc,gpWt) = _edgeIntegrationRule->giveIntegrationPointsAndWeights();
        
        // Determine outward direction from associated domain element
        double b = 0;
        std::vector<Cell*> assocDomCell = analysisModel().domainManager().giveDomainCellsAssociatedWith(targetCell);
        
        // Retrieve shear moduli
        double G;
        if ( assocDomCell.size() == 1 )
        {
            std::vector<Material*> material = this->giveMaterialSetFor(assocDomCell[0]);
            G = material[0]->giveParameter("ShearModulus");
        }
        else
        {
            std::vector<Material*> material = this->giveMaterialSetFor(assocDomCell[0]);
            double G1 = material[0]->giveParameter("ShearModulus");
            material = this->giveMaterialSetFor(assocDomCell[1]);
            double G2 = material[0]->giveParameter("ShearModulus");
            G = G1 - G2;
        }
        
        std::vector<Node*> domCellNode = analysisModel().domainManager().giveNodesOf(assocDomCell[0]);
        
        if ( (domCellNode[0] == node[0] && domCellNode[1] == node[1]) ||
             (domCellNode[1] == node[0] && domCellNode[2] == node[1]) ||
             (domCellNode[2] == node[0] && domCellNode[0] == node[1]) )
            b = 1.;
        else if ( (domCellNode[0] == node[1] && domCellNode[1] == node[0]) ||
             (domCellNode[1] == node[1] && domCellNode[2] == node[0]) ||
             (domCellNode[2] == node[1] && domCellNode[0] == node[0]) )
            b = -1.;
        else
        {
            // Sanity check
            throw std::runtime_error("Error: Cannot determine outward direction from associated domain cell!\nSource: " + _name);
        }
        
        for ( int i = 0; i < nGaussPts; i++)
        {
            // Retrieve shape functions and their derivatives
            RealVector psi = _edgeBasisFunction->giveBasisFunctionsAt(gpLoc[i]);
            std::vector<RealVector> dpsi = _edgeBasisFunction->giveBasisFunctionDerivativesAt(gpLoc[i]);
            
            double x = psi.dot(xNode);
            double y = psi.dot(yNode);
            
            double dxde = dpsi[0].dot(xNode);
            double dyde = dpsi[0].dot(yNode);
            
            // Gauss point contribution to RHS
            rhs += G*b*psi*(x*dxde + y*dyde)*gpWt(i);
        }
        
        // Note that bndCond.targetDof() assumes 1-based notation but first
        // element of nodalDof is stored at index 0.
        int dofNum = _nodalDof[ bndCond.targetDof() - 1 ];
        
        
        rowDof.assign(3, nullptr);
        rowDof[0] = analysisModel().domainManager().giveNodalDof(dofNum, node[0]);
        rowDof[1] = analysisModel().domainManager().giveNodalDof(dofNum, node[1]);
        rowDof[2] = analysisModel().domainManager().giveNodalDof(dofNum, node[2]);
    }
    
    return std::make_tuple(std::move(rowDof), std::move(rhs));
}
// ---------------------------------------------------------------------------
void StVenantTorsion_Fe_Tri6::initializeNumericsAt( Cell* targetCell )
{
    std::vector<RealVector> gpCoor;
    RealVector gpWt;
    
    int nGaussPts = _integrationRule->giveNumberOfIntegrationPoints();
    std::tie(gpCoor, gpWt) = _integrationRule->giveIntegrationPointsAndWeights();
    
    targetCell->numericsStatus = new CellNumericsStatus_StVenantTorsion_Fe_Tri6(nGaussPts);
}
// ----------------------------------------------------------------------------
void StVenantTorsion_Fe_Tri6::readAdditionalDataFrom( FILE* fp )
{
    verifyKeyword(fp, "TwistPerLength", _name);
    _beta = getRealInputFrom(fp, "Failed to read rate of twist from input file!", _name);
    
    verifyKeyword(fp, "DofGroup", _name);
    _dofGrp = getIntegerInputFrom(fp, "Failed to read DOF group from input file!", _name);
    
    // Create DOF associated with Lagrange multiplier
    _az = analysisModel().dofManager().createNumericsDofWithGroup(_dofGrp);
    _kx = analysisModel().dofManager().createNumericsDofWithGroup(_dofGrp);
    _ky = analysisModel().dofManager().createNumericsDofWithGroup(_dofGrp);
    
    analysisModel().dofManager().setStageFor(_az, _stage[0]);
    analysisModel().dofManager().setStageFor(_kx, _stage[0]);
    analysisModel().dofManager().setStageFor(_ky, _stage[0]);
    
    analysisModel().dofManager().setSubsystemFor(_az, _subsystem[0]);
    analysisModel().dofManager().setSubsystemFor(_kx, _subsystem[0]);
    analysisModel().dofManager().setSubsystemFor(_ky, _subsystem[0]);
}
// ----------------------------------------------------------------------------
void StVenantTorsion_Fe_Tri6::setDofStagesAt( Cell* targetCell )
{
    std::vector<Node*> node = analysisModel().domainManager().giveNodesOf(targetCell);
    int nNodes = node.size();
    
    if ( nNodes != 6 )
        throw std::runtime_error("Error: Cell with " + std::to_string(nNodes) + " detected!\n" + _name + " requires 6 nodes per cell.");
    
    std::vector<Dof*> dof = this->giveNodalDofsAt(targetCell);

#pragma GCC ivdep
    for ( int i = 0; i < 6; i++ )
    {
        analysisModel().dofManager().setStageFor(dof[i], _stage[0]);
        analysisModel().dofManager().setSubsystemFor(dof[i], _subsystem[0]);
    }
}

// Private methods
// ----------------------------------------------------------------------------
CellNumericsStatus_StVenantTorsion_Fe_Tri6*
StVenantTorsion_Fe_Tri6::getNumericsStatusAt( Cell* targetCell )
{
    auto cns = dynamic_cast<CellNumericsStatus_StVenantTorsion_Fe_Tri6*>(targetCell->numericsStatus);
    if ( !cns )
        throw std::runtime_error("Error: Unable to retrieve numerics status at cell!\nSource: " + _name);
    
    return cns;
}
// ----------------------------------------------------------------------------
RealMatrix StVenantTorsion_Fe_Tri6::giveBmatAt( Cell* targetCell, const RealVector& natCoor )
{
    RealMatrix dpsi;
    
    RealMatrix jmat = this->giveJacobianMatrixAt(targetCell, natCoor);
    std::vector<RealVector> dpsiNat = _basisFunction->giveBasisFunctionDerivativesAt(natCoor);
    
    RealMatrix dpsiNatMat(2,6);

    for ( int i = 0; i < 6; i++)
    {
        dpsiNatMat(0,i) = dpsiNat[0](i);
        dpsiNatMat(1,i) = dpsiNat[1](i);
    }
    
    dpsi = inv(jmat)*dpsiNatMat;
    
    return dpsi;
}
// ----------------------------------------------------------------------------
RealVector StVenantTorsion_Fe_Tri6::giveNmatAt( const RealVector& natCoor )
{
    RealVector psi = _basisFunction->giveBasisFunctionsAt(natCoor);
    
    return psi;
}
// ----------------------------------------------------------------------------
RealMatrix StVenantTorsion_Fe_Tri6::giveJacobianMatrixAt( Cell* targetCell, const RealVector& natCoor )
{
    std::vector<Node*> node = analysisModel().domainManager().giveNodesOf(targetCell);
#ifndef NDEBUG
    if ( (int)node.size() != 6 )
        throw std::runtime_error("\nError: Basis function requires 6 nodes be specified for calculation of Jacobian matrix!\nSource: " + _name);
#endif
    
    std::vector<RealVector> dpsi = _basisFunction->giveBasisFunctionDerivativesAt(natCoor);
    RealVector x(6), y(6);
    
#pragma GCC ivdep
    for ( int i = 0; i < 6; i++ )
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
RealVector StVenantTorsion_Fe_Tri6::giveLocalValuesAt( Cell* targetCell, ValueType valType )
{
    std::vector<Dof*> dof = this->giveNodalDofsAt(targetCell);
    RealVector u(6);
    
#pragma GCC ivdep
    for ( int i = 0; i < 6; i++ )
        u(i) = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof[i], valType);
    
    return u;
}
// ----------------------------------------------------------------------------
std::vector<Dof*> StVenantTorsion_Fe_Tri6::giveNodalDofsAt( Cell* targetCell )
{
    std::vector<Node*> node = analysisModel().domainManager().giveNodesOf(targetCell);
    std::vector<Dof*> dof(6, nullptr);
    
#pragma GCC ivdep
    for ( int i = 0; i < 6; i++ )
        dof[i] = analysisModel().domainManager().giveNodalDof(_nodalDof[0], node[i]);
    
    return dof;
}