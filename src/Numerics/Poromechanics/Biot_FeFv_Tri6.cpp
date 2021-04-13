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

#include "Biot_FeFv_Tri6.hpp"
#include <cmath>
#include "Core/AnalysisModel.hpp"
#include "Core/DomainManager.hpp"
#include "Core/EvalPoint.hpp"
#include "Core/ObjectFactory.hpp"
#include "Materials/Material.hpp"
#include "Util/RealMatrix.hpp"
#include "Util/linearAlgebra.hpp"
#include "Util/readOperations.hpp"

#include "IntegrationRules/Legendre_1D.hpp"
#include "IntegrationRules/Legendre_2D_Tri.hpp"
#include "BasisFunctions/Line_P2.hpp"
#include "BasisFunctions/Triangle_P1.hpp"
#include "BasisFunctions/Triangle_P2.hpp"

using namespace broomstyx;

registerBroomstyxObject(Numerics, Biot_FeFv_Tri6)

// Cell numerics status
CellNumericsStatus_Biot_FeFv_Tri6::CellNumericsStatus_Biot_FeFv_Tri6( int nGaussPts )
    : _cellCenter(RealVector(3))
    , _area(0.)
    , _head(0.)
    , _centerFlux {0., 0.}
    , _headOnFace {0., 0., 0.}
    , _fluxOnFace {0., 0., 0.}
    , _fluxIsPrescribedOnFace {false, false, false}
    , _headIsPrescribedOnFace {false, false, false}
    , _hasHeadConstraint(false)
    , _hasNotComputedTransmissibilities(true)
    , _transmissibility{0., 0., 0.}
    , _materialStatus(nullptr)
{
    _nGaussPts = nGaussPts;
    _gp.assign(_nGaussPts, EvalPoint());

    for ( int i = 0; i < _nGaussPts; i++ )
        _gp[i].numericsStatus = new EvalPtNumericsStatus_Biot_FeFv_Tri6();
}

CellNumericsStatus_Biot_FeFv_Tri6::~CellNumericsStatus_Biot_FeFv_Tri6()
{
    for ( int i = 0; i < _nGaussPts; i++ )
        delete _gp[i].numericsStatus;
}

// Evaluation point numerics status
EvalPtNumericsStatus_Biot_FeFv_Tri6::EvalPtNumericsStatus_Biot_FeFv_Tri6()
    : _strain(RealVector(4))
    , _stress(RealVector(4))
    , _gradU(RealMatrix(2,2))
    , _materialStatus {nullptr, nullptr}
{}

EvalPtNumericsStatus_Biot_FeFv_Tri6::~EvalPtNumericsStatus_Biot_FeFv_Tri6() {}

// Constructor
Biot_FeFv_Tri6::Biot_FeFv_Tri6()
{
    _dim = 2;
    _dofPerCell = 1;
    _dofPerNode = 2;
    
    _nNodes = 6;
    _nMaterials = 3;
    _nStages = 1;
    _nSubsystems = 1;
    
    _name = "Biot_FeFv_Tri6";
    
    _basisFunction = new Triangle_P2();
    _edgeBasisFunction = new Line_P2();
    _integrationRule = new Legendre_2D_Tri(3);
    _edgeIntegrationRule = new Legendre_1D(3);
    
    // Default is no gravity term
    _sgn = 0;
    _vrtIndex = 2;

    this->formExtrapolationMatrix();
}

// Destructor
Biot_FeFv_Tri6::~Biot_FeFv_Tri6()
{
    delete _basisFunction;
    delete _edgeBasisFunction;
    delete _integrationRule;
    delete _edgeIntegrationRule;
}

// Public methods
// ----------------------------------------------------------------------------
void Biot_FeFv_Tri6::deleteNumericsAt( Cell* targetCell )
{
    auto cns = this->getNumericsStatusAt(targetCell);
    std::vector<Material*> material = this->giveMaterialSetFor(targetCell);

    for ( int i = 0; i < cns->_nGaussPts; i++ )
    {
        auto gpns = this->getNumericsStatusAt(cns->_gp[i]);
        material[0]->destroy(gpns->_materialStatus[0]);
        material[1]->destroy(gpns->_materialStatus[1]);
    }
    
    material[2]->destroy(cns->_materialStatus);

    delete cns;
}
// ----------------------------------------------------------------------------
void Biot_FeFv_Tri6::finalizeDataAt( Cell* targetCell )
{
    // Pointer to numerics status
    auto cns = this->getNumericsStatusAt(targetCell);
    
    // Value of hydraulic head at cell
    Dof* dof1 = analysisModel().domainManager().giveCellDof(_cellDof[0], targetCell);
    cns->_head = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof1, converged_value);
    
    // Permeability tensor for cell
    std::vector<Material*> material = this->giveMaterialSetFor(targetCell);
    RealMatrix kmat;
    kmat = material[2]->giveModulusFrom(RealVector({_rhoF*_gAccel*cns->_head}), cns->_materialStatus);
    
    // Get cell neighbors
    std::vector<Cell*> neighbor = analysisModel().domainManager().giveNeighborsOf(targetCell);
    
    // A. Calculate flux at faces
    std::vector<Node*> node;
    node = analysisModel().domainManager().giveNodesOf(targetCell);
    std::vector<std::vector<Node*> > face = this->giveFaceNodesOf(targetCell);
    
    // Cycle through faces
    for ( int i = 0; i < 3; i++ )
    {
        double length = this->giveLengthOf(face[i]);
        double normalFlux;
        if ( cns->_headIsPrescribedOnFace[i] )
        {
            RealVector nhat = giveOutwardUnitNormalOf(face[i]);
            double d1 = this->giveDistanceToMidpointOf(face[i], cns->_cellCenter);

            RealVector kvec;
            kvec = kmat*nhat;
            double K1 = std::sqrt(kvec.dot(kvec))*_rhoF*_gAccel/_mu;

            normalFlux = K1*(cns->_head - cns->_headOnFace[i])/d1;
        }
        else if ( cns->_fluxIsPrescribedOnFace[i] )
            normalFlux = cns->_fluxOnFace[i];
        else
        {
            // Hydraulic potential at neighbor cell
            Dof* dof2 = analysisModel().domainManager().giveCellDof(_cellDof[0], neighbor[i]);
            double h2 = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof2, converged_value);

            normalFlux = cns->_transmissibility[i]*(cns->_head - h2)/length;
        }
        cns->_fluxOnFace[i] = normalFlux;
    }
    
    // B. RT0 reconstruction of flux at cell center
    std::vector<Node*> vertexNode(3, nullptr);
    vertexNode[0] = node[2];
    vertexNode[1] = node[0];
    vertexNode[2] = node[1];

    RealVector centerFlux(3);
    for ( int i = 0; i < 3; i++ )
    {
        RealVector nhat2d = this->giveOutwardUnitNormalOf(face[i]);
        RealVector nhat({nhat2d(0), nhat2d(1), 0.});
        
        RealVector x = analysisModel().domainManager().giveCoordinatesOf(node[i]);
        RealVector x0 = analysisModel().domainManager().giveCoordinatesOf(vertexNode[i]);
        
        double alpha = cns->_fluxOnFace[i]/nhat.dot(x - x0);
        
        centerFlux = centerFlux + alpha*(cns->_cellCenter - x0);
    }
    cns->_centerFlux[0] = centerFlux(0);
    cns->_centerFlux[1] = centerFlux(1);
    
    // C. Calculation of stresses and strains

    // Retrieve DOF values
    std::vector<Dof*> dof = this->giveNodalDofsAt(targetCell);
    RealVector uVec = this->giveLocalDisplacementsAt(dof, converged_value);
    
    // Displacement values in matrix form
    RealMatrix u(6,2);
    u = {{uVec(0),  uVec(1)},
         {uVec(2),  uVec(3)},
         {uVec(4),  uVec(5)},
         {uVec(6),  uVec(7)},
         {uVec(8),  uVec(9)},
         {uVec(10), uVec(11)}};

    for ( int i = 0; i < cns->_nGaussPts; i++ )
    {
        // Numerics status at Gauss point
        auto gpns = this->getNumericsStatusAt(cns->_gp[i]);
        
        // Shape functions derivatives
        RealMatrix Jmat = this->giveJacobianMatrixAt(targetCell, cns->_gp[i].coordinates);
        
        std::vector<RealVector> dpsiNat = _basisFunction->giveBasisFunctionDerivativesAt(cns->_gp[i].coordinates);
        RealMatrix dpsiNatMat(2,6);
        
        for ( int j = 0; j < 6; j++)
        {
            dpsiNatMat(0,j) = dpsiNat[0](j);
            dpsiNatMat(1,j) = dpsiNat[1](j);
        }
        
        RealMatrix dpsi = inv(Jmat)*dpsiNatMat;

        RealMatrix bmatU;
        RealVector bmatDiv;
        std::tie(bmatU,bmatDiv) = this->giveBmatAt(targetCell, cns->_gp[i].coordinates);
        
        // Strain vector
        gpns->_gradU = dpsi*u;
        gpns->_strain = bmatU*uVec;
        
        // Stress
        material[1]->updateStatusFrom(gpns->_strain, gpns->_materialStatus[1]);
        gpns->_stress = material[1]->giveForceFrom(gpns->_strain, gpns->_materialStatus[1]);
    }
}
// ----------------------------------------------------------------------------
double Biot_FeFv_Tri6::giveCellFieldValueAt( Cell* targetCell, int fieldNum )
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
// ---------------------------------------------------------------------------
RealVector Biot_FeFv_Tri6::giveCellNodeFieldValuesAt( Cell* targetCell, int fieldNum )
{
    RealVector gpVals, wt, nodeVals(6);
    std::string fieldTag;
    
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
        s_xx = _extrapolationMatrix*gpVals;
        std::tie(gpVals,wt) = this->giveFieldOutputAt(targetCell, subTag = "s_yy");
        s_yy = _extrapolationMatrix*gpVals;
        std::tie(gpVals,wt) = this->giveFieldOutputAt(targetCell, subTag = "s_xy");
        s_xy = _extrapolationMatrix*gpVals;
        std::tie(gpVals,wt) = this->giveFieldOutputAt(targetCell, subTag = "ux_x");
        e_xx = _extrapolationMatrix*gpVals;
        std::tie(gpVals,wt) = this->giveFieldOutputAt(targetCell, subTag = "uy_y");
        e_yy = _extrapolationMatrix*gpVals;
        std::tie(gpVals,wt) = this->giveFieldOutputAt(targetCell, subTag = "g_xy");
        g_xy = _extrapolationMatrix*gpVals;
        
        nodeVals(0) = 0.5*(s_xx(0)*e_xx(0) + s_yy(0)*e_yy(0) + s_xy(0)*g_xy(0));
        nodeVals(1) = 0.5*(s_xx(1)*e_xx(1) + s_yy(1)*e_yy(1) + s_xy(1)*g_xy(1));
        nodeVals(2) = 0.5*(s_xx(2)*e_xx(2) + s_yy(2)*e_yy(2) + s_xy(2)*g_xy(2));
        nodeVals(3) = 0.5*(0.25*(s_xx(0) + s_xx(1))*(e_xx(0) + e_xx(1))
                         + 0.25*(s_yy(0) + s_yy(1))*(e_yy(0) + e_yy(1))
                         + 0.25*(s_xy(0) + s_xy(1))*(g_xy(0) + g_xy(1)));
        nodeVals(4) = 0.5*(0.25*(s_xx(1) + s_xx(2))*(e_xx(1) + e_xx(2))
                         + 0.25*(s_yy(1) + s_yy(2))*(e_yy(1) + e_yy(2))
                         + 0.25*(s_xy(1) + s_xy(2))*(g_xy(1) + g_xy(2)));
        nodeVals(5) = 0.5*(0.25*(s_xx(2) + s_xx(0))*(e_xx(2) + e_xx(0))
                         + 0.25*(s_yy(2) + s_yy(0))*(e_yy(2) + e_yy(0))
                         + 0.25*(s_xy(2) + s_xy(0))*(g_xy(2) + g_xy(0)));
    }
    else
    {
        std::tie(gpVals,wt) = this->giveFieldOutputAt(targetCell, fieldTag);
        RealVector linVals = _extrapolationMatrix*gpVals;
        
        nodeVals(0) = linVals(0);
        nodeVals(1) = linVals(1);
        nodeVals(2) = linVals(2);
        nodeVals(3) = 0.5*(linVals(0) + linVals(1));
        nodeVals(4) = 0.5*(linVals(1) + linVals(2));
        nodeVals(5) = 0.5*(linVals(2) + linVals(0));
    }
    
    return nodeVals;
}
// ----------------------------------------------------------------------------
std::vector<RealVector> Biot_FeFv_Tri6::giveEvaluationPointsFor( Cell *targetCell)
{
    std::vector<Node*> node = analysisModel().domainManager().giveNodesOf(targetCell);
    int nNodes = node.size();
    
    if ( nNodes != 6 )
        throw std::runtime_error("Error: Cell with " + std::to_string(nNodes) + " nodes detected!\n" + _name + " requires 6 nodes per cell.");
    
    int nGaussPts = _integrationRule->giveNumberOfIntegrationPoints();
    std::vector<RealVector> coor(nGaussPts, RealVector());
    std::vector<RealVector> gpLoc;
    RealVector gpWt;
    std::tie(gpLoc, gpWt) = _integrationRule->giveIntegrationPointsAndWeights();
    
    RealVector x(6), y(6), z(6);
    for ( int i = 0; i < 6; i++)
    {
        RealVector nodeCoor = analysisModel().domainManager().giveCoordinatesOf(node[i]);
        x(i) = nodeCoor(0);
        y(i) = nodeCoor(1);
        z(i) = nodeCoor(2);
    }
    
    for ( int i = 0; i < nGaussPts; i++)
    {
        coor[i].init(3);
        RealVector psi= _basisFunction->giveBasisFunctionsAt(gpLoc[i]);
        
        coor[i](0) = x.dot(psi);
        coor[i](1) = y.dot(psi);
        coor[i](2) = z.dot(psi);
    }
    
    return coor;
}
// ----------------------------------------------------------------------------
std::tuple< RealVector, RealVector >
Biot_FeFv_Tri6::giveFieldOutputAt( Cell* targetCell, const std::string& fieldTag  )
{
    RealVector fieldVal(3), weight(3);
    
    auto cns = this->getNumericsStatusAt(targetCell);
    std::vector<Material*> material = this->giveMaterialSetFor(targetCell);
    
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
            fieldVal(i) = gpns->_gradU(1,0);
        else if ( fieldTag == "ux_y" )
            fieldVal(i) = gpns->_gradU(0,1);
        else if ( fieldTag == "uy_y" )
            fieldVal(i) = gpns->_gradU(1,1);
        else if ( fieldTag == "g_xy" )
            fieldVal(i) = gpns->_strain(3);
        else if ( fieldTag == "h" )
        {
            Dof* dof = analysisModel().domainManager().giveCellDof(_cellDof[0], targetCell);
            fieldVal(i) = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof, current_value);
        }
        else if (fieldTag == "p" )
        {
            Dof* dof = analysisModel().domainManager().giveCellDof(_cellDof[0], targetCell);
            double head = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof, current_value);
            
            // std::vector<RealVector> coor = this->giveEvaluationPointsFor(targetCell);
            fieldVal(i) = _rhoF*_gAccel*(head + _sgn*cns->_cellCenter(_vrtIndex));
        }
        else if ( fieldTag == "q_x" )
            fieldVal(i) = cns->_centerFlux[0];
        else if ( fieldTag == "q_y" )
            fieldVal(i) = cns->_centerFlux[1];
        else if ( fieldTag == "sp_xx" )
        {
            Dof* dof = analysisModel().domainManager().giveCellDof(_cellDof[0], targetCell);
            double head = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof, current_value);
            
            std::vector<RealVector> coor = this->giveEvaluationPointsFor(targetCell);
            fieldVal(i) = gpns->_stress(0) - _alpha*_rhoF*_gAccel*(head + _sgn*cns->_cellCenter(_vrtIndex));
        }
        else if ( fieldTag == "sp_yy" )
        {
            Dof* dof = analysisModel().domainManager().giveCellDof(_cellDof[0], targetCell);
            double head = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof, current_value);
            
            std::vector<RealVector> coor = this->giveEvaluationPointsFor(targetCell);
            fieldVal(i) = gpns->_stress(1) - _alpha*_rhoF*_gAccel*(head + _sgn*cns->_cellCenter(_vrtIndex));
        }
        else if ( fieldTag == "sp_zz" )
        {
            Dof* dof = analysisModel().domainManager().giveCellDof(_cellDof[0], targetCell);
            double head = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof, current_value);
            
            std::vector<RealVector> coor = this->giveEvaluationPointsFor(targetCell);
            fieldVal(i) = gpns->_stress(2) - _alpha*_rhoF*_gAccel*(head + _sgn*cns->_cellCenter(_vrtIndex));
        }
        else
            throw std::runtime_error("Invalid tag '" + fieldTag + "' encountered in field output request!\nSource: " + _name);
    }
    
    return std::make_tuple(std::move(fieldVal), std::move(weight));
}
// ----------------------------------------------------------------------------
std::tuple< std::vector<Dof*>, std::vector<Dof*>, RealVector >
Biot_FeFv_Tri6::giveStaticCoefficientMatrixAt( Cell*           targetCell
                                             , int             stage
                                             , int             subsys
                                             , const TimeData& time )
{
    std::vector<Dof*> rowDof, colDof;
    RealVector coefVal;
    
    if ( stage == _stage[0] && (subsys == _subsystem[0] || subsys == UNASSIGNED) )
    {
        auto cns = this->getNumericsStatusAt(targetCell);
        
        // Get cell neighbors
        std::vector<Cell*> neighbor = analysisModel().domainManager().giveNeighborsOf(targetCell);
        
        int nCols = 0;
        for ( int i = 0; i < 3; i++ )
            if ( neighbor[i] )
                ++nCols;
        
        rowDof.assign(157 + nCols, nullptr);
        colDof.assign(157 + nCols, nullptr);
        coefVal.init(157 + nCols);
        
        // Retrieve DOFs local to element
        std::vector<Dof*> nodalDof = this->giveNodalDofsAt(targetCell);
        Dof* dof_h = analysisModel().domainManager().giveCellDof(_cellDof[0], targetCell);
        
        // Retrieve material set for element
        std::vector<Material*> material = this->giveMaterialSetFor(targetCell);
        
        // A. Mechanics part

        RealMatrix kmatUU(12,12);
        RealVector kmatUP(12);

        // Loop through Gauss points
        for ( int i = 0; i < cns->_nGaussPts; i++)
        {
            // Jacobian matrix and determinant
            RealMatrix Jmat = this->giveJacobianMatrixAt(targetCell, cns->_gp[i].coordinates);
            double J = Jmat(0,0)*Jmat(1,1) - Jmat(1,0)*Jmat(0,1);
            
            // Bmat
            RealVector bmatDiv;
            RealMatrix bmat;
            std::tie(bmat,bmatDiv) = this->giveBmatAt(targetCell, cns->_gp[i].coordinates);
            
            // Get tangent modulus
            auto gpns = this->getNumericsStatusAt(cns->_gp[i]);
            RealMatrix cmat = material[1]->giveModulusFrom(gpns->_strain, gpns->_materialStatus[1]);
            
            // Add kmat contribution from Gauss point
            kmatUU += trp(bmat)*cmat*(bmat*(J*cns->_gp[i].weight));
            kmatUP += -_alpha*_rhoF*_gAccel*bmatDiv*((J*cns->_gp[i].weight));
        }

        int counter = 0;
        for ( int i = 0; i < 12; i++)
            for ( int j = 0; j < 12; j++)
            {
                rowDof[counter] = nodalDof[i];
                colDof[counter] = nodalDof[j];
                coefVal(counter) = kmatUU(i,j);
                ++counter;
            }

        for ( int i = 0; i < 12; i++ )
        {
            rowDof[counter] = nodalDof[i];
            colDof[counter] = dof_h;
            coefVal(counter) = kmatUP(i);
            ++counter;
        }

        // ---------------------------------------------------
        // B. Flow part
        
        std::vector<std::vector<Node*> > face = this->giveFaceNodesOf(targetCell);

        rowDof[156] = dof_h;
        colDof[156] = dof_h;
        
        // Cycle through faces
        for ( int i = 0; i < 3; i++ )
        {
            if ( cns->_headIsPrescribedOnFace[i] )
            {
                // Unit outward normal of face
                RealVector nhat = this->giveOutwardUnitNormalOf(face[i]);
                double length = this->giveLengthOf(face[i]);
                
                // Distance from cell center to face midpoint
                double d1 = this->giveDistanceToMidpointOf(face[i], cns->_cellCenter);

                // Permeability tensor for cell
                std::vector<Material*> material = this->giveMaterialSetFor(targetCell);
                RealMatrix kmat;
                double head = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof_h, current_value);
                kmat = material[2]->giveModulusFrom(RealVector({_rhoF*_gAccel*head}), cns->_materialStatus);

                // Hydraulic conductivity
                RealVector kvec;
                kvec = kmat*nhat;
                double K1 = std::sqrt(kvec.dot(kvec))*_rhoF*_gAccel/_mu;

                coefVal(156) += K1/d1*length;
            }
            else if ( !cns->_fluxIsPrescribedOnFace[i] )
            {
                // DOF address for neighbouring cell
                ++counter;
                rowDof[counter] = rowDof[156];
                colDof[counter] = analysisModel().domainManager().giveCellDof(_cellDof[0], neighbor[i]);
                
                // Transmissibility coefficient
                if ( cns->_hasNotComputedTransmissibilities )
                {
                    for ( int i = 0; i < 3; i++ )
                        if ( !cns->_headIsPrescribedOnFace[i] && !cns->_fluxIsPrescribedOnFace[i] )
                            cns->_transmissibility[i] = this->giveTransmissibilityCoefficientAt(face[i], targetCell, neighbor[i]);
                    
                    cns->_hasNotComputedTransmissibilities = false;
                }

                coefVal(156) += cns->_transmissibility[i];
                coefVal(counter) = -cns->_transmissibility[i];
            }
        }
    }
    
    return std::make_tuple(std::move(rowDof), std::move(colDof), std::move(coefVal));
}
// ----------------------------------------------------------------------------
std::tuple< std::vector<Dof*>, RealVector >
Biot_FeFv_Tri6::giveStaticLeftHandSideAt( Cell*           targetCell
                                        , int             stage
                                        , int             subsys
                                        , const TimeData& time )
{
    std::vector<Dof*> rowDof;
    RealVector lhs;
    
    if ( stage == _stage[0] && (subsys == _subsystem[0] || subsys == UNASSIGNED) )
    {
        rowDof.assign(13, nullptr);
        lhs.init(13);
        
        auto cns = this->getNumericsStatusAt(targetCell);
        
        // Get coordinate of cell center
        std::vector<RealVector> epCoor = this->giveEvaluationPointsFor(targetCell);
        RealVector cellCoor = epCoor[0];

        // Retrieve material set for cell
        std::vector<Material*> material = this->giveMaterialSetFor(targetCell);

        // A. Calculate normal flux at faces
        RealVector outflow(3);
        
        std::vector<Node*> node;
        node = analysisModel().domainManager().giveNodesOf(targetCell);
        std::vector<std::vector<Node*> > face = this->giveFaceNodesOf(targetCell);

        // Permeability tensor for cell
        Dof* dof1 = analysisModel().domainManager().giveCellDof(_cellDof[0], targetCell);
        double h1 = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof1, current_value);
        RealMatrix kmat;
        kmat = material[2]->giveModulusFrom(RealVector({_rhoF*_gAccel*h1}), cns->_materialStatus);
        
        // Get cell neighbors
        std::vector<Cell*> neighbor = analysisModel().domainManager().giveNeighborsOf(targetCell);

        // Cycle through faces
        for ( int i = 0; i < 3; i++ )
        {
            // Length of face
            double length = this->giveLengthOf(face[i]);

            // Flux calculation
            if ( cns->_headIsPrescribedOnFace[i] )
            {
                // Unit outward normal of face
                RealVector nhat = giveOutwardUnitNormalOf(face[i]);

                // Distance from cell center to face midpoint
                double d1 = this->giveDistanceToMidpointOf(face[i], cellCoor);

                // Hydraulic conductivity
                RealVector kvec;
                kvec = kmat*nhat;
                double K1 = std::sqrt(kvec.dot(kvec))*_rhoF*_gAccel/_mu;

                outflow(i) = length*K1*h1/d1;
            }
            else if ( cns->_fluxIsPrescribedOnFace[i] )
                outflow(i) = length*cns->_fluxOnFace[i];
            else 
            {
                // Hydraulic potential at neighbor cell
                Dof* dof2 = analysisModel().domainManager().giveCellDof(_cellDof[0], neighbor[i]);
                double h2 = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof2, current_value);

                if ( cns->_hasNotComputedTransmissibilities )
                {
                    for ( int i = 0; i < 3; i++ )
                        if ( !cns->_headIsPrescribedOnFace[i] && !cns->_fluxIsPrescribedOnFace[i] )
                            cns->_transmissibility[i] = this->giveTransmissibilityCoefficientAt(face[i], targetCell, neighbor[i]);
                    
                    cns->_hasNotComputedTransmissibilities = false;
                }

                outflow(i) = cns->_transmissibility[i]*(h1 - h2);
            }
        }

        rowDof[12] = analysisModel().domainManager().giveCellDof(_cellDof[0], targetCell);
        lhs(12) = outflow(0) + outflow(1) + outflow(2);
        
        // B. Calculation of stresses and strains
        // Retrieve nodal DOFs local to element
        
        std::vector<Dof*> dof = this->giveNodalDofsAt(targetCell);
        rowDof[0] = dof[0];
        rowDof[1] = dof[1];
        rowDof[2] = dof[2];
        rowDof[3] = dof[3];
        rowDof[4] = dof[4];
        rowDof[5] = dof[5];
        rowDof[6] = dof[6];
        rowDof[7] = dof[7];
        rowDof[8] = dof[8];
        rowDof[9] = dof[9];
        rowDof[10] = dof[10];
        rowDof[11] = dof[11];

        // Retrieve DOF values
        RealVector uVec = this->giveLocalDisplacementsAt(dof, current_value);
        RealVector fmat(12);
        
        // Loop through Gauss points
        for ( int i = 0; i < cns->_nGaussPts; i++ )
        {
            // Jacobian matrix and determinant
            RealMatrix Jmat = this->giveJacobianMatrixAt(targetCell, cns->_gp[i].coordinates);
            double J = Jmat(0,0)*Jmat(1,1) - Jmat(1,0)*Jmat(0,1);
            
            // Compute local strains
            auto gpns = this->getNumericsStatusAt(cns->_gp[i]);
            RealMatrix bmat;
            RealVector bmatDiv;
            std::tie(bmat, bmatDiv) = this->giveBmatAt(targetCell, cns->_gp[i].coordinates);
            gpns->_strain = bmat*uVec;
            
            // Update material state
            material[1]->updateStatusFrom(gpns->_strain, gpns->_materialStatus[1]);
            
            // Compute stress
            gpns->_stress = material[1]->giveForceFrom(gpns->_strain, gpns->_materialStatus[1]);
            
            // Add lhs contribution from Gauss point
            fmat += (trp(bmat)*gpns->_stress - _alpha*_rhoF*_gAccel*h1*bmatDiv)*(J*cns->_gp[i].weight);
        }

        lhs(0) = fmat(0);
        lhs(1) = fmat(1);
        lhs(2) = fmat(2);
        lhs(3) = fmat(3);
        lhs(4) = fmat(4);
        lhs(5) = fmat(5);
        lhs(6) = fmat(6);
        lhs(7) = fmat(7);
        lhs(8) = fmat(8);
        lhs(9) = fmat(9);
        lhs(10) = fmat(10);
        lhs(11) = fmat(11);
    }
    
    return std::make_tuple(std::move(rowDof), std::move(lhs));
}
// ----------------------------------------------------------------------------
std::tuple< std::vector<Dof*>, RealVector >
Biot_FeFv_Tri6::giveStaticRightHandSideAt( Cell*                    targetCell
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

            for ( int i = 0; i < nGaussPts; i++ )
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
            rowDof[2] = analysisModel().domainManager().giveNodalDof(dofNum, node[2]);
        }
        else if ( bndCond.conditionType() == "Flux" || bndCond.conditionType() == "HydraulicHead" )
        {
            // --> Flow boundary conditions
            std::vector<Cell*> domCell = analysisModel().domainManager().giveDomainCellsAssociatedWith(targetCell);
            for ( Cell* curDomCell : domCell )
            {
                Numerics* domCellNumerics = analysisModel().domainManager().giveNumericsFor(curDomCell);
                if ( this == domCellNumerics )
                {
                    // Retrieve nodes of boundary cell
                    std::vector<Node*> bndCellNode = analysisModel().domainManager().giveNodesOf(targetCell);

                    if ( bndCellNode.size() != 3 )
                        throw std::runtime_error("Numerics '" + _name + "'' only accepts boundary cells with exactly three nodes!");

                    std::vector<std::vector<Node*> > face = this->giveFaceNodesOf(curDomCell);

                    // Initialize rowDof and rhs vector
                    rowDof = {analysisModel().domainManager().giveCellDof(_cellDof[0], curDomCell)};
                    // rhs.init(1);

                    // Cycle through faces to determine which one corresponds 
                    // to current boundary cell
                    for ( int i = 0; i < 3; i++ )
                    {
                        if ( (bndCellNode[0] == face[i][0] && bndCellNode[1] == face[i][1]) || 
                             (bndCellNode[1] == face[i][0] && bndCellNode[0] == face[i][1]) )
                        {
                            double length = this->giveLengthOf(bndCellNode);
                            RealVector coor0 = analysisModel().domainManager().giveCoordinatesOf(bndCellNode[0]);
                            RealVector coor1 = analysisModel().domainManager().giveCoordinatesOf(bndCellNode[1]);
                            RealVector midpt;
                            midpt = 0.5*(coor0 + coor1);

                            double bcVal = bndCond.valueAt(midpt, time);

                            auto cns = this->getNumericsStatusAt(curDomCell);
                            if ( bndCond.conditionType() == "HydraulicHead" ) // Dirichlet BC
                            {
                                // Store boundary condition value on face
                                cns->_headOnFace[i] = bcVal;

                                // Hydraulic conductivity component normal to face
                                RealVector nhat = giveOutwardUnitNormalOf(face[i]);
                                // std::vector<RealVector> epCoor = this->giveEvaluationPointsFor(curDomCell);
                                double d = this->giveDistanceToMidpointOf(face[i], cns->_cellCenter);
                                
                                std::vector<Material*> material = this->giveMaterialSetFor(curDomCell);
                                RealMatrix kmat;
                                kmat = material[2]->giveModulusFrom(RealVector({_rhoF*_gAccel*bcVal}), cns->_materialStatus);
                                
                                RealVector kvec;
                                kvec = kmat*nhat;
                                double K = std::sqrt(kvec.dot(kvec))*_rhoF*_gAccel/_mu;

                                rhs = {length*K*(cns->_headOnFace[i])/d};
                            }
                            else
                            {
                                cns->_fluxOnFace[i] = bcVal;
                                rhs = {-bcVal*length};
                            }
                        }
                    }
                }
            }
        }
    }
    
    return std::make_tuple(std::move(rowDof), std::move(rhs));
}
// ----------------------------------------------------------------------------
std::tuple< std::vector<Dof*>, RealVector >
Biot_FeFv_Tri6::giveStaticRightHandSideAt( Cell*                 targetCell
                                         , int                   stage
                                         , int                   subsys
                                         , const FieldCondition& fldCond
                                         , const TimeData&       time )
{
    std::vector<Dof*> rowDof;
    RealVector rhs;
    
    if ( stage == _stage[0] && (subsys == _subsystem[0] || subsys == UNASSIGNED) )
    {
        auto cns = this->getNumericsStatusAt(targetCell);
        
        // Get element density
        int label = analysisModel().domainManager().giveLabelOf(targetCell);
        std::vector<Material*> material = analysisModel().domainManager().giveMaterialSetForDomain(label);
        
        // Retrieve nodal DOFs for element
        std::vector<Dof*> dof = this->giveNodalDofsAt(targetCell);
        
        // Evaluate body forces arising from field conditions
        std::string condType = fldCond.conditionType();
        if ( condType == "SelfWeight" )
        {
            // -=-=-=-=-=-
            // Initialize local right hand side
            rhs.init(6);
            
            for ( int i = 0; i < cns->_nGaussPts; i++ )
            {
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
                for ( int j = 0; j < 6; j++)
                    rhs(j) += _sgn*(rho - _alpha*_rhoF)*_gAccel*(J*cns->_gp[i].weight)*psi(j);
            }
                    
            rowDof.assign(6, nullptr);
                
            if ( _vrtIndex == 1 )
            {
                rowDof[0] = dof[0];
                rowDof[1] = dof[2];
                rowDof[2] = dof[4];
                rowDof[3] = dof[6];
                rowDof[4] = dof[8];
                rowDof[5] = dof[10];
            }
            else if ( _vrtIndex == 2 )
            {
                rowDof[0] = dof[1];
                rowDof[1] = dof[3];
                rowDof[2] = dof[5];
                rowDof[3] = dof[7];
                rowDof[4] = dof[9];
                rowDof[5] = dof[11];
            }
        }        
    }
    
    return std::make_tuple(std::move(rowDof), std::move(rhs));
}
// ----------------------------------------------------------------------------
std::tuple< std::vector<Dof*>, std::vector<Dof*>, RealVector >
Biot_FeFv_Tri6::giveTransientCoefficientMatrixAt( Cell*           targetCell
                                                , int             stage
                                                , int             subsys
                                                , const TimeData& time )
{
    std::vector<Dof*> rowDof, colDof;
    RealVector coefVal;
    
    if ( stage == _stage[0] && (subsys == _subsystem[0] || subsys == UNASSIGNED) )
    {
        rowDof.assign(13, nullptr);
        colDof.assign(13, nullptr);
        coefVal.init(13);
        
        auto cns = this->getNumericsStatusAt(targetCell);
        
        std::vector<Dof*> nodalDof = this->giveNodalDofsAt(targetCell);
        Dof* cellDof = analysisModel().domainManager().giveCellDof(_cellDof[0], targetCell);
        
        RealVector kmatPU(12);

        // Loop through Gauss points
        for ( int i = 0; i < cns->_nGaussPts; i++ )
        {
            // Jacobian matrix and determinant
            RealMatrix Jmat = this->giveJacobianMatrixAt(targetCell, cns->_gp[i].coordinates);
            double J = Jmat(0,0)*Jmat(1,1) - Jmat(1,0)*Jmat(0,1);
            
            // Bmat
            RealVector bmatDiv;
            RealMatrix bmat;
            std::tie(bmat,bmatDiv) = this->giveBmatAt(targetCell, cns->_gp[i].coordinates);
            
            kmatPU += _alpha*bmatDiv*((J*cns->_gp[i].weight));
        }

        for ( int i = 0; i < 12; i++ )
        {
            rowDof[i] = cellDof;
            colDof[i] = nodalDof[i];
            coefVal(i) = kmatPU(i);
        }
        
        rowDof[12] = cellDof;
        colDof[12] = cellDof;
        coefVal(12) = cns->_area*_rhoF*_gAccel/_M;
    }
    
    return std::make_tuple(std::move(rowDof), std::move(colDof), std::move(coefVal));
}
// ----------------------------------------------------------------------------
std::tuple< std::vector<Dof*>, RealVector >
Biot_FeFv_Tri6::giveTransientLeftHandSideAt( Cell*           targetCell
                                           , int             stage
                                           , int             subsys
                                           , const TimeData& time
                                           , ValueType       valType )
{
    std::vector<Dof*> rowDof;
    RealVector lhs;
    
    if ( stage == _stage[0] && (subsys == _subsystem[0] || subsys == UNASSIGNED) )
    {
        std::vector<Dof*> nodalDof = this->giveNodalDofsAt(targetCell);
        Dof* cellDof = analysisModel().domainManager().giveCellDof(_cellDof[0], targetCell);

        rowDof.assign(1, cellDof);
        lhs.init(1);

        auto cns = this->getNumericsStatusAt(targetCell);
        
        // Retrieve nodal displacements and hydraulic head
        RealVector u = this->giveLocalDisplacementsAt(nodalDof, valType);
        double h = analysisModel().dofManager().giveValueOfPrimaryVariableAt(cellDof, valType);

        // Loop through Gauss points
        for ( int i = 0; i < cns->_nGaussPts; i++ )
        {
            // Jacobian matrix and determinant
            RealMatrix Jmat = this->giveJacobianMatrixAt(targetCell, cns->_gp[i].coordinates);
            double J = Jmat(0,0)*Jmat(1,1) - Jmat(1,0)*Jmat(0,1);
            
            // Bmat
            RealVector bmatDiv;
            RealMatrix bmat;
            std::tie(bmat,bmatDiv) = this->giveBmatAt(targetCell, cns->_gp[i].coordinates);
            
            lhs(0) += _alpha*bmatDiv.dot(u)*((J*cns->_gp[i].weight));
        }
        
        lhs(0) += cns->_area*(_rhoF*_gAccel*h/_M);
    }
    
    return std::make_tuple(std::move(rowDof), std::move(lhs));
}
// ----------------------------------------------------------------------------
void Biot_FeFv_Tri6::imposeConstraintAt( Cell*                    targetCell
                                       , int                      stg
                                       , const BoundaryCondition& bndCond
                                       , const TimeData&          time )
{
    // Essential BCs on nodal DOFs
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
    
    // Constraints pertaining to cell DOFS
    if ( bndCond.conditionType() == "CellConstraint" )
    {
        int dofNum = analysisModel().dofManager().giveIndexForCellDof(bndCond.targetDof());
        std::vector<RealVector> coor = this->giveEvaluationPointsFor(targetCell);
        Dof* targetDof = analysisModel().domainManager().giveCellDof(dofNum, targetCell);
        
        double bcVal = bndCond.valueAt(coor[0],time);
        analysisModel().dofManager().setConstraintValueAt(targetDof, bcVal);
        auto cns = this->getNumericsStatusAt(targetCell);
        cns->_hasHeadConstraint = true;
    }
    
    if ( bndCond.conditionType() == "HydraulicHead" || bndCond.conditionType() == "Flux" )
    {
        std::vector<Cell*> domCell = analysisModel().domainManager().giveDomainCellsAssociatedWith(targetCell);
        for ( Cell* curDomCell : domCell )
        {
            Numerics* domCellNumerics = analysisModel().domainManager().giveNumericsFor(curDomCell);
            if ( this == domCellNumerics )
            {
                auto cns = this->getNumericsStatusAt(curDomCell);
                std::vector<Node*> bndCellNode = analysisModel().domainManager().giveNodesOf(targetCell);

                if ( bndCellNode.size() != 3 )
                    throw std::runtime_error("Numerics type '" + _name + "' only accepts boundary cells with exactly three nodes!");

                std::vector< std::vector<Node*> > face = this->giveFaceNodesOf(curDomCell);

                // Cycle through faces to determine which one corresponds to current boundary cell
                for ( int i = 0; i < 3; i++ )
                {
                    if ( (bndCellNode[0] == face[i][0] && bndCellNode[1] == face[i][1]) ||
                         (bndCellNode[1] == face[i][0] && bndCellNode[0] == face[i][1]) )
                    {
                        // Set dirichlet boundary flag for appropriate face
                        if ( bndCond.conditionType() == "HydraulicHead" )
                            cns->_headIsPrescribedOnFace[i] = true;
                        else
                            cns->_fluxIsPrescribedOnFace[i] = true;
                    }
                }
            }
        }
    }
}
// ----------------------------------------------------------------------------
void Biot_FeFv_Tri6::imposeInitialConditionAt( Cell*                   targetCell
                                             , const InitialCondition& initCond )
{
    std::string condType = initCond.conditionType();
    if ( condType == "CellDof" )
    {
        Dof* cellDof = analysisModel().domainManager().giveCellDof(initCond.targetDofNumber(), targetCell);
        
        // Location of cell DOF
        std::vector<RealVector> coor = this->giveEvaluationPointsFor(targetCell);
        double val = initCond.valueAt(coor[0]);
        
        analysisModel().dofManager().updatePrimaryVariableAt(cellDof, val, converged_value);
    }
    else
        throw std::runtime_error("ERROR: Unrecognized initial condition type '" + condType + "' encountered! Source: " + _name);
}
// ----------------------------------------------------------------------------
void Biot_FeFv_Tri6::initializeMaterialsAt( Cell* targetCell )
{
    auto cns = this->getNumericsStatusAt(targetCell);
    std::vector<Material*> material = this->giveMaterialSetFor(targetCell);

    for ( int i = 0; i < cns->_nGaussPts; i++ )
    {
        auto gpns = this->getNumericsStatusAt(cns->_gp[i]);
        gpns->_materialStatus[0] = material[0]->createMaterialStatus();
        gpns->_materialStatus[1] = material[1]->createMaterialStatus();
    }

    cns->_materialStatus = material[2]->createMaterialStatus();
}
// ----------------------------------------------------------------------------
void Biot_FeFv_Tri6::initializeNumericsAt( Cell* targetCell )
{
    std::vector<RealVector> gpCoor;
    RealVector gpWt;

    int nGaussPts = _integrationRule->giveNumberOfIntegrationPoints();
    std::tie(gpCoor, gpWt) = _integrationRule->giveIntegrationPointsAndWeights();

    targetCell->numericsStatus = new CellNumericsStatus_Biot_FeFv_Tri6(nGaussPts);
    auto cns = this->getNumericsStatusAt(targetCell);
    for ( int i = 0; i < nGaussPts; i++ )
    {
        cns->_gp[i].coordinates = gpCoor[i];
        cns->_gp[i].weight = gpWt(i);
    }

    // Calculate and store cell center and area based on vertex nodes of triangle
    std::vector<Node*> node = analysisModel().domainManager().giveNodesOf(targetCell);
    RealVector A, B, C;
    A = analysisModel().domainManager().giveCoordinatesOf(node[0]);
    B = analysisModel().domainManager().giveCoordinatesOf(node[1]);
    C = analysisModel().domainManager().giveCoordinatesOf(node[2]);
    cns->_cellCenter = (A + B + C)/3.;
    cns->_area = A(0)*(B(1) - C(1)) + B(0)*(C(1) - A(1)) + C(0)*(A(1) - B(1));
}
// ----------------------------------------------------------------------------
void Biot_FeFv_Tri6::readAdditionalDataFrom( FILE* fp )
{
    verifyKeyword(fp, "BiotCoefficient", _name);
    _alpha = getRealInputFrom(fp, "Failed to read Biot coefficient from input file!", _name);
    
    verifyKeyword(fp, "BiotModulus", _name);
    _M = getRealInputFrom(fp, "Failed to read Biot modulus from input file!", _name);
    
    verifyKeyword(fp, "FluidDensity", _name);
    _rhoF = getRealInputFrom(fp, "Failed to read fluid density from input file!", _name);
    
    verifyKeyword(fp, "FluidViscosity", _name);
    _mu = getRealInputFrom(fp, "Failed to read fluid viscosity from input file!", _name);
    
    verifyKeyword(fp, "GravitationalAcceleration", _name);
    _gAccel = getRealInputFrom(fp, "Failed to read gravitational acceleration from input file!", _name);
}
// ----------------------------------------------------------------------------
void Biot_FeFv_Tri6::removeConstraintsOn( Cell* targetCell )
{
    auto cns = this->getNumericsStatusAt(targetCell);
    cns->_headIsPrescribedOnFace[0] = false;
    cns->_headIsPrescribedOnFace[1] = false;
    cns->_headIsPrescribedOnFace[2] = false;
    cns->_fluxIsPrescribedOnFace[0] = false;
    cns->_fluxIsPrescribedOnFace[1] = false;
    cns->_fluxIsPrescribedOnFace[2] = false;
    cns->_hasHeadConstraint = false;
    cns->_hasNotComputedTransmissibilities = true;
}
// ----------------------------------------------------------------------------
void Biot_FeFv_Tri6::setDofStagesAt( Cell* targetCell )
{
    // Nodal DOFs
    std::vector<Node*> node = analysisModel().domainManager().giveNodesOf(targetCell);
    for ( int i = 0; i < (int)node.size(); i++)
    {
        Dof *dof_x, *dof_y;
        
        dof_x = analysisModel().domainManager().giveNodalDof(_nodalDof[0],node[i]);
        dof_y = analysisModel().domainManager().giveNodalDof(_nodalDof[1],node[i]);
        
        analysisModel().dofManager().setStageFor(dof_x, _stage[0]);
        analysisModel().dofManager().setStageFor(dof_y, _stage[0]);
        
        analysisModel().dofManager().setSubsystemFor(dof_x, _subsystem[0]);
        analysisModel().dofManager().setSubsystemFor(dof_y, _subsystem[0]);
    }
        
    // Cell DOF
    Dof* dof_h = analysisModel().domainManager().giveCellDof(_cellDof[0], targetCell);
    analysisModel().dofManager().setStageFor(dof_h, _stage[0]);
    analysisModel().dofManager().setSubsystemFor(dof_h, _subsystem[0]);
}

// Private methods
// ---------------------------------------------------------------------------
void Biot_FeFv_Tri6::formExtrapolationMatrix()
{
    std::vector<RealVector> gpLoc;
    RealVector gpWt;
    Triangle_P1 stressInterpolation;

    int nPoints = _integrationRule->giveNumberOfIntegrationPoints();
    std::tie(gpLoc, gpWt) = _integrationRule->giveIntegrationPointsAndWeights();
    
    RealMatrix interpolationMatrix(nPoints, nPoints);
    for ( int i = 0; i < nPoints; i++ )
    {
        RealVector psi = stressInterpolation.giveBasisFunctionsAt(gpLoc[i]);
        for ( int j = 0; j < nPoints; j++ )
            interpolationMatrix(i,j) = psi(j);
    }

    _extrapolationMatrix = inv(interpolationMatrix);
}
// ----------------------------------------------------------------------------
CellNumericsStatus_Biot_FeFv_Tri6*
Biot_FeFv_Tri6::getNumericsStatusAt( Cell* targetCell )
{
    auto cns = dynamic_cast<CellNumericsStatus_Biot_FeFv_Tri6*>(targetCell->numericsStatus);
    if ( !cns )
        throw std::runtime_error("Error: Unable to retrieve cell numerics status at cell!\nSource: " + _name);
    
    return cns;
}
// ---------------------------------------------------------------------------
EvalPtNumericsStatus_Biot_FeFv_Tri6*
Biot_FeFv_Tri6::getNumericsStatusAt( EvalPoint& gp )
{
    auto gpns = dynamic_cast<EvalPtNumericsStatus_Biot_FeFv_Tri6*>(gp.numericsStatus);
    if ( !gpns )
        throw std::runtime_error("Error: Unable to retrieve numerics status at integration point!\nSource: " + _name);
    
    return gpns;
}
// ----------------------------------------------------------------------------
std::tuple< RealMatrix, RealVector > Biot_FeFv_Tri6::giveBmatAt( Cell* targetCell, const RealVector& natCoor )
{
    RealMatrix dpsi, bmat(4,12);
    
    RealMatrix jmat = this->giveJacobianMatrixAt(targetCell, natCoor);
    std::vector<RealVector> dpsiNat = _basisFunction->giveBasisFunctionDerivativesAt(natCoor);
    
    RealMatrix dpsiNatMat(2,6);

    for ( int i = 0; i < 6; i++)
    {
        dpsiNatMat(0,i) = dpsiNat[0](i);
        dpsiNatMat(1,i) = dpsiNat[1](i);
    }
    
    dpsi = inv(jmat)*dpsiNatMat;
    
    bmat = {
        {dpsi(0,0), 0,         dpsi(0,1), 0,         dpsi(0,2), 0,         dpsi(0,3), 0,         dpsi(0,4), 0,         dpsi(0,5), 0},
        {0,         dpsi(1,0), 0,         dpsi(1,1), 0,         dpsi(1,2), 0,         dpsi(1,3), 0,         dpsi(1,4), 0,         dpsi(1,5)},
        {0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0},
        {dpsi(1,0), dpsi(0,0), dpsi(1,1), dpsi(0,1), dpsi(1,2), dpsi(0,2), dpsi(1,3), dpsi(0,3), dpsi(1,4), dpsi(0,4), dpsi(1,5), dpsi(0,5)}
    };

    // Fill entries for BmatDiv
    RealVector bmatDiv({dpsi(0,0), dpsi(1,0), dpsi(0,1), dpsi(1,1), dpsi(0,2), dpsi(1,2), dpsi(0,3), dpsi(1,3), dpsi(0,4), dpsi(1,4), dpsi(0,5), dpsi(1,5)});

    return std::make_tuple(std::move(bmat), std::move(bmatDiv));
}
// ----------------------------------------------------------------------------
std::vector<std::vector<Node*> > Biot_FeFv_Tri6::giveFaceNodesOf( Cell* targetCell )
{
    std::vector<Node*> node = analysisModel().domainManager().giveNodesOf(targetCell);
    std::vector<Node*> faceTemplate(2, nullptr);
    std::vector<std::vector<Node*> > face(3, faceTemplate);
    
    face[0][0] = node[0];
    face[0][1] = node[1];
    face[1][0] = node[1];
    face[1][1] = node[2];
    face[2][0] = node[2];
    face[2][1] = node[0];
    
    return face;
}
// ----------------------------------------------------------------------------
double Biot_FeFv_Tri6::giveDistanceToMidpointOf( std::vector<Node*>& face
                                               , RealVector&         coor )
{
    RealVector coor0, coor1, dx;
    coor0 = analysisModel().domainManager().giveCoordinatesOf(face[0]);
    coor1 = analysisModel().domainManager().giveCoordinatesOf(face[1]);
    dx = 0.5*(coor0 + coor1) - coor;
    
    return std::sqrt(dx.dot(dx));
}
// ----------------------------------------------------------------------------
RealMatrix Biot_FeFv_Tri6::giveJacobianMatrixAt( Cell* targetCell, const RealVector& natCoor )
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
double Biot_FeFv_Tri6::giveLengthOf( std::vector<Node*>& face )
{
    RealVector coor0, coor1, dx;
    coor0 = analysisModel().domainManager().giveCoordinatesOf(face[0]);
    coor1 = analysisModel().domainManager().giveCoordinatesOf(face[1]);
    dx = coor1 - coor0;
    
    return std::sqrt(dx.dot(dx));
}
// ----------------------------------------------------------------------------
RealVector Biot_FeFv_Tri6::giveLocalDisplacementsAt( std::vector<Dof*>& dof, ValueType valType )
{
    RealVector u(12);
    
#pragma GCC ivdep
    for ( int i = 0; i < 12; i++ )
        u(i) = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof[i], valType);

    return u;
}
// ----------------------------------------------------------------------------
std::vector<Dof*> Biot_FeFv_Tri6::giveNodalDofsAt( Cell* targetCell )
{
    std::vector<Node*> node = analysisModel().domainManager().giveNodesOf(targetCell);
    std::vector<Dof*> dof(12, nullptr);

#pragma GCC ivdep
    for ( int i = 0; i < 6; i++ )
    {
        dof[2*i] = analysisModel().domainManager().giveNodalDof(_nodalDof[0], node[i]);
        dof[2*i+1] = analysisModel().domainManager().giveNodalDof(_nodalDof[1], node[i]);
    }
    
    return dof;
}
// ----------------------------------------------------------------------------
RealVector Biot_FeFv_Tri6::giveOutwardUnitNormalOf( std::vector<Node*>& face )
{
    RealVector coor0, coor1, dx, normVec;
    coor0 = analysisModel().domainManager().giveCoordinatesOf(face[0]);
    coor1 = analysisModel().domainManager().giveCoordinatesOf(face[1]);
    dx = coor1 - coor0;
    double dist = std::sqrt(dx.dot(dx));
    
    normVec = {dx(1)/dist, -dx(0)/dist};
    
    return normVec;
}
// ----------------------------------------------------------------------------
double Biot_FeFv_Tri6::giveTransmissibilityCoefficientAt
    ( std::vector<Node*>& face
    , Cell* targetCell
    , Cell* neighborCell )
{
    // Unit outward normal of face
    RealVector nhat;
    nhat = this->giveOutwardUnitNormalOf(face);
    double length = this->giveLengthOf(face);

    // A. Target cell
    
    // Distance from cell center to midpoint of face
    auto cns = this->getNumericsStatusAt(targetCell);
    double d1 = this->giveDistanceToMidpointOf(face, cns->_cellCenter);
    
    // Permeability tensor for cell
    std::vector<Material*> material = this->giveMaterialSetFor(targetCell);
    
    Dof* dof_h = analysisModel().domainManager().giveCellDof(_cellDof[0], targetCell);
    double head = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof_h, current_value);
    
    RealMatrix kmat;
    RealVector kvec;
    kmat = material[2]->giveModulusFrom(RealVector({_rhoF*_gAccel*head}), cns->_materialStatus);
    kvec = kmat*nhat;
    double K1 = std::sqrt(kvec.dot(kvec))*_rhoF*_gAccel/_mu;
    
    // B. Neighbor cell
    
    // Distance from cell center to midpoint of face
    cns = this->getNumericsStatusAt(neighborCell);
    double d2 = this->giveDistanceToMidpointOf(face, cns->_cellCenter);
    
    // Permeability tensor for neigbor cell
    material = this->giveMaterialSetFor(neighborCell);

    dof_h = analysisModel().domainManager().giveCellDof(_cellDof[0], neighborCell);
    head = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof_h, current_value);

    kmat = material[2]->giveModulusFrom(RealVector({_rhoF*_gAccel*head}), cns->_materialStatus);
    kvec = kmat*nhat;
    double K2 = std::sqrt(kvec.dot(kvec))*_rhoF*_gAccel/_mu;
    
    return length/(d1/K1 + d2/K2);
}