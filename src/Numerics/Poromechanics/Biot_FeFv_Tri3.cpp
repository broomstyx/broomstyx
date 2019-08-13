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

#include "Biot_FeFv_Tri3.hpp"
#include <cmath>
#include "../../Core/AnalysisModel.hpp"
#include "../../Core/DomainManager.hpp"
#include "../../Core/ObjectFactory.hpp"
#include "../../Materials/Material.hpp"
#include "../../Util/linearAlgebra.hpp"
#include "../../Util/readOperations.hpp"

using namespace broomstyx;

registerBroomstyxObject(Numerics, Biot_FeFv_Tri3)

// Numerics status
NumericsStatus_Biot_FeFv_Tri3::NumericsStatus_Biot_FeFv_Tri3()
    : _strain(RealVector(4))
    , _stress(RealVector(4))
    , _gradU(RealMatrix(2,2))
    , _head(0.)
    , _centerFlux {0., 0.}
    , _headOnFace {0., 0., 0.}
    , _fluxOnFace {0., 0., 0.}
    , _fluxIsPrescribedOnFace {false, false, false}
    , _headIsPrescribedOnFace {false, false, false}
    , _hasHeadConstraint(false)
    , _hasNotComputedTransmissibilities(true)
    , _transmissibility{0., 0., 0.}
    , _materialStatus {nullptr, nullptr, nullptr}
{}

NumericsStatus_Biot_FeFv_Tri3::~NumericsStatus_Biot_FeFv_Tri3() {}

// Constructor
Biot_FeFv_Tri3::Biot_FeFv_Tri3()
{
    _dim = 2;
    _dofPerCell = 1;
    _dofPerNode = 2;
    
    _nNodes = 3;
    _nMaterials = 3;
    _nStages = 1;
    _nSubsystems = 1;
    
    _name = "Biot_FeFv_Tri3";
    
    // Lone Gauss point coordinate and weight
    _gpNatCoor = {1./3., 1./3.};
    _wt = 0.5;
    
    // Pre-calculate shape functions and derivatives at Gauss point
    _basisFunctionValues = _basisFunction.giveBasisFunctionsAt(_gpNatCoor);
    std::vector<RealVector> dpsiNat = _basisFunction.giveBasisFunctionDerivativesAt(_gpNatCoor);
    _basisFunctionDerivatives.init(2,3);
    for ( int i = 0; i < 3; i++)
    {
        _basisFunctionDerivatives(0,i) = dpsiNat[0](i);
        _basisFunctionDerivatives(1,i) = dpsiNat[1](i);
    }
    
    // Default is no gravity term
    _sgn = 0;
    _vrtIndex = 2;
}

// Destructor
Biot_FeFv_Tri3::~Biot_FeFv_Tri3() {}

// Public methods
// ----------------------------------------------------------------------------
void Biot_FeFv_Tri3::deleteNumericsAt( Cell* targetCell )
{
    auto cns = this->getNumericsStatusAt(targetCell);
    auto material = this->giveMaterialSetFor(targetCell);
    
    material[0]->destroy(cns->_materialStatus[0]);
    material[1]->destroy(cns->_materialStatus[1]);
    material[2]->destroy(cns->_materialStatus[2]);
}
// ----------------------------------------------------------------------------
void Biot_FeFv_Tri3::finalizeDataAt( Cell* targetCell )
{
    // Pointer to numerics status
    auto cns = this->getNumericsStatusAt(targetCell);
    
    // Get coordinate of cell center
    std::vector<RealVector> epCoor = this->giveEvaluationPointsFor(targetCell);
    RealVector cellCoor = epCoor[0];
    
    // Value of hydraulic head at cell
    Dof* dof1 = analysisModel().domainManager().giveCellDof(_cellDof[0], targetCell);
    cns->_head = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof1, converged_value);
    
    // Permeability tensor for cell
    double k_xx, k_yy, k_xy;
    std::vector<Material*> material = this->giveMaterialSetFor(targetCell);
    k_xx = material[2]->giveParameter("Permeability_xx");
    k_yy = material[2]->giveParameter("Permeability_yy");
    k_xy = material[2]->giveParameter("Permeability_xy");
    
    RealMatrix kmat_cell({{k_xx, k_xy},
                          {k_xy, k_yy}});

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
            double d1 = this->giveDistanceToMidpointOf(face[i], cellCoor);

            RealVector kvec;
            kvec = kmat_cell*nhat;
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
        
        centerFlux = centerFlux + alpha*(cellCoor - x0);
    }
    cns->_centerFlux[0] = centerFlux(0);
    cns->_centerFlux[1] = centerFlux(1);
    
    // C. Calculation of stresses and strains
    
    // Retrieve nodal DOFs local to element
    std::vector<Dof*> dof = this->giveNodalDofsAt(targetCell);
    
    // Retrieve DOF values
    RealVector uVec = this->giveLocalDisplacementsAt(dof, converged_value);
    
    // Compute shape functions gradients
    RealMatrix bmatGrad = this->giveGradBmatAt(targetCell);
    
    // Compute displacement gradient
    RealMatrix uMat({{uVec(0), uVec(1)},
                     {uVec(2), uVec(3)},
                     {uVec(4), uVec(5)}});
    RealMatrix gradU = bmatGrad*uMat;
    
    // Construct strain vector
    RealMatrix bmatU;
    RealVector bmatDiv;
    
    std::tie(bmatU,bmatDiv) = this->giveBmatAt(targetCell);
    cns->_strain = bmatU*uVec;
    
    // Get constitutive force
    material[1]->updateStatusFrom(cns->_strain, cns->_materialStatus[1]);
    cns->_stress = material[1]->giveForceFrom(cns->_strain, cns->_materialStatus[1]);
    
    RealVector fmatU;
    fmatU = _wt*cns->_Jdet*(trp(bmatU)*cns->_stress);
    // fmatU = _wt*cns->_Jdet*(trp(bmatU)*cns->_stress - _alpha*_rhoF*_gAccel*cns->_head*bmatDiv);
    
    analysisModel().dofManager().addToSecondaryVariableAt(dof[0], fmatU(0));
    analysisModel().dofManager().addToSecondaryVariableAt(dof[1], fmatU(1));
    analysisModel().dofManager().addToSecondaryVariableAt(dof[2], fmatU(2));
    analysisModel().dofManager().addToSecondaryVariableAt(dof[3], fmatU(3));
    analysisModel().dofManager().addToSecondaryVariableAt(dof[4], fmatU(4));
    analysisModel().dofManager().addToSecondaryVariableAt(dof[5], fmatU(5));
}
// ----------------------------------------------------------------------------
double Biot_FeFv_Tri3::giveCellFieldValueAt( Cell* targetCell, int fieldNum )
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
std::vector<RealVector> Biot_FeFv_Tri3::giveEvaluationPointsFor( Cell *targetCell )
{
    std::vector<Node*> node = analysisModel().domainManager().giveNodesOf(targetCell);
    
    std::vector<RealVector> coor(1, RealVector());
    coor[0].init(3);
    
    for ( int i = 0; i < 3; i++)
    {
        RealVector nodeCoor = analysisModel().domainManager().giveCoordinatesOf(node[i]);
        coor[0](0) += nodeCoor(0)/3.0;
        coor[0](1) += nodeCoor(1)/3.0;
        coor[0](2) += nodeCoor(2)/3.0;
    }
    
    return coor;
}
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
std::tuple< RealVector, RealVector >
Biot_FeFv_Tri3::giveFieldOutputAt( Cell* targetCell, const std::string& fieldTag  )
{
    RealVector fieldVal(1), weight(1);
    auto cns = this->getNumericsStatusAt(targetCell);
    weight(0) = _wt*cns->_Jdet;
    
    if ( fieldTag == "unassigned" )
        fieldVal(0) = 0.;
    else if ( fieldTag == "s_xx" )
        fieldVal(0) = cns->_stress(0);
    else if ( fieldTag == "s_yy" )
        fieldVal(0) = cns->_stress(1);
    else if ( fieldTag == "s_zz" )
        fieldVal(0) = cns->_stress(2);
    else if ( fieldTag == "s_xy" )
        fieldVal(0) = cns->_stress(3);
    else if ( fieldTag == "ux_x" )
        fieldVal(0) = cns->_gradU(0,0);
    else if ( fieldTag == "uy_x" )
        fieldVal(0) = cns->_gradU(1,0);
    else if ( fieldTag == "ux_y" )
        fieldVal(0) = cns->_gradU(0,1);
    else if ( fieldTag == "uy_y" )
        fieldVal(0) = cns->_gradU(1,1);
    else if ( fieldTag == "g_xy" )
        fieldVal(0) = cns->_strain(3);
    else if ( fieldTag == "h" )
    {
        Dof* dof = analysisModel().domainManager().giveCellDof(_cellDof[0], targetCell);
        fieldVal(0) = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof, current_value);
    }
    else if (fieldTag == "p" )
    {
        Dof* dof = analysisModel().domainManager().giveCellDof(_cellDof[0], targetCell);
        double head = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof, current_value);
        
        std::vector<RealVector> coor = this->giveEvaluationPointsFor(targetCell);
        fieldVal(0) = _rhoF*_gAccel*(head + _sgn*coor[0](_vrtIndex));
    }
    else if ( fieldTag == "q_x" )
        fieldVal(0) = cns->_centerFlux[0];
    else if ( fieldTag == "q_y" )
        fieldVal(0) = cns->_centerFlux[1];
    else if ( fieldTag == "sp_xx" )
    {
        Dof* dof = analysisModel().domainManager().giveCellDof(_cellDof[0], targetCell);
        double head = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof, current_value);
        
        std::vector<RealVector> coor = this->giveEvaluationPointsFor(targetCell);
        fieldVal(0) = cns->_stress(0) - _alpha*_rhoF*_gAccel*(head + _sgn*coor[0](_vrtIndex));
    }
    else if ( fieldTag == "sp_yy" )
    {
        Dof* dof = analysisModel().domainManager().giveCellDof(_cellDof[0], targetCell);
        double head = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof, current_value);
        
        std::vector<RealVector> coor = this->giveEvaluationPointsFor(targetCell);
        fieldVal(0) = cns->_stress(1) - _alpha*_rhoF*_gAccel*(head + _sgn*coor[0](_vrtIndex));
    }
    else if ( fieldTag == "sp_zz" )
    {
        Dof* dof = analysisModel().domainManager().giveCellDof(_cellDof[0], targetCell);
        double head = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof, current_value);
        
        std::vector<RealVector> coor = this->giveEvaluationPointsFor(targetCell);
        fieldVal(0) = cns->_stress(2) - _alpha*_rhoF*_gAccel*(head + _sgn*coor[0](_vrtIndex));
    }
    else
        throw std::runtime_error("Invalid tag '" + fieldTag + "' encountered in field output request!\nSource: " + _name);
    
    return std::make_tuple(std::move(fieldVal), std::move(weight));
}
// ----------------------------------------------------------------------------
std::tuple< std::vector<Dof*>, std::vector<Dof*>, RealVector >
Biot_FeFv_Tri3::giveStaticCoefficientMatrixAt( Cell*           targetCell
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
        
        rowDof.assign(43 + nCols, nullptr);
        colDof.assign(43 + nCols, nullptr);
        coefVal.init(43 + nCols);
        
        // Retrieve DOFs local to element
        std::vector<Dof*> nodalDof = this->giveNodalDofsAt(targetCell);
        Dof* cellDof = analysisModel().domainManager().giveCellDof(_cellDof[0], targetCell);
        
        // Retrieve material set for element
        std::vector<Material*> material = this->giveMaterialSetFor(targetCell);
        
        // A. Mechanics part
        // Get tangent modulus
        RealMatrix cmat = material[1]->giveModulusFrom(cns->_strain, cns->_materialStatus[1]);
                
        // Calculate stiffness matrix
        RealMatrix bmatU;
        RealVector bmatDiv;
        std::tie(bmatU,bmatDiv) = this->giveBmatAt(targetCell);
        
        RealMatrix kmatUU;
        RealVector kmatUP;
        kmatUU = _wt*cns->_Jdet*trp(bmatU)*cmat*bmatU;
        kmatUP = -_wt*cns->_Jdet*_alpha*_rhoF*_gAccel*bmatDiv;
        
        int counter = 0;
        for ( int i = 0; i < 6; i++ )
            for ( int j = 0; j < 6; j++ )
            {
                rowDof[counter] = nodalDof[i];
                colDof[counter] = nodalDof[j];
                coefVal(counter) = kmatUU(i,j);
                ++counter;
            }
        
        for ( int i = 0; i < 6; i++ )
        {
            rowDof[counter] = nodalDof[i];
            colDof[counter] = cellDof;
            coefVal(counter) = kmatUP(i);
            ++counter;
        }
        
        // ---------------------------------------------------
        // B. Flow part
        
        std::vector<std::vector<Node*> > face = this->giveFaceNodesOf(targetCell);

        rowDof[42] = cellDof;
        colDof[42] = cellDof;
        
        // Cycle through faces
        for ( int i = 0; i < 3; i++ )
        {
            if ( cns->_headIsPrescribedOnFace[i] )
            {
                // Unit outward normal of face
                RealVector nhat = this->giveOutwardUnitNormalOf(face[i]);
                double length = this->giveLengthOf(face[i]);
                
                // Get coordinate of cell center
                std::vector<RealVector> epCoor;
                epCoor = giveEvaluationPointsFor(targetCell);
                RealVector cellCoor;
                cellCoor = epCoor[0];

                // Distance from cell center to face midpoint
                double d1 = this->giveDistanceToMidpointOf(face[i], cellCoor);

                // Permeability tensor for cell
                // Permeability tensor for cell
                std::vector<Material*> material = this->giveMaterialSetFor(targetCell);
                double k_xx, k_yy, k_xy;
                k_xx = material[2]->giveParameter("Permeability_xx");
                k_yy = material[2]->giveParameter("Permeability_yy");
                k_xy = material[2]->giveParameter("Permeability_xy");

                RealMatrix kmat({{k_xx, k_xy},
                                 {k_xy, k_yy}});
                        
                // Hydraulic conductivity
                RealVector kvec;
                kvec = kmat*nhat;
                double K1 = std::sqrt(kvec.dot(kvec))*_rhoF*_gAccel/_mu;

                coefVal(42) += K1/d1*length;
            }
            else if ( !cns->_fluxIsPrescribedOnFace[i] )
            {
                // DOF address for neighbouring cell
                ++counter;
                rowDof[counter] = rowDof[42];
                colDof[counter] = analysisModel().domainManager().giveCellDof(_cellDof[0], neighbor[i]);
                
                // Transmissibility coefficient
                if ( cns->_hasNotComputedTransmissibilities )
                {
                    for ( int i = 0; i < 3; i++ )
                        if ( !cns->_headIsPrescribedOnFace[i] && !cns->_fluxIsPrescribedOnFace[i] )
                            cns->_transmissibility[i] = this->giveTransmissibilityCoefficientAt(face[i], targetCell, neighbor[i]);
                    
                    cns->_hasNotComputedTransmissibilities = false;
                }

                coefVal(42) += cns->_transmissibility[i];
                coefVal(counter) = -cns->_transmissibility[i];
            }
        }
    }
    
    return std::make_tuple(std::move(rowDof), std::move(colDof), std::move(coefVal));
}
// ----------------------------------------------------------------------------
std::tuple< std::vector<Dof*>, RealVector >
Biot_FeFv_Tri3::giveStaticLeftHandSideAt( Cell*           targetCell
                                        , int             stage
                                        , int             subsys
                                        , const TimeData& time )
{
    std::vector<Dof*> rowDof;
    RealVector lhs;
    
    if ( stage == _stage[0] && (subsys == _subsystem[0] || subsys == UNASSIGNED) )
    {
        rowDof.assign(7, nullptr);
        lhs.init(7);
        
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
        // Cell permeability tensor
        double k_xx, k_yy, k_xy;
        k_xx = material[2]->giveParameter("Permeability_xx");
        k_yy = material[2]->giveParameter("Permeability_yy");
        k_xy = material[2]->giveParameter("Permeability_xy");

        RealMatrix kmat_cell({{k_xx, k_xy},
                              {k_xy, k_yy}});

        // Hydraulic potential at cell
        Dof* dof1 = analysisModel().domainManager().giveCellDof(_cellDof[0], targetCell);
        double h1 = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof1, current_value);

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
                kvec = kmat_cell*nhat;
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

        rowDof[6] = analysisModel().domainManager().giveCellDof(_cellDof[0], targetCell);
        lhs(6) = outflow(0) + outflow(1) + outflow(2);
        
        // B. Calculation of stresses and strains
        // Retrieve nodal DOFs local to element
        
        std::vector<Dof*> dof = this->giveNodalDofsAt(targetCell);
        rowDof[0] = dof[0];
        rowDof[1] = dof[1];
        rowDof[2] = dof[2];
        rowDof[3] = dof[3];
        rowDof[4] = dof[4];
        rowDof[5] = dof[5];

        // Retrieve DOF values
        RealVector uVec = this->giveLocalDisplacementsAt(dof, current_value);

        // strain
        RealMatrix bmatU;
        RealVector bmatDiv;
        
        std::tie(bmatU,bmatDiv) = this->giveBmatAt(targetCell);
        
        cns->_strain = bmatU*uVec;
        cns->_stress = material[1]->giveForceFrom(cns->_strain, cns->_materialStatus[1]);
        RealVector fmat;
        fmat = _wt*cns->_Jdet*(trp(bmatU)*cns->_stress - _alpha*_rhoF*_gAccel*h1*bmatDiv);

        lhs(0) = fmat(0);
        lhs(1) = fmat(1);
        lhs(2) = fmat(2);
        lhs(3) = fmat(3);
        lhs(4) = fmat(4);
        lhs(5) = fmat(5);
    }
    
    return std::make_tuple(std::move(rowDof), std::move(lhs));
}
// ----------------------------------------------------------------------------
std::tuple< std::vector<Dof*>, RealVector >
Biot_FeFv_Tri3::giveStaticRightHandSideAt( Cell*                    targetCell
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

            if ( (int)node.size() != 2 && (int)node.size() != 1 )
                throw std::runtime_error("Incompatible boundary cell encountered for numerics 'CPEF3VBiot'!");

            if ( (int)node.size() == 2 ) // Boundary conditions for 2-node line
            {
                RealVector coor0, coor1, dx;
                coor0 = analysisModel().domainManager().giveCoordinatesOf(node[0]);
                coor1 = analysisModel().domainManager().giveCoordinatesOf(node[1]);
                dx = coor1 - coor0;
                double length = std::sqrt(dx.dot(dx));

                // BC evaluated at midpoint of boundary cell
                RealVector midpt;
                midpt = 0.5*(coor0 + coor1);

                double bcVal = bndCond.valueAt(midpt, time);

                // Construct local RHS vector (only for relevant DOFs)
                rhs.init(2);
                rhs(0) = 0.5*length*bcVal;
                rhs(1) = 0.5*length*bcVal;

                // Construct global address vector
                rowDof.assign(2, nullptr);

                // Note that bndCond.targetDof() assumes 1-based notation but
                // first element of nodalDof is stored at index 0!
                // So you need to make adjustments
                int dofNum = _nodalDof[bndCond.targetDof() - 1];

                rowDof[0] = analysisModel().domainManager().giveNodalDof(dofNum, node[0]);
                rowDof[1] = analysisModel().domainManager().giveNodalDof(dofNum, node[1]);
            }
            // else // Boundary conditions for 1-node point
            // {
            //     // Determine proper value of boundary condition
            //     double bcVal;
            //     if ( bndCond.specification() == "Constant" )
            //         bcVal = bndCond.startingValue();
            //     else if ( bndCond.specification() == "Linear" )
            //     {
            //         // BC is constant over element edge, but linear in time
            //         double startVal = bndCond.startingValue();
            //         double endVal = bndCond.endingValue();
            //         bcVal = startVal + (endVal - startVal)*(time.target - time.start)/(time.end - time.start);
            //     }
            //     else
            //     {
            //         RealVector coor = analysisModel().domainManager().giveCoordinatesOf(node[0]);
            //         bcVal = giveUserFunctionValueAt(coor, time, bndCond.userFunctionLabel());
            //     }

            //     // Construct local RHS vector (only for relevant DOFs)
            //     rhs.init(1);
            //     rhs(0) = bcVal;

            //     // Construct global address vector
            //     rowDof.assign(1, nullptr);

            //     // Note that bndCond.targetDof() assumes 1-based notation but
            //     // first element of nodalDof is stored at index 0!
            //     // So you need to make adjustments
            //     int dofNum = _nodalDof[bndCond.targetDof() - 1];

            //     rowDof[0] = analysisModel().domainManager().giveNodalDof(dofNum, node[0]);
            // }
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

                    if ( bndCellNode.size() != 2 )
                        throw std::runtime_error("Numerics '" + _name + "'' only accepts boundary cells with exactly two nodes!");

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
                                std::vector<RealVector> epCoor = this->giveEvaluationPointsFor(curDomCell);
                                double d = this->giveDistanceToMidpointOf(face[i], epCoor[0]);
                                
                                std::vector<Material*> material = this->giveMaterialSetFor(curDomCell);
                                double k_xx, k_yy, k_xy;
                                k_xx = material[2]->giveParameter("Permeability_xx");
                                k_yy = material[2]->giveParameter("Permeability_yy");
                                k_xy = material[2]->giveParameter("Permeability_xy");

                                RealMatrix kmat_cell({{k_xx, k_xy},
                                                      {k_xy, k_yy}});
                                                      
                                RealVector kvec;
                                kvec = kmat_cell*nhat;
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
Biot_FeFv_Tri3::giveStaticRightHandSideAt( Cell*                 targetCell
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
        double rho = material[0]->giveParameter("Density");
        
        // Retrieve nodal DOFs for element
        std::vector<Dof*> dof = this->giveNodalDofsAt(targetCell);
        
        // Evaluate body forces arising from field conditions
        std::string condType = fldCond.conditionType();
        if ( condType == "SelfWeight" )
        {
            double force = _sgn*(rho - _alpha*_rhoF)*_gAccel*_wt*cns->_Jdet/3.;
            
            rhs.init(3);
            rhs(0) = force;
            rhs(1) = force;
            rhs(2) = force;
                    
            rowDof.assign(3, nullptr);
                
            if ( _vrtIndex == 1 )
            {
                rowDof[0] = dof[0];
                rowDof[1] = dof[2];
                rowDof[2] = dof[4];
            }
            else if ( _vrtIndex == 2 )
            {
                rowDof[0] = dof[1];
                rowDof[1] = dof[3];
                rowDof[2] = dof[5];
            }
        }        
    }
    
    return std::make_tuple(std::move(rowDof), std::move(rhs));
}
// ----------------------------------------------------------------------------
std::tuple< std::vector<Dof*>, std::vector<Dof*>, RealVector >
Biot_FeFv_Tri3::giveTransientCoefficientMatrixAt( Cell*           targetCell
                                                , int             stage
                                                , int             subsys
                                                , const TimeData& time )
{
    std::vector<Dof*> rowDof, colDof;
    RealVector coefVal;
    
    if ( stage == _stage[0] && (subsys == _subsystem[0] || subsys == UNASSIGNED) )
    {
        rowDof.assign(7, nullptr);
        colDof.assign(7, nullptr);
        coefVal.init(7);
        
        auto cns = this->getNumericsStatusAt(targetCell);
        
        std::vector<Dof*> nodalDof = this->giveNodalDofsAt(targetCell);
        Dof* cellDof = analysisModel().domainManager().giveCellDof(_cellDof[0], targetCell);
        
        RealMatrix bmatU;
        RealVector bmatDiv;
        std::tie(bmatU,bmatDiv) = this->giveBmatAt(targetCell);
        
        RealVector kmatPU;
        kmatPU = _wt*cns->_Jdet*_alpha*bmatDiv;
        
        for (int i = 0; i < 6; i++)
        {
            rowDof[i] = cellDof;
            colDof[i] = nodalDof[i];
            coefVal(i) = kmatPU(i);
        }
        
        rowDof[6] = cellDof;
        colDof[6] = cellDof;
        coefVal(6) = _wt*cns->_Jdet*_rhoF*_gAccel/_M;
    }
    
    return std::make_tuple(std::move(rowDof), std::move(colDof), std::move(coefVal));
}
// ----------------------------------------------------------------------------
std::tuple< std::vector<Dof*>, RealVector >
Biot_FeFv_Tri3::giveTransientLeftHandSideAt( Cell*           targetCell
                                           , int             stage
                                           , int             subsys
                                           , const TimeData& time
                                           , ValueType       valType )
{
    std::vector<Dof*> rowDof;
    RealVector lhs;
    
    if ( stage == _stage[0] && (subsys == _subsystem[0] || subsys == UNASSIGNED) )
    {
        auto cns = this->getNumericsStatusAt(targetCell);
        
        std::vector<Dof*> nodalDof = this->giveNodalDofsAt(targetCell);
        Dof* cellDof = analysisModel().domainManager().giveCellDof(_cellDof[0], targetCell);
        
        RealVector bmatDiv;
        std::tie(std::ignore,bmatDiv) = this->giveBmatAt(targetCell);
        
        RealVector u = this->giveLocalDisplacementsAt(nodalDof, valType);
        double h = analysisModel().dofManager().giveValueOfPrimaryVariableAt(cellDof, valType);
        
        rowDof.assign(1, cellDof);
        lhs.init(1);
        
        lhs(0) = _wt*cns->_Jdet*(_alpha*bmatDiv.dot(u) + _rhoF*_gAccel*h/_M);
    }
    
    return std::make_tuple(std::move(rowDof), std::move(lhs));
}
// ----------------------------------------------------------------------------
void Biot_FeFv_Tri3::imposeConstraintAt( Cell*                    targetCell
                                       , int                      stg
                                       , const BoundaryCondition& bndCond
                                       , const TimeData&          time )
{
    // Essential BCs on nodal DOFs
    if ( bndCond.conditionType() == "NodalConstraint" )
    {
        // Retrieve nodes of boundary element
        std::vector<Node*> node = analysisModel().domainManager().giveNodesOf(targetCell);
        int targetDofNum = _nodalDof[bndCond.targetDof() - 1];

        for ( int j = 0; j < (int)node.size(); j++)
        {
            Dof* targetDof = analysisModel().domainManager().giveNodalDof(targetDofNum, node[j]);
            RealVector coor = analysisModel().domainManager().giveCoordinatesOf(node[j]);
            double bcVal = bndCond.valueAt(coor, time);
            analysisModel().dofManager().setConstraintValueAt(targetDof, bcVal);
        }
    }
    
    // Constraints pertaining to cell DOFS
    if ( bndCond.conditionType() == "CellConstraint" )
    {
        int targetDofNum = _cellDof[bndCond.targetDof() - 1];
        std::vector<RealVector> coor = this->giveEvaluationPointsFor(targetCell);
        Dof* targetDof = analysisModel().domainManager().giveCellDof(targetDofNum, targetCell);
        
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

                if ( bndCellNode.size() != 2 )
                    throw std::runtime_error("Numerics type 'CPEF3VBiot' only accepts boundary cells with exactly two nodes!");

                std::vector< std::vector<Node*> > face = this->giveFaceNodesOf(curDomCell);

                // Cycle through faces to determine which one corresponds to current boundary cell
                for ( int i = 0; i < 3; i++)
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
void Biot_FeFv_Tri3::imposeInitialConditionAt( Cell*                   targetCell
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
        throw std::runtime_error("ERROR: Unrecognized initial condition type '" + condType + "' encountered! Source: CPEF3VBiot");
}
// ----------------------------------------------------------------------------
void Biot_FeFv_Tri3::initializeMaterialsAt( Cell* targetCell )
{
    auto cns = this->getNumericsStatusAt(targetCell);
    std::vector<Material*> material = this->giveMaterialSetFor(targetCell);
    
    cns->_materialStatus[0] = material[0]->createMaterialStatus();
    cns->_materialStatus[1] = material[1]->createMaterialStatus();
    cns->_materialStatus[2] = material[2]->createMaterialStatus();
}
// ----------------------------------------------------------------------------
void Biot_FeFv_Tri3::initializeNumericsAt( Cell* targetCell )
{
    targetCell->numericsStatus = new NumericsStatus_Biot_FeFv_Tri3();
    auto cns = this->getNumericsStatusAt(targetCell);
    
    // Pre-calculate det(J) and inv(J);
    RealMatrix Jmat = this->giveJacobianMatrixAt(targetCell, _gpNatCoor);
    cns->_Jdet = Jmat(0,0)*Jmat(1,1) - Jmat(1,0)*Jmat(0,1);
    
    // Sanity check
    if ( cns->_Jdet <= 0 )
        throw std::runtime_error("Calculation of negative area detected!\nSource: " + _name);
    
    cns->_JmatInv = inv(Jmat);
}
// ----------------------------------------------------------------------------
void Biot_FeFv_Tri3::readAdditionalDataFrom( FILE* fp )
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
void Biot_FeFv_Tri3::removeConstraintsOn( Cell* targetCell )
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
void Biot_FeFv_Tri3::setDofStagesAt( Cell* targetCell )
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
    Dof* dof_h = analysisModel().domainManager().giveCellDof(1, targetCell);
    analysisModel().dofManager().setStageFor(dof_h, _stage[0]);
    analysisModel().dofManager().setSubsystemFor(dof_h, _subsystem[0]);
}

// Private methods
// ----------------------------------------------------------------------------
NumericsStatus_Biot_FeFv_Tri3*
Biot_FeFv_Tri3::getNumericsStatusAt( Cell* targetCell )
{
    auto cns = dynamic_cast<NumericsStatus_Biot_FeFv_Tri3*>(targetCell->numericsStatus);
    if ( !cns )
        throw std::runtime_error("Error: Unable to retrieve numerics status at cell!\nSource: " + _name);
    
    return cns;
}
// ----------------------------------------------------------------------------
std::tuple< RealMatrix, RealVector >
Biot_FeFv_Tri3::giveBmatAt( Cell* targetCell )
{
    auto cns = this->getNumericsStatusAt(targetCell);
    RealMatrix dpsi = cns->_JmatInv*_basisFunctionDerivatives;
    
    // Fill entries for Bmat
    RealMatrix bmat({{dpsi(0,0), 0.,        dpsi(0,1), 0.,        dpsi(0,2), 0.},
                     {0.,        dpsi(1,0), 0.,        dpsi(1,1), 0.,        dpsi(1,2)},
                     {0.,        0.,        0.,        0.,        0.,        0.},
                     {dpsi(1,0), dpsi(0,0), dpsi(1,1), dpsi(0,1), dpsi(1,2), dpsi(0,2)}});
    
    // Fill entries for BmatDiv
    RealVector bmatDiv({dpsi(0,0), dpsi(1,0), dpsi(0,1), dpsi(1,1), dpsi(0,2), dpsi(1,2)});

    return std::make_tuple(std::move(bmat), std::move(bmatDiv));
}
// ----------------------------------------------------------------------------
std::vector<std::vector<Node*> > Biot_FeFv_Tri3::giveFaceNodesOf( Cell* targetCell )
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
RealMatrix Biot_FeFv_Tri3::giveGradBmatAt( Cell* targetCell )
{
    auto cns = this->getNumericsStatusAt(targetCell);
    return cns->_JmatInv*_basisFunctionDerivatives;
}
// ----------------------------------------------------------------------------
double Biot_FeFv_Tri3::giveDistanceToMidpointOf( std::vector<Node*>& face
                                               , RealVector&         coor )
{
    RealVector coor0, coor1, dx;
    coor0 = analysisModel().domainManager().giveCoordinatesOf(face[0]);
    coor1 = analysisModel().domainManager().giveCoordinatesOf(face[1]);
    dx = 0.5*(coor0 + coor1) - coor;
    
    return std::sqrt(dx.dot(dx));
}
// ----------------------------------------------------------------------------
RealMatrix Biot_FeFv_Tri3::giveJacobianMatrixAt( Cell* targetCell, const RealVector& natCoor )
{
    std::vector<Node*> node = analysisModel().domainManager().giveNodesOf(targetCell);
#ifndef NDEBUG
    if ( (int)node.size() != 3 )
        throw std::runtime_error("\nError: Basis function requires 3 nodes be specified for calculation of Jacobian matrix!\nSource: " + _name);
#endif
    
    RealMatrix coorMat(3,2);
    
#pragma GCC ivdep
    for ( int i = 0; i < 3; i++ )
    {
        RealVector nodeCoor = analysisModel().domainManager().giveCoordinatesOf(node[i]);
        coorMat(i,0) = nodeCoor(0);
        coorMat(i,1) = nodeCoor(1);
    }
    
   return _basisFunctionDerivatives*coorMat;
}
// ----------------------------------------------------------------------------
double Biot_FeFv_Tri3::giveLengthOf( std::vector<Node*>& face )
{
    RealVector coor0, coor1, dx;
    coor0 = analysisModel().domainManager().giveCoordinatesOf(face[0]);
    coor1 = analysisModel().domainManager().giveCoordinatesOf(face[1]);
    dx = coor1 - coor0;
    
    return std::sqrt(dx.dot(dx));
}
// ----------------------------------------------------------------------------
RealVector Biot_FeFv_Tri3::giveLocalDisplacementsAt( std::vector<Dof*>& dof, ValueType valType )
{
    // Displacements
    RealVector u(6);
    u(0) = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof[0], valType);
    u(1) = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof[1], valType);
    u(2) = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof[2], valType);
    u(3) = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof[3], valType);
    u(4) = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof[4], valType);
    u(5) = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof[5], valType);
    
    return u;
}
// ----------------------------------------------------------------------------
std::vector<Dof*> Biot_FeFv_Tri3::giveNodalDofsAt( Cell* targetCell )
{
    std::vector<Node*> node = analysisModel().domainManager().giveNodesOf(targetCell);
    
    std::vector<Dof*> dof(6, nullptr);
    dof[0] = analysisModel().domainManager().giveNodalDof(_nodalDof[0], node[0]);
    dof[1] = analysisModel().domainManager().giveNodalDof(_nodalDof[1], node[0]);
    dof[2] = analysisModel().domainManager().giveNodalDof(_nodalDof[0], node[1]);
    dof[3] = analysisModel().domainManager().giveNodalDof(_nodalDof[1], node[1]);
    dof[4] = analysisModel().domainManager().giveNodalDof(_nodalDof[0], node[2]);
    dof[5] = analysisModel().domainManager().giveNodalDof(_nodalDof[1], node[2]);
    
    return dof;
}
// ----------------------------------------------------------------------------
RealVector Biot_FeFv_Tri3::giveOutwardUnitNormalOf( std::vector<Node*>& face )
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
double Biot_FeFv_Tri3::giveTransmissibilityCoefficientAt
    ( std::vector<Node*>& face
    , Cell* targetCell
    , Cell* neighborCell )
{
    // Unit outward normal of face
    RealVector nhat;
    nhat = this->giveOutwardUnitNormalOf(face);
    
    double length = this->giveLengthOf(face);

    // A. Target cell
    // Coordinates of cell center
    std::vector<RealVector> ep = this->giveEvaluationPointsFor(targetCell);
    RealVector coor1 = ep[0];
    
    // Distance from cell center to midpoint of face
    double d1 = this->giveDistanceToMidpointOf(face, coor1);
    
    // Permeability tensor for cell
    std::vector<Material*> material = this->giveMaterialSetFor(targetCell);
    double k_xx, k_yy, k_xy;
    k_xx = material[2]->giveParameter("Permeability_xx");
    k_yy = material[2]->giveParameter("Permeability_yy");
    k_xy = material[2]->giveParameter("Permeability_xy");

    RealMatrix kmat({{k_xx, k_xy},
                     {k_xy, k_yy}});
                                 
    RealVector kvec;
    kvec = kmat*nhat;
    double k1 = std::sqrt(kvec.dot(kvec));
    
    // B. Neighbor cell
    // Coordinates of cell center
    ep = this->giveEvaluationPointsFor(neighborCell);
    RealVector coor2 = ep[0];
    
    // Distance from cell center to midpoint of face
    double d2 = this->giveDistanceToMidpointOf(face, coor2);
    
    // Permeability tensor for neigbor cell
    material = this->giveMaterialSetFor(neighborCell);
    k_xx = material[2]->giveParameter("Permeability_xx");
    k_yy = material[2]->giveParameter("Permeability_yy");
    k_xy = material[2]->giveParameter("Permeability_xy");

    kmat = {{k_xx, k_xy},
            {k_xy, k_yy}};
                     
    kvec = kmat*nhat;
    double k2 = std::sqrt(kvec.dot(kvec));
    
    return length/(_mu*d1/k1 + _mu*d2/k2);
}