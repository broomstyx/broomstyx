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

#include "PhaseFieldFracture_FeFv_Tri3.hpp"
#include <cmath>
#include <string>
#include "../../Core/AnalysisModel.hpp"
#include "../../Core/ObjectFactory.hpp"
#include "../../Core/DofManager.hpp"
#include "../../Core/DomainManager.hpp"
#include "../../Materials/Material.hpp"
#include "../../Util/linearAlgebra.hpp"
#include "../../Util/readOperations.hpp"

using namespace broomstyx;

registerBroomstyxObject(Numerics, PhaseFieldFracture_FeFv_Tri3)

// Numerics status
NumericsStatus_PhaseFieldFracture_FeFv_Tri3::NumericsStatus_PhaseFieldFracture_FeFv_Tri3()
    : _strain(RealVector(4))
    , _stress(RealVector(4))
    , _gradU(RealMatrix(2,2))
    , _surfEgy(0)
    , _bulkEgy(0)
    , _gradPhi(RealVector(2))
    , _hasPhsFldConstraint(false)
    , _hasPhsFldPrescribedOnFace {false, false, false}
    , _hasPhsFldGradientPrescribedOnFace {false, false, false}
    , _valueOnFace {0., 0., 0.}
    , _hasNotComputedTransmissibilities(true)
    , _transmissibility{0., 0., 0.}
    , _materialStatus {nullptr, nullptr}
{}
    
NumericsStatus_PhaseFieldFracture_FeFv_Tri3::~NumericsStatus_PhaseFieldFracture_FeFv_Tri3() {}

// Constructor
PhaseFieldFracture_FeFv_Tri3::PhaseFieldFracture_FeFv_Tri3()
{
    _dim = 2;
    _dofPerCell = 1;
    _dofPerNode = 2;
    
    _nNodes = 3;
    _nMaterials = 2;
    _nStages = 1;
    _nSubsystems = 2;
    
    _name = "PhaseFieldFracture_FeFv_Tri3";
    
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
}

// Destructor
PhaseFieldFracture_FeFv_Tri3::~PhaseFieldFracture_FeFv_Tri3() {}

// Public methods
// ----------------------------------------------------------------------------
void PhaseFieldFracture_FeFv_Tri3::deleteNumericsAt( Cell* targetCell )
{
    auto cns = this->getNumericsStatusAt(targetCell);
    auto material = this->giveMaterialSetFor(targetCell);
    
    material[0]->destroy(cns->_materialStatus[0]);
    material[1]->destroy(cns->_materialStatus[1]);
}
// ----------------------------------------------------------------------------
void PhaseFieldFracture_FeFv_Tri3::finalizeDataAt( Cell* targetCell )
{
    // Pointer to numerics status
    auto cns = this->getNumericsStatusAt(targetCell);
    
    // Get coordinate of cell center
    std::vector<RealVector> epCoor = this->giveEvaluationPointsFor(targetCell);
    RealVector cellCoor = epCoor[0];
    
    // Phase-field value at cell
    Dof* dof1 = analysisModel().domainManager().giveCellDof(_cellDof[0], targetCell);
    cns->_phi = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof1, converged_value);
    
    // A. Calculate phase-field flux at faces
    std::vector<Node*> node;
    node = analysisModel().domainManager().giveNodesOf(targetCell);
    
    std::vector<std::vector<Node*> > face = this->giveFaceNodesOf(targetCell);
    
    // Get cell neighbors
    std::vector<Cell*> neighbor = analysisModel().domainManager().giveNeighborsOf(targetCell);

    // Cycle through faces
    double flux_sum = 0;
    for ( int i = 0; i < 3; i++ )
    {
        // Length of face
        double length = this->giveLengthOf(face[i]);
        
        // Flux calculation
        double d_pf_n;
        
        if ( cns->_hasPhsFldPrescribedOnFace[i] )
        {
            // Distance from cell center to face midpoint
            double d1 = this->giveDistanceToMidpointOf(face[i], cellCoor);

            d_pf_n = ( cns->_valueOnFace[i] - cns->_phi )/d1;
            flux_sum += _Gc*_l*d_pf_n*length;
        }
        else if ( cns->_hasPhsFldGradientPrescribedOnFace[i] )
        {
            d_pf_n = cns->_valueOnFace[i];
            flux_sum += d_pf_n*length;
        }
        else 
        {
            // phase-field at neighbor cell
            Dof* dof2 = analysisModel().domainManager().giveCellDof(_cellDof[0], neighbor[i]);
            double phi2 = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof2, converged_value);
            
            d_pf_n = cns->_transmissibility[i]*(phi2 - cns->_phi)/(length*_Gc*_l);
            flux_sum += cns->_transmissibility[i]*(phi2 - cns->_phi);
        }
        
        cns->_valueOnFace[i] = d_pf_n;
    }
    
    // B. RT0 reconstruction of flux at cell center
    std::vector<Node*> vertexNode(3, nullptr);
    vertexNode[0] = node[2];
    vertexNode[1] = node[0];
    vertexNode[2] = node[1];
    
    RealVector d_pf(3);

    for ( int i = 0; i < 3; i++ )
    {
        RealVector temp = this->giveOutwardUnitNormalOf(face[i]);
        RealVector nhat(3); // 3rd dimension of nhat is zero
        nhat(0) = temp(0);
        nhat(1) = temp(1);
        
        RealVector x = analysisModel().domainManager().giveCoordinatesOf(node[i]);
        RealVector x0 = analysisModel().domainManager().giveCoordinatesOf(vertexNode[i]);
        
        double alpha = cns->_valueOnFace[i]/nhat.dot(x-x0);
        
        d_pf = d_pf + alpha*(cellCoor - x0);
    }
    
    cns->_gradPhi = {d_pf(0), d_pf(1)};
    
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
    RealMatrix bmatU = this->giveBmatAt(targetCell);
    cns->_strain = bmatU*uVec;
    
    // Retrieve material set for element
    std::vector<Material*> material = this->giveMaterialSetFor(targetCell);
    
    // Get constitutive force and elastic bulk energy
    RealVector conState({cns->_strain(0),
                         cns->_strain(1),
                         cns->_strain(2),
                         cns->_strain(3),
                         cns->_phi});

    material[1]->updateStatusFrom(conState, cns->_materialStatus[1], "FinalizeHistoryField");
    cns->_stress = material[1]->giveForceFrom(conState, cns->_materialStatus[1], "Mechanics");
    cns->_bulkEgy = material[1]->givePotentialFrom(conState, cns->_materialStatus[1]);
    cns->_surfEgy = _Gc/2.0*(_l*cns->_gradPhi.dot(cns->_gradPhi) + cns->_phi*cns->_phi/_l);
    
    RealVector fmatU;
    fmatU = _wt*cns->_Jdet*(trp(bmatU)*cns->_stress);
    
    analysisModel().dofManager().addToSecondaryVariableAt(dof[0], fmatU(0));
    analysisModel().dofManager().addToSecondaryVariableAt(dof[1], fmatU(1));
    analysisModel().dofManager().addToSecondaryVariableAt(dof[2], fmatU(2));
    analysisModel().dofManager().addToSecondaryVariableAt(dof[3], fmatU(3));
    analysisModel().dofManager().addToSecondaryVariableAt(dof[4], fmatU(4));
    analysisModel().dofManager().addToSecondaryVariableAt(dof[5], fmatU(5));
}
// ----------------------------------------------------------------------------
double PhaseFieldFracture_FeFv_Tri3::giveCellFieldValueAt( Cell* targetCell, int fieldNum )
{
    RealVector val, wt;
    std::string fieldTag;
    
    try
    {
        fieldTag = _cellFieldOutput.at(fieldNum);
    }
    catch (std::exception& e)
    {
        fieldTag = "unassigned";
    }
    
    std::tie(val,wt) = this->giveFieldOutputAt(targetCell, fieldTag);
    
    return val(0);
}
// ----------------------------------------------------------------------------
RealVector PhaseFieldFracture_FeFv_Tri3::giveCellNodeFieldValuesAt( Cell* targetCell, int fieldNum )
{
    double val = this->giveCellFieldValueAt(targetCell, fieldNum);
    RealVector cellNodeValues({val, val, val});
    
    return cellNodeValues;
}
// ----------------------------------------------------------------------------
std::vector<RealVector> PhaseFieldFracture_FeFv_Tri3::giveEvaluationPointsFor( Cell *targetCell )
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
std::tuple< RealVector, RealVector >
PhaseFieldFracture_FeFv_Tri3::giveFieldOutputAt( Cell* targetCell, const std::string& fieldTag )
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
        fieldVal(0) = cns->_gradU(0,1);
    else if ( fieldTag == "ux_y" )
        fieldVal(0) = cns->_gradU(1,0);
    else if ( fieldTag == "uy_y" )
        fieldVal(0) = cns->_gradU(1,1);
    else if ( fieldTag == "g_xy" )
        fieldVal(0) = cns->_strain(3);
    else if ( fieldTag == "pf" )
        fieldVal(0) = cns->_phi;
    else if ( fieldTag == "ene_s" )
        fieldVal(0) = cns->_surfEgy;
    else if ( fieldTag == "ene_b" )
        fieldVal(0) = cns->_bulkEgy;
    else if ( fieldTag == "cr_len" )
        fieldVal(0) = cns->_surfEgy/_Gc;
    else if ( fieldTag == "pf_x" )
        fieldVal(0) = cns->_gradPhi(0);
    else if ( fieldTag == "pf_y" )
        fieldVal(0) = cns->_gradPhi(1);
    else
        throw std::runtime_error("Invalid tag '" + fieldTag + "' supplied in field output request made to numerics 'CPEF3F'!");
    
    return std::make_tuple(std::move(fieldVal), std::move(weight));
}
// ----------------------------------------------------------------------------
RealVector PhaseFieldFracture_FeFv_Tri3::giveNumericsParameter( const std::string& paramTag )
{
    RealVector paramVec;
    
    if ( paramTag == "Gc_ell")
    {
        paramVec.init(2);
        paramVec(0) = _Gc;
        paramVec(1) = _l;
    }
    else if ( paramTag == "Gc" )
    {
        paramVec.init(1);
        paramVec(0) = _Gc;
    }
    else if ( paramTag == "ell" )
    {
        paramVec.init(1);
        paramVec(0) = _l;
    }
    else
        throw std::runtime_error("nERROR: Unknown parameter '" + paramTag + "' requested from numerics '" + _name + "'!\n");
    
    return paramVec;
}
// ----------------------------------------------------------------------------
std::tuple< std::vector<Dof*>, std::vector<Dof*>, RealVector >
PhaseFieldFracture_FeFv_Tri3::giveStaticCoefficientMatrixAt
        ( Cell*           targetCell
        , int             stage
        , int             subsys
        , const TimeData& time )
{
    std::vector<Dof*> rowDof, colDof;
    RealVector coefVal;
    
    if ( stage == _stage[0] )
    {
        auto cns = this->getNumericsStatusAt(targetCell);
        
        // Get cell neighbors
        std::vector<Cell*> neighbor = analysisModel().domainManager().giveNeighborsOf(targetCell);
        
        int nCols = 0;
        for ( int i = 0; i < 3; i++ )
            if ( !cns->_hasPhsFldPrescribedOnFace[i] )
                ++nCols;
        
        int vecLength = 37 + nCols; // Default output lengths for unassigned subsystems
        
        if ( subsys == _subsystem[0] )
            vecLength = 36;
        else if ( subsys == _subsystem[1] )
            vecLength = 1 + nCols;
        
        rowDof.assign(vecLength, nullptr);
        colDof.assign(vecLength, nullptr);
        coefVal.init(vecLength);
        
        // Retrieve nodal DOFs local to element
        std::vector<Dof*> dof = giveNodalDofsAt(targetCell);
        Dof* dof_phi = analysisModel().domainManager().giveCellDof(_cellDof[0], targetCell);
        
        // Retrieve material set for element
        std::vector<Material*> material = this->giveMaterialSetFor(targetCell);
        
        // Assemble constitutive state
        RealVector conState({cns->_strain(0),
                             cns->_strain(1),
                             cns->_strain(2),
                             cns->_strain(3),
                             cns->_phi});
        
        int counter = 0;
        if ( subsys == _subsystem[0] || subsys == UNASSIGNED )
        {
            // Get tangent modulus (linear momentum eq.)
            RealMatrix cmat = material[1]->giveModulusFrom(conState, cns->_materialStatus[1], "Mechanics");

            // Calculate stiffness KmatUU
            RealMatrix bmatU = this->giveBmatAt(targetCell);
            RealMatrix kmatUU = _wt*cns->_Jdet*trp(bmatU)*cmat*bmatU;

            for ( int i = 0; i < 6; i++ )
                for ( int j = 0; j < 6; j++ )
                {
                    rowDof[counter] = dof[i];
                    colDof[counter] = dof[j];
                    coefVal(counter) = kmatUU(i,j);
                    ++counter;
                }
        }
        
        if ( subsys == _subsystem[1] || subsys == UNASSIGNED )
        {
            int startIdx = counter;
            
            // Compute g''(phi)*elasticEnergy
            RealMatrix ddgPhi_Psi0 = material[1]->giveModulusFrom(conState, cns->_materialStatus[1], "PhaseField");
            
            double coef = _wt*cns->_Jdet*(_Gc/_l + ddgPhi_Psi0(0,0));
            
            rowDof[startIdx] = dof_phi;
            colDof[startIdx] = dof_phi;
            coefVal(startIdx) = coef;
            
            std::vector<std::vector<Node*> > face = this->giveFaceNodesOf(targetCell);
            
            // Cycle through faces
            for (int i = 0; i < 3; i++ )
            {
                double length = this->giveLengthOf(face[i]);
                if ( cns->_hasPhsFldPrescribedOnFace[i] )
                {
                    // Get coordinate of cell center
                    std::vector<RealVector> epCoor;
                    epCoor = this->giveEvaluationPointsFor(targetCell);
                    RealVector cellCoor = epCoor[0];
                    
                    // Distance from cell to face midpoint
                    double d1 = this->giveDistanceToMidpointOf(face[i], cellCoor);
                    
                    coefVal(startIdx) += _Gc*_l/d1*length;
                }
                else if ( !cns->_hasPhsFldGradientPrescribedOnFace[i] )
                {
                    // DOF address for neighbor cell
                    ++counter;
                    rowDof[counter] = dof_phi;
                    colDof[counter] = analysisModel().domainManager().giveCellDof(_cellDof[0], neighbor[i]);
                    
                    coefVal(startIdx) += cns->_transmissibility[i];
                    coefVal(counter) = -cns->_transmissibility[i];
                }
            }
        }
    }
    
    return std::make_tuple(std::move(rowDof), std::move(colDof), std::move(coefVal));
}
// ----------------------------------------------------------------------------
std::tuple< std::vector<Dof*>, RealVector >
PhaseFieldFracture_FeFv_Tri3::giveStaticLeftHandSideAt
        ( Cell*           targetCell
        , int             stage
        , int             subsys
        , const TimeData& time )
{
    std::vector<Dof*> rowDof;
    RealVector lhs;
    
    if ( stage == _stage[0] )
    {
        auto cns = this->getNumericsStatusAt(targetCell);
        
        // Retrieve nodal DOFs local to element
        std::vector<Dof*> dof = this->giveNodalDofsAt(targetCell);

        int vecLength = 9; // Default output lengths for unassigned subsystems
        
        if ( subsys == _subsystem[0] )
            vecLength = 6;
        else if ( subsys == _subsystem[1] )
            vecLength = 3;
        
        rowDof.assign(vecLength, nullptr);
        lhs.init(vecLength);
        
        // Retrieve current value of local displacements
        RealVector uVec = this->giveLocalDisplacementsAt(dof, current_value);
        
        // Retrieve phase-field
        Dof* dof_phi = analysisModel().domainManager().giveCellDof(_cellDof[0], targetCell);
        cns->_phi = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof_phi, current_value);
        
        // Compute local strains
        RealMatrix bmatU = this->giveBmatAt(targetCell);
        cns->_strain = bmatU*uVec;
        
        // Assemble constitutive state vector
        RealVector conState({cns->_strain(0), 
                             cns->_strain(1),
                             cns->_strain(2),
                             cns->_strain(3),
                             cns->_phi});
        
        // Retrieve material set for element
        std::vector<Material*> material = this->giveMaterialSetFor(targetCell);
        
        // Update material state
        material[1]->updateStatusFrom(conState, cns->_materialStatus[1]);
        
        // Compute stress and store in gpData for use in calculating tangent stiffness
        cns->_stress = material[1]->giveForceFrom(conState, cns->_materialStatus[1], "Mechanics");
        
        int offset = 0;
        if ( subsys == _subsystem[0] || subsys == UNASSIGNED )
        {
            rowDof[0] = dof[0];
            rowDof[1] = dof[1];
            rowDof[2] = dof[2];
            rowDof[3] = dof[3];
            rowDof[4] = dof[4];
            rowDof[5] = dof[5];
            offset = 6;
            
            // Calculate FmatU
            RealVector fmatU;
            fmatU = _wt*cns->_Jdet*trp(bmatU)*cns->_stress;
            
            lhs(0) = fmatU(0);
            lhs(1) = fmatU(1);
            lhs(2) = fmatU(2);
            lhs(3) = fmatU(3);
            lhs(4) = fmatU(4);
            lhs(5) = fmatU(5);
        }
        
        if ( subsys == _subsystem[1] || subsys == UNASSIGNED )
        {
            // Get coordinate of cell center and cell area
            std::vector<RealVector> epCoor = this->giveEvaluationPointsFor(targetCell);
            RealVector cellCoor = epCoor[0];

            // normal phase-field flux on edges
            RealVector normFlux(3);

            std::vector<Node*> node;
            node = analysisModel().domainManager().giveNodesOf(targetCell);

            std::vector<std::vector<Node*> > face = this->giveFaceNodesOf(targetCell);

            // Get cell neighbors
            std::vector<Cell*> neighbor = analysisModel().domainManager().giveNeighborsOf(targetCell);

            // Cycle through faces
            if ( cns->_hasNotComputedTransmissibilities )
            {
                for ( int i = 0; i < 3; i++ )
                    if ( !cns->_hasPhsFldPrescribedOnFace[i] && !cns->_hasPhsFldGradientPrescribedOnFace[i] )
                        cns->_transmissibility[i] = this->giveTransmissibilityCoefficientAt(face[i], targetCell, neighbor[i]);
                
                cns->_hasNotComputedTransmissibilities = false;
            }
            
            for ( int i = 0; i < 3; i++ )
            {
                if ( cns->_hasPhsFldPrescribedOnFace[i] ) // Cell face is a boundary where phase-field is specified
                {
                    double length = this->giveLengthOf(face[i]);
                    double d1 = this->giveDistanceToMidpointOf(face[i], cellCoor);

                    normFlux(i) = length*_Gc*_l*cns->_phi/d1;
                }
                else if ( !cns->_hasPhsFldGradientPrescribedOnFace[i] )
                {
                    // Phase-field at neighbor cell
                    Dof* dof_phi2 = analysisModel().domainManager().giveCellDof(_cellDof[0], neighbor[i]);
                    double phi2 = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof_phi2, current_value);

                    normFlux(i) = cns->_transmissibility[i]*(cns->_phi - phi2);
                }
            }
            
            // Update material state and calculate g'(phi)*elasticEnergy
            RealVector dgPhi_Psi0 = material[1]->giveForceFrom(conState, cns->_materialStatus[1], "PhaseField");
            
            rowDof[offset] = dof_phi;
            rowDof[offset + 1] = dof_phi;
            rowDof[offset + 2] = dof_phi;
            
            // Calculate FmatPhi
            lhs(offset) = _wt*cns->_Jdet*_Gc/_l*cns->_phi;
            lhs(offset + 1) = normFlux(0) + normFlux(1) + normFlux(2);
            lhs(offset + 2) = _wt*cns->_Jdet*dgPhi_Psi0(0);
        }
    }
    
    return std::make_tuple(std::move(rowDof), std::move(lhs));
}
// ----------------------------------------------------------------------------
std::tuple< std::vector<Dof*>, RealVector >
PhaseFieldFracture_FeFv_Tri3::giveStaticRightHandSideAt
        ( Cell*                    targetCell
        , int                      stage
        , int                      subsys
        , const BoundaryCondition& bndCond
        , const TimeData&          time )
{
    std::vector<Dof*> rowDof;
    RealVector rhs;
    
    if ( stage == _stage[0] )
    {
        if ( bndCond.conditionType() == "Traction" )
        {
            if ( subsys == _subsystem[0] || subsys == UNASSIGNED )
            {
                // ---------------------------------------------
                // Note: Boundary element must be a 2-node line!
                // ---------------------------------------------

                // Retrieve nodes of boundary element
                std::vector<Node*> node = analysisModel().domainManager().giveNodesOf(targetCell);

                if ( (int)node.size() != 2 )
                    throw std::runtime_error("Incompatible boundary cell encountered for '" + _name + "'!");

                // Compute length of boundary element
                RealVector coor0, coor1, dx;
                coor0 = analysisModel().domainManager().giveCoordinatesOf(node[0]);
                coor1 = analysisModel().domainManager().giveCoordinatesOf(node[1]);
                dx = coor1 - coor0;
                double length = std::sqrt(dx.dot(dx));

                // BC evaluated at midpoint of boundary cell
                RealVector midpt = 0.5*(coor0 + coor1);
                
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
        }
        else if ( bndCond.conditionType() == "PhaseField" )
        {
            std::vector<Cell*> domCell = analysisModel().domainManager().giveDomainCellsAssociatedWith(targetCell);
            
            for ( Cell* curDomCell : domCell )
            {
                Numerics* domCellNumerics = analysisModel().domainManager().giveNumericsFor(curDomCell);
                if ( this == domCellNumerics )
                {
                    auto cns = this->getNumericsStatusAt(curDomCell);
                    
                    std::vector<RealVector> cellCoor = this->giveEvaluationPointsFor(curDomCell);
                    std::vector<Node*> bndCellNode;

                    // Retrieve nodes of boundary cell
                    bndCellNode = analysisModel().domainManager().giveNodesOf(targetCell);

                    if ( bndCellNode.size() != 2 )
                        throw std::runtime_error("Numerics type '" + _name + "' only accepts boundary cells with two nodes for phase-field BCs!");

                    std::vector<std::vector<Node*> > face = this->giveFaceNodesOf(curDomCell);

                    // Initialize rowDof and rhs vector
                    rowDof.assign(1, nullptr);
                    rowDof[0] = analysisModel().domainManager().giveCellDof(_cellDof[0], curDomCell);
                    rhs.init(1);

                    // Cycle through faces to determine which one corresponds to current boundary cell
                    for ( int i = 0; i < 3; i++ )
                    {
                        if ( (bndCellNode[0] == face[i][0] && bndCellNode[1] == face[i][1]) || 
                             (bndCellNode[1] == face[i][0] && bndCellNode[0] == face[i][1]) )
                        {
                            // BC is evaluate at midpoint of boundary cell
                            RealVector coor0, coor1, midpt;
                            coor0 = analysisModel().domainManager().giveCoordinatesOf(bndCellNode[0]);
                            coor1 = analysisModel().domainManager().giveCoordinatesOf(bndCellNode[1]);
                            midpt = 0.5*(coor0 + coor1);
                            double bcVal = bndCond.valueAt(midpt, time);
                            
                            cns->_hasPhsFldPrescribedOnFace[i] = true;
                            cns->_valueOnFace[i] = bcVal;
                                
                            double length = this->giveLengthOf(face[i]);
                            double d = this->giveDistanceToMidpointOf(face[i], cellCoor[0]);
                            rhs(0) = length*_Gc*_l*bcVal/d;
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
PhaseFieldFracture_FeFv_Tri3::giveStaticRightHandSideAt
        ( Cell*                 targetCell
        , int                   stage
        , int                   subsys
        , const FieldCondition& fldCond
        , const TimeData&       time )
{
    std::vector<Dof*> rowDof;
    RealVector rhs;
    
    if ( stage == _stage[0] )
        if ( subsys == _subsystem[0] || subsys == UNASSIGNED )
        {
            auto cns = this->getNumericsStatusAt(targetCell);
            
            // Evaluate body forces arising from field conditions
            std::string name = fldCond.conditionType();
            if ( name == "Acceleration_X" || name == "Acceleration_Y" )
            {
                // Get element density
                int label = analysisModel().domainManager().giveLabelOf(targetCell);
                std::vector<Material*> material = analysisModel().domainManager().giveMaterialSetForDomain(label);
                double rho = material[0]->giveMaterialVariable("Density", cns->_materialStatus[0]);

                // Retrieve nodal DOFs for element
                std::vector<Dof*> dof = giveNodalDofsAt(targetCell);
                
                std::vector<RealVector> ep = this->giveEvaluationPointsFor(targetCell);
                double accel = fldCond.valueAt(ep[0], time);
                double force = accel*rho*_wt*cns->_Jdet/3.0;
                
                rhs = {force, force, force};
                
                rowDof.assign(3, nullptr);
                if ( name == "Acceleration_X" )
                {
                    rowDof[0] = dof[0];
                    rowDof[1] = dof[3];
                    rowDof[2] = dof[6];
                }
                else
                {
                    rowDof[0] = dof[1];
                    rowDof[1] = dof[4];
                    rowDof[2] = dof[7];
                }
            }
        }
    
    return std::make_tuple(std::move(rowDof), std::move(rhs));
}
// ----------------------------------------------------------------------------
void PhaseFieldFracture_FeFv_Tri3::imposeConstraintAt
        ( Cell*                    targetCell
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

        for ( int j = 0; j < (int)node.size(); j++ )
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
        auto cns = this->getNumericsStatusAt(targetCell);
        
        int targetDofNum = _cellDof[bndCond.targetDof() - 1];
        Dof* targetDof = analysisModel().domainManager().giveCellDof(targetDofNum, targetCell);
        
        std::vector<RealVector> ep = this->giveEvaluationPointsFor(targetCell);
        double bcVal = bndCond.valueAt(ep[0], time);
        analysisModel().dofManager().setConstraintValueAt(targetDof, bcVal);
        
        cns->_hasPhsFldConstraint = true;
    }
    
    if ( bndCond.conditionType() == "PhaseField" || bndCond.conditionType() == "PhaseFieldGradient" )
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
                    throw std::runtime_error("Numerics type 'PhsFieldBrittleFracture_HybFeFv_Tri3' only accepts boundary cells with exactly two nodes for phase-field BCs!");

                std::vector< std::vector<Node*> > face = this->giveFaceNodesOf(curDomCell);

                // Cycle through faces to determine which one corresponds to current boundary cell
                for ( int i = 0; i < 3; i++)
                {
                    if ( (bndCellNode[0] == face[i][0] && bndCellNode[1] == face[i][1]) ||
                         (bndCellNode[1] == face[i][0] && bndCellNode[0] == face[i][1]) )
                    {
                        // Set boundary flag for appropriate face
                        if ( bndCond.conditionType() == "PhaseField" )
                            cns->_hasPhsFldPrescribedOnFace[i] = true;
                        else
                            cns->_hasPhsFldGradientPrescribedOnFace[i] = true;
                    }
                }
            }
        }
    }
}
// ----------------------------------------------------------------------------
void PhaseFieldFracture_FeFv_Tri3::initializeMaterialsAt( Cell* targetCell )
{
    auto cns = this->getNumericsStatusAt(targetCell);
    std::vector<Material*> material = this->giveMaterialSetFor(targetCell);
    
    cns->_materialStatus[0] = material[0]->createMaterialStatus();
    cns->_materialStatus[1] = material[1]->createMaterialStatus();
}
// ----------------------------------------------------------------------------
void PhaseFieldFracture_FeFv_Tri3::initializeNumericsAt( Cell* targetCell )
{
    targetCell->numericsStatus = new NumericsStatus_PhaseFieldFracture_FeFv_Tri3();
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
void PhaseFieldFracture_FeFv_Tri3::performPostprocessingAt( Cell* targetCell, std::string tag )
{
    if ( tag == "InitializeHistoryField" )
    {
        // Pointer to numerics status
        auto cns = this->getNumericsStatusAt(targetCell);

        // Get coordinate of cell center
        std::vector<RealVector> epCoor = this->giveEvaluationPointsFor(targetCell);
        RealVector cellCoor = epCoor[0];

        // Phase-field value at cell
        Dof* dof1 = analysisModel().domainManager().giveCellDof(_cellDof[0], targetCell);
        double phi = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof1, converged_value);

        // A. Calculate phase-field flux at faces
        std::vector<Node*> node;
        node = analysisModel().domainManager().giveNodesOf(targetCell);

        std::vector<std::vector<Node*> > face = this->giveFaceNodesOf(targetCell);

        // Get cell neighbors
        std::vector<Cell*> neighbor = analysisModel().domainManager().giveNeighborsOf(targetCell);

        // Cycle through faces
        double flux_sum = 0;
        for ( int i = 0; i < 3; i++ )
        {
            // Length of face
            double length = this->giveLengthOf(face[i]);

            // Flux calculation
            double d_pf_n;
            if ( cns->_hasPhsFldPrescribedOnFace[i] )
            {
                // Distance from cell center to face midpoint
                double d1 = this->giveDistanceToMidpointOf(face[i], cellCoor);

                d_pf_n = ( cns->_valueOnFace[i] - cns->_phi )/d1;
                flux_sum += _Gc*_l*d_pf_n*length;
            }
            else if ( cns->_hasPhsFldGradientPrescribedOnFace[i] )
            {
                d_pf_n = cns->_valueOnFace[i];
                flux_sum += d_pf_n*length;
            }
            else 
            {
                Dof* dof2 = analysisModel().domainManager().giveCellDof(_cellDof[0], neighbor[i]);
                double phi2 = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof2, converged_value);

                flux_sum += cns->_transmissibility[i]*(phi2 - cns->_phi);
            }
        }

        // Required history field value (note sign in front of flux_sum due to conventions used above)
        double gprimeH = -(-flux_sum/(_wt*cns->_Jdet) + _Gc/_l*phi);
        
        RealVector conState({phi, gprimeH});
        
        // Initialize history field
        std::vector<Material*> material = this->giveMaterialSetFor(targetCell);
        material[1]->updateStatusFrom(conState, cns->_materialStatus[1], "InitializeHistoryField");
    }
    else
        throw std::runtime_error("Unrecognized post-processing tag '" + tag + "'!\nSource: " + _name);
}
// ----------------------------------------------------------------------------
void PhaseFieldFracture_FeFv_Tri3::printPostIterationMessage( int stage )
{
    // Outputs the number of cells where phase-field is larger than 
    // irreversibility threshold
    
    int nCells = analysisModel().domainManager().giveNumberOfDomainCells();
    int crackCellCount = 0;
    
#ifdef _OPENMP
#pragma omp parallel for reduction (+:crackCellCount)
#endif
    for ( int i = 0; i < nCells; i++ )
    {
        Cell* curCell = analysisModel().domainManager().giveDomainCell(i);
        Numerics* numerics = analysisModel().domainManager().giveNumericsFor(curCell);
        
        if ( numerics == this )
        {
            Dof* phiDof = analysisModel().domainManager().giveCellDof(_cellDof[0], curCell);
            double phi = analysisModel().dofManager().giveValueOfPrimaryVariableAt(phiDof, current_value);
            
            if ( phi >= 0.9 )
                crackCellCount += 1;
        }
    }
    
    std::printf("\n    Phase-field threshold (0.9) exceeded in %d cells.\n", crackCellCount);
}
// ----------------------------------------------------------------------------
void PhaseFieldFracture_FeFv_Tri3::readAdditionalDataFrom( FILE* fp )
{
    std::string key, src = "PhsFieldBrittleFracture_HybFeFv_Tri3 (Numerics)";
    
    verifyKeyword(fp, key = "CharacteristicLength", src);
    _l = getRealInputFrom(fp, "Failed to read characteristic length from input file!", src);
    
    verifyKeyword(fp, key = "CriticalEnergyReleaseRate", src);
    _Gc = getRealInputFrom(fp, "Failed to read critical energy release rate from input file!", src);
}
// ----------------------------------------------------------------------------
void PhaseFieldFracture_FeFv_Tri3::removeConstraintsOn( Cell* targetCell )
{
    auto cns = this->getNumericsStatusAt(targetCell);
    cns->_hasPhsFldConstraint = false;
    cns->_hasPhsFldPrescribedOnFace[0] = false;
    cns->_hasPhsFldPrescribedOnFace[1] = false;
    cns->_hasPhsFldPrescribedOnFace[2] = false;
    cns->_hasPhsFldGradientPrescribedOnFace[0] = false;
    cns->_hasPhsFldGradientPrescribedOnFace[1] = false;
    cns->_hasPhsFldGradientPrescribedOnFace[2] = false;
    cns->_hasNotComputedTransmissibilities = true;
}
// ----------------------------------------------------------------------------
void PhaseFieldFracture_FeFv_Tri3::setDofStagesAt( Cell* targetCell )
{
    // A. Nodal DOFs
    std::vector<Node*> node = analysisModel().domainManager().giveNodesOf(targetCell);
    for ( int i = 0; i < (int)node.size(); i++ )
    {
        Dof *dof_x, *dof_y;
        
        dof_x = analysisModel().domainManager().giveNodalDof(_nodalDof[0], node[i]);
        dof_y = analysisModel().domainManager().giveNodalDof(_nodalDof[1], node[i]);
        
        // Set DOF stage numbers
        analysisModel().dofManager().setStageFor(dof_x, _stage[0]);
        analysisModel().dofManager().setStageFor(dof_y, _stage[0]);
        
        // Set DOF subsystem numbers
        analysisModel().dofManager().setSubsystemFor(dof_x, _subsystem[0]);
        analysisModel().dofManager().setSubsystemFor(dof_y, _subsystem[0]);
    }
        
    // B. Element DOFs
    Dof* dof_phi = analysisModel().domainManager().giveCellDof(_cellDof[0], targetCell);
    analysisModel().dofManager().setStageFor(dof_phi, _stage[0]);
    analysisModel().dofManager().setSubsystemFor(dof_phi, _subsystem[1]);
}

// Private methods
// ----------------------------------------------------------------------------
NumericsStatus_PhaseFieldFracture_FeFv_Tri3*
PhaseFieldFracture_FeFv_Tri3::getNumericsStatusAt( Cell* targetCell )
{
    auto cns = dynamic_cast<NumericsStatus_PhaseFieldFracture_FeFv_Tri3*>(targetCell->numericsStatus);
    if ( !cns )
        throw std::runtime_error("Error: Unable to retrieve numerics status at cell!\nSource: " + _name);
    
    return cns;
}
// ----------------------------------------------------------------------------
RealMatrix PhaseFieldFracture_FeFv_Tri3::giveBmatAt( Cell* targetCell )
{
    auto cns = this->getNumericsStatusAt(targetCell);
    RealMatrix dpsi = cns->_JmatInv*_basisFunctionDerivatives;
    
    // Fill entries for Bmat
    RealMatrix bmat({{dpsi(0,0), 0.,        dpsi(0,1), 0.,        dpsi(0,2), 0.},
                     {0.,        dpsi(1,0), 0.,        dpsi(1,1), 0.,        dpsi(1,2)},
                     {0.,        0.,        0.,        0.,        0.,        0.},
                     {dpsi(1,0), dpsi(0,0), dpsi(1,1), dpsi(0,1), dpsi(1,2), dpsi(0,2)}});
    return bmat;
}
// ----------------------------------------------------------------------------
std::vector<std::vector<Node*> > PhaseFieldFracture_FeFv_Tri3::giveFaceNodesOf( Cell* targetCell )
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
RealMatrix PhaseFieldFracture_FeFv_Tri3::giveGradBmatAt( Cell* targetCell )
{
    auto cns = this->getNumericsStatusAt(targetCell);
    return cns->_JmatInv*_basisFunctionDerivatives;
}
// ----------------------------------------------------------------------------
double PhaseFieldFracture_FeFv_Tri3::giveDistanceToMidpointOf( std::vector<Node*>& face
                                                               , RealVector&         coor )
{
    RealVector coor0, coor1, dx;
    coor0 = analysisModel().domainManager().giveCoordinatesOf(face[0]);
    coor1 = analysisModel().domainManager().giveCoordinatesOf(face[1]);
    dx = 0.5*(coor0 + coor1) - coor;
    
    return std::sqrt(dx.dot(dx));
}
// ----------------------------------------------------------------------------
RealMatrix PhaseFieldFracture_FeFv_Tri3::giveJacobianMatrixAt( Cell* targetCell, const RealVector& natCoor )
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
double PhaseFieldFracture_FeFv_Tri3::giveLengthOf( std::vector<Node*>& face )
{
    RealVector coor0, coor1, dx;
    coor0 = analysisModel().domainManager().giveCoordinatesOf(face[0]);
    coor1 = analysisModel().domainManager().giveCoordinatesOf(face[1]);
    dx = coor1 - coor0;
    
    return std::sqrt(dx.dot(dx));
}
// ----------------------------------------------------------------------------
RealVector PhaseFieldFracture_FeFv_Tri3::giveLocalDisplacementsAt( std::vector<Dof*>& dof, ValueType valType )
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
std::vector<Dof*> PhaseFieldFracture_FeFv_Tri3::giveNodalDofsAt( Cell* targetCell )
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
RealVector PhaseFieldFracture_FeFv_Tri3::giveOutwardUnitNormalOf( std::vector<Node*>& face )
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
double PhaseFieldFracture_FeFv_Tri3::giveTransmissibilityCoefficientAt
        ( std::vector<Node*>& face
        , Cell*               targetCell
        , Cell*               neighborCell )
{
    // A. Target cell
    // Coordinates of cell center
    std::vector<RealVector> ep = this->giveEvaluationPointsFor(targetCell);
    RealVector coor1 = ep[0];
    
    // length of face
    double length = this->giveLengthOf(face);
    
    // Distance from cell center to midpoint of face
    double d1 = this->giveDistanceToMidpointOf(face, coor1);
    
    // Gc_l
    double GcEll_1 = _Gc*_l;
    
    // B. Neighbor cell
    // Coordinates of cell center
    ep = this->giveEvaluationPointsFor(neighborCell);
    RealVector coor2 = ep[0];
    
    // Distance from cell center to midpoint of face
    double d2 = this->giveDistanceToMidpointOf(face, coor2);
    
    // Gc_l
    int label = analysisModel().domainManager().giveLabelOf(neighborCell);
    Numerics* numerics = analysisModel().domainManager().giveNumericsForDomain(label);
    RealVector param = numerics->giveNumericsParameter("Gc_ell");
    double GcEll_2 = param(0)*param(1);
    
    return length/(d1/GcEll_1 + d2/GcEll_2);
}