#include "HeatTransfer_Fe_Tri3.hpp"
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

registerBroomstyxObject(Numerics, HeatTransfer_Fe_Tri3)

// Numerics status
NumericsStatus_HeatTransfer_Fe_Tri3::NumericsStatus_HeatTransfer_Fe_Tri3()
    : _T(RealVector(1))
    , _gradT(RealVector(2))
    , _materialStatus {nullptr}
{}

NumericsStatus_HeatTransfer_Fe_Tri3::~NumericsStatus_HeatTransfer_Fe_Tri3()
{}

// Constructor
HeatTransfer_Fe_Tri3::HeatTransfer_Fe_Tri3()
{    
    _dim = 2;
    _dofPerNode = 1;
    _nNodes = 3;
    _nMaterials = 3;
    _nStages = 1;
    _nSubsystems = 1;
    _name = "HeatTransfer_Fe_Tri3";
    
    // Lone Gauss point coordinate is at natural coordinate (1/3, 1/3), weight = 0.5
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

    // Pre-define mass matrix
    _massMatrix = {{1./6.,  1./12., 1./12.},
                   {1./12., 1./6.,  1./12.},
                   {1./12., 1./12., 1./6.}};
}

// Destructor
HeatTransfer_Fe_Tri3::~HeatTransfer_Fe_Tri3()
{}

// Private methods
// ----------------------------------------------------------------------------
void HeatTransfer_Fe_Tri3::deleteNumericsAt(Cell* targetCell)
{
    auto cns = this->getNumericsStatusAt(targetCell);
    delete cns;
}
// ----------------------------------------------------------------------------
void HeatTransfer_Fe_Tri3::finalizeDataAt( Cell* targetCell )
{
    // Retrieve numerics status at cell
    auto cns = this->getNumericsStatusAt(targetCell);
    
    // Retrieve nodal DOFs local to element
    std::vector<Dof*> dof = this->giveNodalDofsAt(targetCell);
    
    // Retrieve DOF values
    RealVector TVec = this->giveLocalValuesAt(dof, converged_value);

    // Calculate temperature at midpoint (Gauss point location)
    cns->_T(0) = (TVec(0) + TVec(1) + TVec(2))/3.;
    
    // Construct gradient
    RealMatrix bmat = this->giveBmatAt(targetCell);
    cns->_gradT = bmat*TVec;
}
// ----------------------------------------------------------------------------
double HeatTransfer_Fe_Tri3::giveCellFieldValueAt( Cell* targetCell, int fieldNum )
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
RealVector HeatTransfer_Fe_Tri3::giveCellNodeFieldValuesAt( Cell* targetCell, int fieldNum )
{
    double val = this->giveCellFieldValueAt(targetCell, fieldNum);
    
    RealVector cellNodeValues(3);
    cellNodeValues(0) = val;
    cellNodeValues(1) = val;
    cellNodeValues(2) = val;
    
    return cellNodeValues;
}
// ----------------------------------------------------------------------------
std::vector<RealVector> HeatTransfer_Fe_Tri3::giveEvaluationPointsFor( Cell *targetCell)
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
HeatTransfer_Fe_Tri3::giveFieldOutputAt( Cell* targetCell, const std::string& fieldTag )
{
    RealVector fieldVal(1), weight(1);
    
    auto cns = this->getNumericsStatusAt(targetCell);
    weight(0) = _wt*cns->_Jdet;
    
    std::vector<Material*> material = this->giveMaterialSetFor(targetCell);
    
    if ( fieldTag == "unassigned" )
        fieldVal(0) = 0.;
    else if ( fieldTag == "T_x" )
        fieldVal(0) = cns->_gradT(0);
    else if ( fieldTag == "T_y" )
        fieldVal(0) = cns->_gradT(1);
    else
        throw std::runtime_error("Invalid tag '" + fieldTag + "' supplied in field output request made to numerics '" + _name + "'!");
    
    return std::make_tuple(std::move(fieldVal), std::move(weight));
}
// ----------------------------------------------------------------------------
std::tuple< std::vector<Dof*>
          , std::vector<Dof*>
          , RealVector >
HeatTransfer_Fe_Tri3::giveStaticCoefficientMatrixAt( Cell*           targetCell
                                              , int             stage
                                              , int             subsys
                                              , const TimeData& time )
{
    std::vector<Dof*> rowDof, colDof;
    RealVector coefVal;
    
    if ( stage == _stage[0] && (subsys == _subsystem[0] || subsys == UNASSIGNED) )
    {
        auto cns = this->getNumericsStatusAt(targetCell);

        // Retrieve nodal DOFs local to element
        std::vector<Dof*> dof = this->giveNodalDofsAt(targetCell);

        // Retrieve material set for element
        std::vector<Material*> material = this->giveMaterialSetFor(targetCell);

        // Get tangent modulus
        RealMatrix cmat = material[2]->giveModulusFrom(cns->_T, cns->_materialStatus);
        
        // Calculate stiffness matrix
        RealMatrix bmat = this->giveBmatAt(targetCell);
        RealMatrix kmat = _wt*cns->_Jdet*trp(bmat)*cmat*bmat;
        
        rowDof.assign(9, nullptr);
        colDof.assign(9, nullptr);
        coefVal.init(9);
        
        int counter = 0;
        for ( int i = 0; i < 3; i++)
            for ( int j = 0; j < 3; j++)
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
HeatTransfer_Fe_Tri3::giveStaticLeftHandSideAt(Cell*            targetCell
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
        
        // Element area and local displacements
        rowDof = this->giveNodalDofsAt(targetCell);
        RealVector T = this->giveLocalValuesAt(rowDof, current_value);
            
        // Compute local temperature and gradients
        RealMatrix bmat = giveBmatAt(targetCell);
        cns->_T(0) = (T(0) + T(1) + T(2))/3.;
        cns->_gradT = bmat*T;

        // Retrieve material set for element
        std::vector<Material*> material = this->giveMaterialSetFor(targetCell);

        // Get tangent modulus
        RealMatrix cmat = material[2]->giveModulusFrom(cns->_T, cns->_materialStatus);
        
        // Calculate lhs
        lhs = _wt*cns->_Jdet*trp(bmat)*cmat*cns->_gradT;
    }
    
    return std::make_tuple(std::move(rowDof), std::move(lhs));
}
// ----------------------------------------------------------------------------
std::tuple< std::vector<Dof*>
          , RealVector >
HeatTransfer_Fe_Tri3::giveStaticRightHandSideAt( Cell*                    targetCell
                                          , int                      stage
                                          , int                      subsys
                                          , const BoundaryCondition& bndCond
                                          , const TimeData&          time )
{   
    std::vector<Dof*> rowDof;
    RealVector rhs;
    
    if ( stage == _stage[0] && (subsys == _subsystem[0] || subsys == UNASSIGNED) )
    {
        std::string cndType = bndCond.conditionType();
        if ( cndType == "NormalHeatFlux" )
        {
            // Retrieve nodes of boundary element
            std::vector<Node*> node = analysisModel().domainManager().giveNodesOf(targetCell);
            
            // ----------------------------------------------------------------
            // Note: Boundary element must be a 2-node line
            // ----------------------------------------------------------------
            
            if ( (int)node.size() != 2 )
                throw std::runtime_error("Error: Neumann boundary condition for '" + _name + "' requires 2-node boundary elements!");
            
            RealVector coor0, coor1, dx;
            coor0 = analysisModel().domainManager().giveCoordinatesOf(node[0]);
            coor1 = analysisModel().domainManager().giveCoordinatesOf(node[1]);
            dx = coor1 - coor0;
            double length = std::sqrt(dx.dot(dx));

            // Determine proper value of boundary condition
            RealVector midpt = 0.5*(coor0 + coor1);
            double bcVal = bndCond.valueAt(midpt, time);

            // Construct local RHS vector (only for relevant DOFs)
            rhs.init(2);
            rhs(0) = 0.5*length*bcVal;
            rhs(1) = 0.5*length*bcVal;

            // Construct global address vector
            rowDof.assign(2, nullptr);

            int dofNum = analysisModel().dofManager().giveIndexForNodalDof(bndCond.targetDof());
            rowDof[0] = analysisModel().domainManager().giveNodalDof(dofNum, node[0]);
            rowDof[1] = analysisModel().domainManager().giveNodalDof(dofNum, node[1]);
        }
        else if ( cndType == "PointSource" )
        {
            // Retrieve nodes of boundary element
            std::vector<Node*> node = analysisModel().domainManager().giveNodesOf(targetCell);

            // ----------------------------------------------------------------
            // Note: Boundary element must be a 1-node point
            // ----------------------------------------------------------------

            if ( (int)node.size() != 1 )
                throw std::runtime_error("Error: Point source definition for '" + _name + "' requires 1-node elements!");

            RealVector coor0;
            coor0 = analysisModel().domainManager().giveCoordinatesOf(node[0]);

            // Determine proper value of boundary condition
            double bcVal = bndCond.valueAt(coor0, time);

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
HeatTransfer_Fe_Tri3::giveStaticRightHandSideAt( Cell*                 targetCell
                                          , int                   stage
                                          , int                   subsys
                                          , const FieldCondition& fldCond
                                          , const TimeData&       time )
{
    std::vector<Dof*> rowDof;
    RealVector rhs;
    
    if ( stage == _stage[0] && (subsys == _subsystem[0] || subsys == UNASSIGNED) )
    {
        // Evaluate source term defined over whole domain
        std::string cndType = fldCond.conditionType();
        if ( cndType == "SourceTerm" )
        {
            auto cns = this->getNumericsStatusAt(targetCell);
            
            // Retrieve nodal DOFs for element
            rowDof = this->giveNodalDofsAt(targetCell);
            
            // Note that field condition is assumed to be piecewise linear over each cell
            std::vector<RealVector> coor = this->giveEvaluationPointsFor(targetCell);
            double val = fldCond.valueAt(coor[0], time)*_wt*cns->_Jdet/3.0;
            
            rhs.init(3);
            rhs(0) = val;
            rhs(1) = val;
            rhs(2) = val;
        }
    }
    
    return std::make_tuple(std::move(rowDof), std::move(rhs));
}
// ----------------------------------------------------------------------------
std::tuple< std::vector<Dof*>, std::vector<Dof*>, RealVector >
HeatTransfer_Fe_Tri3::giveTransientCoefficientMatrixAt( Cell*           targetCell
                                                      , int             stage
                                                      , int             subsys
                                                      , const TimeData& time )
{
    std::vector<Dof*> rowDof, colDof;
    RealVector coefVal;

    if ( stage == _stage[0] && (subsys == _subsystem[0] || subsys == UNASSIGNED) )
    {
        rowDof.assign(9, nullptr);
        colDof.assign(9, nullptr);
        coefVal.init(9);

        auto cns = this->getNumericsStatusAt(targetCell);

        std::vector<Dof*> nodalDof = this->giveNodalDofsAt(targetCell);

        // Retrieve material set for element
        std::vector<Material*> material = this->giveMaterialSetFor(targetCell);

        double rho = material[0]->giveParameter("Density");
        double cp = material[1]->giveParameter("HeatCapacity");

        RealMatrix Mmat = _wt*cns->_Jdet*rho*cp*_massMatrix;

        rowDof[0] = nodalDof[0];    colDof[0] = nodalDof[0];
        rowDof[1] = nodalDof[0];    colDof[1] = nodalDof[1];
        rowDof[2] = nodalDof[0];    colDof[2] = nodalDof[2];
        rowDof[3] = nodalDof[1];    colDof[3] = nodalDof[0];
        rowDof[4] = nodalDof[1];    colDof[4] = nodalDof[1];
        rowDof[5] = nodalDof[1];    colDof[5] = nodalDof[2];
        rowDof[6] = nodalDof[2];    colDof[6] = nodalDof[0];
        rowDof[7] = nodalDof[2];    colDof[7] = nodalDof[1];
        rowDof[8] = nodalDof[2];    colDof[8] = nodalDof[2];

        coefVal = {Mmat(0,0), Mmat(0,1), Mmat(0,2),
                   Mmat(1,0), Mmat(1,1), Mmat(1,2),
                   Mmat(2,0), Mmat(2,1), Mmat(2,2)};
    }

    return std::make_tuple(std::move(rowDof), std::move(colDof), std::move(coefVal));
}
// ----------------------------------------------------------------------------
std::tuple< std::vector<Dof*>, RealVector >
HeatTransfer_Fe_Tri3::giveTransientLeftHandSideAt( Cell*           targetCell
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

        rowDof = this->giveNodalDofsAt(targetCell);
        RealVector T = this->giveLocalValuesAt(rowDof, valType);

        // Retrieve material set for element
        std::vector<Material*> material = this->giveMaterialSetFor(targetCell);

        double rho = material[0]->giveParameter("Density");
        double cp = material[1]->giveParameter("HeatCapacity");

        lhs = _wt*cns->_Jdet*rho*cp*_massMatrix*T;
    }

    return std::make_tuple(std::move(rowDof), std::move(lhs));
}
// ----------------------------------------------------------------------------
void HeatTransfer_Fe_Tri3::imposeConstraintAt( Cell*                    targetCell
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
void HeatTransfer_Fe_Tri3::initializeNumericsAt( Cell* targetCell )
{
    targetCell->numericsStatus = new NumericsStatus_HeatTransfer_Fe_Tri3();
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
void HeatTransfer_Fe_Tri3::setDofStagesAt( Cell* targetCell )
{
    // Element uses only one stage i.e. stage[0], and all nodal degrees of 
    // freedom get assigned to this stage
    std::vector<Node*> node = analysisModel().domainManager().giveNodesOf(targetCell);
    
    for ( int i = 0; i < (int)node.size(); i++)
    {
        Dof* dof = analysisModel().domainManager().giveNodalDof(_nodalDof[0],node[i]);
        analysisModel().dofManager().setStageFor(dof, _stage[0]);
    }
}

// Private methods
// ---------------------------------------------------------------------------
NumericsStatus_HeatTransfer_Fe_Tri3*
HeatTransfer_Fe_Tri3::getNumericsStatusAt( Cell* targetCell )
{
    auto cns = dynamic_cast<NumericsStatus_HeatTransfer_Fe_Tri3*>(targetCell->numericsStatus);
    if ( !cns )
        throw std::runtime_error("Error: Unable to retrieve numerics status at cell!\nSource: " + _name);
    
    return cns;
}
// ----------------------------------------------------------------------------
RealMatrix HeatTransfer_Fe_Tri3::giveBmatAt( Cell* targetCell )
{
    auto cns = this->getNumericsStatusAt(targetCell);
    return cns->_JmatInv*_basisFunctionDerivatives;
}
// ----------------------------------------------------------------------------
RealMatrix HeatTransfer_Fe_Tri3::giveJacobianMatrixAt( Cell* targetCell, const RealVector& natCoor )
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
RealVector HeatTransfer_Fe_Tri3::giveLocalValuesAt( std::vector<Dof*>& dof, ValueType valType)
{    
    RealVector u(3);
    u(0) = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof[0], valType);
    u(1) = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof[1], valType);
    u(2) = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof[2], valType);

    return u;
}
// ----------------------------------------------------------------------------
std::vector<Dof*> HeatTransfer_Fe_Tri3::giveNodalDofsAt( Cell* targetCell )
{
    std::vector<Node*> node = analysisModel().domainManager().giveNodesOf(targetCell);
    std::vector<Dof*> dof(3, nullptr);

    dof[0] = analysisModel().domainManager().giveNodalDof(_nodalDof[0], node[0]);
    dof[1] = analysisModel().domainManager().giveNodalDof(_nodalDof[0], node[1]);
    dof[2] = analysisModel().domainManager().giveNodalDof(_nodalDof[0], node[2]);

    return dof;
}
