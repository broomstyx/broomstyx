#include "DarcyFlow_2D_1Phase_Fv_Tri.hpp"
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

registerBroomstyxObject(Numerics, DarcyFlow_2D_1Phase_Fv_Tri)

// Numerics status constructor
NumericsStatus_DarcyFlow_2D_1Phase_Fv_Tri::NumericsStatus_DarcyFlow_2D_1Phase_Fv_Tri()
    : _centerFlux{0, 0}
    , _headOnFace{0, 0, 0}
    , _fluxOnFace{0, 0, 0}
    , _fluxIsPrescribedOnFace{false, false, false}
    , _headIsPrescribedOnFace{false, false, false}
    , _hasHeadConstraint(false)
    , _materialStatus(nullptr)
{}

NumericsStatus_DarcyFlow_2D_1Phase_Fv_Tri::~NumericsStatus_DarcyFlow_2D_1Phase_Fv_Tri() {}

// Constructor
DarcyFlow_2D_1Phase_Fv_Tri::DarcyFlow_2D_1Phase_Fv_Tri()
{
    _dim = 2;
    _dofPerCell = 1;
    _nNodes = 3;
    _nMaterials = 1;
    _nStages = 1;
    _nSubsystems = 1;
    
    _name = "DarcyFlow_2D_1Phase_Fv_Tri";
}

// Destructor
DarcyFlow_2D_1Phase_Fv_Tri::~DarcyFlow_2D_1Phase_Fv_Tri() {}

// Public methods
// -------------------------------------------------------------------------
void DarcyFlow_2D_1Phase_Fv_Tri::deleteNumericsAt( Cell* targetCell )
{
    auto cns = this->getNumericsStatusAt(targetCell);
    auto material = this->giveMaterialSetFor(targetCell);
    
    material[0]->destroy(cns->_materialStatus);
    delete cns;
}
// -------------------------------------------------------------------------
void DarcyFlow_2D_1Phase_Fv_Tri::finalizeDataAt( Cell* targetCell )
{
//    double area = this->giveAreaOf(targetCell);
    auto cns = this->getNumericsStatusAt(targetCell);
    
    // Get coordinate of cell center
    std::vector<RealVector> epCoor = this->giveEvaluationPointsFor(targetCell);
    RealVector cellCoor = epCoor[0];
    
    // Calculate vertical coordinate at cell center
//    double z = this->giveVerticalCoordinateAt(targetCell);
    
    // Retrieve material set for cell
    std::vector<Material*> material = this->giveMaterialSetFor(targetCell);
    
    // Get cell vertex nodes and faces
    std::vector<Node*> node = analysisModel().domainManager().giveNodesOf(targetCell);
    std::vector< std::vector<Node*> > face = this->giveFaceNodesOf(targetCell);
    
    // Get cell neighbors
    std::vector<Cell*> neighbor = analysisModel().domainManager().giveNeighborsOf(targetCell);
    
    // A. Calculate normal flux at faces
    // Permeability tensor for cell
    double k_xx, k_yy, k_xy;
    k_xx = material[0]->giveParameter("Permeability_xx");
    k_yy = material[0]->giveParameter("Permeability_yy");
    k_xy = material[0]->giveParameter("Permeability_xy");
    
    RealMatrix kmat_cell({{k_xx, k_xy},
                          {k_xy, k_yy}});
    
    // Hydraulic potential at cell
    Dof* dof_h1 = analysisModel().domainManager().giveCellDof(_cellDof[0], targetCell);
    double h1 = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof_h1, converged_value);
    
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

            normalFlux = K1*(h1 - cns->_headOnFace[i])/d1;
        }
        else if ( cns->_fluxIsPrescribedOnFace[i] )
            normalFlux = cns->_fluxOnFace[i];
        else 
        {
            // Hydraulic potential at neighbor cell
            Dof* dof2 = analysisModel().domainManager().giveCellDof(_cellDof[0], neighbor[i]);
            double h2 = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof2, converged_value);

            double transCoef = this->giveTransmissibilityCoefficientAt(face[i], targetCell, neighbor[i]);
            normalFlux = transCoef*(h1 - h2)/length;
        }
        cns->_fluxOnFace[i] = normalFlux;
    }
    
    // RT0 reconstruction of flux at cell center
    std::vector<Node*> vertexNode(3, nullptr);
    vertexNode[0] = node[2];
    vertexNode[1] = node[0];
    vertexNode[2] = node[1];
    
    RealVector centerFlux(3);
    for ( int i = 0; i < 3; i++ )
    {
        RealVector temp = this->giveOutwardUnitNormalOf(face[i]);
        RealVector nhat({temp(0), temp(1), 0.});
        
        RealVector x = analysisModel().domainManager().giveCoordinatesOf(node[i]);
        RealVector x0 = analysisModel().domainManager().giveCoordinatesOf(vertexNode[i]);
        
        double alpha = cns->_fluxOnFace[i]/nhat.dot(x-x0);
        
        centerFlux = centerFlux + alpha*(cellCoor - x0);
    }
    cns->_centerFlux[0] = centerFlux(0);
    cns->_centerFlux[1] = centerFlux(1);
}
//-----------------------------------------------------------------------------
double DarcyFlow_2D_1Phase_Fv_Tri::giveCellFieldValueAt( Cell* targetCell, int fieldNum )
{
    RealVector val, wt;
    std::string fieldTag;
    
    try
    {
        fieldTag = _cellFieldOutput.at(fieldNum);
    }
    catch ( std::exception& e )
    {
        fieldTag = "unassigned";
    }
    std::tie(val,wt) = this->giveFieldOutputAt(targetCell, fieldTag);
    
    return val(0);
}
// ----------------------------------------------------------------------------
std::vector<RealVector>
DarcyFlow_2D_1Phase_Fv_Tri::giveEvaluationPointsFor( Cell* targetCell )
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
DarcyFlow_2D_1Phase_Fv_Tri::giveFieldOutputAt( Cell* targetCell, const std::string& fieldTag )
{
    RealVector fieldVal(1), weight(1);
    weight(0) = this->giveAreaOf(targetCell);
    
    if ( fieldTag == "unassigned" )
        fieldVal(0) = 0.;
    else if ( fieldTag == "h" )
    {
        Dof* dof = analysisModel().domainManager().giveCellDof(_cellDof[0], targetCell);
        fieldVal(0) = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof, converged_value);
    }
    else if ( fieldTag == "p" )
    {
        Dof* dof = analysisModel().domainManager().giveCellDof(_cellDof[0], targetCell);
        double head = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof, current_value);
        
        std::vector<RealVector> coor = this->giveEvaluationPointsFor(targetCell);
        fieldVal(0) = _rhoF*_gAccel*(head + _sgn*coor[0](_vrtIndex));
    }
    else if ( fieldTag == "q_x" )
    {
        auto cns = this->getNumericsStatusAt(targetCell);
        fieldVal(0) = cns->_centerFlux[0];
    }
    else if ( fieldTag == "q_y" )
    {
        auto cns = this->getNumericsStatusAt(targetCell);
        fieldVal(0) = cns->_centerFlux[1];
    }
    else
        throw std::runtime_error("Invalid tag '" + fieldTag + "' supplied in field output request made to numerics '" + _name + "'!");
    
    return std::make_tuple(std::move(fieldVal), std::move(weight));
}
// ----------------------------------------------------------------------------
std::tuple< std::vector<Dof*>
          , std::vector<Dof*>
          , RealVector >
DarcyFlow_2D_1Phase_Fv_Tri::giveStaticCoefficientMatrixAt
    ( Cell*           targetCell
    , int             stage
    , int             subsys
    , const TimeData& time )
{
    std::vector<Dof*> rowDof, colDof;
    RealVector coefVal;
    
    if ( stage == _stage[0] && ( subsys == _subsystem[0] || subsys == UNASSIGNED ) )
    {
        std::vector< std::vector<Node*> > face = this->giveFaceNodesOf(targetCell);
        auto cns = this->getNumericsStatusAt(targetCell);
        
        // Get cell neighbors
        std::vector<Cell*> neighbor = analysisModel().domainManager().giveNeighborsOf(targetCell);
        
        int nCols = 1;
        for ( int i = 0; i < 3; i++)
            if ( neighbor[i] )
                ++nCols;
        
        // Row numbers
        rowDof.assign(nCols, nullptr);
        colDof.assign(nCols, nullptr);
        rowDof[0] = analysisModel().domainManager().giveCellDof(_cellDof[0], targetCell);
        colDof[0] = rowDof[0];

        // Initialize local coefficient matrix
        coefVal.init(nCols);
        
        // Cycle through faces
        int counter = 0;
        for ( int i = 0; i < 3; i++)
        {
            if ( !neighbor[i] && cns->_headIsPrescribedOnFace[i] )
            {
                double length = this->giveLengthOf(face[i]);
                RealVector nhat = this->giveOutwardUnitNormalOf(face[i]);
                
                // Get coordinate of cell center
                std::vector<RealVector> epCoor;
                epCoor = giveEvaluationPointsFor(targetCell);
                RealVector cellCoor({epCoor[0](0), epCoor[0](1), 0.});
                
                double d1 = this->giveDistanceToMidpointOf(face[i], cellCoor);

                // Permeability tensor for cell
                std::vector<Material*> material = this->giveMaterialSetFor(targetCell);
                double k_xx, k_yy, k_xy;
                k_xx = material[0]->giveParameter("Permeability_xx");
                k_yy = material[0]->giveParameter("Permeability_yy");
                k_xy = material[0]->giveParameter("Permeability_xy");

                RealMatrix kmat({{k_xx, k_xy},
                                 {k_xy, k_yy}});
                          
                RealVector kvec;
                kvec = kmat*nhat;
                double k1 = std::sqrt(kvec.dot(kvec));
                coefVal(0) += k1/(_mu*d1)*length;
            }
            else if ( neighbor[i] )
            {
                ++counter;
                rowDof[counter] = rowDof[0];
                colDof[counter] = analysisModel().domainManager().giveCellDof(_cellDof[0], neighbor[i]);
                
                double transCoef = this->giveTransmissibilityCoefficientAt(face[i], targetCell, neighbor[i]);
                coefVal(0) += transCoef;
                coefVal(counter) = -transCoef;
            }
        }
    }
    
    return std::make_tuple(std::move(rowDof), std::move(colDof), std::move(coefVal));
}
// ----------------------------------------------------------------------------
std::tuple< std::vector<Dof*>
          , RealVector >
DarcyFlow_2D_1Phase_Fv_Tri::giveStaticLeftHandSideAt
            ( Cell*           targetCell
            , int             stage
            , int             subsys
            , const TimeData& time )
{
    std::vector<Dof*> rowDof;
    RealVector lhs;
    
    if ( stage == _stage[0] && ( subsys == _subsystem[0] || subsys == UNASSIGNED ) )
    {
        rowDof.assign(1, nullptr);
        lhs.init(1);
        
        // Material set
        std::vector<Material*> material = this->giveMaterialSetFor(targetCell);
        
        // Cell faces and neighbors
        std::vector< std::vector<Node*> > face = this->giveFaceNodesOf(targetCell);
        std::vector<Cell*> neighbor = analysisModel().domainManager().giveNeighborsOf(targetCell);
        
        // Cell area, Cell center coordinates
        std::vector<RealVector> epCoor = this->giveEvaluationPointsFor(targetCell);
        
        // Retrieve hydraulic head
        Dof* dof_h = analysisModel().domainManager().giveCellDof(_cellDof[0], targetCell);
        double h = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof_h, current_value);
        
        if ( subsys == _subsystem[0] || subsys == UNASSIGNED )
        {
            auto cns = this->getNumericsStatusAt(targetCell);
            
            // Cell permeability tensor
            double k_xx, k_yy, k_xy;
            k_xx = material[0]->giveParameter("Permeability_xx");
            k_yy = material[0]->giveParameter("Permeability_yy");
            k_xy = material[0]->giveParameter("Permeability_xy");

            RealMatrix kmat_cell({{k_xx, k_xy},
                                  {k_xy, k_yy}});
            
            // Calculate normal flux at faces
            RealVector outflow(3);
            
            // Cycle through faces
            for ( int i = 0; i < 3; i++ )
            {
                // Flux calculation
                if ( cns->_headIsPrescribedOnFace[i] ) // Hydraulic head is specified at cell face
                {
                    // Length of face
                    double length = this->giveLengthOf(face[i]);
                    
                    // Hydraulic conductivity component normal to face
                    RealVector nhat = giveOutwardUnitNormalOf(face[i]);
                    double d1 = this->giveDistanceToMidpointOf(face[i], epCoor[0]);
                    RealVector kvec;
                    kvec = kmat_cell*nhat;
                    double K1 = std::sqrt(kvec.dot(kvec))*_rhoF*_gAccel/_mu;
                    
                    outflow(i) = length*K1*h/d1;
                }
                else if ( !cns->_fluxIsPrescribedOnFace[i] )
                {
                    // Hydraulic potential at neighbor cell
                    Dof* dof_h2 = analysisModel().domainManager().giveCellDof(_cellDof[0], neighbor[i]);
                    double h2 = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof_h2, current_value);
                    
                    double transCoef = this->giveTransmissibilityCoefficientAt(face[i], targetCell, neighbor[i]);
                    outflow(i) = transCoef*(h - h2);
                }
            }
            
            rowDof[6] = dof_h;
            lhs(6) = outflow(0) + outflow(1) + outflow(2);
        }
    }
    
    return std::make_tuple(std::move(rowDof), std::move(lhs));
}
// ----------------------------------------------------------------------------
std::tuple< std::vector<Dof*>, RealVector >
DarcyFlow_2D_1Phase_Fv_Tri::giveStaticRightHandSideAt
    ( Cell*                    targetCell
    , int                      stage
    , int                      subsys
    , const BoundaryCondition& bndCond
    , const TimeData&          time )
{
    std::vector<Dof*> rowDof;
    RealVector rhs;
    
    if ( stage == _stage[0] && ( subsys == _subsystem[0] || subsys == UNASSIGNED ) )
    {
        if ( bndCond.conditionType() == "Flux" || bndCond.conditionType() == "HydraulicHead" )
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
                            
                            if ( bndCond.conditionType() == "HydraulicHead" ) // Hydraulic head is specified at cell face
                            {
                                cns->_headOnFace[i] = bcVal;
                                
                                // Hydraulic conductivity component normal to face
                                RealVector nhat = giveOutwardUnitNormalOf(face[i]);
                                std::vector<RealVector> epCoor = this->giveEvaluationPointsFor(curDomCell);
                                double d = this->giveDistanceToMidpointOf(face[i], epCoor[0]);
                                
                                std::vector<Material*> material = this->giveMaterialSetFor(curDomCell);
                                double k_xx, k_yy, k_xy;
                                k_xx = material[0]->giveParameter("Permeability_xx");
                                k_yy = material[0]->giveParameter("Permeability_yy");
                                k_xy = material[0]->giveParameter("Permeability_xy");

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
DarcyFlow_2D_1Phase_Fv_Tri::giveStaticRightHandSideAt
    ( Cell*                 targetCell
    , int                   stage
    , int                   subsys
    , const FieldCondition& fldCond
    , const TimeData&       time )
{
    std::vector<Dof*> rowDof;
    RealVector rhs;
    
    if ( stage == _stage[0] && ( subsys == _subsystem[0] || subsys == UNASSIGNED ) )
    {
        // Evaluate body forces arising from field conditions
        std::string name = fldCond.conditionType();
        if ( name == "DistributedSource" )
        {
            auto coor = this->giveEvaluationPointsFor(targetCell);
            double source = fldCond.valueAt(coor[0], time);
            
            double area = this->giveAreaOf(targetCell);
            rhs = {source*area};

            Dof* cellDof = analysisModel().domainManager().giveCellDof(_cellDof[0], targetCell);
            rowDof.assign(1, cellDof);
        }
        else if ( name == "PointSource" )
        {
            auto coor = this->giveEvaluationPointsFor(targetCell);
            double source = fldCond.valueAt(coor[0], time);
            
            rhs = {source};
            
            Dof* cellDof = analysisModel().domainManager().giveCellDof(_cellDof[0], targetCell);
//            rowDof.assign(1, cellDof);
            rowDof = {cellDof};
        }
    }
    
    return std::make_tuple(std::move(rowDof), std::move(rhs));
}
// ----------------------------------------------------------------------------
void DarcyFlow_2D_1Phase_Fv_Tri::imposeConstraintAt
    ( Cell*                    targetCell
    , int                      stg
    , const BoundaryCondition& bndCond
    , const TimeData&          time )
{   
    // Constraints pertaining to cell DOFS
    if ( bndCond.conditionType() == "CellConstraint" )
    {
        if ( bndCond.targetDof() == _cellDof[0] )
        {
            Dof* targetDof = analysisModel().domainManager().giveCellDof(_cellDof[0], targetCell);

            std::vector<RealVector> coor = this->giveEvaluationPointsFor(targetCell);
            double bcVal = bndCond.valueAt(coor[0], time);
            analysisModel().dofManager().setConstraintValueAt(targetDof, bcVal);
            
            auto cns = this->getNumericsStatusAt(targetCell);
            cns->_hasHeadConstraint = true;
        }
        else
            throw std::runtime_error("ERROR! Constraint imposed at invalid DOF.\nSource: " + _name);
    }
    
    if ( bndCond.conditionType() == "HydraulicHead" || bndCond.conditionType() == "Flux" )
    {
        std::vector<Node*> bndCellNode = analysisModel().domainManager().giveNodesOf(targetCell);
        
        if ( bndCellNode.size() != 2 )
            throw std::runtime_error("Numerics '" + _name + "' only accepts boundary cells with exactly two nodes!");
        
        std::vector<Cell*> domCell = analysisModel().domainManager().giveDomainCellsAssociatedWith(targetCell);
        for ( Cell* curDomCell : domCell )
        {
            Numerics* domCellNumerics = analysisModel().domainManager().giveNumericsFor(curDomCell);
            if ( this == domCellNumerics )
            {
                std::vector< std::vector<Node*> > face = this->giveFaceNodesOf(curDomCell);
                auto cns = this->getNumericsStatusAt(curDomCell);

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
void DarcyFlow_2D_1Phase_Fv_Tri::initializeMaterialsAt( Cell* targetCell )
{
    std::vector<Material*> material = this->giveMaterialSetFor(targetCell);
    auto cns = this->getNumericsStatusAt(targetCell);
    cns->_materialStatus = material[0]->createMaterialStatus();
}
// ----------------------------------------------------------------------------
void DarcyFlow_2D_1Phase_Fv_Tri::initializeNumericsAt( Cell* targetCell )
{
    targetCell->numericsStatus = new NumericsStatus_DarcyFlow_2D_1Phase_Fv_Tri();
}
// ----------------------------------------------------------------------------
void DarcyFlow_2D_1Phase_Fv_Tri::readAdditionalDataFrom( FILE* fp )
{
    std::string key;
    verifyKeyword(fp, key = "SpecificStorage", _name);
    _S = getRealInputFrom(fp, "Failed to read fluid specific storage from input file!", _name);
    
    verifyKeyword(fp, key = "FluidDensity", _name);
    _rhoF = getRealInputFrom(fp, "Failed to read fluid density from input file!", _name);
    
    verifyKeyword(fp, key = "FluidViscosity", _name);
    _mu = getRealInputFrom(fp, "Failed to read fluid viscosity from input file!", _name);
    
    verifyKeyword(fp, key = "GravitationalAcceleration", _name);
    key = getStringInputFrom(fp, "Failed to read direction of gravitational acceleration from input file!", _name);
    if ( key == "-x" )
    {
        _sgn = -1;
        _vrtIndex = 0;
    }
    else if ( key == "+x" )
    {
        _sgn = 1;
        _vrtIndex = 0;
    }
    else if ( key == "-y" )
    {
        _sgn = -1;
        _vrtIndex = 1;
    }
    else if ( key == "-z" )
    {
        _sgn = 1;
        _vrtIndex = 2;
    }
    else if ( key == "+z" )
    {
        _sgn = 1;
        _vrtIndex = 2;
    }
    else
        throw std::runtime_error("Invalid option '" + key + "' specified for direction of gravity acceleration in input file!\nSource: " + _name);
    
    _gAccel = getRealInputFrom(fp, "Failed to read value of gravitational acceleration from input file!", _name);
}
// ----------------------------------------------------------------------------
void DarcyFlow_2D_1Phase_Fv_Tri::removeConstraintsOn( Cell* targetCell )
{
    auto cns = this->getNumericsStatusAt(targetCell);
    cns->_headIsPrescribedOnFace[0] = false;
    cns->_headIsPrescribedOnFace[1] = false;
    cns->_headIsPrescribedOnFace[2] = false;
    cns->_fluxIsPrescribedOnFace[0] = false;
    cns->_fluxIsPrescribedOnFace[1] = false;
    cns->_fluxIsPrescribedOnFace[2] = false;
    cns->_hasHeadConstraint = false;
}
// ----------------------------------------------------------------------------
void DarcyFlow_2D_1Phase_Fv_Tri::setDofStagesAt( Cell* targetCell )
{
    Dof* dof = analysisModel().domainManager().giveCellDof(1, targetCell);
    analysisModel().dofManager().setStageFor(dof, _stage[0]);
}

// Private methods
// ---------------------------------------------------------------------------
NumericsStatus_DarcyFlow_2D_1Phase_Fv_Tri*
DarcyFlow_2D_1Phase_Fv_Tri::getNumericsStatusAt( Cell* targetCell )
{
    auto cns = dynamic_cast<NumericsStatus_DarcyFlow_2D_1Phase_Fv_Tri*>(targetCell->numericsStatus);
    if ( !cns )
        throw std::runtime_error("Error: Unable to retrieve numerics status at cell!\nSource: " + _name);
    
    return cns;
}
// ---------------------------------------------------------------------------
double DarcyFlow_2D_1Phase_Fv_Tri::giveAreaOf( Cell* targetCell )
{
    double area;
    
    RealVector coor_0, coor_1, coor_2;
    std::vector<Node*> node = analysisModel().domainManager().giveNodesOf(targetCell);
    
    coor_0 = analysisModel().domainManager().giveCoordinatesOf(node[0]);
    coor_1 = analysisModel().domainManager().giveCoordinatesOf(node[1]);
    coor_2 = analysisModel().domainManager().giveCoordinatesOf(node[2]);
    
    area = 0.5*(coor_0(0)*(coor_1(1) - coor_2(1)) + 
                coor_1(0)*(coor_2(1) - coor_0(1)) +
                coor_2(0)*(coor_0(1) - coor_1(1)));
    
    if ( area <= 0 )
        throw std::runtime_error("Calculation of negative cell area detected!\nSource: " + _name);
    
    return area;
}
// ----------------------------------------------------------------------------
double DarcyFlow_2D_1Phase_Fv_Tri::giveDistanceToMidpointOf( const std::vector<Node*>& face
                                                           , const RealVector&   coor )
{
    RealVector coor0, coor1;
    coor0 = analysisModel().domainManager().giveCoordinatesOf(face[0]);
    coor1 = analysisModel().domainManager().giveCoordinatesOf(face[1]);
    
    RealVector dx;
    dx = 0.5*(coor0 + coor1) - coor;
    
    return std::sqrt(dx.dot(dx));
}
// ----------------------------------------------------------------------------
std::vector<std::vector<Node*> >
DarcyFlow_2D_1Phase_Fv_Tri::giveFaceNodesOf( Cell* targetCell )
{
    std::vector<Node*> node = analysisModel().domainManager().giveNodesOf(targetCell);
    std::vector<std::vector<Node*> > face(3, std::vector<Node*>(2, nullptr));
    
    face[0][0] = node[0];
    face[0][1] = node[1];
    face[1][0] = node[1];
    face[1][1] = node[2];
    face[2][0] = node[2];
    face[2][1] = node[0];
    
    return face;
}
// ----------------------------------------------------------------------------
double DarcyFlow_2D_1Phase_Fv_Tri::giveLengthOf( const std::vector<Node*>& face )
{
    RealVector coor0, coor1;
    coor0 = analysisModel().domainManager().giveCoordinatesOf(face[0]);
    coor1 = analysisModel().domainManager().giveCoordinatesOf(face[1]);
    
    RealVector dx = coor1 - coor0;
    
    return std::sqrt(dx.dot(dx));
}
// ----------------------------------------------------------------------------
RealVector DarcyFlow_2D_1Phase_Fv_Tri::giveOutwardUnitNormalOf( const std::vector<Node*>& face )
{
    RealVector coor0, coor1;
    coor0 = analysisModel().domainManager().giveCoordinatesOf(face[0]);
    coor1 = analysisModel().domainManager().giveCoordinatesOf(face[1]);
    
    RealVector dx;
    dx = coor1 - coor0;
    double dist = std::sqrt(dx.dot(dx));
    
    RealVector normVec({dx(1)/dist, -dx(0)/dist});
    
    return normVec;
}
// ----------------------------------------------------------------------------
double DarcyFlow_2D_1Phase_Fv_Tri::giveTransmissibilityCoefficientAt
    ( const std::vector<Node*>& face
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
    k_xx = material[0]->giveParameter("Permeability_xx");
    k_yy = material[0]->giveParameter("Permeability_yy");
    k_xy = material[0]->giveParameter("Permeability_xy");

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
    k_xx = material[0]->giveParameter("Permeability_xx");
    k_yy = material[0]->giveParameter("Permeability_yy");
    k_xy = material[0]->giveParameter("Permeability_xy");

    kmat = {{k_xx, k_xy},
            {k_xy, k_yy}};
                     
    kvec = kmat*nhat;
    double k2 = std::sqrt(kvec.dot(kvec));
    
    return length/(_mu*d1/k1 + _mu*d2/k2);
}
// ----------------------------------------------------------------------------
double DarcyFlow_2D_1Phase_Fv_Tri::giveVerticalCoordinateAt( Cell* targetCell )
{
    // Compute vertical coordinate at cell center
    std::vector<Node*> node = analysisModel().domainManager().giveNodesOf(targetCell);
    RealVector p0, p1, p2;
    p0 = analysisModel().domainManager().giveCoordinatesOf(node[0]);
    p1 = analysisModel().domainManager().giveCoordinatesOf(node[1]);
    p2 = analysisModel().domainManager().giveCoordinatesOf(node[2]);
    
    return _sgn*(p0(_vrtIndex) + p1(_vrtIndex) + p2(_vrtIndex))/3.;
}