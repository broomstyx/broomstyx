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

#include "PhaseFieldFracture_Fe_Tri3.hpp"
#include <cmath>
#include <string>
#include <vector>

#include "Core/AnalysisModel.hpp"
#include "Core/DomainManager.hpp"
#include "Core/ObjectFactory.hpp"
#include "Materials/Material.hpp"
#include "Util/linearAlgebra.hpp"
#include "Util/readOperations.hpp"

using namespace broomstyx;

registerBroomstyxObject( Numerics, PhaseFieldFracture_Fe_Tri3 )

#define CRACKLENGTH_CUTOFF 0.99

// Numerics status
NumericsStatus_PhaseFieldFracture_Fe_Tri3::NumericsStatus_PhaseFieldFracture_Fe_Tri3()
    : _area( 0. )
	, _phi( 0. )
    , _phiOld( 0. )
	, _strain( RealVector( 4 ) )
    , _stress( RealVector( 4 ) )
    , _gradU( RealMatrix( 2,2 ) )
    , _surfEgy( 0. )
    , _bulkEgy( 0. )
    , _Gc( 0. )
    , _histFld( 0. )
    , _materialStatus { nullptr, nullptr, nullptr }
{}

NumericsStatus_PhaseFieldFracture_Fe_Tri3::~NumericsStatus_PhaseFieldFracture_Fe_Tri3() = default;

// Constructor
PhaseFieldFracture_Fe_Tri3::PhaseFieldFracture_Fe_Tri3()
{
    _dim = 2;
    _dofPerNode = 3;
    _nNodes = 3;
    _nMaterials = 3;
    _nStages = 1;
    _nSubsystems = 2;
    _name = "PhaseFieldFracture_Fe_Tri3";
    
    // Initialize private variables
    _lc = 0.;
    _phiIrrev = 0.;

    // Lone Gauss point coordinate and weight
    _gpNatCoor = { 1./3., 1./3. };
    _wt = 0.5;

    _analysisMode = Unset;

    // Pre-calculate shape functions and derivatives at Gauss point
    _basisFunctionValues = _basisFunction.giveBasisFunctionsAt( _gpNatCoor );
    std::vector<RealVector> dpsiNat = _basisFunction.giveBasisFunctionDerivativesAt( _gpNatCoor );

    _basisFunctionDerivatives = { { dpsiNat[ 0 ]( 0 ), dpsiNat[ 0 ]( 1 ), dpsiNat[ 0 ]( 2 ) },
                                  { dpsiNat[ 1 ]( 0 ), dpsiNat[ 1 ]( 1 ), dpsiNat[ 1 ]( 2 ) } };

    // Pre-define mass matrix
    _massMatrix = { { 1./6.,  1./12., 1./12. },
                    { 1./12., 1./6.,  1./12. },
                    { 1./12., 1./12., 1./6. } };
}

// Destructor
PhaseFieldFracture_Fe_Tri3::~PhaseFieldFracture_Fe_Tri3() = default;

// Public methods
// ----------------------------------------------------------------------------
void PhaseFieldFracture_Fe_Tri3::deleteNumericsAt( Cell* targetCell )
{
    auto cns = this->getNumericsStatusAt( targetCell );
    auto material = this->giveMaterialSetFor( targetCell );
    
    material[ 0 ]->destroy( cns->_materialStatus[ 0 ] );
    material[ 1 ]->destroy( cns->_materialStatus[ 1 ] );
    material[ 2 ]->destroy( cns->_materialStatus[ 2 ] );
}
// ----------------------------------------------------------------------------
void PhaseFieldFracture_Fe_Tri3::finalizeDataAt( Cell* targetCell, const TimeData& time )
{
    auto cns = this->getNumericsStatusAt( targetCell );
    
    // Retrieve nodal DOFs local to element
    std::vector<Dof*> dof = this->giveNodalDofsAt( targetCell );
    
    // Retrieve DOF values
    RealVector uVec, phiVec;
    std::tie( uVec,phiVec ) = giveLocalVariablesAt( dof, converged_value );
    
    // Compute displacement gradient
    RealMatrix uMat( { { uVec(0), uVec(1) },
                       { uVec(2), uVec(3) },
                       { uVec(4), uVec(5) } } );
    
    cns->_gradU = cns->_dPsi * uMat;
    
    // Construct strain vector
    RealMatrix bmatU = this->giveBmatUAt( targetCell );
    cns->_strain = bmatU * uVec;

    // Construct phase-field and its gradient
    cns->_phi = _basisFunctionValues.dot( phiVec );
    RealVector dphi = cns->_dPsi * phiVec;

    // Update old value of phase-field
    if ( cns->_phi > cns->_phiOld )
        cns->_phiOld = cns->_phi;
    
    // Retrieve material set for element
    std::vector<Material*> material = this->giveMaterialSetFor( targetCell );
    
    // Get constitutive force and elastic bulk energy
    RealVector conState( { cns->_strain( 0 ),
                           cns->_strain( 1 ),
                           cns->_strain( 2 ),
                           cns->_strain( 3 ),
                           cns->_phi } );

    cns->_stress = material[ 1 ]->giveForceFrom( conState, cns->_materialStatus[ 1 ], "Mechanics" );
    cns->_bulkEgy = material[ 1 ]->givePotentialFrom( conState, cns->_materialStatus[ 1 ] );

    RealVector conForce = material[ 1 ]->giveForceFrom( conState, cns->_materialStatus[ 1 ], "PhaseField" );
    double d_degFcn = conForce( 0 );
    double drvForce = conForce( 1 );

    RealVector conState2( { cns->_phi,
                            -_lc * d_degFcn * drvForce,
                            time.giveCurrentTime(),
                            time.giveTargetTime() } );

    // Finalize history field
    if ( drvForce > cns->_histFld )
        cns->_histFld = drvForce;

    RealMatrix Gc = material[ 2 ]->giveModulusFrom( conState2, cns->_materialStatus[ 2 ] );
    cns->_Gc = Gc( 0,0 );
    cns->_surfEgy = cns->_Gc / 2.0 * ( _lc * dphi.dot( dphi ) + phiVec.dot( _massMatrix * phiVec ) / _lc );
}
// ----------------------------------------------------------------------------
double PhaseFieldFracture_Fe_Tri3::giveCellFieldValueAt( Cell* targetCell, int fieldNum )
{
    RealVector val, wt;
    std::string fieldTag;
    
    try
    {
        fieldTag = _cellFieldOutput.at( fieldNum );
    }
    catch ( std::exception& e )
    {
        fieldTag = "unassigned";
    }
    
    std::tie( val,wt ) = this->giveFieldOutputAt( targetCell, fieldTag );
    
    return val( 0 );
}
// ----------------------------------------------------------------------------
RealVector PhaseFieldFracture_Fe_Tri3::giveCellNodeFieldValuesAt( Cell* targetCell, int fieldNum )
{
    double val = this->giveCellFieldValueAt( targetCell, fieldNum );
    
    RealVector cellNodeValues( 3 );
    cellNodeValues( 0 ) = val;
    cellNodeValues( 1 ) = val;
    cellNodeValues( 2 ) = val;
    
    return cellNodeValues;
}
// ----------------------------------------------------------------------------
std::vector<RealVector> PhaseFieldFracture_Fe_Tri3::giveEvaluationPointsFor( Cell *targetCell )
{
    std::vector<Node*> node = analysisModel().domainManager().giveNodesOf( targetCell );
    
    std::vector<RealVector> coor( 1, RealVector() );
    coor[ 0 ].init( 3 );
    
    for ( int i = 0; i < 3; i++ )
    {
        RealVector nodeCoor = analysisModel().domainManager().giveCoordinatesOf( node[ i ] );
        coor[ 0 ]( 0 ) += nodeCoor( 0 ) / 3.0;
        coor[ 0 ]( 1 ) += nodeCoor( 1 ) / 3.0;
        coor[ 0 ]( 2 ) += nodeCoor( 2 ) / 3.0;
    }
    
    return coor;
}
// ----------------------------------------------------------------------------
std::tuple< RealVector,RealVector > 
PhaseFieldFracture_Fe_Tri3::giveFieldOutputAt( Cell* targetCell, const std::string& fieldTag )
{
    RealVector fieldVal(1), weight(1);
    
    auto cns = this->getNumericsStatusAt( targetCell );
    std::vector<Material*> material = this->giveMaterialSetFor( targetCell );

    weight( 0 ) = cns->_area;
    
    if ( fieldTag == "unassigned" )
        fieldVal( 0 ) = 0.;
    else if ( fieldTag == "s_xx" )
        fieldVal( 0 ) = cns->_stress( 0 );
    else if ( fieldTag == "s_yy" )
        fieldVal( 0 ) = cns->_stress( 1 );
    else if ( fieldTag == "s_zz" )
    {
        if ( _analysisMode == PlaneStrain )
            fieldVal( 0 ) = cns->_stress( 2 );
        else
            fieldVal( 0 ) = 0.;
    }
    else if ( fieldTag == "s_xy" )
        fieldVal( 0 ) = cns->_stress( 3 );
    else if ( fieldTag == "ux_x" || fieldTag == "e_xx" )
        fieldVal( 0 ) = cns->_gradU( 0,0 );
    else if ( fieldTag == "uy_x" )
        fieldVal( 0 ) = cns->_gradU( 0,1 );
    else if ( fieldTag == "ux_y" )
        fieldVal( 0 ) = cns->_gradU( 1,0 );
    else if ( fieldTag == "uy_y" || fieldTag == "e_yy" )
        fieldVal( 0 ) = cns->_gradU( 1,1 );
    else if ( fieldTag == "g_xy" )
        fieldVal( 0 ) = cns->_strain( 3 );
    else if ( fieldTag == "pf" )
        fieldVal( 0 ) = cns->_phi;
    else if ( fieldTag == "ene_s" )
        fieldVal( 0 ) = cns->_surfEgy;
    else if ( fieldTag == "ene_b" )
        fieldVal( 0 ) = cns->_bulkEgy;
    else if ( fieldTag == "cr_len" )
    {
        if ( cns->_phi < CRACKLENGTH_CUTOFF )
            fieldVal( 0 ) = cns->_surfEgy / cns->_Gc;
        else
            fieldVal( 0 ) = 0.;
    }
    else if ( fieldTag == "Gc" )
        fieldVal( 0 ) = cns->_Gc;
        
    else if ( fieldTag == "eigVec1_1" || fieldTag == "eigVec1_2" || fieldTag == "eigVec2_1" || fieldTag == "eigVec2_2" )
        fieldVal( 0 ) = material[ 1 ]->giveMaterialVariable( fieldTag, cns->_materialStatus[ 1 ] );
    else
        throw std::runtime_error( "Invalid tag '" + fieldTag + "' supplied in field output request made to numerics '" + _name + "'!" );
    
    return std::make_tuple( std::move( fieldVal ), std::move( weight ) );
}
// ----------------------------------------------------------------------------
std::tuple< std::vector<Dof*>
          , std::vector<Dof*>
          , RealVector >
PhaseFieldFracture_Fe_Tri3::giveStaticCoefficientMatrixAt( Cell*           targetCell
                                                         , int             stage
                                                         , int             subsys
                                                         , const TimeData& time )
{
    std::vector<Dof*> rowDof, colDof;
    RealVector coefVal;
    
    if ( stage == _stage[ 0 ] )
    {
        // Retrieve nodal DOFs local to element
        std::vector<Dof*> dof = giveNodalDofsAt(targetCell);
        
        int vecLength = 45; // Default output lengths for unassigned subsystems
        if ( subsys == _subsystem[ 0 ] )
            vecLength = 36;
        else if ( subsys == _subsystem[ 1 ] )
            vecLength = 9;
        
        rowDof.assign( vecLength, nullptr );
        colDof.assign( vecLength, nullptr );
        coefVal.init( vecLength );
        
        // Retrieve numerics status
        auto cns = this->getNumericsStatusAt( targetCell );
        
        // Retrieve material set for element
        std::vector<Material*> material = this->giveMaterialSetFor( targetCell );
        
        // Assemble constitutive state
        RealVector conState( { cns->_strain( 0 ),
                               cns->_strain( 1 ),
                               cns->_strain( 2 ),
                               cns->_strain( 3 ),
                               cns->_phi } );
        
        int counter = 0;
        if ( subsys == _subsystem[ 0 ] || subsys == UNASSIGNED )
        {
            // Get tangent modulus (linear momentum eq.)
            RealMatrix cmat = material[ 1 ]->giveModulusFrom( conState, cns->_materialStatus[ 1 ], "Mechanics" );
            
            // Calculate stiffness KmatUU
            RealMatrix bmatU = this->giveBmatUAt( targetCell );
            RealMatrix kmatUU = cns->_area * trp( bmatU ) * cmat * bmatU;
            
            for ( int i = 0; i < 6; i++ )
                for ( int j = 0; j < 6; j++ )
                {
                    rowDof[ counter ] = dof[ i ];
                    colDof[ counter ] = dof[ j ];
                    coefVal( counter ) = kmatUU( i,j );
                    ++counter;
                }
        }
        if ( subsys == _subsystem[ 1 ] || subsys == UNASSIGNED )
        {
            // Compute g''(phi)*elasticEnergy
            RealMatrix conMod = material[ 1 ]->giveModulusFrom( conState, cns->_materialStatus[ 1 ], "PhaseField" );
            double dd_degFcn = conMod( 0,0 );
            double drvForce = conMod( 1,1 );

            if ( cns->_phi > _phiIrrev && drvForce < cns->_histFld )
                drvForce = cns->_histFld;
            
            // Calculate KmatPhiPhi
            RealMatrix kmatPhiPhi = cns->_area * cns->_Gc * _lc * trp( cns->_dPsi ) * cns->_dPsi
                                  + cns->_area * ( cns->_Gc / _lc + dd_degFcn * drvForce ) * _massMatrix;

            for ( int i = 0; i < 3; i++ )
                for ( int j = 0; j < 3; j++ )
                {
                    rowDof[ counter ] = dof[ 6+i ];
                    colDof[ counter ] = dof[ 6+j ];
                    coefVal( counter ) = kmatPhiPhi( i,j );
                    ++counter;
                }
        }
    }
    
    return std::make_tuple(std::move(rowDof), std::move(colDof), std::move(coefVal));
}
// ----------------------------------------------------------------------------
std::tuple< std::vector<Dof*>
          , RealVector >
PhaseFieldFracture_Fe_Tri3::giveStaticLeftHandSideAt( Cell*           targetCell
                                                    , int             stage
                                                    , int             subsys
                                                    , const TimeData& time )
{
    std::vector<Dof*> rowDof;
    RealVector lhs;
    
    if ( stage == _stage[ 0 ] )
    {
        // Retrieve nodal DOFs local to element
        std::vector<Dof*> dof = this->giveNodalDofsAt( targetCell );

        int vecLength = 15; // Default output lengths for unassigned subsystems
        if ( subsys == _subsystem[ 0 ] )
            vecLength = 6;
        else if ( subsys == _subsystem[ 1 ] )
            vecLength = 9;
        
        rowDof.assign( vecLength, nullptr );
        lhs.init( vecLength );

        // Retrieve current value of local displacements and phase-field
        RealVector uVec, phiVec;
        std::tie( uVec,phiVec ) = giveLocalVariablesAt( dof, current_value );
        
        // Retrieve numerics status
        auto cns = this->getNumericsStatusAt( targetCell );
        
        // Calculate phase-field
        cns->_phi = _basisFunctionValues.dot( phiVec );
        
        // Compute local strains
        RealMatrix bmatU = this->giveBmatUAt( targetCell );
        cns->_strain = bmatU * uVec;
        
        // Assemble constitutive state vector
        RealVector conState( { cns->_strain( 0 ),
                               cns->_strain( 1 ),
                               cns->_strain( 2 ),
                               cns->_strain( 3 ),
                               cns->_phi } );
                             
        // Retrieve material set for element
        std::vector<Material*> material = this->giveMaterialSetFor( targetCell );
        
        material[ 1 ]->updateStatusFrom( conState, cns->_materialStatus[ 1 ] );

        int offset = 0;
        if ( subsys == _subsystem[ 0 ] || subsys == UNASSIGNED )
        {
            rowDof[ 0 ] = dof[ 0 ];
            rowDof[ 1 ] = dof[ 1 ];
            rowDof[ 2 ] = dof[ 2 ];
            rowDof[ 3 ] = dof[ 3 ];
            rowDof[ 4 ] = dof[ 4 ];
            rowDof[ 5 ] = dof[ 5 ];
            offset = 6;
            
            // Update material state and compute damage reduced stress
            RealVector stress = material[ 1 ]->giveForceFrom( conState, cns->_materialStatus[ 1 ], "Mechanics" );
            
            // Calculate FmatU
            RealVector fmatU = cns->_area * trp( bmatU ) * stress;
            
            lhs( 0 ) = fmatU( 0 );
            lhs( 1 ) = fmatU( 1 );
            lhs( 2 ) = fmatU( 2 );
            lhs( 3 ) = fmatU( 3 );
            lhs( 4 ) = fmatU( 4 );
            lhs( 5 ) = fmatU( 5 );
        }
        
        if ( subsys == _subsystem[ 1 ] || subsys == UNASSIGNED )
        {
            rowDof[ offset + 0 ] = dof[ 6 ];
            rowDof[ offset + 1 ] = dof[ 7 ];
            rowDof[ offset + 2 ] = dof[ 8 ];
            rowDof[ offset + 3 ] = dof[ 6 ];
            rowDof[ offset + 4 ] = dof[ 7 ];
            rowDof[ offset + 5 ] = dof[ 8 ];
            rowDof[ offset + 6 ] = dof[ 6 ];
            rowDof[ offset + 7 ] = dof[ 7 ];
            rowDof[ offset + 8 ] = dof[ 8 ];

            RealVector dphi = cns->_dPsi * phiVec;
            
            // Update material state and calculate g'(phi)*elasticEnergy
            RealVector conForce = material[ 1 ]->giveForceFrom( conState, cns->_materialStatus[ 1 ], "PhaseField" );
            double d_degFcn = conForce( 0 );
            double drvForce = conForce( 1 );

            RealVector conState2( { cns->_phi,
                                    -_lc * d_degFcn * drvForce,
                                    time.giveCurrentTime(),
                                    time.giveTargetTime() } );

            material[ 2 ]->updateStatusFrom( conState2, cns->_materialStatus[ 2 ] );
            RealMatrix Gc = material[ 2 ]->giveModulusFrom( conState2, cns->_materialStatus[ 2 ] );
            cns->_Gc = Gc( 0,0 );

            // Enforce irreversibility
            if ( cns->_phi > _phiIrrev && drvForce < cns->_histFld )
                drvForce = cns->_histFld;

            // Calculate FmatPhi
            RealVector fmatPhi1, fmatPhi2, fmatPhi3;
            fmatPhi1 = cns->_area * ( cns->_Gc / _lc * cns->_phi ) * _basisFunctionValues;
            fmatPhi2 = cns->_area * cns->_Gc * _lc * trp( cns->_dPsi ) * dphi;
            fmatPhi3 = cns->_area * d_degFcn * drvForce * _basisFunctionValues;

            lhs( offset + 0 ) = fmatPhi1( 0 );
            lhs( offset + 1 ) = fmatPhi1( 1 );
            lhs( offset + 2 ) = fmatPhi1( 2 );
            lhs( offset + 3 ) = fmatPhi2( 0 );
            lhs( offset + 4 ) = fmatPhi2( 1 );
            lhs( offset + 5 ) = fmatPhi2( 2 );
            lhs( offset + 6 ) = fmatPhi3( 0 );
            lhs( offset + 7 ) = fmatPhi3( 1 );
            lhs( offset + 8 ) = fmatPhi3( 2 );
        }
    }
    
    return std::make_tuple( std::move( rowDof ), std::move( lhs ) );
}
// ----------------------------------------------------------------------------
std::tuple< std::vector<Dof*>
          , RealVector >
PhaseFieldFracture_Fe_Tri3::giveStaticRightHandSideAt( Cell*                    targetCell
                                                     , int                      stage
                                                     , int                      subsys
                                                     , const BoundaryCondition& bndCond
                                                     , const TimeData&          time )
{   
    std::vector<Dof*> rowDof;
    RealVector rhs;
    
    if ( stage == _stage[ 0 ] && ( subsys == _subsystem[ 0 ] || subsys == UNASSIGNED ) )
    {
        if ( bndCond.conditionType() == "Traction" )
        {
            // Retrieve nodes of boundary element
            std::vector<Node*> node = analysisModel().domainManager().giveNodesOf( targetCell );
            
            // ----------------------------------------------------------------
            // Note: Boundary element must be a 2-node line
            // ----------------------------------------------------------------
            
            if ( (int)node.size() != 2 )
                throw std::runtime_error( "Error: Traction boundary condition for '" + _name + "' requires 2-node boundary elements!" );
            
            RealVector coor0, coor1, dx;
            coor0 = analysisModel().domainManager().giveCoordinatesOf( node[ 0 ] );
            coor1 = analysisModel().domainManager().giveCoordinatesOf( node[ 1 ] );
            dx = coor1 - coor0;
            double length = std::sqrt( dx.dot( dx ) );

            // Determine proper value of boundary condition
            RealVector midpt = 0.5 * ( coor0 + coor1 );
            double bcVal = bndCond.valueAt( midpt, time );

            // Construct local RHS vector (only for relevant DOFs)
            rhs.init( 2 );
            rhs( 0 ) = 0.5 * length * bcVal;
            rhs( 1 ) = 0.5 * length * bcVal;

            // Construct global address vector
            rowDof.assign( 2, nullptr );

            int dofNum = analysisModel().dofManager().giveIndexForNodalDof( bndCond.targetDof() );
            rowDof[ 0 ] = analysisModel().domainManager().giveNodalDof( dofNum, node[ 0 ] );
            rowDof[ 1 ] = analysisModel().domainManager().giveNodalDof( dofNum, node[ 1 ] );
        }
        else if ( bndCond.conditionType() == "ConcentratedForce" ) // Boundary conditions for 1-node point
        {
            // Retrieve nodes of boundary element
            std::vector<Node*> node = analysisModel().domainManager().giveNodesOf( targetCell );
            
            // ----------------------------------------------------------------
            // Note: Boundary element must be a 1-node point
            // ----------------------------------------------------------------
            
            if ( (int)node.size() != 1 )
                throw std::runtime_error( "Error: Concentrated force boundary condition for '" + _name + "' requires 1-node boundary elements!" );
            
            // Determine proper value of boundary condition
            RealVector coor = analysisModel().domainManager().giveCoordinatesOf( node[ 0 ] );
            double bcVal = bndCond.valueAt(coor, time);

            // Construct local RHS vector (only for relevant DOFs)
            rhs.init( 1 );
            rhs( 0 ) = bcVal;

            // Construct global address vector
            rowDof.assign( 1, nullptr );

            int dofNum = analysisModel().dofManager().giveIndexForNodalDof( bndCond.targetDof() );
            rowDof[ 0 ] = analysisModel().domainManager().giveNodalDof( dofNum, node[ 0 ] );
        }
    }
    
    return std::make_tuple( std::move( rowDof ), std::move( rhs ) );
}
// ----------------------------------------------------------------------------
std::tuple< std::vector<Dof*>
          , RealVector >
PhaseFieldFracture_Fe_Tri3::giveStaticRightHandSideAt( Cell*                 targetCell
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
        if ( cndType == "Acceleration_X" || cndType == "Acceleration_Y" )
        {
            auto cns = this->getNumericsStatusAt( targetCell );
            
            // Get element density
            std::vector<Material*> material = this->giveMaterialSetFor( targetCell );
            double rho = material[ 0 ]->giveMaterialVariable( "Density", cns->_materialStatus[ 0 ] );

            // Retrieve nodal DOFs for element
            std::vector<Dof*> dof = this->giveNodalDofsAt( targetCell );
            
            // Note that field condition is assumed to be piecewise linear over each cell
            std::vector<RealVector> coor = this->giveEvaluationPointsFor( targetCell );
            double accel = fldCond.valueAt( coor[ 0 ], time );
            double force = accel * rho * cns->_area / 3.0;
            
            rhs = { force, force, force };
            
            rowDof.assign( 3, nullptr );
            if ( cndType == "Acceleration_X" )
            {
                rowDof[ 0 ] = dof[ 0 ];
                rowDof[ 1 ] = dof[ 2 ];
                rowDof[ 2 ] = dof[ 4 ];
            }
            else
            {
                rowDof[ 0 ] = dof[ 1 ];
                rowDof[ 1 ] = dof[ 3 ];
                rowDof[ 2 ] = dof[ 5 ];
            }
        }        
    }
    
    return std::make_tuple(std::move(rowDof), std::move(rhs));
}
// ----------------------------------------------------------------------------
void PhaseFieldFracture_Fe_Tri3::imposeConstraintAt( Cell*                    targetCell
                                                   , int                      stage
                                                   , const BoundaryCondition& bndCond
                                                   , const TimeData&          time )
{
    // Only for essential BCs on nodal DOFs
    if ( bndCond.conditionType() == "NodalConstraint" )
    {
        // Retrieve nodes of boundary element
        std::vector<Node*> node = analysisModel().domainManager().giveNodesOf( targetCell );
        int dofNum = analysisModel().dofManager().giveIndexForNodalDof( bndCond.targetDof() );

        for ( auto& curNode : node)
        {
            Dof* targetDof = analysisModel().domainManager().giveNodalDof( dofNum, curNode );
            RealVector coor = analysisModel().domainManager().giveCoordinatesOf( curNode );
            double bcVal = bndCond.valueAt( coor, time );
            DofManager::setConstraintValueAt( targetDof, bcVal );
        }
    }
}
// ----------------------------------------------------------------------------
void PhaseFieldFracture_Fe_Tri3::initializeMaterialsAt( Cell* targetCell )
{
    auto cns = this->getNumericsStatusAt( targetCell );
    std::vector<Material*> material = this->giveMaterialSetFor( targetCell );
    
    cns->_materialStatus[ 0 ] = material[ 0 ]->createMaterialStatus();
    cns->_materialStatus[ 1 ] = material[ 1 ]->createMaterialStatus();
    cns->_materialStatus[ 2 ] = material[ 2 ]->createMaterialStatus();
}
// ----------------------------------------------------------------------------
void PhaseFieldFracture_Fe_Tri3::initializeNumericsAt( Cell* targetCell )
{
    targetCell->numericsStatus = new NumericsStatus_PhaseFieldFracture_Fe_Tri3();
    auto cns = this->getNumericsStatusAt( targetCell );
    
    // Pre-calculate det(J) and inv(J);
    RealMatrix Jmat = this->giveJacobianMatrixAt( targetCell );
    double Jdet = Jmat( 0,0 ) * Jmat( 1,1 ) - Jmat( 1,0 ) * Jmat( 0,1 );
    cns->_area = _wt * Jdet;

    // Sanity check
    if ( Jdet <= 0 )
        throw std::runtime_error( "Calculation of negative area detected!\nSource: " + _name );
    
    // Calculate shape function derivatives in the actual space
    RealMatrix JmatInv = inv( Jmat );
    cns->_dPsi = JmatInv * _basisFunctionDerivatives;
}
// ----------------------------------------------------------------------------
void PhaseFieldFracture_Fe_Tri3::printPostIterationMessage( int stage )
{
    if ( stage == _stage[ 0 ] )
    {
        // Calculate crack length
        double totalCrackLength = 0;
        
        int nCells = analysisModel().domainManager().giveNumberOfDomainCells();
#ifdef _OPENMP
#pragma omp parallel for reduction (+:totalCrackLength)
#endif
        for ( int i = 0; i < nCells; i++ )
        {
            Cell* curCell = analysisModel().domainManager().giveDomainCell( i );
            
            Numerics* numerics = analysisModel().domainManager().giveNumericsFor( curCell );
            
            if ( numerics == this )
            {
                auto cns = this->getNumericsStatusAt( curCell );
                std::vector<Dof*> dof = this->giveNodalDofsAt( curCell );
                RealVector uVec, phiVec;
                std::tie( uVec,phiVec ) = giveLocalVariablesAt( dof, current_value );
                RealVector dphi = cns->_dPsi * phiVec;
                
                totalCrackLength += 0.5 * cns->_area * ( _lc * dphi.dot( dphi ) + phiVec.dot( _massMatrix * phiVec ) / _lc );
            }
        }
        
        std::printf( "\n    Total crack length = %E\n", totalCrackLength );
    }
}
// ----------------------------------------------------------------------------
void PhaseFieldFracture_Fe_Tri3::readAdditionalDataFrom( FILE *fp )
{
    verifyKeyword( fp, "AnalysisMode", _name );
    std::string mode = getStringInputFrom( fp, "Failed to read analysis mode from input file!", _name );
    if ( mode == "PlaneStress" )
        _analysisMode = PlaneStress;
    else if ( mode == "PlaneStrain" )
        _analysisMode = PlaneStrain;
    else
        throw std::runtime_error( "ERROR: Invalid analysis mode '" + mode + "' encountered in input file!\nSource: " + _name );

    verifyKeyword( fp, "CharacteristicLength", _name );
    _lc = getRealInputFrom( fp, "Failed to read characteristic length from input file!", _name );
    
    verifyKeyword( fp, "IrreversibilityThreshold", _name );
    _phiIrrev = getRealInputFrom( fp, "Failed to read irreversibility threshold from input file!", _name );
}
// ----------------------------------------------------------------------------
void PhaseFieldFracture_Fe_Tri3::setDofStagesAt( Cell* targetCell )
{    
    // A. Nodal DOFs
    // Numerics uses only one stage (_stage[0])
    std::vector<Node*> node = analysisModel().domainManager().giveNodesOf( targetCell );
    
    for ( auto& curNode : node)
    {
        Dof *dof_x, *dof_y, *dof_phi;
        
        dof_x = analysisModel().domainManager().giveNodalDof( _nodalDof[ 0 ], curNode );
        dof_y = analysisModel().domainManager().giveNodalDof( _nodalDof[ 1 ], curNode );
        dof_phi = analysisModel().domainManager().giveNodalDof( _nodalDof[ 2 ], curNode );

        DofManager::setStageFor( dof_x, _stage[ 0 ] );
        DofManager::setStageFor( dof_y, _stage[ 0 ] );
        DofManager::setStageFor( dof_phi, _stage[ 0 ] );
        
        // Mechanics is assigned to _subsytem[0], phase-field equation to _subsystem[1]
        DofManager::setSubsystemFor( dof_x, _subsystem[ 0 ] );
        DofManager::setSubsystemFor( dof_y, _subsystem[ 0 ] );
        DofManager::setSubsystemFor( dof_phi, _subsystem[ 1 ] );
    }
}

// Private methods
// ----------------------------------------------------------------------------
NumericsStatus_PhaseFieldFracture_Fe_Tri3*
PhaseFieldFracture_Fe_Tri3::getNumericsStatusAt( Cell* targetCell )
{
    auto cns = dynamic_cast<NumericsStatus_PhaseFieldFracture_Fe_Tri3*>( targetCell->numericsStatus );
    if ( !cns )
        throw std::runtime_error( "Error: Unable to retrieve numerics status at cell!\nSource: " + _name );
    
    return cns;
}
// ----------------------------------------------------------------------------
RealMatrix PhaseFieldFracture_Fe_Tri3::giveBmatUAt( Cell* targetCell )
{
    auto cns = this->getNumericsStatusAt( targetCell );
    RealMatrix& dpsi = cns->_dPsi;
    
    // Fill entries for BmatU
    RealMatrix bmatU;
    if ( _analysisMode == PlaneStrain )
    {
        bmatU = { { dpsi( 0,0 ), 0.,          dpsi( 0,1 ), 0.,          dpsi( 0,2 ), 0. },
                  { 0.,          dpsi( 1,0 ), 0.,          dpsi( 1,1 ), 0.,          dpsi( 1,2 ) },
                  { 0.,          0.,          0.,          0.,          0.,          0. },
                  { dpsi( 1,0 ), dpsi( 0,0 ), dpsi( 1,1 ), dpsi( 0,1 ), dpsi( 1,2 ), dpsi( 0,2 ) } };
    }
    else if ( _analysisMode == PlaneStress )
    {
        bmatU = { { dpsi( 0,0 ), 0.,          dpsi( 0,1 ), 0.,          dpsi( 0,2 ), 0. },
                  { 0.,          dpsi( 1,0 ), 0.,          dpsi( 1,1 ), 0.,          dpsi( 1,2 ) },
                  { dpsi( 1,0 ), dpsi( 0,0 ), dpsi( 1,1 ), dpsi( 0,1 ), dpsi( 1,2 ), dpsi( 0,2 ) } };
    }

    return bmatU;
}
// ----------------------------------------------------------------------------
RealMatrix PhaseFieldFracture_Fe_Tri3::giveJacobianMatrixAt( Cell* targetCell )
{
    std::vector< Node* > node = analysisModel().domainManager().giveNodesOf( targetCell );
#ifndef NDEBUG
    if ( (int)node.size() != 3 )
        throw std::runtime_error( "\nError: Basis function requires 3 nodes be specified for calculation of Jacobian matrix!\nSource: " + _name );
#endif
    
    RealVector nc0 = analysisModel().domainManager().giveCoordinatesOf( node[ 0 ] );
    RealVector nc1 = analysisModel().domainManager().giveCoordinatesOf( node[ 1 ] );
    RealVector nc2 = analysisModel().domainManager().giveCoordinatesOf( node[ 2 ] );

    RealMatrix coorMat( { { nc0( 0 ), nc0( 1 ) },
                          { nc1( 0 ), nc1( 1 ) },
                          { nc2( 0 ), nc2( 1 ) } } );

    return _basisFunctionDerivatives * coorMat;
}
// ----------------------------------------------------------------------------
std::tuple< RealVector
          , RealVector >
PhaseFieldFracture_Fe_Tri3::giveLocalVariablesAt( std::vector<Dof*>& dof, ValueType valType )
{
    // Displacements
    RealVector u( 6 );
    u( 0 ) = DofManager::giveValueOfPrimaryVariableAt( dof[ 0 ], valType );
    u( 1 ) = DofManager::giveValueOfPrimaryVariableAt( dof[ 1 ], valType );
    u( 2 ) = DofManager::giveValueOfPrimaryVariableAt( dof[ 2 ], valType );
    u( 3 ) = DofManager::giveValueOfPrimaryVariableAt( dof[ 3 ], valType );
    u( 4 ) = DofManager::giveValueOfPrimaryVariableAt( dof[ 4 ], valType );
    u( 5 ) = DofManager::giveValueOfPrimaryVariableAt( dof[ 5 ], valType );
    
    // Phase-field
    RealVector phi( 3 );
    phi( 0 ) = DofManager::giveValueOfPrimaryVariableAt( dof[ 6 ], valType );
    phi( 1 ) = DofManager::giveValueOfPrimaryVariableAt( dof[ 7 ], valType );
    phi( 2 ) = DofManager::giveValueOfPrimaryVariableAt( dof[ 8 ], valType );
    
    return std::make_tuple( std::move( u ), std::move( phi ) );
}
// ----------------------------------------------------------------------------
std::vector< Dof* > PhaseFieldFracture_Fe_Tri3::giveNodalDofsAt( Cell* targetCell )
{
    std::vector< Node* > node = analysisModel().domainManager().giveNodesOf( targetCell );
    std::vector< Dof* > dof( 9, nullptr );
    
    dof[ 0 ] = analysisModel().domainManager().giveNodalDof( _nodalDof[ 0 ], node[ 0 ] );
    dof[ 1 ] = analysisModel().domainManager().giveNodalDof( _nodalDof[ 1 ], node[ 0 ] );
    dof[ 2 ] = analysisModel().domainManager().giveNodalDof( _nodalDof[ 0 ], node[ 1 ] );
    dof[ 3 ] = analysisModel().domainManager().giveNodalDof( _nodalDof[ 1 ], node[ 1 ] );
    dof[ 4 ] = analysisModel().domainManager().giveNodalDof( _nodalDof[ 0 ], node[ 2 ] );
    dof[ 5 ] = analysisModel().domainManager().giveNodalDof( _nodalDof[ 1 ], node[ 2 ] );
    dof[ 6 ] = analysisModel().domainManager().giveNodalDof( _nodalDof[ 2 ], node[ 0 ] );
    dof[ 7 ] = analysisModel().domainManager().giveNodalDof( _nodalDof[ 2 ], node[ 1 ] );
    dof[ 8 ] = analysisModel().domainManager().giveNodalDof( _nodalDof[ 2 ], node[ 2 ] );

    return dof;
}
