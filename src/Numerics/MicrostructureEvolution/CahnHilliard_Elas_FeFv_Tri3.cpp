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

#include "CahnHilliard_Elas_FeFv_Tri3.hpp"
#include "Core/AnalysisModel.hpp"
#include "Core/ObjectFactory.hpp"
#include "Materials/Material.hpp"
#include "Util/linearAlgebra.hpp"
#include "Util/readOperations.hpp"

using namespace broomstyx;

registerBroomstyxObject( Numerics, CahnHilliard_Elas_FeFv_Tri3 )

// Constructor for cell numerics status
CahnHilliard_Elas_FeFv_Tri3::CellNumericsStatus::CellNumericsStatus()
    : _area( 0. )
    , _c( 0. )
    , _Psi( 0. )
    , _strain( RealVector( 3 ) )
    , _stress( RealVector( 3 ) )
    , _Egy_chem( 0. )
    , _Egy_elas( 0. )
    , _ShapeFuncDeriv( RealMatrix( 2,3 ) )
    , _Cmat( RealMatrix( 3,3 ) )
    , _hasConcentrationGradientPrescribedOnFace { false, false, false }
    , _hasPsiGradientPrescribedOnFace{ false, false, false }
    , _hasNotComputedTransmissibilities( true )
    , _transmissibility{ 0., 0., 0. }
    , _materialStatus { nullptr, nullptr, nullptr, nullptr }
{}

// Constructor
CahnHilliard_Elas_FeFv_Tri3::CahnHilliard_Elas_FeFv_Tri3()
{
    _dim = 2;
    _dofPerCell = 2;
    _dofPerNode = 2;

    _nNodes = 3;
    _nMaterials = 4;
    _nStages = 1;
    _nSubsystems = 2;

    _name = "CahnHilliard_Elas_FeFv_Tri3";

    // Lone Gauss point coordinate and weight
    _gpNatCoor = { 1./3., 1./3. };
    _wt = 0.5;

    _M = 0.;
    _kappa = 0.;
    _A = 0.;
    _eps_m = 0.;
    _Nv = 0.;

    _I = { 1., 1., 0. };

    // Pre-calculate shape functions and derivatives at Gauss point
    _basisFunctionValues = _basisFunction.giveBasisFunctionsAt( _gpNatCoor );
    std::vector<RealVector> dpsiNat = _basisFunction.giveBasisFunctionDerivativesAt( _gpNatCoor );
    _basisFunctionDerivatives.init( 2,3 );
    for ( int i = 0; i < 3; i++)
    {
        _basisFunctionDerivatives( 0,i ) = dpsiNat[ 0 ]( i );
        _basisFunctionDerivatives( 1,i ) = dpsiNat[ 1 ]( i );
    }
}

// Public methods
// ----------------------------------------------------------------------------
void CahnHilliard_Elas_FeFv_Tri3::deleteNumericsAt( Cell* targetCell )
{
    auto cns = this->getNumericsStatusAt( targetCell );
    auto material = this->giveMaterialSetFor( targetCell );

    material[ 0 ]->destroy( cns->_materialStatus[ 0 ] );
    material[ 1 ]->destroy( cns->_materialStatus[ 1 ] );
    material[ 2 ]->destroy( cns->_materialStatus[ 2 ] );
    material[ 3 ]->destroy( cns->_materialStatus[ 3 ] );

    delete cns;
}
// ----------------------------------------------------------------------------
void CahnHilliard_Elas_FeFv_Tri3::finalizeDataAt( Cell* targetCell, const TimeData& time )
{
    // Pointer to numerics status
    auto cns = this->getNumericsStatusAt( targetCell );

    // Get coordinate of cell center

    // Finalize values at cells
    Dof* dof_c = analysisModel().domainManager().giveCellDof( _cellDof[ 0 ], targetCell );
    Dof* dof_Psi = analysisModel().domainManager().giveCellDof( _cellDof[ 1 ], targetCell );

    cns->_c = DofManager::giveValueOfPrimaryVariableAt( dof_c, converged_value );
    cns->_Psi = DofManager::giveValueOfPrimaryVariableAt( dof_Psi, converged_value );

    // Retrieve nodal DOFs local to element
    std::vector<Dof*> dof = this->giveNodalDofsAt( targetCell );

    // Retrieve DOF values
    RealVector uVec = giveLocalDisplacementsAt( dof, converged_value );

    // Retrieve material set for element
    std::vector<Material*> material = this->giveMaterialSetFor( targetCell );

    // Container for constitutive state
    RealVector conState;

    // Update internal variables for concentration functions
    conState = { cns->_c };
    material[ 2 ]->updateStatusFrom( conState, cns->_materialStatus[ 2 ] );
    material[ 3 ]->updateStatusFrom( conState, cns->_materialStatus[ 3 ] );

    // Calculate interpolation for misfit strain
    double beta = material[ 2 ]->givePotentialFrom( conState, cns->_materialStatus[ 2 ] );

    // Calculate elastic strain
    RealMatrix bmatU = this->giveBmatAt( targetCell );
    cns->_strain = bmatU * uVec - beta * _eps_m * _I;

    // Update internaal variables for elasticity tensors
    material[ 0 ]->updateStatusFrom( cns->_strain, cns->_materialStatus[ 0 ] );
    material[ 1 ]->updateStatusFrom( cns->_strain, cns->_materialStatus[ 1 ] );

    // Calculate elasticity tensors
    RealMatrix Cmat_m = material[ 0 ]->giveModulusFrom( cns->_strain, cns->_materialStatus[ 0 ] );
    RealMatrix Cmat_p = material[ 1 ]->giveModulusFrom( cns->_strain, cns->_materialStatus[ 1 ] );

    // Assemble elastiity tensor based on concentration value
    double alpha = material[ 3 ]->givePotentialFrom( conState, cns->_materialStatus[ 3 ] );
    cns->_Cmat = Cmat_m + alpha * ( Cmat_p - Cmat_m );

    // Calculate stress
    cns->_stress = cns->_Cmat * cns->_strain;

    // A. Calculate concentration flux at faces
    std::vector<Node*> node;
    node = analysisModel().domainManager().giveNodesOf( targetCell );
    std::vector<std::vector<Node*> > face = giveFaceNodesOf( targetCell );

    // Get cell neighbors
    std::vector<Cell*> neighbor = analysisModel().domainManager().giveNeighborsOf( targetCell );

    // Cycle through faces
    double c_normalGradient[ 3 ];
    for ( int i = 0; i < 3; i++ )
    {
        // Length of face
        double length = giveLengthOf( face[ i ] );

        if ( cns->_hasConcentrationGradientPrescribedOnFace[ i ] )
            c_normalGradient[ i ] = 0.;
        else
        {
            // phase-field at neighbor cell
            Dof* dof_c2 = analysisModel().domainManager().giveCellDof( _cellDof[ 0 ], neighbor[ i ] );
            double c2 = DofManager::giveValueOfPrimaryVariableAt( dof_c2, converged_value );

            c_normalGradient[ i ] = cns->_transmissibility[ i ] * ( c2 - cns->_c ) / length;
        }
    }

    // B. RT0 reconstruction of flux at cell center
    std::vector<Node*> vertexNode( 3, nullptr );
    vertexNode[ 0 ] = node[ 2 ];
    vertexNode[ 1 ] = node[ 0 ];
    vertexNode[ 2 ] = node[ 1 ];

    RealVector d_c( 3 );

    // Get coordinate of cell center
    std::vector<RealVector> epCoor = this->giveEvaluationPointsFor( targetCell );
    RealVector cellCoor = epCoor[ 0 ];

    for ( int i = 0; i < 3; i++ )
    {
        double length = giveLengthOf( face[ i ] );
        RealVector vtxNodeCoor = analysisModel().domainManager().giveCoordinatesOf( vertexNode[ i ] );
        d_c = d_c + c_normalGradient[ i ] * length / ( 2. * cns->_area ) * ( cellCoor - vtxNodeCoor );
    }
    RealVector grad_c( { d_c( 0 ), d_c( 1 ) } );

    // Calculate chemical energy density
    cns->_Egy_chem = _Nv * ( _A * cns->_c * cns->_c * ( 1. - cns->_c ) * ( 1. - cns->_c ) + _kappa * grad_c.dot( grad_c ) );
    cns->_Egy_elas = _Nv * 0.5 * cns->_stress.dot( cns->_strain );
}
// ----------------------------------------------------------------------------
double CahnHilliard_Elas_FeFv_Tri3::giveCellFieldValueAt( Cell* targetCell, int fieldNum )
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
RealVector CahnHilliard_Elas_FeFv_Tri3::giveCellNodeFieldValuesAt( Cell* targetCell, int fieldNum )
{
    double val = this->giveCellFieldValueAt( targetCell, fieldNum );
    RealVector cellNodeValues( { val, val, val } );

    return cellNodeValues;
}
// ----------------------------------------------------------------------------
std::vector< RealVector > CahnHilliard_Elas_FeFv_Tri3::giveEvaluationPointsFor( Cell *targetCell )
{
    std::vector<Node*> node = analysisModel().domainManager().giveNodesOf( targetCell );

    std::vector<RealVector> coor( 1, RealVector() );
    coor[ 0 ].init( 3 );

    for ( int i = 0; i < 3; i++ )
    {
        RealVector nodeCoor = analysisModel().domainManager().giveCoordinatesOf(node[i]);
        coor[ 0 ]( 0 ) += nodeCoor( 0 ) / 3.0;
        coor[ 0 ]( 1 ) += nodeCoor( 1 ) / 3.0;
        coor[ 0 ]( 2 ) += nodeCoor( 2 ) / 3.0;
    }

    return coor;
}
// ----------------------------------------------------------------------------
std::tuple< RealVector, RealVector >
CahnHilliard_Elas_FeFv_Tri3::giveFieldOutputAt( Cell* targetCell, const std::string& fieldTag )
{
    RealVector fieldVal( 1 ), weight( 1 );

    auto cns = this->getNumericsStatusAt( targetCell );
    weight( 0 ) = cns->_area;

    if ( fieldTag == "unassigned" )
        fieldVal( 0 ) = 0.;
    else if ( fieldTag == "s_xx" )
        fieldVal( 0 ) = cns->_stress( 0 );
    else if ( fieldTag == "s_yy" )
        fieldVal( 0 ) = cns->_stress( 1 );
    else if ( fieldTag == "s_xy" )
        fieldVal( 0 ) = cns->_stress( 2 );
    else if ( fieldTag == "e_xx" )
        fieldVal( 0 ) = cns->_strain( 0 );
    else if ( fieldTag == "e_yy" )
        fieldVal( 0 ) = cns->_strain( 1 );
    else if ( fieldTag == "g_xy" )
        fieldVal( 0 ) = cns->_strain( 2 );
    else if ( fieldTag == "c" )
    {
        Dof* dof_c = analysisModel().domainManager().giveCellDof( _cellDof[ 0 ], targetCell );
        cns->_c = DofManager::giveValueOfPrimaryVariableAt( dof_c, converged_value );
        fieldVal(0) = cns->_c;
    }
    else if ( fieldTag == "Psi" )
    {
        Dof* dof_Psi = analysisModel().domainManager().giveCellDof( _cellDof[ 1 ], targetCell );
        cns->_Psi = DofManager::giveValueOfPrimaryVariableAt( dof_Psi, converged_value );
        fieldVal(0) = cns->_Psi;
    }
    else if ( fieldTag == "egy_chem" )
        fieldVal( 0 ) = cns->_Egy_chem;
    else if ( fieldTag == "egy_elas" )
        fieldVal( 0 ) = cns->_Egy_elas;
    else
        throw std::runtime_error( "Invalid tag '" + fieldTag + "' supplied in field output request!\nSource: " + _name );

    return std::make_tuple(std::move( fieldVal ), std::move( weight ) );
}
// ----------------------------------------------------------------------------
std::tuple< std::vector< Dof* >, std::vector< Dof* >, RealVector >
CahnHilliard_Elas_FeFv_Tri3::giveStaticCoefficientMatrixAt( Cell*           targetCell
                                                          , int             stage
                                                          , int             subsys
                                                          , const TimeData& time )
{
    std::vector<Dof*> rowDof, colDof;
    RealVector coefVal;

    if ( stage == _stage[ 0 ] )
    {
        auto cns = this->getNumericsStatusAt( targetCell );

        // Get cell neighbors
        std::vector<Cell*> neighbor = analysisModel().domainManager().giveNeighborsOf( targetCell );

        int nCols = 0;
        for ( int i = 0; i < 3; i++ )
        {
            if ( !cns->_hasConcentrationGradientPrescribedOnFace[ i ] )
                ++nCols;
            if ( !cns->_hasPsiGradientPrescribedOnFace[ i ] )
                ++nCols;
        }

        int vecLength = 39 + nCols; // Default output lengths for unassigned subsystems

        if ( subsys == _subsystem[ 0 ] )
            vecLength = 36;
        else if ( subsys == _subsystem[ 1 ] )
            vecLength = 3 + nCols;

        rowDof.assign( vecLength, nullptr );
        colDof.assign( vecLength, nullptr );
        coefVal.init( vecLength );

        // Retrieve nodal DOFs local to element
        std::vector< Dof* > dof = giveNodalDofsAt( targetCell );

        int counter = 0;
        if ( subsys == _subsystem[ 0 ] || subsys == UNASSIGNED )
        {
            // Cmat is already stored when calculating LHS, so no need to recompute
            // Just need to assemble local stiffness matrix
            RealMatrix bmatU = this->giveBmatAt( targetCell );
            RealMatrix kmatUU = cns->_area * trp( bmatU ) * cns->_Cmat * bmatU;

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
            int startIdx = counter;

            Dof* dof_c = analysisModel().domainManager().giveCellDof( _cellDof[ 0 ], targetCell );
            Dof* dof_Psi = analysisModel().domainManager().giveCellDof( _cellDof[ 1 ], targetCell );

            // Retrieve material set for element
            std::vector<Material*> material = this->giveMaterialSetFor( targetCell );

            rowDof[ startIdx ] = dof_Psi;
            colDof[ startIdx ] = dof_Psi;
            coefVal( startIdx ) = 0.;

            // Tangent matrix related to r_Psi

            // Cycle through faces
            for ( int i = 0; i < 3; i++ )
            {
                if ( !cns->_hasPsiGradientPrescribedOnFace[ i ] )
                {
                    // DOF address for neighbor cell
                    ++counter;
                    rowDof[ counter ] = dof_Psi;
                    colDof[ counter ] = analysisModel().domainManager().giveCellDof( _cellDof[ 1 ], neighbor[ i ] );

                    coefVal( startIdx ) += cns->_transmissibility[ i ];
                    coefVal( counter ) = -cns->_transmissibility[ i ];
                }
            }

            ++counter;
            startIdx = counter;

            // Tangent matrix related to r_c

            // Calculate d/dc ( delta_Felas_delta_c )
            RealMatrix Cmat_m = material[ 0 ]->giveModulusFrom( cns->_strain, cns->_materialStatus[ 0 ] );
            RealMatrix Cmat_p = material[ 1 ]->giveModulusFrom( cns->_strain, cns->_materialStatus[ 1 ] );
            RealMatrix Delta_Cmat = Cmat_p - Cmat_m;

            RealVector conState( { cns->_c } );
            RealVector d_beta = material[ 2 ]->giveForceFrom( conState, cns->_materialStatus[ 2 ] );
            RealVector d_alpha = material[ 3 ]->giveForceFrom( conState, cns->_materialStatus[ 3 ] );
            RealMatrix dd_beta = material[ 2 ]->giveModulusFrom( conState, cns->_materialStatus[ 2 ] );
            RealMatrix dd_alpha = material[ 3 ]->giveModulusFrom( conState, cns->_materialStatus[ 3 ] );

            double d_dc_delt_Fel = _Nv / 2. * (
                    cns->_strain.dot( dd_alpha( 0,0 ) * Delta_Cmat * cns->_strain )
                    - 4. * cns->_strain.dot( d_alpha( 0 ) * Delta_Cmat * ( -d_beta( 0 ) * _eps_m * _I ) )
                    + 2. * d_beta( 0 ) * _eps_m * _I.dot( cns->_Cmat * ( d_beta( 0 ) * _eps_m * _I ) )
                    - 2. * cns->_strain.dot( cns->_Cmat * ( dd_beta( 0,0 ) * _eps_m * _I ) )
            );

            rowDof[ startIdx ] = dof_c;
            colDof[ startIdx ] = dof_c;
            coefVal( startIdx ) = cns->_area * ( _A * ( 6. * cns->_c * cns->_c - 6. * cns->_c + 1. ) + 0.5 / ( _Nv * _kappa ) * d_dc_delt_Fel );

            // Cycle through faces
            for ( int i = 0; i < 3; i++ )
            {
                if ( !cns->_hasConcentrationGradientPrescribedOnFace[ i ] )
                {
                    // DOF address for neighbor cell
                    ++counter;
                    rowDof[ counter ] = dof_c;
                    colDof[ counter ] = analysisModel().domainManager().giveCellDof( _cellDof[ 0 ], neighbor[ i ] );

                    coefVal( startIdx ) += cns->_transmissibility[ i ];
                    coefVal( counter ) = -cns->_transmissibility[ i ];
                }
            }

            ++counter;
            rowDof[ counter ] = dof_c;
            colDof[ counter ] = dof_Psi;
            coefVal( counter ) = -0.5 * cns->_area / ( _Nv * _kappa );
        }
    }

    return std::make_tuple( std::move( rowDof ), std::move( colDof ), std::move( coefVal ) );
}
// ----------------------------------------------------------------------------
std::tuple< std::vector< Dof* >, RealVector >
CahnHilliard_Elas_FeFv_Tri3::giveStaticLeftHandSideAt( Cell*           targetCell
                                                     , int             stage
                                                     , int             subsys
                                                     , const TimeData& time )
{
    std::vector<Dof*> rowDof;
    RealVector lhs;

    if ( stage == _stage[ 0 ] )
    {
        auto cns = this->getNumericsStatusAt( targetCell );
        int vecLength = 13; // Default output lengths for unassigned subsystems

        if ( subsys == _subsystem[ 0 ] )
            vecLength = 6;
        else if ( subsys == _subsystem[ 1 ] )
            vecLength = 7;

        rowDof.assign( vecLength, nullptr );
        lhs.init( vecLength );

        // Retrieve nodal DOFs local to element
        std::vector<Dof*> dof = this->giveNodalDofsAt( targetCell );

        // Retrieve current value of local displacements
        RealVector uVec = giveLocalDisplacementsAt( dof, current_value );

        // Retrieve concentration
        Dof* dof_c = analysisModel().domainManager().giveCellDof( _cellDof[ 0 ], targetCell );
        cns->_c = DofManager::giveValueOfPrimaryVariableAt( dof_c, current_value );

        // Retrieve material set for element
        std::vector<Material*> material = this->giveMaterialSetFor( targetCell );

        // Container for constitutive state
        RealVector conState;

        // Update internal variables for concentration functions
        conState = { cns->_c };
        material[ 2 ]->updateStatusFrom( conState, cns->_materialStatus[ 2 ] );
        material[ 3 ]->updateStatusFrom( conState, cns->_materialStatus[ 3 ] );

        // Calculate interpolation for misfit strain
        double beta = material[ 2 ]->givePotentialFrom( conState, cns->_materialStatus[ 2 ] );

        // Calculate elastic strain
        RealMatrix bmatU = this->giveBmatAt( targetCell );
        cns->_strain = bmatU * uVec - beta * _eps_m * _I;

        // Update internaal variables for elasticity tensors
        material[ 0 ]->updateStatusFrom( cns->_strain, cns->_materialStatus[ 0 ] );
        material[ 1 ]->updateStatusFrom( cns->_strain, cns->_materialStatus[ 1 ] );

        // Calculate elasticity tensors
        RealMatrix Cmat_m = material[ 0 ]->giveModulusFrom( cns->_strain, cns->_materialStatus[ 0 ] );
        RealMatrix Cmat_p = material[ 1 ]->giveModulusFrom( cns->_strain, cns->_materialStatus[ 1 ] );

        // Assemble elastiity tensor based on concentration value
        double alpha = material[ 3 ]->givePotentialFrom( conState, cns->_materialStatus[ 3 ] );
        cns->_Cmat = Cmat_m + alpha * ( Cmat_p - Cmat_m );

        int offset = 0;
        if ( subsys == _subsystem[ 0 ] || subsys == UNASSIGNED )
        {
            // Stress equilibrium residual
            rowDof[ 0 ] = dof[ 0 ];
            rowDof[ 1 ] = dof[ 1 ];
            rowDof[ 2 ] = dof[ 2 ];
            rowDof[ 3 ] = dof[ 3 ];
            rowDof[ 4 ] = dof[ 4 ];
            rowDof[ 5 ] = dof[ 5 ];
            offset = 6;

            // Calculate stress
            cns->_stress = cns->_Cmat * cns->_strain;

            // Calculate FmatU
            RealVector fmatU;
            fmatU = cns->_area * trp( bmatU ) * cns->_stress;

            lhs( 0 ) = fmatU( 0 );
            lhs( 1 ) = fmatU( 1 );
            lhs( 2 ) = fmatU( 2 );
            lhs( 3 ) = fmatU( 3 );
            lhs( 4 ) = fmatU( 4 );
            lhs( 5 ) = fmatU( 5 );
        }

        if ( subsys == _subsystem[ 1 ] || subsys == UNASSIGNED )
        {
            // Get faces of target cell
            std::vector<std::vector<Node*> > face = giveFaceNodesOf( targetCell );

            // Get cell neighbors
            std::vector<Cell*> neighbor = analysisModel().domainManager().giveNeighborsOf( targetCell );

            // Cycle through faces and compute transmissibilities if not done yet
            if ( cns->_hasNotComputedTransmissibilities )
            {
                for ( int i = 0; i < 3; i++ )
                    if ( !cns->_hasConcentrationGradientPrescribedOnFace[ i ] || !cns->_hasPsiGradientPrescribedOnFace[ i ] )
                        cns->_transmissibility[ i ] = this->giveTransmissibilityCoefficientAt( face[ i ], targetCell, neighbor[ i ] );

                cns->_hasNotComputedTransmissibilities = false;
            }

            // Retrieve DOF for Psi and its value
            Dof* dof_Psi = analysisModel().domainManager().giveCellDof( _cellDof[ 1 ], targetCell );
            cns->_Psi = DofManager::giveValueOfPrimaryVariableAt( dof_Psi, current_value );

            // Compute components of r_Psi
            RealVector Psi_normFlux( 3 );
            for ( int i = 0; i < 3; i++ )
            {
                if ( !cns->_hasConcentrationGradientPrescribedOnFace[ i ] )
                {
                    // Phase-field at neighbor cell
                    Dof* dof_Psi2 = analysisModel().domainManager().giveCellDof( _cellDof[ 1 ], neighbor[ i ] );
                    double Psi2 = DofManager::giveValueOfPrimaryVariableAt( dof_Psi2, current_value );

                    Psi_normFlux( i ) = cns->_transmissibility[ i ] * ( cns->_Psi - Psi2 );
                }
            }

            // Compute components of r_c
            // delta_Felas/delta_c
            RealVector d_beta = material[ 2 ]->giveForceFrom( conState, cns->_materialStatus[ 2 ] );
            RealVector d_alpha = material[ 3 ]->giveForceFrom( conState, cns->_materialStatus[ 3 ] );
            double delta_Felas_delta_c = ( _Nv / 2. ) *
                    ( d_alpha( 0 ) * cns->_strain.dot( ( Cmat_p - Cmat_m ) * cns->_strain )
                    - 2. * d_beta( 0 ) * _eps_m * cns->_strain.dot( cns->_Cmat * _I ) );

            // Concentration dependent terms
            double conc_term = _A * ( 2. * std::pow( cns->_c, 3. ) - 3. * std::pow( cns->_c, 2. ) + cns->_c );
            double c_sourceTerm = cns->_area * ( conc_term + ( delta_Felas_delta_c - cns->_Psi ) / ( 2. * _Nv * _kappa ) );

            RealVector c_normFlux( 3 );
            for ( int i = 0; i < 3; i++ )
            {
                if ( !cns->_hasConcentrationGradientPrescribedOnFace[ i ] )
                {
                    // Phase-field at neighbor cell
                    Dof* dof_c2 = analysisModel().domainManager().giveCellDof( _cellDof[ 0 ], neighbor[ i ] );
                    double c2 = DofManager::giveValueOfPrimaryVariableAt( dof_c2, current_value );

                    c_normFlux( i ) = cns->_transmissibility[ i ] * ( cns->_c - c2 );
                }
            }

            // r_psi components
            rowDof[ offset ] = dof_Psi;
            rowDof[ offset + 1 ] = dof_Psi;
            rowDof[ offset + 2 ] = dof_Psi;

            lhs( offset ) = Psi_normFlux( 0 );
            lhs( offset + 1 ) = Psi_normFlux( 1 );
            lhs( offset + 2 ) = Psi_normFlux( 2 );

            // r_c components
            rowDof[ offset + 3 ] = dof_c;
            rowDof[ offset + 4 ] = dof_c;
            rowDof[ offset + 5 ] = dof_c;
            rowDof[ offset + 6 ] = dof_c;

            // Calculate FmatPhi
            lhs( offset + 3 ) = c_normFlux( 0 );
            lhs( offset + 4 ) = c_normFlux( 1 );
            lhs( offset + 5 ) = c_normFlux( 2 );
            lhs( offset + 6 ) = c_sourceTerm;
        }
    }

    return std::make_tuple( std::move( rowDof ), std::move( lhs ) );
}
// ----------------------------------------------------------------------------
std::tuple< std::vector< Dof* >, RealVector >
CahnHilliard_Elas_FeFv_Tri3::giveStaticRightHandSideAt( Cell*                    targetCell
                                                      , int                      stage
                                                      , int                      subsys
                                                      , const BoundaryCondition& bndCond
                                                      , const TimeData&          time )
{
    std::vector<Dof*> rowDof;
    RealVector rhs;

    if ( stage == _stage[ 0 ] )
    {
        if ( bndCond.conditionType() == "Traction" )
        {
            if ( subsys == _subsystem[ 0 ] || subsys == UNASSIGNED )
            {
                // ---------------------------------------------
                // Note: Boundary element must be a 2-node line!
                // ---------------------------------------------

                // Retrieve nodes of boundary element
                std::vector<Node*> node = analysisModel().domainManager().giveNodesOf( targetCell );

                if ( (int)node.size() != 2 )
                    throw std::runtime_error( "Incompatible boundary cell encountered for '" + _name + "'!" );

                // Compute length of boundary element
                RealVector coor0, coor1, dx;
                coor0 = analysisModel().domainManager().giveCoordinatesOf( node[ 0 ] );
                coor1 = analysisModel().domainManager().giveCoordinatesOf( node[ 1 ] );
                dx = coor1 - coor0;
                double length = std::sqrt( dx.dot( dx ) );

                // BC evaluated at midpoint of boundary cell
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
            double bcVal = bndCond.valueAt( coor, time );

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
std::tuple< std::vector< Dof* >, std::vector< Dof* >, RealVector >
CahnHilliard_Elas_FeFv_Tri3::giveTransientCoefficientMatrixAt( Cell*           targetCell
                                                             , int             stage
                                                             , int             subsys
                                                             , const TimeData& time )
{
    std::vector<Dof*> rowDof, colDof;
    RealVector coefVal;

    if ( stage == _stage[ 0 ] && (subsys == _subsystem[ 1 ] || subsys == UNASSIGNED ) )
    {
        auto cns = this->getNumericsStatusAt( targetCell );

        Dof* dof_c = analysisModel().domainManager().giveCellDof( _cellDof[ 0 ], targetCell );
        Dof* dof_Psi = analysisModel().domainManager().giveCellDof( _cellDof[ 1 ], targetCell );

        rowDof.assign( 1, dof_Psi );
        colDof.assign( 1, dof_c );
        coefVal = { cns->_area * _M };
    }

    return std::make_tuple( std::move( rowDof ), std::move( colDof ), std::move( coefVal ) );
}
// ----------------------------------------------------------------------------
std::tuple< std::vector< Dof* >, RealVector >
CahnHilliard_Elas_FeFv_Tri3::giveTransientLeftHandSideAt( Cell*           targetCell
                                                        , int             stage
                                                        , int             subsys
                                                        , const TimeData& time
                                                        , ValueType       valType )
{
    std::vector<Dof*> rowDof;
    RealVector lhs;

    if ( stage == _stage[ 0 ] && ( subsys == _subsystem[ 1 ] || subsys == UNASSIGNED ) )
    {
        auto cns = this->getNumericsStatusAt( targetCell );

        Dof* dof_c = analysisModel().domainManager().giveCellDof( _cellDof[ 0 ], targetCell );
        double c = DofManager::giveValueOfPrimaryVariableAt( dof_c, valType );

        Dof* dof_Psi = analysisModel().domainManager().giveCellDof( _cellDof[ 1 ], targetCell );

        rowDof.assign( 1, dof_Psi );
        lhs = { cns->_area * _M * c };
    }

    return std::make_tuple( std::move( rowDof ), std::move( lhs ) );
}
// ----------------------------------------------------------------------------
void CahnHilliard_Elas_FeFv_Tri3::imposeConstraintAt( Cell*                    targetCell
                                                    , int                      stg
                                                    , const BoundaryCondition& bndCond
                                                    , const TimeData&          time )
{
    // Essential BCs on nodal DOFs
    if ( bndCond.conditionType() == "NodalConstraint" )
    {
        // Retrieve nodes of boundary element
        std::vector< Node* > node = analysisModel().domainManager().giveNodesOf( targetCell );
        int dofNum = analysisModel().dofManager().giveIndexForNodalDof( bndCond.targetDof() );

        for ( auto& curNode : node)
        {
            Dof* targetDof = analysisModel().domainManager().giveNodalDof( dofNum, curNode );

            RealVector coor = analysisModel().domainManager().giveCoordinatesOf( curNode );
            double bcVal = bndCond.valueAt( coor, time );

            DofManager::setConstraintValueAt( targetDof, bcVal );
        }
    }

    // Constraints pertaining to cell DOFS
    if ( bndCond.conditionType() == "ConcentrationGradient" || bndCond.conditionType() == "PsiGradient" )
    {
        std::vector<Cell*> domCell = analysisModel().domainManager().giveDomainCellsAssociatedWith( targetCell );
        for ( Cell* curDomCell : domCell )
        {
            Numerics* domCellNumerics = analysisModel().domainManager().giveNumericsFor( curDomCell );
            if ( this == domCellNumerics )
            {
                auto cns = this->getNumericsStatusAt( curDomCell );
                std::vector<Node*> bndCellNode = analysisModel().domainManager().giveNodesOf( targetCell );

                // Sanity check on mesh
                if ( bndCellNode.size() != 2 )
                    throw std::runtime_error( "Numerics type '" + _name + "' only accepts boundary cells with exactly two nodes!" );

                std::vector< std::vector<Node*> > face = giveFaceNodesOf( curDomCell );

                // Cycle through faces to determine which one corresponds to current boundary cell
                for ( int i = 0; i < 3; i++ )
                {
                    if ( (bndCellNode[ 0 ] == face[ i ][ 0 ] && bndCellNode[ 1 ] == face[ i ][ 1 ] ) ||
                         (bndCellNode[ 1 ] == face[ i ][ 0 ] && bndCellNode[ 0 ] == face[ i ][ 1 ] ) )
                    {
                        // Set boundary flag for appropriate face
                        if ( bndCond.conditionType() == "ConcentrationGradient" )
                            cns->_hasConcentrationGradientPrescribedOnFace[ i ] = true;
                        else
                            cns->_hasPsiGradientPrescribedOnFace[ i ] = true;
                    }
                }
            }
        }
    }
}
// ----------------------------------------------------------------------------
void CahnHilliard_Elas_FeFv_Tri3::imposeInitialConditionAt( Cell* targetCell, const InitialCondition& initCond )
{
    if ( _cellDof[ 0 ] != initCond.targetDofNumber() )
        throw std::runtime_error( "Error imposing initial condition! Specified DOF in initial condition \
            does not correspond to the cell phase-field DOF.\nSource: " + _name );
    else
    {
        std::vector<RealVector> epCoor = this->giveEvaluationPointsFor( targetCell );
        double initVal = initCond.valueAt( epCoor[ 0 ] );

        Dof* pfDof = analysisModel().domainManager().giveCellDof( _cellDof[ 0 ], targetCell );
        DofManager::updatePrimaryVariableAt( pfDof, initVal, converged_value );
    }
}
// ----------------------------------------------------------------------------
void CahnHilliard_Elas_FeFv_Tri3::initializeMaterialsAt( Cell* targetCell )
{
    auto cns = this->getNumericsStatusAt( targetCell );
    std::vector<Material*> material = this->giveMaterialSetFor( targetCell );

    cns->_materialStatus[ 0 ] = material[ 0 ]->createMaterialStatus();
    cns->_materialStatus[ 1 ] = material[ 1 ]->createMaterialStatus();
    cns->_materialStatus[ 2 ] = material[ 2 ]->createMaterialStatus();
    cns->_materialStatus[ 3 ] = material[ 3 ]->createMaterialStatus();
}
// ----------------------------------------------------------------------------
void CahnHilliard_Elas_FeFv_Tri3::initializeNumericsAt( Cell* targetCell )
{
    targetCell->numericsStatus = new CellNumericsStatus();
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
    cns->_ShapeFuncDeriv = JmatInv * _basisFunctionDerivatives;
}
// ----------------------------------------------------------------------------
void CahnHilliard_Elas_FeFv_Tri3::readAdditionalDataFrom( FILE* fp )
{
    verifyKeyword( fp, "Mobility", _name );
    _M = getRealInputFrom( fp, "Failed to read mobility from input file!", _name );
    verifyKeyword( fp, "GradientEnergyCoef", _name );
    _kappa = getRealInputFrom( fp, "Failed to read gradient energy coefficient from input file!", _name );
    verifyKeyword( fp, "EnergyBarrier", _name );
    _A = getRealInputFrom( fp, "Failed to read energy barrier from input file!", _name );
    verifyKeyword( fp, "LatticeMisfit", _name );
    _eps_m = getRealInputFrom( fp, "Failed to read lattice misfit from input file!", _name );
    verifyKeyword( fp, "AtomsPerUnitVolume", _name );
    _Nv = getRealInputFrom( fp, "Failed to read number of atoms per unit volume from input file!", _name );
}
// ----------------------------------------------------------------------------
void CahnHilliard_Elas_FeFv_Tri3::removeConstraintsOn( Cell* targetCell )
{
    auto cns = this->getNumericsStatusAt(targetCell);
    cns->_hasConcentrationGradientPrescribedOnFace[ 0 ] = false;
    cns->_hasConcentrationGradientPrescribedOnFace[ 1 ] = false;
    cns->_hasConcentrationGradientPrescribedOnFace[ 2 ] = false;
    cns->_hasPsiGradientPrescribedOnFace[ 0 ] = false;
    cns->_hasPsiGradientPrescribedOnFace[ 1 ] = false;
    cns->_hasPsiGradientPrescribedOnFace[ 2 ] = false;
    cns->_hasNotComputedTransmissibilities = true;
}
// ----------------------------------------------------------------------------
void CahnHilliard_Elas_FeFv_Tri3::setDofStagesAt( Cell* targetCell )
{
    // A. Nodal DOFs
    std::vector<Node*> node = analysisModel().domainManager().giveNodesOf( targetCell );
    for ( auto& curNode : node )
    {
        Dof *dof_x, *dof_y;

        dof_x = analysisModel().domainManager().giveNodalDof( _nodalDof[ 0 ], curNode );
        dof_y = analysisModel().domainManager().giveNodalDof( _nodalDof[ 1 ], curNode );

        // Set DOF stage numbers
        DofManager::setStageFor( dof_x, _stage[ 0 ] );
        DofManager::setStageFor( dof_y, _stage[ 0 ] );

        // Set DOF subsystem numbers
        DofManager::setSubsystemFor( dof_x, _subsystem[ 0 ] );
        DofManager::setSubsystemFor( dof_y, _subsystem[ 0 ] );
    }

    // B. Cell DOFs
    Dof* dof_c = analysisModel().domainManager().giveCellDof( _cellDof[ 0 ], targetCell );
    Dof* dof_psi = analysisModel().domainManager().giveCellDof( _cellDof[ 1 ], targetCell );

    DofManager::setStageFor( dof_c, _stage[ 0 ] );
    DofManager::setSubsystemFor( dof_c, _subsystem[ 1 ] );
    DofManager::setStageFor( dof_psi, _stage[ 0 ] );
    DofManager::setSubsystemFor( dof_psi, _subsystem[ 1 ] );

}

// Private methods
// ----------------------------------------------------------------------------
CahnHilliard_Elas_FeFv_Tri3::CellNumericsStatus*
CahnHilliard_Elas_FeFv_Tri3::getNumericsStatusAt( Cell* targetCell )
{
    auto cns = dynamic_cast< CellNumericsStatus* >( targetCell->numericsStatus );
    if ( !cns )
        throw std::runtime_error( "Error: Unable to retrieve numerics status at cell!\nSource: " + _name );

    return cns;
}
// ----------------------------------------------------------------------------
RealMatrix CahnHilliard_Elas_FeFv_Tri3::giveBmatAt( Cell* targetCell )
{
    auto cns = this->getNumericsStatusAt( targetCell );

    // Fill entries for Bmat
    RealMatrix& dpsi = cns->_ShapeFuncDeriv;
    RealMatrix bmat( { { dpsi( 0,0 ), 0.,          dpsi( 0,1 ), 0.,          dpsi( 0,2 ), 0. },
                       { 0.,          dpsi( 1,0 ), 0.,          dpsi( 1,1 ), 0.,          dpsi( 1,2 ) },
                       { dpsi( 1,0 ), dpsi( 0,0 ), dpsi( 1,1 ), dpsi( 0,1 ), dpsi( 1,2 ), dpsi( 0,2 ) } } );
    return bmat;
}
// ----------------------------------------------------------------------------
std::vector<std::vector<Node*> > CahnHilliard_Elas_FeFv_Tri3::giveFaceNodesOf( Cell* targetCell )
{
    std::vector<Node*> node = analysisModel().domainManager().giveNodesOf( targetCell );
    std::vector<Node*> faceTemplate( 2, nullptr );
    std::vector<std::vector<Node*> > face( 3, faceTemplate );

    face[ 0 ][ 0 ] = node[ 0 ];
    face[ 0 ][ 1 ] = node[ 1 ];
    face[ 1 ][ 0 ] = node[ 1 ];
    face[ 1 ][ 1 ] = node[ 2 ];
    face[ 2 ][ 0 ] = node[ 2 ];
    face[ 2 ][ 1 ] = node[ 0 ];

    return face;
}
// ----------------------------------------------------------------------------
double CahnHilliard_Elas_FeFv_Tri3::giveDistanceToMidpointOf( std::vector<Node*>& face
                                                            , RealVector&         coor )
{
    RealVector coor0, coor1, dx;
    coor0 = analysisModel().domainManager().giveCoordinatesOf( face[ 0 ] );
    coor1 = analysisModel().domainManager().giveCoordinatesOf( face[ 1 ] );
    dx = 0.5 * ( coor0 + coor1 ) - coor;

    return std::sqrt( dx.dot( dx ) );
}
// ----------------------------------------------------------------------------
RealMatrix CahnHilliard_Elas_FeFv_Tri3::giveJacobianMatrixAt( Cell* targetCell )
{
    std::vector<Node*> node = analysisModel().domainManager().giveNodesOf( targetCell );
    RealMatrix coorMat( 3,2 );

    for ( int i = 0; i < 3; i++ )
    {
        RealVector nodeCoor = analysisModel().domainManager().giveCoordinatesOf( node[ i ] );
        coorMat( i,0 ) = nodeCoor( 0 );
        coorMat( i,1 ) = nodeCoor( 1 );
    }

    return _basisFunctionDerivatives * coorMat;
}
// ----------------------------------------------------------------------------
double CahnHilliard_Elas_FeFv_Tri3::giveLengthOf( std::vector<Node*>& face )
{
    RealVector coor0, coor1, dx;
    coor0 = analysisModel().domainManager().giveCoordinatesOf( face[ 0 ] );
    coor1 = analysisModel().domainManager().giveCoordinatesOf( face[ 1 ] );
    dx = coor1 - coor0;

    return std::sqrt( dx.dot( dx ) );
}
// ----------------------------------------------------------------------------
RealVector CahnHilliard_Elas_FeFv_Tri3::giveLocalDisplacementsAt( std::vector<Dof*>& dof, ValueType valType )
{
    // Displacements
    RealVector u( 6 );
    u( 0 ) = DofManager::giveValueOfPrimaryVariableAt( dof[ 0 ], valType );
    u( 1 ) = DofManager::giveValueOfPrimaryVariableAt( dof[ 1 ], valType );
    u( 2 ) = DofManager::giveValueOfPrimaryVariableAt( dof[ 2 ], valType );
    u( 3 ) = DofManager::giveValueOfPrimaryVariableAt( dof[ 3 ], valType );
    u( 4 ) = DofManager::giveValueOfPrimaryVariableAt( dof[ 4 ], valType );
    u( 5 ) = DofManager::giveValueOfPrimaryVariableAt( dof[ 5 ], valType );

    return u;
}
// ----------------------------------------------------------------------------
std::vector<Dof*> CahnHilliard_Elas_FeFv_Tri3::giveNodalDofsAt( Cell* targetCell )
{
    std::vector<Node*> node = analysisModel().domainManager().giveNodesOf( targetCell );

    std::vector<Dof*> dof( 6, nullptr );
    dof[ 0 ] = analysisModel().domainManager().giveNodalDof( _nodalDof[ 0 ], node[ 0 ] );
    dof[ 1 ] = analysisModel().domainManager().giveNodalDof( _nodalDof[ 1 ], node[ 0 ] );
    dof[ 2 ] = analysisModel().domainManager().giveNodalDof( _nodalDof[ 0 ], node[ 1 ] );
    dof[ 3 ] = analysisModel().domainManager().giveNodalDof( _nodalDof[ 1 ], node[ 1 ] );
    dof[ 4 ] = analysisModel().domainManager().giveNodalDof( _nodalDof[ 0 ], node[ 2 ] );
    dof[ 5 ] = analysisModel().domainManager().giveNodalDof( _nodalDof[ 1 ], node[ 2 ] );

    return dof;
}
// ----------------------------------------------------------------------------
double CahnHilliard_Elas_FeFv_Tri3::giveTransmissibilityCoefficientAt( std::vector<Node*>& face
        , Cell*               targetCell
        , Cell*               neighborCell )
{
    // A. Target cell
    // Coordinates of cell center
    std::vector<RealVector> ep = this->giveEvaluationPointsFor( targetCell );
    RealVector coor1 = ep[ 0 ];

    // length of face
    double length = giveLengthOf( face );

    // Distance from cell center to midpoint of face
    double d1 = giveDistanceToMidpointOf( face, coor1 );

    // B. Neighbor cell
    // Coordinates of cell center
    ep = this->giveEvaluationPointsFor( neighborCell );
    RealVector coor2 = ep[ 0 ];

    // Distance from cell center to midpoint of face
    double d2 = giveDistanceToMidpointOf( face, coor2 );

    return length / ( d1 + d2 );
}