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

// Input file declaration:
// 
// *NUMERICS <n_numerics>
//    ...
//    <tag> PhaseFieldFracture_FeFv_Tri3
//            NodalDof <dof1> <dof2>
//            CellDof <dof3>
//            Stage <stg>
//            NodalFieldOutput <n1>
//              <nodalField_1>   <dataIndex_1>
//              <nodalField_2>   <dataIndex_2>
//                   ...             ...
//              <nodalField_n1>  <dataIndex_n1>
//            CellFieldOutput <n2>
//              <cellField_1>    <cellField_1>
//              <cellField_2>    <cellField_2>
//                   ...             ...
//              <cellField_n2>   <cellField_n2>
//            CharacteristicLength <lc>
//            CriticalEnergyReleaseRate <Gc>
//            IrreversibilityThreshold <phi_c>
//
// ----------------------------------------------------------------------------
// Material set must contain exactly two materials:
//  material 1 --> solid density
//  material 2 --> damage constitutive law
//
// ----------------------------------------------------------------------------
// Field Data labels:
// 
//  s_xx  : reduced stress along xx
//  s_yy  : reduced stress along yy
//  s_zz  : reduced stress along zz
//  s_xy  : reduced stress along xy
//  ux_x  : du1/dx1
//  uy_x  : du2/dx1
//  ux_y  : du1/dx2
//  uy_y  : du2/dx2
//  g_xy  : total strain along xy
//    pf  : phase-field
// ene_s  : surface energy density
// ene_b  : elastic strain energy density accounting for damage
// cr_len : total crack length
//  pf_x  : x-component of phase-field gradient
//  pf_y  : y-component of phase-field gradient
//
// ----------------------------------------------------------------------------

#ifndef PHASEFIELDFRACTURE_FEFV_TRI3_HPP
#define PHASEFIELDFRACTURE_FEFV_TRI3_HPP

#include "Numerics/Numerics.hpp"
#include "Core/DofManager.hpp"
#include "Core/NumericsManager.hpp"
#include "BasisFunctions/Triangle_P1.hpp"

namespace broomstyx
{
    class MaterialStatus;
    
    class NumericsStatus_PhaseFieldFracture_FeFv_Tri3 final : public NumericsStatus
    {
        friend class PhaseFieldFracture_FeFv_Tri3;
        
    public:
        NumericsStatus_PhaseFieldFracture_FeFv_Tri3();
        virtual ~NumericsStatus_PhaseFieldFracture_FeFv_Tri3();
        
    private:
        double     _area;
        double     _phi;
        RealVector _strain;
        RealVector _stress;
        RealMatrix _gradU;
        double     _surfEgy;
        double     _bulkEgy;
        RealVector _gradPhi;
        RealMatrix _dPsi;
        
        bool   _hasPhsFldConstraint;
        bool   _hasPhsFldPrescribedOnFace[3];
        bool   _hasPhsFldGradientPrescribedOnFace[3];
        double _valueOnFace[3];
        bool   _hasNotComputedTransmissibilities;
        double _transmissibility[3];
        
        MaterialStatus* _materialStatus[2];
    };
    
    class PhaseFieldFracture_FeFv_Tri3 final : public Numerics
    {
    public:
        PhaseFieldFracture_FeFv_Tri3();
        virtual ~PhaseFieldFracture_FeFv_Tri3();

        void   deleteNumericsAt( Cell* targetCell ) override;
        void   finalizeDataAt( Cell* targetCell, const TimeData& time ) override;
        double giveCellFieldValueAt( Cell* targetCell, int fieldNum ) override;
        
        RealVector 
            giveCellNodeFieldValuesAt( Cell* targetCell, int fieldNum ) override;
        
        std::vector<RealVector> 
            giveEvaluationPointsFor( Cell* targetCell ) override;
        
        std::tuple< RealVector,RealVector > 
            giveFieldOutputAt( Cell* targetCell, const std::string& fieldTag ) override;
        
        RealVector
            giveNumericsParameter( const std::string& paramTag ) override;
        
        std::tuple< std::vector<Dof*>
                  , std::vector<Dof*>
                  , RealVector >
            giveStaticCoefficientMatrixAt( Cell*           targetCell
                                         , int             stage
                                         , int             subsys
                                         , const TimeData& time ) override;
        
        std::tuple< std::vector<Dof*>, RealVector >
            giveStaticLeftHandSideAt( Cell*           targetCell
                                    , int             stage
                                    , int             subsys
                                    , const TimeData& time ) override;
        
        std::tuple< std::vector<Dof*>, RealVector >
            giveStaticRightHandSideAt( Cell*                    targetCell
                                     , int                      stage
                                     , int                      subsys
                                     , const BoundaryCondition& bndCond
                                     , const TimeData&          time ) override;

        std::tuple< std::vector<Dof*>, RealVector >
            giveStaticRightHandSideAt( Cell*                 targetCell
                                     , int                   stage
                                     , int                   subsys
                                     , const FieldCondition& fldCond
                                     , const TimeData&       time ) override;
        
        void imposeConstraintAt( Cell*                    targetCell
                               , int                      stage
                               , const BoundaryCondition& bndCond
                               , const TimeData&          time ) override;
        
        void initializeMaterialsAt( Cell* targetCell ) override;
        void initializeNumericsAt( Cell* targetCell ) override;
        void performPostprocessingAt( Cell* targetCell, std::string tag ) override;
        void printPostIterationMessage( int stage ) override;
        void readAdditionalDataFrom( FILE* fp ) override;
        void removeConstraintsOn( Cell* targetCell ) override;
        void setDofStagesAt( Cell* targetCell ) override;

    private:
        Triangle_P1 _basisFunction;
        RealVector  _basisFunctionValues;
        RealMatrix  _basisFunctionDerivatives;
        
        RealVector  _gpNatCoor;
        double      _wt;
        
        double _l;
        double _Gc;
        
        NumericsStatus_PhaseFieldFracture_FeFv_Tri3*
                   getNumericsStatusAt( Cell* targetCell );
        RealMatrix giveBmatAt( Cell* targetCell );
        std::vector< std::vector<Node*> >
                   giveFaceNodesOf( Cell* targetCell );
        double     giveDistanceToMidpointOf( std::vector<Node*>& face, RealVector& coor);
        RealMatrix giveJacobianMatrixAt( Cell* targetCell, const RealVector& natCoor );
        double     giveLengthOf( std::vector<Node*>& face );
        RealVector giveLocalDisplacementsAt( std::vector<Dof*>& dof, ValueType valType );
        std::vector<Dof*> 
                   giveNodalDofsAt( Cell* targetCell );
        RealVector giveOutwardUnitNormalOf( std::vector<Node*>& face );
        double     giveTransmissibilityCoefficientAt( std::vector<Node*>& face
                                                    , Cell*               targetCell
                                                    , Cell*               neighborCell );
    };
}

#endif /* PHASEFIELDFRACTURE_FEFV_TRI3_HPP */
