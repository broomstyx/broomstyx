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

/****************************************************************************
 * 
 * Element type 'DarcyFlow_2D_1Phase_Fv_Tri' is a 2-dimensional continuum 
 * element with 3 faces/vertices.
 * 
 * The element implements the control volume equation for Darcy flow using
 * a cell-centered finite volume scheme and the two-point flux approximation.
 * Each cell has one pressure degree of freedom located at the cell centroid.
 * 
 * --------------------------------------------------------------------------
 * Input file declaration:
 * 
 * *NUMERICS <n_numerics>
 *    ...
 *    <tag> DarcyFlow_2D_1Phase_Fv_Tri
 *            CellDof <dof1>
 *            Stage <stg>
 *            Subsystem <subsys>
 *            CellFieldOutput <n2>
 *                   ...             ...
 *                   ...             ...
 *            
 * 
 * Field data labels:
 *  h : hydraulic head at cell center
 *  p : pressure at cell center
 *  q_x : x-component of reconstructed flux at cell center
 *  q_y : y-component of reconstructed flux at cell center
 *
 * Material set:
 *  material 1 --> permeability of the porous medium
 *  
 */

#ifndef DARCYFLOW_2D_1PHASE_FV_TRI_HPP
#define	DARCYFLOW_2D_1PHASE_FV_TRI_HPP

#include "../Numerics.hpp"

namespace broomstyx
{
    class MaterialStatus;
    
    class NumericsStatus_DarcyFlow_2D_1Phase_Fv_Tri : public NumericsStatus
    {
        friend class DarcyFlow_2D_1Phase_Fv_Tri;
        
    public:
        NumericsStatus_DarcyFlow_2D_1Phase_Fv_Tri();
        virtual ~NumericsStatus_DarcyFlow_2D_1Phase_Fv_Tri();
        
    private:
        double _centerFlux[2];
        double _headOnFace[3];
        double _fluxOnFace[3];
        bool _fluxIsPrescribedOnFace[3];
        bool _headIsPrescribedOnFace[3];
        bool _hasHeadConstraint;
        MaterialStatus* _materialStatus;
    };
    
    class DarcyFlow_2D_1Phase_Fv_Tri : public Numerics
    {
    public:
        DarcyFlow_2D_1Phase_Fv_Tri();
        virtual ~DarcyFlow_2D_1Phase_Fv_Tri();

        void   deleteNumericsAt( Cell* targetCell ) override;
        void   finalizeDataAt( Cell* targetCell ) override;
        double giveCellFieldValueAt( Cell* targetCell, int fieldNum ) override;
        
        std::vector<RealVector> giveEvaluationPointsFor( Cell* targetCell ) override;
        
        std::tuple< RealVector, RealVector >
            giveFieldOutputAt( Cell* targetCell, const std::string& fieldTag ) override;
        
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
        void
            imposeConstraintAt( Cell*                    targetCell
                              , int                      stage
                              , const BoundaryCondition& bndCond
                              , const TimeData&          time ) override;
        
        void initializeMaterialsAt( Cell* targetCell ) override;
        void initializeNumericsAt( Cell* targetCell ) override;
        void readAdditionalDataFrom( FILE* fp ) override;        
        void removeConstraintsOn( Cell* targetCell ) override;
        void setDofStagesAt( Cell* targetCell ) override;
        
    private:
        double _S;
        double _rhoF;
        double _mu;
        double _gAccel;
        double _sgn;
        int    _vrtIndex;

        NumericsStatus_DarcyFlow_2D_1Phase_Fv_Tri*
                   getNumericsStatusAt( Cell* targetCell );        
        double     giveAreaOf( Cell* targetCell );
        double     giveDistanceToMidpointOf( const std::vector<Node*>& face, const RealVector& coor );
        std::vector< std::vector<Node*> > 
                   giveFaceNodesOf( Cell* targetCell );
        double     giveLengthOf( const std::vector<Node*>& face );
        RealVector giveOutwardUnitNormalOf( const std::vector<Node*>& face );
        double     giveTransmissibilityCoefficientAt( const std::vector<Node*>& face
                                                , Cell* targetCell
                                                , Cell* neighborCell );
        double     giveVerticalCoordinateAt( Cell* targetCell );
    };
}

#endif	/* DARCYFLOW_2D_1PHASE_FV_TRI_HPP */