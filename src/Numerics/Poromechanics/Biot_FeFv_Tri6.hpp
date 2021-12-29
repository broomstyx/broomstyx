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

#ifndef BIOT_FEFV_TRI6_HPP
#define	BIOT_FEFV_TRI6_HPP

#include "Numerics/Numerics.hpp"

namespace broomstyx
{
    class EvalPoint;
    class IntegrationRule;
    class MaterialStatus;
    class ScalarBasisFunction;
    
    // Cell numerics status
    class CellNumericsStatus_Biot_FeFv_Tri6 : public NumericsStatus
    {
        friend class Biot_FeFv_Tri6;
        
    public:
        CellNumericsStatus_Biot_FeFv_Tri6( int nGaussPts );
        virtual ~CellNumericsStatus_Biot_FeFv_Tri6();

    private:
        std::vector<EvalPoint> _gp;
        int _nGaussPts;

        RealVector _cellCenter;

        double _area;
        double _head;
        double _centerFlux[2];
        double _headOnFace[3];
        double _fluxOnFace[3];
        bool   _fluxIsPrescribedOnFace[3];
        bool   _headIsPrescribedOnFace[3];
        bool   _hasHeadConstraint;
        bool   _hasNotComputedTransmissibilities;
        double _transmissibility[3];

        MaterialStatus* _materialStatus;
    };

    // Integration point numerics status
    class EvalPtNumericsStatus_Biot_FeFv_Tri6 final : public NumericsStatus
    {  
        friend class Biot_FeFv_Tri6;

    public:
        EvalPtNumericsStatus_Biot_FeFv_Tri6();
        ~EvalPtNumericsStatus_Biot_FeFv_Tri6();
        
    private:
        RealVector _strain;
        RealVector _stress;
        RealMatrix _gradU;
        RealMatrix _dPsi;
        
        MaterialStatus* _materialStatus[2];
    };
    
    class Biot_FeFv_Tri6 final : public Numerics
    {
    public:
        Biot_FeFv_Tri6();
        virtual ~Biot_FeFv_Tri6();

        void   deleteNumericsAt( Cell* targetCell ) override;
        void   finalizeDataAt( Cell* targetCell ) override;
        double giveCellFieldValueAt( Cell* targetCell, int fieldNum ) override;
        RealVector giveCellNodeFieldValuesAt( Cell* targetCell, int fieldNum );
        
        std::vector<RealVector> giveEvaluationPointsFor( Cell* targetCell ) override;
        
        std::tuple< RealVector,RealVector > 
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
        
        std::tuple< std::vector<Dof*>, std::vector<Dof*>, RealVector >
            giveTransientCoefficientMatrixAt( Cell*           targetCell
                                            , int             stage
                                            , int             subsys
                                            , const TimeData& time ) override;
        
        virtual std::tuple< std::vector<Dof*>, RealVector >
            giveTransientLeftHandSideAt( Cell*           targetCell
                                       , int             stage
                                       , int             subsys
                                       , const TimeData& time
                                       , ValueType       valType ) override;
        
        void imposeConstraintAt( Cell*                    targetCell
                               , int                      stage
                               , const BoundaryCondition& bndCond
                               , const TimeData&          time ) override;
        
        void imposeInitialConditionAt( Cell*                   targetCell
                                     , const InitialCondition& initCond ) override;
        
        void initializeMaterialsAt( Cell* targetCell ) override;
        void initializeNumericsAt( Cell* targetCell ) override;
        void readAdditionalDataFrom( FILE* fp ) override;
        void removeConstraintsOn( Cell* targetCell ) override;
        void setDofStagesAt( Cell* targetCell ) override;
        
    private:
        double _alpha;
        double _S;
        double _gAccel;
        double _rhoF;
        double _mu;
        
        double _sgn;
        int    _vrtIndex;

        ScalarBasisFunction*   _basisFunction;
        ScalarBasisFunction*   _edgeBasisFunction;
        
        IntegrationRule* _integrationRule;
        IntegrationRule* _edgeIntegrationRule;

        RealMatrix _extrapolationMatrix;
        
        void formExtrapolationMatrix();
        CellNumericsStatus_Biot_FeFv_Tri6*
                   getNumericsStatusAt( Cell* targetCell );
        EvalPtNumericsStatus_Biot_FeFv_Tri6*
                   getNumericsStatusAt( EvalPoint& gp );
        std::tuple< RealMatrix, RealVector >
                   giveBmatAt( Cell* targetCell, const RealVector& natCoor );
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

#endif	/* BIOT_FEFV_TRI6_HPP */