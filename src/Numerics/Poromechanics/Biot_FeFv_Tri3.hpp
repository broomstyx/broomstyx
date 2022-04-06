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

#ifndef BIOT_FEFV_TRI3_HPP
#define	BIOT_FEFV_TRI3_HPP

#include "Numerics/Numerics.hpp"
#include "BasisFunctions/Triangle_P1.hpp"

namespace broomstyx
{
    class MaterialStatus;
    
    class NumericsStatus_Biot_FeFv_Tri3 : public NumericsStatus
    {
        friend class Biot_FeFv_Tri3;
        
    public:
        NumericsStatus_Biot_FeFv_Tri3();
        virtual ~NumericsStatus_Biot_FeFv_Tri3();
        
    private:
        RealVector _strain;
        RealVector _stress;
        RealMatrix _gradU;
        RealMatrix _dPsi;
        
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
        
        MaterialStatus* _materialStatus[3];
    };
    
    class Biot_FeFv_Tri3 final : public Numerics
    {
    public:
        Biot_FeFv_Tri3();
        virtual ~Biot_FeFv_Tri3();

        void   deleteNumericsAt( Cell* targetCell ) override;
        void   finalizeDataAt( Cell* targetCell, const TimeData& time ) override;
        double giveCellFieldValueAt( Cell* targetCell, int fieldNum ) override;
        
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
        
        std::tuple< std::vector<Dof*>, RealVector >
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
        
        Triangle_P1 _basisFunction;
        RealVector  _basisFunctionValues;
        RealMatrix  _basisFunctionDerivatives;
        
        RealVector  _gpNatCoor;
        double      _wt;
        
        NumericsStatus_Biot_FeFv_Tri3*
                   getNumericsStatusAt( Cell* targetCell );
        std::tuple< RealMatrix, RealVector >
                   giveBmatAt( Cell* targetCell );
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

#endif	/* HPEF3VBIOT_HPP */
