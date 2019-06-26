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

/* --------------------------------------------------------------------------
 * Input file declaration:
 * 
 * *NUMERICS <n_numerics>
 *    ...
 *    <tag> PlaneStrain_Fe_Tri3
 *            NodalDof <dof1> <dof2>
 *            Stage <stg>
 *            NodalFieldOutput <n1>
 *                   ...             ...
 *                   ...             ...
 *            CellFieldOutput <n2>
 *                   ...             ...
 *                   ...             ...
 * 
 * Field data tags:
 * 
 *  s_xx : stress along xx
 *  s_yy : stress along yy
 *  s_zz : stress along zz
 *  s_xy : stress along xy
 *  ux_x : du1/dx1
 *  uy_x : du2/dx1
 *  ux_y : du1/dx2
 *  uy_y : du2/dx2
 *  g_xy : total engineering strain along xy
 * ep_xx : plastic strain along xx
 * ep_yy : plastic strain along yy
 * ep_zz : plastic strain along zz
 * ep_xy : plastic strain along xy
 *   ene : elastic strain energy
 * 
 * Material set:
 *  material 1 --> solid density
 *  material 2 --> constitutive law for the solid (elasticity)
 * 
 ****************************************************************************/

#ifndef PLANESTRAIN_FE_TRI3_HPP
#define	PLANESTRAIN_FE_TRI3_HPP

#include "../Numerics.hpp"
#include "../../Core/DofManager.hpp"
#include "../../Core/NumericsManager.hpp"
#include "../../BasisFunctions/Triangle_P1.hpp"

namespace broomstyx
{
    class MaterialStatus;
        
    class NumericsStatus_PlaneStrain_Fe_Tri3 final : public NumericsStatus
    {
        friend class PlaneStrain_Fe_Tri3;
        
    public:
        NumericsStatus_PlaneStrain_Fe_Tri3();
        virtual ~NumericsStatus_PlaneStrain_Fe_Tri3();
        
    private:
        RealVector _strain;
        RealVector _stress;
        RealMatrix _gradU;
        RealMatrix _JmatInv;
        double     _Jdet;
        MaterialStatus* _materialStatus[2];
    };
    
    class PlaneStrain_Fe_Tri3 final : public Numerics 
    {
    public:
        PlaneStrain_Fe_Tri3();
        virtual ~PlaneStrain_Fe_Tri3();

        void   deleteNumericsAt( Cell* targetCell ) override;
        void   finalizeDataAt( Cell* targetCell ) override;
        double giveCellFieldValueAt( Cell* targetCell, int fieldNum ) override;
        
        RealVector 
            giveCellNodeFieldValuesAt( Cell* targetCell, int fieldNum ) override;
        
        std::vector<RealVector> 
            giveEvaluationPointsFor( Cell* targetCell ) override;
        
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
        
        void imposeConstraintAt( Cell*                    targetCell
                               , int                      stage
                               , const BoundaryCondition& bndCond
                               , const TimeData&          time ) override;
        
        void initializeMaterialsAt( Cell* targetCell ) override;
        void initializeNumericsAt( Cell* targetCell ) override;
        void setDofStagesAt( Cell* targetCell ) override;

    private:
        Triangle_P1 _basisFunction;
        RealVector  _basisFunctionValues;
        RealMatrix  _basisFunctionDerivatives;
        
        RealVector  _gpNatCoor;
        double      _wt;
        
        NumericsStatus_PlaneStrain_Fe_Tri3*
                   getNumericsStatusAt( Cell* targetCell );
        RealMatrix giveBmatAt( Cell* targetCell );
        RealMatrix giveGradBmatAt( Cell* targetCell );
        RealMatrix giveJacobianMatrixAt( Cell* targetCell, const RealVector& natCoor );
        RealVector giveLocalDisplacementsAt( std::vector<Dof*>& dof, ValueType valType );
        std::vector<Dof*> giveNodalDofsAt( Cell* targetCell );
    };
}

#endif	/* PLANESTRAIN_FE_TRI3_HPP */