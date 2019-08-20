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
 *  s_yz : stress along yz
 *  s_xz : stress along xz
 *  s_xy : stress along xy
 *  ux_x : du1/dx1
 *  uy_x : du2/dx1
 *  uz_x : du3/dx1
 *  ux_y : du1/dx2
 *  uy_y : du2/dx2
 *  uz_y : du3/dx2
 *  ux_z : du1/dx3
 *  uy_z : du2/dx3
 *  uz_z : du3/dx3
 *  g_yz : total engineering strain along yz
 *  g_xz : total engineering strain along xz
 *  g_xy : total engineering strain along xy
 *   ene : elastic strain energy
 * 
 * Material set:
 *  material 1 --> solid density
 *  material 2 --> constitutive law for the solid (elasticity)
 * 
 ****************************************************************************/

#ifndef MECH_FE_TET4_HPP
#define	MECH_FE_TET4_HPP

#include "../Numerics.hpp"
#include "../../Core/DofManager.hpp"
#include "../../Core/NumericsManager.hpp"
#include "../../BasisFunctions/Tetrahedron_P1.hpp"

namespace broomstyx
{
    class MaterialStatus;
        
    class NumericsStatus_Mech_Fe_Tet4 final : public NumericsStatus
    {
        friend class Mech_Fe_Tet4;
        
    public:
        NumericsStatus_Mech_Fe_Tet4();
        virtual ~NumericsStatus_Mech_Fe_Tet4();
        
    private:
        RealVector _strain;
        RealVector _stress;
        RealMatrix _gradU;
        RealMatrix _JmatInv;
        double     _Jdet;
        MaterialStatus* _materialStatus[2];
    };
    
    class Mech_Fe_Tet4 final : public Numerics 
    {
    public:
        Mech_Fe_Tet4();
        virtual ~Mech_Fe_Tet4();

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
        Tetrahedron_P1 _basisFunction;
        RealVector  _basisFunctionValues;
        RealMatrix  _basisFunctionDerivatives;
        
        RealVector  _gpNatCoor;
        double      _wt;
        
        NumericsStatus_Mech_Fe_Tet4*
                   getNumericsStatusAt( Cell* targetCell );
        RealMatrix giveBmatAt( Cell* targetCell );
        RealMatrix giveGradBmatAt( Cell* targetCell );
        RealMatrix giveJacobianMatrixAt( Cell* targetCell, const RealVector& natCoor );
        RealVector giveLocalDisplacementsAt( std::vector<Dof*>& dof, ValueType valType );
        std::vector<Dof*> giveNodalDofsAt( Cell* targetCell );
    };
}

#endif	/* MECH_FE_TET4_HPP */