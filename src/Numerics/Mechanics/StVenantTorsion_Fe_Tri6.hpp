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

#ifndef STVENANTTORSION_FE_TRI6_HPP
#define STVENANTTORSION_FE_TRI6_HPP

#include "../Numerics.hpp"
#include "../../BasisFunctions/Line_P2.hpp"
#include "../../BasisFunctions/Triangle_P2.hpp"
#include "../../IntegrationRules/Legendre_1D.hpp"
#include "../../IntegrationRules/Legendre_2D_Tri.hpp"

namespace broomstyx
{
    class Dof;
    class EvalPoint;
    class MaterialStatus;

    // Cell numerics status
    class CellNumericsStatus_StVenantTorsion_Fe_Tri6 final : public NumericsStatus
    {
        friend class StVenantTorsion_Fe_Tri6;

    public:
        CellNumericsStatus_StVenantTorsion_Fe_Tri6( int nGaussPts );
        virtual ~CellNumericsStatus_StVenantTorsion_Fe_Tri6();

    private:
        RealVector _nodalStrain_zx;
        RealVector _nodalStrain_zy;
        RealVector _nodalStrain_mag;
        RealVector _nodalStress_zx;
        RealVector _nodalStress_zy;
        RealVector _nodalStress_mag;
        RealVector _zDisp;
        double     _torqueContrib;
    };

    class StVenantTorsion_Fe_Tri6 final : public Numerics
    {
    public:
        StVenantTorsion_Fe_Tri6();
        virtual ~StVenantTorsion_Fe_Tri6();

        void deleteNumericsAt( Cell* targetCell ) override;
        void finalizeDataAt( Cell* targetCell ) override;

        RealVector giveCellNodeFieldValuesAt( Cell* targetCell, int fieldNum ) override;

        std::tuple< std::vector<Dof*>
                  , std::vector<Dof*>
                  , RealVector >
             giveStaticCoefficientMatrixAt( Cell*           targetCell
                                          , int             stage
                                          , int             subsys
                                          , const TimeData& time ) override;

        std::tuple< std::vector<Dof*>, RealVector >
             giveStaticRightHandSideAt( Cell*                    targetCell
                                      , int                      stage
                                      , int                      subsys
                                      , const BoundaryCondition& bndCond
                                      , const TimeData&          time ) override;

        void initializeNumericsAt( Cell* targetCell ) override;
        void readAdditionalDataFrom( FILE* fp ) override;
        void setDofStagesAt( Cell* targetCell ) override;

    private:
        Triangle_P2 _basisFunction;
        Line_P2     _edgeBasisFunction;

        Legendre_2D_Tri _integrationRule;
        Legendre_1D     _edgeIntegrationRule;
        Legendre_2D_Tri _torqueIntegrationRule;

        double _beta;

        Dof* _az;
        Dof* _kx;
        Dof* _ky;
        int  _dofGrp;

        CellNumericsStatus_StVenantTorsion_Fe_Tri6*
                   getNumericsStatusAt( Cell* targetCell );
        RealMatrix giveBmatAt( Cell* targetCell, const RealVector& natCoor );
        RealVector giveNmatAt( const RealVector& natCoor );
        RealMatrix giveJacobianMatrixAt( Cell* targetCell, const RealVector& natCoor );
        RealVector giveLocalValuesAt( Cell* targetCell, ValueType valType );
        std::vector<Dof*> giveNodalDofsAt( Cell* targetCell );
    };
}

#endif /* STVENANTTORSION_FE_TRI6_HPP */
