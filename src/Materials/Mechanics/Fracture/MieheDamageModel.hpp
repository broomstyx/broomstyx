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

#ifndef MIEHEDAMAGEMODEL_HPP
#define MIEHEDAMAGEMODEL_HPP

#include "Materials/Material.hpp"
#include <tuple>

namespace broomstyx
{
    class MaterialStatus_MieheDamageModel final : public MaterialStatus
    {
        friend class MieheDamageModel;
        
    public:
        MaterialStatus_MieheDamageModel();
        ~MaterialStatus_MieheDamageModel() override = default;
        
    private:
        MaterialStatus* _materialStatus;
        RealVector      _prinStrain;
        RealMatrix      _M;
        double          _drivingForce;
        double          _eta;
        double          _eps3;
    };
    
    class MieheDamageModel final : public Material
    {
    public:
        MieheDamageModel();
        ~MieheDamageModel() override = default;
        
        MaterialStatus* createMaterialStatus() override;
        void       destroy( MaterialStatus*& matStatus ) override;
        double     givePotentialFrom( const RealVector& conState, const MaterialStatus* matStatus ) override;
        RealVector giveForceFrom( const RealVector& conState, const MaterialStatus* matStatus, const std::string& label ) override;
        double     giveMaterialVariable( const std::string& label, const MaterialStatus* matStatus ) override;
        RealMatrix giveModulusFrom( const RealVector& conState, const MaterialStatus* matStatus, const std::string& label ) override;
        void       readParamatersFrom( FILE* fp ) override;
        void       updateStatusFrom( const RealVector& conState, MaterialStatus* matStatus ) override;

    private:
        enum AnalysisMode
        {
            Unset,
            PlaneStress,
            PlaneStrain
        };

        AnalysisMode _analysisMode;
        double    _E;
        double    _nu;
        double    _lambda;
        double    _G;
        Material* _degradationFunction;
        int       _maxIdx;

        MaterialStatus_MieheDamageModel* accessMaterialStatus( MaterialStatus* matStatus );
        const MaterialStatus_MieheDamageModel* accessConstMaterialStatus( const MaterialStatus* matStatus );
        std::tuple< RealVector, RealVector > getStrainAndPhaseFieldFrom( const RealVector& conState );
    };
}

#endif /* MIEHEDAMAGEMODEL_HPP */