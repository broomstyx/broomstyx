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

#ifndef AMORDAMAGEMODEL_HPP
#define AMORDAMAGEMODEL_HPP

#include "Materials/Material.hpp"
#include <tuple>

namespace broomstyx
{
    class MaterialStatus_AmorDamageModel final : public MaterialStatus
    {
        friend class AmorDamageModel;
        
    public:
        MaterialStatus_AmorDamageModel();
        ~MaterialStatus_AmorDamageModel() override = default;
        
    private:
        MaterialStatus* _materialStatus[ 2 ];
        double _volElasticEnergy;
        double _devElasticEnergy;
    };
    
    class AmorDamageModel final : public Material
    {
    public:
        AmorDamageModel();
        ~AmorDamageModel() override = default;
        
        MaterialStatus* createMaterialStatus() override;
        void       destroy( MaterialStatus*& matStatus ) override;
        double     givePotentialFrom( const RealVector& conState, const MaterialStatus* matStatus ) override;
        double     givePotentialFrom( const RealVector& conState, const MaterialStatus* matStatus, const std::string& label ) override;
        RealVector giveForceFrom( const RealVector& conState, const MaterialStatus* matStatus, const std::string& label ) override;
        RealMatrix giveModulusFrom( const RealVector& conState, const MaterialStatus* matStatus, const std::string& label ) override;
        void       readParamatersFrom( FILE* fp ) override;
        void       updateStatusFrom( const RealVector& conState, MaterialStatus* matStatus ) override;
        void       updateStatusFrom( const RealVector& conState, MaterialStatus* matStatus, const std::string& label ) override;
        
    private:
        enum AnalysisMode
        {
            TwoDimensional,
            PlaneStress,
            PlaneStrain,
            Axisymmetric,
            ThreeDimensional
        };

        AnalysisMode _analysisMode;
        Material*    _elasticityModel;
        Material*    _degradationFunction;
        double       _kT;
        double       _kC;
        RealMatrix   _I;
        RealMatrix   _Pvol;

        MaterialStatus_AmorDamageModel* 
                   accessMaterialStatus( MaterialStatus* matStatus );
        const MaterialStatus_AmorDamageModel* 
                   accessConstMaterialStatus( const MaterialStatus* matStatus );

        std::tuple< RealVector, RealVector > retrieveStrainAndPhasefieldFrom( const RealVector& conState );
    };
}

#endif /* AMORDAMAGEMODEL_HPP */