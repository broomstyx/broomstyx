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
        virtual ~MaterialStatus_MieheDamageModel();
        
    private:
        MaterialStatus* _materialStatus;
        RealVector      _prinStrain;
        RealMatrix      _M;
        double          _historyField;
        double          _drivingForce;
    };
    
    class MieheDamageModel final : public Material
    {
    public:
        MieheDamageModel();
        virtual ~MieheDamageModel();
        
        MaterialStatus* createMaterialStatus() override;
        void       destroy( MaterialStatus*& matStatus ) override;
        double     givePotentialFrom( const RealVector& conState, const MaterialStatus* matStatus ) override;
        RealVector giveForceFrom( const RealVector& conState, const MaterialStatus* matStatus, const std::string& label ) override;
        double     giveMaterialVariable( const std::string& label, const MaterialStatus* matStatus ) override;
        RealMatrix giveModulusFrom( const RealVector& conState, const MaterialStatus* matStatus, const std::string& label ) override;
        void       readParamatersFrom( FILE* fp ) override;
        void       updateStatusFrom( const RealVector& conState, MaterialStatus* matStatus ) override;
        void       updateStatusFrom( const RealVector& conState, MaterialStatus* matStatus, const std::string& label ) override;
        
    private:
        int       _analysisMode;
        double    _E;
        double    _nu;
        double    _lambda;
        double    _G;
        Material* _degradationFunction;
        double    _phiIrrev;
        
        MaterialStatus_MieheDamageModel* accessMaterialStatus( MaterialStatus* matStatus )
        {
            auto mst = dynamic_cast<MaterialStatus_MieheDamageModel*>(matStatus);
            if ( !mst )
                throw std::runtime_error("Error: Unable to access material status!\nSource: " + _name);
            
            return mst;
        }

        const MaterialStatus_MieheDamageModel* accessConstMaterialStatus( const MaterialStatus* matStatus )
        {
            auto mst = dynamic_cast<const MaterialStatus_MieheDamageModel*>(matStatus);
            if ( !mst )
                throw std::runtime_error("Error: Unable to access material status!\nSource: " + _name);
            
            return mst;
        }
        
        std::tuple< RealVector, RealVector > getStrainAndPhaseFieldFrom( const RealVector& conState )
        {
            RealVector strain;
            RealVector phi(1);

            if ( _analysisMode == 2 /* Plane strain */)
            {
                strain = {conState(0), conState(1), conState(2), conState(3)};
                phi(0) = conState(4);
            }
            else
                throw std::runtime_error("Error: Cannot resolve constitutive state for current analysis mode!\nSource: " + _name);

            return std::make_tuple(std::move(strain), std::move(phi));
        }
    };
}

#endif /* MIEHEDAMAGEMODEL_HPP */