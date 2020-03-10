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

#ifndef MATERIAL_HPP
#define MATERIAL_HPP

#include <cstdio>
#include "Util/RealMatrix.hpp"
#include "Util/RealVector.hpp"

namespace broomstyx
{
    class Cell;
    
    class MaterialStatus
    {
    public:
        MaterialStatus();
        virtual ~MaterialStatus();
    };

    class Material
    {
    public:
        Material();
        virtual ~Material();
        
        // Disable copy constructor and assignment operator
        Material( const Material& ) = delete;
        Material& operator=( const Material& ) = delete;

        virtual MaterialStatus* createMaterialStatus();
        virtual void destroy( MaterialStatus*& matStatus );
        virtual void initialize();
        virtual void readParamatersFrom( FILE* fp );
        virtual void updateStatusFrom( const RealVector& conState, MaterialStatus* matStatus );
        virtual void updateStatusFrom( const RealVector& conState, MaterialStatus* matStatus, const std::string& label );
        
        // Error-generating virtual methods
        virtual double
            givePotentialFrom( const RealVector& conState, const MaterialStatus* matStatus );
        
        virtual RealVector 
            giveForceFrom( const RealVector& conState, const MaterialStatus* matStatus );
        
        virtual RealVector 
            giveForceFrom( const RealVector& conState, const MaterialStatus* matStatus, const std::string& label );
        
        virtual RealMatrix
            giveModulusFrom( const RealVector& conState, const MaterialStatus* matStatus );
        
        virtual RealMatrix
            giveModulusFrom( const RealVector& conState, const MaterialStatus* matStatus, const std::string& label );
        
        virtual double 
            giveMaterialVariable( const std::string& label, const MaterialStatus* matStatus );
        
        virtual double
            giveParameter( const std::string& label );
        
    protected:
        std::string _name;
        
        void error_unimplemented( const std::string& method );
    };
}

#endif /* MATERIAL_HPP */