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

#ifndef OBJECTFACTORY_HPP
#define OBJECTFACTORY_HPP

#include <map>

#define registerBroomstyxObject(baseClass, derivedClass) static bool dummy_ ## derivedClass __attribute__((unused)) = objectFactory().register ## baseClass( # derivedClass, createDerivedObject<baseClass,derivedClass> );

namespace broomstyx
{
    class ConvergenceCriterion;
    class LinearSolver;
    class Material;
    class MeshReader;
    class Numerics;
    class OutputQuantity;
    class OutputWriter;
    class SolutionMethod;
    class SparseMatrix;
    class UserFunction;
    
    template<typename Base, typename Derived> Base* createDerivedObject()
    {
        return new Derived();
    }

    class ObjectFactory
    {
    public:
        ObjectFactory();
        virtual ~ObjectFactory();
        
        // Disable copy constructor and assignment operator
        ObjectFactory( const ObjectFactory& ) = delete;
        ObjectFactory& operator=( const ObjectFactory& ) = delete;
        
        bool hasError();
        
        ConvergenceCriterion* instantiateConvergenceCriterion( std::string name );
        LinearSolver*         instantiateLinearSolver( std::string name );
        Material*             instantiateMaterial( std::string name );
        MeshReader*           instantiateMeshReader( std::string name );
        Numerics*             instantiateNumerics( std::string name );
        OutputQuantity*       instantiateOutputQuantity( std::string name );
        OutputWriter*         instantiateOutputWriter( std::string name );
        SolutionMethod*       instantiateSolutionMethod( std::string name );
        SparseMatrix*         instantiateSparseMatrix( std::string name );
        UserFunction*         instantiateUserFunction( std::string name );
        
        bool registerConvergenceCriterion( std::string name, ConvergenceCriterion* (*)() );
        bool registerLinearSolver( std::string name, LinearSolver* (*)() );
        bool registerMaterial( std::string name, Material* (*)() );
        bool registerMeshReader( std::string name, MeshReader* (*)() );
        bool registerNumerics( std::string name, Numerics* (*)() );
        bool registerOutputQuantity( std::string name, OutputQuantity* (*)() );
        bool registerOutputWriter( std::string name, OutputWriter* (*)() );
        bool registerSolutionMethod( std::string name, SolutionMethod* (*)() );
        bool registerSparseMatrix( std::string name, SparseMatrix* (*)() );
        bool registerUserFunction( std::string name, UserFunction* (*)() );

    private:
        std::map< std::string, ConvergenceCriterion* (*)() > _convergenceCriterion;
        std::map< std::string, LinearSolver* (*)() >         _linearSolver;
        std::map< std::string, Material* (*)() >             _material;
        std::map< std::string, MeshReader* (*)() >           _meshReader;
        std::map< std::string, Numerics* (*)() >             _numerics;
        std::map< std::string, OutputQuantity* (*)() >       _outputQuantity;
        std::map< std::string, OutputWriter* (*)() >         _outputWriter;
        std::map< std::string, SolutionMethod* (*)() >       _solutionMethod;
        std::map< std::string, SparseMatrix* (*)() >         _sparseMatrix;
        std::map< std::string, UserFunction* (*)() >         _userFunction;
        
        bool _errorInRegistration;
    };

    ObjectFactory& objectFactory();
}

#endif /* OBJECTFACTORY_HPP */