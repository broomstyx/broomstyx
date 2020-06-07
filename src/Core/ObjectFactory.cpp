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

#include "ObjectFactory.hpp"

#define VERBOSE_REGISTRATION false

using namespace broomstyx;

ObjectFactory::ObjectFactory() : _errorInRegistration(false) {}

ObjectFactory::~ObjectFactory() {}

// Object Instantiation
// ----------------------------------------------------------------------------
bool ObjectFactory::hasError()
{
    return _errorInRegistration;
}
// ----------------------------------------------------------------------------
ConvergenceCriterion* ObjectFactory::instantiateConvergenceCriterion( std::string name )
{
    auto it = _convergenceCriterion.find(name);
    if ( it != _convergenceCriterion.end() )
        return (it->second)();
    else
        throw std::runtime_error("LinearSolver '" + name + "' has not been registered.\n");
}
// ----------------------------------------------------------------------------
LinearSolver* ObjectFactory::instantiateLinearSolver( std::string name )
{
    auto it = _linearSolver.find(name);
    if ( it != _linearSolver.end() )
        return (it->second)();
    else
        throw std::runtime_error("LinearSolver '" + name + "' has not been registered.\n");
}
// ----------------------------------------------------------------------------
Material* ObjectFactory::instantiateMaterial( std::string name )
{
    auto it = _material.find(name);
    if ( it != _material.end() )
        return (it->second)();
    else
        throw std::runtime_error("Material '" + name + "' has not been registered.\n");
}
// ----------------------------------------------------------------------------
MeshReader* ObjectFactory::instantiateMeshReader( std::string name )
{
    auto it = _meshReader.find(name);
    if ( it != _meshReader.end() )
        return (it->second)();
    else
        throw std::runtime_error("MeshReader '" + name + "' has not been registered.\n");
}
// ----------------------------------------------------------------------------
Numerics* ObjectFactory::instantiateNumerics( std::string name )
{
    auto it = _numerics.find(name);
    if ( it != _numerics.end() )
        return (it->second)();
    else
        throw std::runtime_error("Numerics '" + name + "' has not been registered.\n");
}
// ----------------------------------------------------------------------------
OutputQuantity* ObjectFactory::instantiateOutputQuantity( std::string name )
{
    auto it = _outputQuantity.find(name);
    if ( it != _outputQuantity.end() )
        return (it->second)();
    else
        throw std::runtime_error("OutputQuantity '" + name + "' has not been registered.\n");
}
// ----------------------------------------------------------------------------
OutputWriter* ObjectFactory::instantiateOutputWriter( std::string name )
{
    auto it = _outputWriter.find(name);
    if ( it != _outputWriter.end() )
        return (it->second)();
    else
        throw std::runtime_error("OutputWriter '" + name + "' has not been registered.\n");
}
// ----------------------------------------------------------------------------
SolutionMethod* ObjectFactory::instantiateSolutionMethod( std::string name )
{
    auto it = _solutionMethod.find(name);
    if ( it != _solutionMethod.end() )
        return (it->second)();
    else
        throw std::runtime_error("SolutionMethod '" + name + "' has not been registered.\n");
}
// ----------------------------------------------------------------------------
SparseMatrix* ObjectFactory::instantiateSparseMatrix( std::string name )
{
    auto it = _sparseMatrix.find(name);
    if ( it != _sparseMatrix.end() )
        return (it->second)();
    else
        throw std::runtime_error("SparseMatrix '" + name + "' has not been registered.\n");
}
// ----------------------------------------------------------------------------
UserFunction* ObjectFactory::instantiateUserFunction( std::string name )
{
    auto it = _userFunction.find(name);
    if ( it != _userFunction.end() )
        return (it->second)();
    else
        throw std::runtime_error("UserFunction '" + name + "' has not been registered.\n");
}

// Class Registration
// ----------------------------------------------------------------------------
bool ObjectFactory::registerConvergenceCriterion( std::string name, ConvergenceCriterion* (*creator)() )
{
    auto it = _convergenceCriterion.find(name);
    if ( it != _convergenceCriterion.end() )
    {
        std::string errmsg = "\nERROR: Multiple registrations of ConvergenceCriterion '" + name + "' detected!\n";
        std::printf("%s", errmsg.c_str());
        _errorInRegistration = true;
    }
    else
    {
        _convergenceCriterion[name] = creator;
        if ( VERBOSE_REGISTRATION )
            std::printf("Registered ConvergenceCriterion class '%s'\n", name.c_str());
    }
    return _errorInRegistration;
}
// ----------------------------------------------------------------------------
bool ObjectFactory::registerLinearSolver( std::string name, LinearSolver* (*creator)() )
{
    auto it = _linearSolver.find(name);
    if ( it != _linearSolver.end() )
    {
        std::string errmsg = "\nERROR: Multiple registrations of LinearSolver '" + name + "' detected!\n";
        std::printf("%s", errmsg.c_str());
        _errorInRegistration = true;
    }
    else
    {
        _linearSolver[name] = creator;
        if ( VERBOSE_REGISTRATION )
            std::printf("Registered LinearSolver class '%s'\n", name.c_str());
    }
    return _errorInRegistration;
}
// ----------------------------------------------------------------------------
bool ObjectFactory::registerMaterial( std::string name, Material* (*creator)() )
{
    auto it = _material.find(name);
    if ( it != _material.end() )
    {
        std::string errmsg = "\nERROR: Multiple registrations of Material '" + name + "' detected!\n";
        std::printf("%s", errmsg.c_str());
        _errorInRegistration = true;
    }
    else
    {
        _material[name] = creator;
        if ( VERBOSE_REGISTRATION )
            std::printf("Registered Material '%s'\n", name.c_str());
    }
    return _errorInRegistration;
}
// ----------------------------------------------------------------------------
bool ObjectFactory::registerMeshReader( std::string name, MeshReader* (*creator)() )
{
    auto it = _meshReader.find(name);
    if ( it != _meshReader.end() )
    {
        std::string errmsg = "\nERROR: Multiple registrations of MeshReader '" + name + "' detected!\n";
        std::printf("%s", errmsg.c_str());
        _errorInRegistration = true;
    }
    else
    {
        _meshReader[name] = creator;
        if ( VERBOSE_REGISTRATION )
            std::printf("Registered MeshReader '%s'\n", name.c_str());
    }
    return _errorInRegistration;
}
// ----------------------------------------------------------------------------
bool ObjectFactory::registerNumerics( std::string name, Numerics* (*creator)() )
{
    auto it = _numerics.find(name);
    if ( it != _numerics.end() )
    {
        std::string errmsg = "\nERROR: Multiple registrations of Numerics '" + name + "' detected!\n";
        std::printf("%s", errmsg.c_str());
        _errorInRegistration = true;
    }
    else
    {
        _numerics[name] = creator;
        if ( VERBOSE_REGISTRATION )
            std::printf("Registered Numerics '%s'\n", name.c_str());
    }
    return _errorInRegistration;
}
// ----------------------------------------------------------------------------
bool ObjectFactory::registerOutputQuantity( std::string name, OutputQuantity* (*creator)() )
{
    auto it = _outputQuantity.find(name);
    if ( it != _outputQuantity.end() )
    {
        std::string errmsg = "\nERROR: Multiple registrations of OutputQuantity '" + name + "' detected!\n";
        std::printf("%s", errmsg.c_str());
        _errorInRegistration = true;
    }
    else
    {
        _outputQuantity[name] = creator;
        if ( VERBOSE_REGISTRATION )
            std::printf("Registered OutputQuantity '%s'\n", name.c_str());
    }
    return _errorInRegistration;
}
// ----------------------------------------------------------------------------
bool ObjectFactory::registerOutputWriter( std::string name, OutputWriter* (*creator)() )
{
    auto it = _outputWriter.find(name);
    if ( it != _outputWriter.end() )
    {
        std::string errmsg = "\nERROR: Multiple registrations of OutputQuantity '" + name + "' detected!\n";
        std::printf("%s", errmsg.c_str());
        _errorInRegistration = true;
    }
    else
    {
        _outputWriter[name] = creator;
        if ( VERBOSE_REGISTRATION )
            std::printf("Registered OutputWriter '%s'\n", name.c_str());
    }
    return _errorInRegistration;
}
// ----------------------------------------------------------------------------
bool ObjectFactory::registerSolutionMethod( std::string name, SolutionMethod* (*creator)() )
{
    auto it = _solutionMethod.find(name);
    if ( it != _solutionMethod.end() )
    {
        std::string errmsg = "\nERROR: Multiple registrations of SolutionMethod '" + name + "' detected!\n";
        std::printf("%s", errmsg.c_str());
        _errorInRegistration = true;
    }
    else
    {
        _solutionMethod[name] = creator;
        if ( VERBOSE_REGISTRATION )
            std::printf("Registered SolutionMethod '%s'\n", name.c_str());
    }
    return _errorInRegistration;
}
// ----------------------------------------------------------------------------
bool ObjectFactory::registerSparseMatrix( std::string name, SparseMatrix* (*creator)() )
{
    auto it = _sparseMatrix.find(name);
    if ( it != _sparseMatrix.end() )
    {
        std::string errmsg = "\nERROR: Multiple registrations of SparseMatrix '" + name + "' detected!\n";
        std::printf("%s", errmsg.c_str());
        _errorInRegistration = true;
    }
    else
    {
        _sparseMatrix[name] = creator;
        if ( VERBOSE_REGISTRATION )
            std::printf("Registered SparseMatrix '%s'\n", name.c_str());
    }
    return _errorInRegistration;
}
// ----------------------------------------------------------------------------
bool ObjectFactory::registerUserFunction( std::string name, UserFunction* (*creator)() )
{
    auto it = _userFunction.find(name);
    if ( it != _userFunction.end() )
    {
        std::string errmsg = "\nERROR: Multiple registrations of SparseMatrix '" + name + "' detected!\n";
        std::printf("%s", errmsg.c_str());
        _errorInRegistration = true;
    }
    else
    {
        _userFunction[name] = creator;
        if ( VERBOSE_REGISTRATION )
            std::printf("Registered SparseMatrix '%s'\n", name.c_str());
    }
    return _errorInRegistration;
}

// Meyers singleton implementation
// ----------------------------------------------------------------------------
ObjectFactory& broomstyx::objectFactory()
{
    static ObjectFactory cf;
    return cf;
}