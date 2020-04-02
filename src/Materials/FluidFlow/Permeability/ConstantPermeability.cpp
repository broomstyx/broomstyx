#include "ConstantPermeability.hpp"
#include <string>
#include "Core/ObjectFactory.hpp"
#include "Util/linearAlgebra.hpp"
#include "Util/readOperations.hpp"

using namespace broomstyx;

registerBroomstyxObject(Material, ConstantPermeability)

// Constructor
ConstantPermeability::ConstantPermeability()
{
    _name = "ConstantPermeability";
}

// Destructor
ConstantPermeability::~ConstantPermeability() {}

// Public methods
// ---------------------------------------------------------------------------
RealVector ConstantPermeability::giveForceFrom( const RealVector& conState, const MaterialStatus* matStatus )
{
    if ( conState.dim() != _analysisMode )
        throw std::runtime_error("Invalid size of vector 'conState' detected!\nSource: " + _name);
    
    return _k*conState;;
}

RealMatrix ConstantPermeability::giveModulusFrom( const RealVector& conState, const MaterialStatus* matStatus )
{
    return _k;
}

double ConstantPermeability::giveParameter( const std::string& str )
{
    if ( _analysisMode == 2 )
    {
        if ( str == "Permeability_xx" )
            return _k(0,0);
        else if ( str == "Permeability_yy" )
            return _k(1,1);
        else if ( str == "Permeability_xy" )
            return _k(0,1);
        else
            throw std::runtime_error("Request for unrecognized parameter '" + str + "' made to material class " + _name);
    }
    else
    {
        if ( str == "Permeability_xx" )
            return _k(0,0);
        else if ( str == "Permeability_yy" )
            return _k(1,1);
        else if ( str == "Permeability_zz" )
            return _k(2,2);
        else if ( str == "Permeability_xy" )
            return _k(0,1);
        else if ( str == "Permeability_xz" )
            return _k(0,2);
        else if ( str == "Permeability_yz" )
            return _k(1,2);
        else
            throw std::runtime_error("Request for unrecognized parameter '" + str + "' made to material class " + _name);
    }
}

void ConstantPermeability::readParamatersFrom( FILE* fp )
{
    std::string str = getStringInputFrom(fp, "Failed to read analysis mode from input file!", _name);
    
    if ( str == "2D" )
    {
        _analysisMode = 2;
        _k.init(2,2);
        _k(0,0) = getRealInputFrom(fp, "Failed to read permeability component xx from input file!", _name);
        _k(1,1) = getRealInputFrom(fp, "Failed to read permeability component yy from input file!", _name);
        _k(0,1) = getRealInputFrom(fp, "Failed to read permeability component xy from input file!", _name);
        _k(1,0) = _k(0,1);
    }
    else if ( str == "3D" )
    {
        _analysisMode = 3;
        _k.init(3,3);
        _k(0,0) = getRealInputFrom(fp, "Failed to read permeability component xx from input file!", _name);
        _k(1,1) = getRealInputFrom(fp, "Failed to read permeability component yy from input file!", _name);
        _k(2,2) = getRealInputFrom(fp, "Failed to read permeability component zz from input file!", _name);
        _k(1,2) = getRealInputFrom(fp, "Failed to read permeability component yz from input file!", _name);
        _k(0,2) = getRealInputFrom(fp, "Failed to read permeability component xz from input file!", _name);
        _k(0,1) = getRealInputFrom(fp, "Failed to read permeability component xy from input file!", _name);
        
        _k(1,0) = _k(0,1);
        _k(2,0) = _k(0,2);
        _k(2,1) = _k(1,2);
    }
    else
        throw std::runtime_error("Unrecognized analysis mode '" + str + "' encountered!\nSource: " + _name);
}