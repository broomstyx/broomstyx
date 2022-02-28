#include "IsotropicThermalConductivity.hpp"
#include <string>
#include "Core/ObjectFactory.hpp"
#include "Util/linearAlgebra.hpp"
#include "Util/readOperations.hpp"

using namespace broomstyx;

registerBroomstyxObject(Material, IsotropicThermalConductivity)

// Constructor
IsotropicThermalConductivity::IsotropicThermalConductivity()
{
    _name = "IsotropicThermalConductivity";
}

// Destructor
IsotropicThermalConductivity::~IsotropicThermalConductivity() {}

// Public methods
// ---------------------------------------------------------------------------
RealMatrix IsotropicThermalConductivity::giveModulusFrom( const RealVector& conState, const MaterialStatus* matStatus )
{
    return _k;
}

double IsotropicThermalConductivity::giveParameter( const std::string& str )
{
    if ( _analysisMode == 2 )
    {
        if ( str == "ThermalConductivity_xx" )
            return _k(0,0);
        else if ( str == "ThermalConductivity_yy" )
            return _k(1,1);
        else if ( str == "ThermalConductivity_xy" )
            return _k(0,1);
        else
            throw std::runtime_error("Request for unrecognized parameter '" + str + "' made to material class " + _name);
    }
    else
    {
        if ( str == "ThermalConductivity_xx" )
            return _k(0,0);
        else if ( str == "ThermalConductivity_yy" )
            return _k(1,1);
        else if ( str == "ThermalConductivity_zz" )
            return _k(2,2);
        else if ( str == "ThermalConductivity_xy" )
            return _k(0,1);
        else if ( str == "ThermalConductivity_xz" )
            return _k(0,2);
        else if ( str == "ThermalConductivity_yz" )
            return _k(1,2);
        else
            throw std::runtime_error("Request for unrecognized parameter '" + str + "' made to material class " + _name);
    }
}

void IsotropicThermalConductivity::readParamatersFrom( FILE* fp )
{
    std::string str = getStringInputFrom(fp, "Failed to read analysis mode from input file!", _name);
    double kval = getRealInputFrom(fp, "Failed to read thermal conductivity from input file!", _name);
    
    if ( str == "2D" )
    {
        _analysisMode = 2;
        _k = {{kval, 0.},
              {0., kval}};
    }
    else if ( str == "3D" )
    {
        _analysisMode = 3;
        _k = {{kval, 0., 0.},
              {0., kval, 0.},
              {0., 0., kval}};
    }
    else
        throw std::runtime_error("Unrecognized analysis mode '" + str + "' encountered!\nSource: " + _name);
}
