#include "Constant_CERR.hpp"
#include "Core/ObjectFactory.hpp"
#include "Util/readOperations.hpp"
#include <stdexcept>

using namespace broomstyx;

registerBroomstyxObject(Material, Constant_CERR)

// Constructor
Constant_CERR::Constant_CERR()
{
    _name = "Constant_CERR";
}

// Destructor
Constant_CERR::~Constant_CERR() {}

// Public methods
// ----------------------------------------------------------------------------
RealMatrix Constant_CERR::giveModulusFrom( const RealVector& conState, const MaterialStatus* matStatus)
{
    return RealMatrix({{_Gc}});
}
// ----------------------------------------------------------------------------
void Constant_CERR::readParamatersFrom( FILE* fp )
{
    _Gc = getRealInputFrom(fp, "Failed to read critical energy release rate from input file!", _name);
}
