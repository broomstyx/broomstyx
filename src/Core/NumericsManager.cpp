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

#include "NumericsManager.hpp"
#include <stdexcept>
#include <set>
#include "AnalysisModel.hpp"
#include "../Core/ObjectFactory.hpp"
#include "../MeshReaders/MeshReader.hpp"
#include "../Util/readOperations.hpp"
#include "../Numerics/Numerics.hpp"

using namespace broomstyx;

// Constructor
NumericsManager::NumericsManager() {}

// Destructor
NumericsManager::~NumericsManager()
{
#ifdef VERBOSE_DESTRUCTION
    std::printf("\n  Destroying NumericsManager... ");
    std::fflush(stdout);
#endif

    for ( Numerics* curNumerics : _numerics )
        delete curNumerics;

#ifdef VERBOSE_DESTRUCTION
    std::printf("done.");
    std::fflush(stdout);
#endif

}

// Public methods
// ----------------------------------------------------------------------------
std::vector<Numerics*> NumericsManager::giveAllNumerics()
{
    return _numerics;
}
// ----------------------------------------------------------------------------
Numerics* NumericsManager::giveNumerics( int label )
{
    Numerics* target = nullptr;
    
    for ( Numerics* curNumerics : _numerics )
        if ( curNumerics->_label == label )
            target = curNumerics;
    
    if ( !target )
        throw std::runtime_error("Error: There is no numerics corresponding to label '" + std::to_string(label) + "'!\nSource: NumericsManager");
    
    return target;
}
// ----------------------------------------------------------------------------
void NumericsManager::readNumericsFrom( FILE *fp )
{
    std::string key, src = "NumericsManager";
    
    int nNumerics = getIntegerInputFrom(fp, "Failed to read number of numerics from input file!", src);
    _numerics.assign(nNumerics, nullptr);
    
    // Track defined labels to make sure no numerics label is repeated
    std::set<int> alreadyDefined;
    
    for ( int i = 0; i < nNumerics; i++)
    {
        // Tag for numerics type
        int label = getIntegerInputFrom(fp, "Failed to read domain name from input file!", src);
        
        auto it = alreadyDefined.find(label);
        if ( it != alreadyDefined.end() )
            throw std::runtime_error("Error: multiple numerics detected for label '" + std::to_string(label) + "'!\nSource: NumericsManager\n");
        
        // Name of numerics type
        std::string typeName = getStringInputFrom(fp, "Failed to read numerics type from input file!", src);
        
        _numerics[i] = objectFactory().instantiateNumerics(typeName);
        _numerics[i]->_label = label;
        
        alreadyDefined.insert(label);
        _numerics[i]->readDataFrom(fp);
    }
}