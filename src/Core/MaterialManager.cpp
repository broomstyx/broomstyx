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

#include "MaterialManager.hpp"
#include <omp.h>
#include <stdexcept>
#include <chrono>
#include "AnalysisModel.hpp"
#include "ObjectFactory.hpp"
#include "Materials/Material.hpp"
#include "Util/readOperations.hpp"

using namespace broomstyx;

MaterialManager::MaterialManager() {}

MaterialManager::~MaterialManager()
{
#ifdef VERBOSE_DESTRUCTION
    std::printf("\n  Destroying MaterialManager... ");
    std::fflush(stdout);
#endif
    
    for (int i = 0; i < (int)_material.size(); i++)
        if ( _material[i] )
            delete _material[i];
    
#ifdef VERBOSE_DESTRUCTION
    std::printf("done.");
    std::fflush(stdout);
#endif
}

// Public methods
// ----------------------------------------------------------------------------
Material* MaterialManager::giveMaterial( int label )
{
    Material* target;
    std::map<int,Material*>::iterator it = _material.find(label);
    if ( it != _material.end() )
        target = (*it).second;
    else
    {
        std::string errmsg = "Error: There is no material corresponding to" +
                std::string(" label '") + std::to_string(label)
                + std::string("'!\nSource: MaterialManager");
        throw std::runtime_error(errmsg);
    }
    
    return target;
}
// ----------------------------------------------------------------------------
void MaterialManager::initializeMaterials()
{
    std::chrono::time_point<std::chrono::system_clock> tic, toc;
    std::chrono::duration<double> tictoc;
    
    std::printf("  %-40s", "Initializing materials ...");
    tic = std::chrono::high_resolution_clock::now();
    
    std::map<int,Material*>::iterator it;
    for (it = _material.begin(); it != _material.end(); it++)
        (*it).second->initialize();
    
    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    std::printf("done (time = %f sec.)\n", tictoc.count());
}
// ----------------------------------------------------------------------------
void MaterialManager::readMaterialsFrom( FILE* fp )
{
    std::string errmsg, src = "MaterialsLibrary";
    int nMat = getIntegerInputFrom(fp,
        errmsg = "Failed to read number of materials from input file!", src);
    
    for ( int i = 0; i < nMat; i++)
    {
        int label = getIntegerInputFrom(fp,
                errmsg = "Failed reading material label from input file!", src);
        std::string materialName = getStringInputFrom(fp,
                errmsg = "Failed reading material type from input file!", src);
        
        Material* newMat = objectFactory().instantiateMaterial(materialName);
        
        std::pair< std::map<int,Material*>::iterator, bool> entry;
        entry = _material.insert(std::pair<int,Material*>(label, newMat));
        if ( !entry.second )
            throw std::runtime_error(
                    errmsg = "Multiple declaration of materials for label '"
                    + std::to_string(label) 
                    + std::string("' detected in input file!\nSource: ") + src);
                
        _material[label]->readParamatersFrom(fp);
    }
}