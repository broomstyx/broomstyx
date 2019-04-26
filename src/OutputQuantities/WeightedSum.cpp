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

#include "WeightedSum.hpp"
#include <string>
#include "../Core/AnalysisModel.hpp"
#include "../Core/ObjectFactory.hpp"
#include "../Core/OutputManager.hpp"
#include "OutputQuantity.hpp"

using namespace broomstyx;

registerBroomstyxObject(OutputQuantity, WeightedSum)

// Constructor
WeightedSum::WeightedSum()
{
    _name = "WeightedSum";
}

// Destructor
WeightedSum::~WeightedSum()
{
    for ( int i = 0; i < (int)_term.size(); i++)
        if ( _term[i] )
            delete _term[i];
}

// Public methods
// ----------------------------------------------------------------------------
double WeightedSum::computeOutput()
{
    double result = 0.;
    
    for ( int i = 0; i < _nTerms; i++)
        result += _factor(i)*_term[i]->computeOutput();
    
    return result;
}
// ----------------------------------------------------------------------------
void WeightedSum::initialize()
{
    for ( int i = 0; i < _nTerms; i++ )
        _term[i]->initialize();
}
// ----------------------------------------------------------------------------
void WeightedSum::readDataFrom( FILE *fp )
{
    std::string key;
    
    // Read number of terms
    _nTerms = getIntegerInputFrom(fp, "Failed to read number of terms from input file!", _name);

    // Initialize vectors
    _term.assign(_nTerms, nullptr);
    _factor.init(_nTerms);
    
    for ( int i = 0; i < _nTerms; i++)
    {
        // Read factor for ith term
        _factor(i) = getRealInputFrom(fp, "Failed to read factor for term # " + std::to_string(i+1) + " from input file!", _name);
        
        // Read output quantity for ith term
        std::string oqName = getStringInputFrom(fp, "Failed to read output quantity for term #" + std::to_string(i+1) + " from input file!", _name);
        _term[i] = objectFactory().instantiateOutputQuantity(oqName);
        
        _term[i]->readDataFrom(fp);
    }
}