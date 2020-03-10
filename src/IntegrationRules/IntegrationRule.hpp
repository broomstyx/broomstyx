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

#ifndef INTEGRATIONRULE_HPP
#define INTEGRATIONRULE_HPP

#include <tuple>
#include <vector>
#include "Util/RealVector.hpp"

namespace broomstyx
{
    class IntegrationRule
    {
    public:
        IntegrationRule( int nPoints );
        virtual ~IntegrationRule();
        
        // Disable copy constructor and assignment operator
        IntegrationRule( const IntegrationRule& ) = delete;
        IntegrationRule& operator=( const IntegrationRule& ) = delete;
    
        int giveNumberOfIntegrationPoints();
        virtual std::tuple< std::vector<RealVector>, RealVector > giveIntegrationPointsAndWeights() = 0;

    protected:
        int _nIntegrationPoints;
    };
}

#endif /* INTEGRATIONRULE_HPP */

