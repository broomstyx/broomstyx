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

#ifndef INITIALCONDITION_HPP
#define INITIALCONDITION_HPP

#include <cstdio>
#include <string>
#include "../Util/RealVector.hpp"

namespace broomstyx
{
    class UserFunction;
    
    class InitialCondition final
    {
    public:
        InitialCondition();
        virtual ~InitialCondition();
        
        std::string domainLabel() const;
        void        readDataFrom( FILE* fp );
        std::string conditionType() const;
        int         targetDofNumber() const;
        double      valueAt( const RealVector& coor ) const;

    private:
        std::string _domainLabel;
        std::string _conditionType;
        int         _dofNumber;
        std::string _specType;
        double      _val;
        
        UserFunction* _usrFcn;
    };
}

#endif /* INITIALCONDITION_HPP */

