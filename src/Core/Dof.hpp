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

#ifndef DOF_HPP
#define	DOF_HPP

#include "DofManager.hpp"

namespace broomstyx
{
    
    class Dof final
    {
        friend class DofManager;

    public:
        Dof( int grp )
            : _group(grp)
            , _stage(UNASSIGNED)
            , _subsystem(UNASSIGNED)
            , _eqNo(UNASSIGNED)
            , _isConstrained(false)
            , _primVarConverged(0.)
            , _primVarCurrent(0.)
            , _secVar(0.)
            , _residual(0.)
            , _isSlave(false)
        {}

        virtual ~Dof() {}
        
        // Disable copy constructor and assignment operator
        Dof( const Dof& ) = delete;
        Dof& operator=( const Dof& ) = delete;

    private:
        int _group;
        int _stage;
        int _subsystem;

        int    _eqNo;
        double _constraintValue;
        bool   _isConstrained;
        double _primVarConverged;
        double _primVarCurrent;
        double _primVarCorrection;
        double _secVar;
        double _residual;

        bool   _isSlave;
        Dof*   _masterDof;
    };

}

#endif	/* DOF_HPP */