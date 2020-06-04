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

#include "ConvergenceChecker.hpp"
#include "Util/readOperations.hpp"

using namespace broomstyx;

ConvergenceChecker::ConvergenceChecker() {}

ConvergenceChecker::~ConvergenceChecker() {}

void ConvergenceChecker::initialize( int nDofGroups )
{
    _nDofGroups = nDofGroups;
    _relTolCor.init(_nDofGroups);
    _relTolRes.init(_nDofGroups);
    _absTolCor.init(_nDofGroups);
    _absTolRes.init(_nDofGroups);

    _normCor.init(_nDofGroups);
    _normRes.init(_nDofGroups);
    _critCor.init(_nDofGroups);
    _critRes.init(_nDofGroups);
}

void ConvergenceChecker::readDataFromFile( FILE* fp )
{
    for ( int i = 0; i < _nDofGroups; i++ )
    {
        // Read DOF group number
        _dofGrpNum[i] = getIntegerInputFrom(fp, "Failed to read DOF group number from input file!", _name);

        // Read DOF group convergence parameters
        verifyKeyword(fp, "Parameters", _name);

        // Correction tolerance
        _relTolCor(i) = getRealInputFrom(fp, "Failed to read correction tolerance from input file!", _name);

        // Residual tolerance
        _relTolRes(i) = getRealInputFrom(fp, "Failed to read residual tolerance from input file!", _name);

        // Absolute tolerance for corrections
        _absTolCor(i) = getRealInputFrom(fp, "Failed to read absolute correction tolerance from input file!", _name);

        // Absolute tolerance for residuals
        _absTolRes(i) = getRealInputFrom(fp, "Failed to read absolute residual tolerance from input file!", _name);
    }
}

void ConvergenceChecker::reportConvergenceStatus()
{
    // Report status
    std::printf("\n\n    DOF Grp   L2-Norm         Criterion");
    std::printf("\n   -------------------------------------------------------");
    
    for ( int i = 0; i < _nDofGroups; i++ )
    {
        // Convergence of corrections
        double critCor = ( _critCor(i)*_relTolCor(i) > _absTolCor(i) ) ? _critCor(i)*_relTolCor(i) : _absTolCor(i);
        std::printf("\n     C  %-6d%e   %e ", _dofGrpNum[i], _normCor(i), critCor);
        
        if ( _normCor(i) <= critCor )
            std::printf(" << CONVERGED");
        
        // Convergence of residuals
        double critRes = ( _critRes(i)*_relTolRes(i) > _absTolRes(i) ) ? _critRes(i)*_relTolRes(i) : _absTolRes(i);
        std::printf("\n     R  %-6d%e   %e ", _dofGrpNum[i], _normRes(i), critRes);
        
        if ( _normRes(i) <= critRes )
            std::printf(" << CONVERGED");
    }
}