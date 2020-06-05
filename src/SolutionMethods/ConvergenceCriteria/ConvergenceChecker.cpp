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

// Constructor
ConvergenceChecker::ConvergenceChecker() {}

// Destructor
ConvergenceChecker::~ConvergenceChecker() {}

// Public methods
// ----------------------------------------------------------------------------------------
RealMatrix ConvergenceChecker::giveConvergenceData()
{
    RealMatrix convData(2*_nDofGroups, 2);

    for ( int i = 0; i < _nDofGroups; i++ )
    {
        convData(2*i, 0) = _corrNorm(i);
        convData(2*i, 1) = _corrCrit(i);
        convData(2*i+1, 0) = _residNorm(i);
        convData(2*i+1, 1) = _residCrit(i);
    }

    return convData;
}
// ----------------------------------------------------------------------------------------
std::vector<int> ConvergenceChecker::giveDofGroupNumbers()
{
    return _dofGrpNum;
}
// ----------------------------------------------------------------------------------------
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