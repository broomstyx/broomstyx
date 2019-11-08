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

#include "Diagnostics.hpp"
#include <cstdio>

#define ZEROTIME_TOL 1.0e-15

using namespace broomstyx;

Diagnostics::Diagnostics()
    : _nCoefMatAssembly(0)
    , _nConvergenceChecks(0)
    , _nLhsAssembly(0)
    , _nRhsAssembly(0)
    , _nSolves(0)
    , _nUpdates(0)
    , _coefMatAssemblyTime(0.)
    , _lhsAssemblyTime(0.)
    , _outputWriteTime(0.)
    , _rhsAssemblyTime(0.)
    , _setupTime(0.)
    , _solveTime(0.)
    , _updateTime(0.)
{}

Diagnostics::~Diagnostics()
{}

// Public methods
// ----------------------------------------------------------------------------
void Diagnostics::addCoefMatAssemblyTime( double duration )
{
    ++_nCoefMatAssembly;
    _coefMatAssemblyTime += duration;
}

void Diagnostics::addConvergenceCheckTime( double duration )
{
    ++_nConvergenceChecks;
    _convergenceCheckTime += duration;
}

void Diagnostics::addLhsAssemblyTime( double duration )
{
    ++_nLhsAssembly;
    _lhsAssemblyTime += duration;
}

void Diagnostics::addOutputWriteTime( double duration )
{
    _outputWriteTime += duration;
}

void Diagnostics::addPostprocessingTime( double duration )
{
    _postprocessingTime += duration;
}

void Diagnostics::addRhsAssemblyTime( double duration )
{
    ++_nRhsAssembly;
    _rhsAssemblyTime += duration;
}

void Diagnostics::addSetupTime( double duration )
{
    _setupTime += duration;
}

void Diagnostics::addSolveTime( double duration )
{
    ++_nSolves;
    _solveTime += duration;
}

void Diagnostics::addUpdateTime( double duration )
{
    ++_nUpdates;
    _updateTime += duration;
}

void Diagnostics::outputDiagnostics()
{
    std::printf("\n================= SIMULATION DIAGNOSTICS =================\n\n");
    std::printf("                   Number    Total time (seconds)\n\n");
    std::printf("%-20s          %f\n", "Problem setup", _setupTime);
    std::printf("%-20s          %f\n", "System Assembly", _coefMatAssemblyTime + _lhsAssemblyTime + _rhsAssemblyTime);
    if ( _coefMatAssemblyTime > ZEROTIME_TOL)
        std::printf("%-20s%-10d(%f)\n", "  Coef. Matrix", _nCoefMatAssembly, _coefMatAssemblyTime);
    if ( _lhsAssemblyTime > ZEROTIME_TOL )
        std::printf("%-20s%-10d(%f)\n", "  Left hand side", _nLhsAssembly, _lhsAssemblyTime);
    if ( _rhsAssemblyTime > ZEROTIME_TOL )
        std::printf("%-20s%-10d(%f)\n", "  Right hand side", _nRhsAssembly, _rhsAssemblyTime);
    std::printf("%-20s%-10d%f\n", "Linear Solve", _nSolves, _solveTime);
    if ( _convergenceCheckTime > ZEROTIME_TOL )
        std::printf("%-20s%-10d%f\n", "Convergence checks", _nConvergenceChecks, _convergenceCheckTime);
    std::printf("%-20s          %f\n", "Update", _updateTime);
    std::printf("%-20s          %f\n", "Postprocessing", _postprocessingTime);
    std::printf("%-20s          %f\n", "Output write", _outputWriteTime);
    std::printf("\n");
    std::printf("%-20s          %f\n", "Sum", _setupTime
                                            + _coefMatAssemblyTime 
                                            + _lhsAssemblyTime 
                                            + _rhsAssemblyTime 
                                            + _solveTime
                                            + _convergenceCheckTime 
                                            + _updateTime 
                                            + _postprocessingTime 
                                            + _outputWriteTime);
    std::printf("\n==========================================================\n\n");
}

// ----------------------------------------------------------------------------
Diagnostics& broomstyx::diagnostics()
{
    static Diagnostics diagnostics;
    return diagnostics;
}