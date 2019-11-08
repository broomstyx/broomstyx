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

#ifndef DIAGNOSTICS_HPP
#define DIAGNOSTICS_HPP

namespace broomstyx
{
    class Diagnostics final
    {
        friend Diagnostics& diagnostics();

    public:
        void addCoefMatAssemblyTime( double duration );
        void addConvergenceCheckTime( double duration );
        void addLhsAssemblyTime( double duration );
        void addOutputWriteTime( double duration );
        void addPostprocessingTime( double duration );
        void addRhsAssemblyTime( double duration );
        void addSetupTime( double duration );
        void addSolveTime( double duration );
        void addUpdateTime( double duration );
        void outputDiagnostics();

    private:
        int _nCoefMatAssembly;
        int _nConvergenceChecks;
        int _nLhsAssembly;
        int _nRhsAssembly;
        int _nSolves;
        int _nUpdates;

        double _coefMatAssemblyTime;
        double _convergenceCheckTime;
        double _lhsAssemblyTime;
        double _outputWriteTime;
        double _postprocessingTime;
        double _rhsAssemblyTime;
        double _setupTime;
        double _solveTime;
        double _updateTime;

        Diagnostics();
        virtual ~Diagnostics();
    };

    Diagnostics& diagnostics();
}

#endif /* DIAGNOSTICS_HPP */