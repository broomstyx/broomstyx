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

/**************************************************************************
 * Input file declaration:
 * 
 * CSV_OUTPUT <n>
 * ...
 *     JIntegral_Elastic <bndryLabel> <OutputLabel>
 *         ReferencePoint <ref_x> <ref_y>
 *         CrackTangent <n_x> <n_y>
 *         StressField <sxx> <syy> <sxy>
 *         DisplacementGradient <du1_dx1> <du2_dx1> <du1_dx2> <du2_dx2>
 *  
 */

#ifndef JINTEGRAL_ELASTIC_HPP
#define	JINTEGRAL_ELASTIC_HPP

#include "OutputQuantity.hpp"
#include "Util/RealVector.hpp"

namespace broomstyx
{
    class Cell;

    class JIntegral_Elastic final : public OutputQuantity {
    public:
        JIntegral_Elastic();
        virtual ~JIntegral_Elastic();

        double computeOutput() override;
        void   initialize() override;
        void   readDataFrom(FILE* fp) override;        

    private:
        int _boundaryLabel;

        RealVector _p_ref;
        RealVector _ct;

        int _stressLabel[3];
        int _gradULabel[4];

        double JIntegral_2node_line( Cell* targetCell );
    };
}
    
#endif	/* J_INTEGRAL_ELASTIC_HPP */