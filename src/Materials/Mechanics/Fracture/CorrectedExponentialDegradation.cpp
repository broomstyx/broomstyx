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

#include "CorrectedExponentialDegradation.hpp"
#include <cmath>
#include "../../../Core/ObjectFactory.hpp"
#include "../../../Util/readOperations.hpp"

using namespace broomstyx;

registerBroomstyxObject(Material, CorrectedExponentialDegradation)

// Constructor
CorrectedExponentialDegradation::CorrectedExponentialDegradation()
{
    _name = "CorrectedExponentialDegradation";
}

// Destructor
CorrectedExponentialDegradation::~CorrectedExponentialDegradation() {}

// Public methods
// ----------------------------------------------------------------------------
double CorrectedExponentialDegradation::givePotentialFrom( const RealVector& conState, const MaterialStatus* matStatus )
{
    double phi = conState(0);
    
    double f0, f1, f2;
    f0 = (1.0 - std::exp(-_k));
    f1 = std::fabs(1.0 - phi);
    f2 = std::exp(-_k*std::pow(f1, _n));
    
    double val;
    if ( phi < 1. )
            val = (1.-_w)*(1-f2)/f0 + _w*(_a2*f1*f1 + _a3*f1*f1*f1);
    else
        val = 0.;
    
    return val;
}
// ----------------------------------------------------------------------------
RealVector CorrectedExponentialDegradation::giveForceFrom( const RealVector& conState, const MaterialStatus* matStatus )
{
    double phi = conState(0);
    
    double f0, f1, f2;
    f0 = (1.0 - std::exp(-_k));
    f1 = std::fabs(1.0 - phi);
    f2 = std::exp(-_k*std::pow(f1, _n));
    
    RealVector force(1);
    if ( phi <= 1. )
        force(0) = -(1.-_w)*_n*_k/f0*std::pow(f1, _n-1.)*f2 - _w*(2.*_a2*f1 + 3.*_a3*f1*f1);
    else
        force(0) = (1.-_w)*_n*_k/f0*std::pow(f1, _n-1.)*f2 + _w*(2.*_a2*f1 + 3.*_a3*f1*f1);
    
    return force;
}
// ----------------------------------------------------------------------------
RealMatrix CorrectedExponentialDegradation::giveModulusFrom( const RealVector& conState, const MaterialStatus* matStatus )
{
    double phi = conState(0);
    
    double f0, f1, f2;
    f0 = (1.0 - std::exp(-_k));
    f1 = std::fabs(1.0 - phi);
    f2 = std::exp(-_k*std::pow(f1, _n));
    
    RealMatrix modulus(1,1);
    
//    // True 2nd derivative
//    modulus(0,0) = (1.-qu)*(n*k/_w*y*((n-1.)*std::pow(x,n-2.) - n*k*std::pow(x,2*n-2.))) + 2.*qu; // val = n*k/_w*y*((n-1.)*std::pow(x,n-2.) - n*k*std::pow(x,2*n-2.));

    // Use constructed 2nd derivative to avoid blow-up
    modulus(0,0) = (1.-_w)*_n*_k/f0*std::pow(f1, _n-2.)*f2 + _w*(2.*_a2 + 3.*_a3*f1);
    
    return modulus;
}

void CorrectedExponentialDegradation::readParamatersFrom( FILE* fp )
{
    _n = getRealInputFrom(fp, "Failed to read parameter 'n' from input file!", _name);
    _w = getRealInputFrom(fp, "Failed to read parameter 'w' from input file!", _name);
    
    // Check that parameters is valid
    if ( _n < 2. )
        throw std::runtime_error("Invalid input detected: parameter 'n' must be greater than\nor equal to 2.0!\n\nSource: " + _name);
    
    // Calculate parameter 'k'
    double phiStar = 1./3.;
    if ( _n > 2.0 )
        phiStar = (-_n - 1. + std::sqrt(5.*_n*_n - 6*_n + 1))/(2.*_n*(_n - 2.));
    
    _k = ((_n - 2.)*phiStar + 1.)/(_n*phiStar*std::pow(1. - phiStar, _n));
    
    // Calculate parameters 'a2' and 'a3' for corrector function
    _a2 = (3.*phiStar*phiStar - 3.)/(3.*phiStar*phiStar - 1.);
    _a3 = 2./(3.*phiStar*phiStar - 1.);
}