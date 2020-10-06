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

#include "Quadrilateral_P2_4.hpp"
#include <cmath>
#include "../Core/AnalysisModel.hpp"
#include "../Core/DomainManager.hpp"
//#include "../Util/linearAlgebra.hpp"


using namespace broomstyx;

Quadrilateral_P2_4::Quadrilateral_P2_4() {}

Quadrilateral_P2_4::~Quadrilateral_P2_4() {}

RealVector Quadrilateral_P2_4::giveBasisFunctionsAt( const RealVector& coor )
{
    RealVector xi({coor(0)});
    RealVector eta({coor(1)});

    RealVector ld_psi_xi  = _ld_basis.giveBasisFunctionsAt(xi);
    RealVector ld_psi_eta = _ld_basis.giveBasisFunctionsAt(eta);

    RealVector psi({ld_psi_xi(0)*ld_psi_eta(0),
                    ld_psi_xi(1)*ld_psi_eta(0),
                    ld_psi_xi(1)*ld_psi_eta(1),
                    ld_psi_xi(0)*ld_psi_eta(1),
                    ld_psi_xi(2)*ld_psi_eta(0),
                    ld_psi_xi(1)*ld_psi_eta(2),
                    ld_psi_xi(2)*ld_psi_eta(1),
                    ld_psi_xi(0)*ld_psi_eta(2),
                    ld_psi_xi(2)*ld_psi_eta(2)});

    return psi;
}
// ----------------------------------------------------------------------------
std::vector<RealVector> Quadrilateral_P2_4::giveBasisFunctionDerivativesAt( const RealVector& coor )
{
    std::vector<RealVector> dpsi;
    dpsi.assign(2, RealVector());

    RealVector xi({coor(0)});
    RealVector eta({coor(1)});

    RealVector ld_psi_xi  = _ld_basis.giveBasisFunctionsAt(xi);
    RealVector ld_psi_eta = _ld_basis.giveBasisFunctionsAt(eta);
    std::vector<RealVector> ld_dpsi_xi  = _ld_basis.giveBasisFunctionDerivativesAt(xi);
    std::vector<RealVector> ld_dpsi_eta = _ld_basis.giveBasisFunctionDerivativesAt(eta);
    
    dpsi[0] = {ld_dpsi_xi[0](0)*ld_psi_eta(0),
               ld_dpsi_xi[0](1)*ld_psi_eta(0),
               ld_dpsi_xi[0](1)*ld_psi_eta(1),
               ld_dpsi_xi[0](0)*ld_psi_eta(1),
               ld_dpsi_xi[0](2)*ld_psi_eta(0),
               ld_dpsi_xi[0](1)*ld_psi_eta(2),
               ld_dpsi_xi[0](2)*ld_psi_eta(1),
               ld_dpsi_xi[0](0)*ld_psi_eta(2),
               ld_dpsi_xi[0](2)*ld_psi_eta(2)};

    dpsi[1] = {ld_psi_xi(0)*ld_dpsi_eta[0](0),
               ld_psi_xi(1)*ld_dpsi_eta[0](0),
               ld_psi_xi(1)*ld_dpsi_eta[0](1),
               ld_psi_xi(0)*ld_dpsi_eta[0](1),
               ld_psi_xi(2)*ld_dpsi_eta[0](0),
               ld_psi_xi(1)*ld_dpsi_eta[0](2),
               ld_psi_xi(2)*ld_dpsi_eta[0](1),
               ld_psi_xi(0)*ld_dpsi_eta[0](2),
               ld_psi_xi(2)*ld_dpsi_eta[0](2)};
    
    return dpsi;
}