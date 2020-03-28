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

#include "StVenantTorsion_ErrNorm_WarpFcn.hpp"
#include "Core/AnalysisModel.hpp"
#include "Core/Dof.hpp"
#include "Core/DofManager.hpp"
#include "Core/DomainManager.hpp"
#include "Core/NumericsManager.hpp"
#include "Core/ObjectFactory.hpp"
#include "User/UserFunction.hpp"
#include "Util/readOperations.hpp"
#include "IntegrationRules/Legendre_2D_Tri.hpp"
#include "BasisFunctions/Triangle_P2.hpp"

using namespace broomstyx;

registerBroomstyxObject(OutputQuantity, StVenantTorsion_ErrNorm_WarpFcn)

StVenantTorsion_ErrNorm_WarpFcn::StVenantTorsion_ErrNorm_WarpFcn()
{
    _name = "StVenantTorsion_ErrNorm_WarpFcn";
}

StVenantTorsion_ErrNorm_WarpFcn::~StVenantTorsion_ErrNorm_WarpFcn()
{
    delete _analyticSoln;
}

double StVenantTorsion_ErrNorm_WarpFcn::computeOutput()
{
    double errNorm = 0.;
    
    Numerics* numerics = analysisModel().numericsManager().giveNumerics(_numericsTag);
    int nCells = analysisModel().domainManager().giveNumberOfDomainCells();

    // Solution is assumed to be 2nd order so L2-norm integrand should be 4th order
    Triangle_P2 basisFcn;
    Legendre_2D_Tri integRule(6);
    
    std::vector<RealVector> gpLoc;
    RealVector gpWt;
    std::tie(gpLoc, gpWt) = integRule.giveIntegrationPointsAndWeights();
    
#pragma omp parallel for reduction ( + : errNorm )
    for ( int i = 0; i < nCells; i++ )
    {
        Cell* curCell = analysisModel().domainManager().giveDomainCell(i);
        Numerics* cellNumerics = analysisModel().domainManager().giveNumericsFor(curCell);
        if ( cellNumerics == numerics )
        {
            // Sanity check
            std::vector<Node*> node = analysisModel().domainManager().giveNodesOf(curCell);
            int nNodes = (int)node.size();
            if ( nNodes != 6 )
                throw std::runtime_error("Error: Invalid number of nodes (=" + std::to_string(nNodes) + ") found for domain cell!\nSource: " + _name);
            
            // Get nodal coordinates and solution
            std::vector<RealVector> coor(6, RealVector());
            RealVector xNode(6), yNode(6), uVec(6);
            for ( int j = 0; j < 6; j++ )
            {
                coor[j] = analysisModel().domainManager().giveCoordinatesOf(node[j]);
                xNode(j) = coor[j](0);
                yNode(j) = coor[j](1);
                Dof* dof = analysisModel().domainManager().giveNodalDof(_dofNum, node[j]);
                uVec(j) = analysisModel().dofManager().giveValueOfPrimaryVariableAt(dof, converged_value);
            }
            
            // Cycle over Gauss points
            for ( int j = 0; j < 6; j++ )
            {
                // Cartesian coordinates of Gauss point
                RealVector nmat = basisFcn.giveBasisFunctionsAt(gpLoc[j]);
                double x = nmat.dot(xNode);
                double y = nmat.dot(yNode);
                double u = nmat.dot(uVec);
                
                std::vector<RealVector> dpsi = basisFcn.giveBasisFunctionDerivativesAt(gpLoc[j]);
                RealMatrix Jmat(2,2);
                Jmat(0,0) = dpsi[0].dot(xNode);
                Jmat(0,1) = dpsi[0].dot(yNode);
                Jmat(1,0) = dpsi[1].dot(xNode);
                Jmat(1,1) = dpsi[1].dot(yNode);
                double J = Jmat(0,0)*Jmat(1,1) - Jmat(0,1)*Jmat(1,0);
                
                RealVector cartCoor({x, y, 0.});
                TimeData dummyTime;                
                double u_exact = _analyticSoln->at(cartCoor, dummyTime);
                
                errNorm += (u - u_exact)*(u - u_exact)*J*gpWt(j);
            }
        }
    }
    
    return errNorm;
}

void StVenantTorsion_ErrNorm_WarpFcn::initialize() {}

void StVenantTorsion_ErrNorm_WarpFcn::readDataFrom( FILE* fp )
{
    _numericsTag = getIntegerInputFrom(fp, "Failed to read numerics label from input file!", _name);
    _domainTag = getStringInputFrom(fp, "Failed to read physical label for domain from input file!", _name);
    
    verifyKeyword(fp, "NodalDof", _name);
    std::string name = getStringInputFrom(fp, "Failed to read DOF name from input file!", _name);
    _dofNum = analysisModel().dofManager().giveIndexForNodalDof(name);
    _fcnName = getStringInputFrom(fp, "Failed to read user function name for analytical solution from input file!", _name);
    
    _analyticSoln = objectFactory().instantiateUserFunction(_fcnName);
}