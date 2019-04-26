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

#include "JIntegral_Elastic.hpp"
#include <cmath>
#include <stdexcept>
#include <cstring>
#include "../Core/AnalysisModel.hpp"
#include "../Core/ObjectFactory.hpp"
#include "../Core/DomainManager.hpp"
#include "../Util/RealVector.hpp"
#include "../Util/RealMatrix.hpp"
#include "../Util/linearAlgebra.hpp"
#include "../Util/readOperations.hpp"

using namespace broomstyx;

registerBroomstyxObject(OutputQuantity, JIntegral_Elastic)

JIntegral_Elastic::JIntegral_Elastic()
{
    _name = "J_Integral_Elastic";
}

JIntegral_Elastic::~JIntegral_Elastic() {}

double JIntegral_Elastic::computeOutput()
{
    double result = 0;
    
    int nCells = analysisModel().domainManager().giveNumberOfBoundaryCells();
    for (int i = 0; i < nCells; i++)
    {
        Cell* curCell = analysisModel().domainManager().giveBoundaryCell(i);
        if (analysisModel().domainManager().giveLabelOf(curCell) == _boundaryLabel)
        {
            double dJ;
            int nNodes = analysisModel().domainManager().giveNumberOfNodesOf(curCell);
            if ( nNodes == 2 )
                dJ = this->JIntegral_2node_line(curCell);
            else
                throw std::runtime_error("Output quantity 'JIntegral_Elastic' only handles 2-node line elements!");
            
            result += dJ;
        }
    }
    
    return result;
}

void JIntegral_Elastic::initialize() {}

void JIntegral_Elastic::readDataFrom( FILE* fp )
{
    _boundaryLabel = getIntegerInputFrom(fp, "Failed reading boundary label from input file.", _name);

    // Read reference point
    verifyKeyword(fp, "ReferencePoint", _name);
    _p_ref.init(2);
    _p_ref(0) = getRealInputFrom(fp, "Failed reading reference point x-coordinate from input file.", _name);
    _p_ref(1) = getRealInputFrom(fp, "Failed reading reference point y-coordinate from input file.", _name);
    
    // Tangential direction for crack
    verifyKeyword(fp, "CrackTangent", _name);
    _ct.init(2);
    _ct(0) = getRealInputFrom(fp, "Failed reading x-component of crack tangent vector from input file.", _name);
    _ct(1) = getRealInputFrom(fp, "Failed reading y-component of crack tangent vector from input file.", _name);
    
    // Normalize crack tangent vector
    double crkLen = std::sqrt(_ct.dot(_ct));
    _ct(0) /= crkLen;
    _ct(1) /= crkLen;
    
    verifyKeyword(fp, "StressField", _name);
    _stressLabel[0] = getIntegerInputFrom(fp, "Failed reading nodal field number for xx-stress from input file.", _name);
    _stressLabel[1] = getIntegerInputFrom(fp, "Failed reading nodal field number for yy-stress from input file.", _name);
    _stressLabel[2] = getIntegerInputFrom(fp, "Failed reading nodal field number for xy-stress from input file.", _name);
    
    verifyKeyword(fp, "DisplacementGradient", _name);
    _gradULabel[0] = getIntegerInputFrom(fp, "Failed reading nodal field for disp. grad. (xx) from input file.", _name);
    _gradULabel[1] = getIntegerInputFrom(fp, "Failed reading nodal field for disp. grad. (xy) from input file.", _name);
    _gradULabel[2] = getIntegerInputFrom(fp, "Failed reading nodal field for disp. grad. (yx) from input file.", _name);
    _gradULabel[3] = getIntegerInputFrom(fp, "Failed reading nodal field for disp. grad. (yy) from input file.", _name);
}

double JIntegral_Elastic::JIntegral_2node_line( Cell* targetCell )
{    
    // Get nodes of boundary cell
    std::vector<Node*> node;
    node = analysisModel().domainManager().giveNodesOf(targetCell);
    
    // Length of boundary cell
    RealVector coor1, coor2, dx;
    coor1 = analysisModel().domainManager().giveCoordinatesOf(node[0]);
    coor2 = analysisModel().domainManager().giveCoordinatesOf(node[1]);
    
    dx = coor2 - coor1;
    double length = std::sqrt(dx.dot(dx));
    
    // Check orientation of unit vectors for segment
    // Note: Integration direction is clockwise w.r.t. reference point
    RealVector dirVec(2), dirVecP(2);
    dirVec(0) = coor1(0) - _p_ref(0);
    dirVec(1) = coor1(1) - _p_ref(1);
    
    dirVecP(0) = -dirVec(1);
    dirVecP(1) = dirVec(0);
    
    if ( dirVecP(0)*dx(0) + dirVecP(1)*dx(1) < 0 )
        dx = -1.0*dx;
    
    // outward unit normal vectors to segment
    RealVector sn(2);
    sn(0) = dx(1)/length;
    sn(1) = -dx(0)/length;
        
    // Stresses
    // A. Endpoint 1
    RealMatrix sigma_1(2,2);
    sigma_1(0,0) = analysisModel().domainManager().giveFieldValueAt(node[0], _stressLabel[0]);
    sigma_1(1,1) = analysisModel().domainManager().giveFieldValueAt(node[0], _stressLabel[1]);
    sigma_1(0,1) = analysisModel().domainManager().giveFieldValueAt(node[0], _stressLabel[2]);
    sigma_1(1,0) = sigma_1(0,1);
    
    // B. Endpoint 2
    RealMatrix sigma_2(2,2);
    sigma_2(0,0) = analysisModel().domainManager().giveFieldValueAt(node[1], _stressLabel[0]);
    sigma_2(1,1) = analysisModel().domainManager().giveFieldValueAt(node[1], _stressLabel[1]);
    sigma_2(0,1) = analysisModel().domainManager().giveFieldValueAt(node[1], _stressLabel[2]);
    sigma_2(1,0) = sigma_2(0,1);
    
    // C. Midpoint
    RealMatrix sigma_midpt = 0.5*(sigma_1 + sigma_2);
    
    // Displacement gradient and strain
    // A. Endpoint 1
    RealMatrix nablaU_1(2,2);
    nablaU_1(0,0) = analysisModel().domainManager().giveFieldValueAt(node[0], _gradULabel[0]);
    nablaU_1(0,1) = analysisModel().domainManager().giveFieldValueAt(node[0], _gradULabel[1]);
    nablaU_1(1,0) = analysisModel().domainManager().giveFieldValueAt(node[0], _gradULabel[2]);
    nablaU_1(1,1) = analysisModel().domainManager().giveFieldValueAt(node[0], _gradULabel[3]);
    
    // B. Endpoint 2
    RealMatrix nablaU_2(2,2);
    nablaU_2(0,0) = analysisModel().domainManager().giveFieldValueAt(node[1], _gradULabel[0]);
    nablaU_2(0,1) = analysisModel().domainManager().giveFieldValueAt(node[1], _gradULabel[1]);
    nablaU_2(1,0) = analysisModel().domainManager().giveFieldValueAt(node[1], _gradULabel[2]);
    nablaU_2(1,1) = analysisModel().domainManager().giveFieldValueAt(node[1], _gradULabel[3]);
    
    // C. Midpoint
    RealMatrix nablaU_midpt = 0.5*(nablaU_1 + nablaU_2);
    
    // Strain
    RealMatrix epsilon_1 = 0.5*(trp(nablaU_1) + nablaU_1);
    RealMatrix epsilon_2 = 0.5*(trp(nablaU_2) + nablaU_2);
    RealMatrix epsilon_midpt = 0.5*(epsilon_1 + epsilon_2);
    
    // Calculate energy density at integration points
    double W_1, W_2, W_midpt;
    W_1 = 0.5*(sigma_1(0,0)*epsilon_1(0,0) + sigma_1(0,1)*epsilon_1(0,1)
            + sigma_1(1,0)*epsilon_1(1,0) + sigma_1(1,1)*epsilon_1(1,1));
    W_2 = 0.5*(sigma_2(0,0)*epsilon_2(0,0) + sigma_2(0,1)*epsilon_2(0,1)
            + sigma_2(1,0)*epsilon_2(1,0) + sigma_2(1,1)*epsilon_2(1,1));
    W_midpt = 0.5*(sigma_midpt(0,0)*epsilon_midpt(0,0) 
            + sigma_midpt(0,1)*epsilon_midpt(0,1)
            + sigma_midpt(1,0)*epsilon_midpt(1,0)
            + sigma_midpt(1,1)*epsilon_midpt(1,1));
    
    // J-integral
    double Q_1, Q_2, Q_midpt;
    Q_1 = W_1*_ct.dot(sn) - (sigma_1*sn).dot(trp(nablaU_1)*_ct);
    Q_2 = W_2*_ct.dot(sn) - (sigma_2*sn).dot(trp(nablaU_2)*_ct);
    Q_midpt = W_midpt*_ct.dot(sn)-(sigma_midpt*sn).dot(trp(nablaU_midpt)*_ct);
    
    return length/6.*(Q_1 + 4.*Q_midpt + Q_2);
}