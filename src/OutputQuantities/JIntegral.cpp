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

#include "JIntegral.hpp"
#include <cmath>
#include <stdexcept>
#include <cstring>
#include "Core/AnalysisModel.hpp"
#include "Core/ObjectFactory.hpp"
#include "Core/DomainManager.hpp"
#include "MeshReaders/MeshReader.hpp"
#include "Numerics/Numerics.hpp"
#include "Util/RealVector.hpp"
#include "Util/RealMatrix.hpp"
#include "Util/linearAlgebra.hpp"
#include "Util/readOperations.hpp"

using namespace broomstyx;

registerBroomstyxObject(OutputQuantity, JIntegral)

JIntegral::JIntegral()
{
    _name = "JIntegral";
}

// -----------------------------------------------------------------------------
JIntegral::~JIntegral() {}

// -----------------------------------------------------------------------------
double JIntegral::computeOutput()
{
    double result = 0.;

    int crackTipPhysNum = analysisModel().domainManager().givePhysicalEntityNumberFor(_crackTipLabel);
    
    Node* crackTipNode = nullptr;
    bool crackTipFound = false;
    int nBndCells = analysisModel().domainManager().giveNumberOfBoundaryCells();
    for ( int i = 0; i < nBndCells; i++ )
    {
        Cell* curBndCell = analysisModel().domainManager().giveBoundaryCell(i);
        if ( analysisModel().domainManager().giveLabelOf(curBndCell) == crackTipPhysNum )
        {
            if ( crackTipFound == true )
                throw std::runtime_error("ERROR: Multiple entities found for specified crack tip label!\nSource: " + _name);
            else
            {
                crackTipFound = true;
                std::vector<Node*> bndCellNode = analysisModel().domainManager().giveNodesOf(curBndCell);
                int nNodes = (int)bndCellNode.size();
                if ( nNodes != 1 )
                    throw std::runtime_error("ERROR: Specified crack tip has more than one node!\nSource: " + _name);
                else
                    crackTipNode = bndCellNode[0];
            }
        }
    }
    if ( !crackTipFound )
        throw std::runtime_error("ERROR: Failed to find node corresponding to crack tip!\nSource: " + _name);
    
    RealVector coor = analysisModel().domainManager().giveCoordinatesOf(crackTipNode);
    RealVector p_ref ({coor(0), coor(1)});
    
    for (  int i = 0; i < (int)_path.size(); i++ )
    {
        // Length and orientation of segment
        RealVector dx = this->giveOrientationOf(_path[i].pathCell);
        double length = std::sqrt(dx.dot(dx));

        // We need a node on the path cell
        std::vector<Node*> node = analysisModel().domainManager().giveNodesOf(_path[i].pathCell);
        RealVector coor0 = analysisModel().domainManager().giveCoordinatesOf(node[0]);

        // Check orientation of unit vectors for segment
        // Note: Integration direction is clockwise w.r.t. reference point
        RealVector dirVec(2), dirVecP(2);
        dirVec(0) = coor0(0) - p_ref(0);
        dirVec(1) = coor0(1) - p_ref(1);

        
        dirVecP(0) = -dirVec(1);
        dirVecP(1) = dirVec(0);

        if ( dirVecP(0)*dx(0) + dirVecP(1)*dx(1) < 0 )
            dx = -1.0*dx;

        // outward unit normal vectors to segment
        RealVector sn(2);
        sn(0) = dx(1)/length;
        sn(1) = -dx(0)/length;
        
        // Nodal weights for integration
        RealVector wt = this->giveNodalWeightsFor(_path[i].pathCell);
        
        for ( int j = 0; j < (int)node.size(); j++)
        {
            // Node-averaged value of stress
            RealMatrix sigma(2,2);
            sigma(0,0) = this->computeNodalAverageOf(_stressLabel[0], _path[i].nodeInfo[j]);
            sigma(0,1) = this->computeNodalAverageOf(_stressLabel[2], _path[i].nodeInfo[j]);
            sigma(1,0) = sigma(0,1);
            sigma(1,1) = this->computeNodalAverageOf(_stressLabel[1], _path[i].nodeInfo[j]);
            
            // Node-averaged value of displacement gradient
            RealMatrix gradU(2,2);
            gradU(0,0) = this->computeNodalAverageOf(_gradULabel[0], _path[i].nodeInfo[j]);
            gradU(0,1) = this->computeNodalAverageOf(_gradULabel[1], _path[i].nodeInfo[j]);
            gradU(1,0) = this->computeNodalAverageOf(_gradULabel[2], _path[i].nodeInfo[j]);
            gradU(1,1) = this->computeNodalAverageOf(_gradULabel[3], _path[i].nodeInfo[j]);
            
            // Node-averaged elastic strain energy
            double ene = this->computeNodalAverageOf(_strainEnergyLabel, _path[i].nodeInfo[j]);
            
            // Calculate surface traction
            RealVector trac = sigma*sn;

            // Calculate gradU components along tangential
            RealVector gradU_t = trp(gradU)*_ct;

            // Calculate contribution to J-integral
            result += length*wt(j)*(ene*_ct.dot(sn) - trac.dot(gradU_t));
        }
    }
    
    return result;
}
// -----------------------------------------------------------------------------
void JIntegral::initialize()
{
    // Determine number of boundary cells constituting calculation path
    int physNum = analysisModel().domainManager().givePhysicalEntityNumberFor(_boundaryLabel);
    int nCells = analysisModel().domainManager().giveNumberOfBoundaryCells();
    int count = 0;
    for (  int iCell = 0; iCell < nCells; iCell++ )
    {
        Cell* curCell = analysisModel().domainManager().giveBoundaryCell(iCell);
        if ( analysisModel().domainManager().giveLabelOf(curCell) == physNum )
            ++count;
    }
    int nPath = count;
    _path.assign(nPath, PathInfo());
    
    // Setup calculation data
    count = 0;
    for (  int iCell = 0; iCell < nCells; iCell++ )
    {
        Cell* curCell = analysisModel().domainManager().giveBoundaryCell(iCell);
        int cellLabel = analysisModel().domainManager().giveLabelOf(curCell);
        
        if ( cellLabel  == physNum )
        {
            _path[count].pathCell = curCell;
            
            std::vector<Node*> cellNode = analysisModel().domainManager().giveNodesOf(curCell);
            int nCellNodes = cellNode.size();
            
            _path[count].nodeInfo.assign(nCellNodes, NodeInfo());
            
            for ( int i = 0; i < nCellNodes; i++)
            {
                std::set<Cell*> cellSet = analysisModel().domainManager().giveAttachedDomainCellsOf(cellNode[i]);
                int nAdjCells = cellSet.size();
                
                _path[count].nodeInfo[i].adjCell.assign(nAdjCells, nullptr);
                _path[count].nodeInfo[i].adjCellNodeOrder.assign(nAdjCells, -1);
                
                int nodeCount = 0;
                for (std::set<Cell*>::iterator it = cellSet.begin(); it != cellSet.end(); ++it)
                {
                    _path[count].nodeInfo[i].adjCell[nodeCount] = *it;
                    std::vector<Node*> adjCellNode = analysisModel().domainManager().giveNodesOf(*it);
                    for ( int j = 0; j < (int)adjCellNode.size(); j++)
                        if ( adjCellNode[j] == cellNode[i] )
                        {
                            _path[count].nodeInfo[i].adjCellNodeOrder[nodeCount] = j;
                            break;
                        }
                    
                    ++nodeCount;
                }
            }
            
            ++count;
        }
    }
}
// -----------------------------------------------------------------------------
void JIntegral::readDataFrom( FILE* fp )
{
    _boundaryLabel = getStringInputFrom(fp, "Failed reading boundary label from input file.", _name);
    
    // Read reference point
    verifyKeyword(fp, "CrackTip", _name);
    _crackTipLabel = getStringInputFrom(fp, "Failed to read crack tip label from input file.", _name);
    
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

    verifyKeyword(fp, "StrainEnergy", _name);
    _strainEnergyLabel = getIntegerInputFrom(fp, "Failed reading nodal field for strain energy from input file.", _name);
}
// -----------------------------------------------------------------------------
double JIntegral::computeNodalAverageOf( int fieldNum, NodeInfo& nodeInfo )
{
    double fieldVal = 0.;
    int nAdjCells = nodeInfo.adjCell.size();
    
    for ( int i = 0; i < nAdjCells; i++)
    {
        int label = analysisModel().domainManager().giveLabelOf(nodeInfo.adjCell[i]);
        Numerics* numerics = analysisModel().domainManager().giveNumericsForDomain(label);
        
        RealVector nodalValues = numerics->giveCellNodeFieldValuesAt(nodeInfo.adjCell[i], fieldNum);
        fieldVal += nodalValues(nodeInfo.adjCellNodeOrder[i]);
    }
    
    fieldVal /= nAdjCells;
    
    return fieldVal;
}
// ----------------------------------------------------------------------------
RealVector JIntegral::giveNodalWeightsFor( Cell* targetCell )
{
    int nNodes = analysisModel().domainManager().giveNumberOfNodesOf(targetCell);
    
    RealVector wt;
    switch ( nNodes )
    {
        case 2:
            wt.init(2);
            wt(0) = 0.5;
            wt(1) = 0.5;
            break;
            
        case 3:
            wt.init(3);
            wt(0) = 1./6.;
            wt(1) = 1./6.;
            wt(2) = 2./3.;
            break;
            
        default:
            throw std::runtime_error("Error: Output quantity 'JIntegral' cannot handle path cells with " 
                    + std::to_string(nNodes) + " nodes!");
    }
    
    return wt;
}
// ----------------------------------------------------------------------------
RealVector JIntegral::giveOrientationOf( Cell* targetCell )
{
    std::vector<Node*> cellNode = analysisModel().domainManager().giveNodesOf(targetCell);
    int nNodes = cellNode.size();
    
    RealVector coorA, coorB;
    switch ( nNodes )
    {
        case 2:
        case 3:
            coorA = analysisModel().domainManager().giveCoordinatesOf(cellNode[0]);
            coorB = analysisModel().domainManager().giveCoordinatesOf(cellNode[1]);
            break;
            
        default:
            throw std::runtime_error("Error: Output quantity 'JIntegral' cannot handle path cells with " 
                    + std::to_string(nNodes) + " nodes!");
    }
    
    return coorB - coorA;
}