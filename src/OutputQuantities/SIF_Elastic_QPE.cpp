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

#include "SIF_Elastic_QPE.hpp"
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

#define PI 3.14159265358979323846

using namespace broomstyx;

registerBroomstyxObject(OutputQuantity, SIF_Elastic_QPE)

SIF_Elastic_QPE::SIF_Elastic_QPE()
{
    _name = "SIF_Elastic_QPE";
}

// -----------------------------------------------------------------------------
SIF_Elastic_QPE::~SIF_Elastic_QPE() {}

// -----------------------------------------------------------------------------
double SIF_Elastic_QPE::computeOutput()
{
    RealVector crackTipLoc, nodePairLoc, dX, crackNormal, crackTangent;
    RealVector posDisp, negDisp;
    Dof* uDof[2];
    double r;

    double k;
    if ( _analysisMode == "PlaneStress" )
        k = (3. - _nu)/(1. + _nu);
    else
        k = 3. - 4.*_nu;
    
    crackTipLoc = analysisModel().domainManager().giveCoordinatesOf(_crackTipNode);
    
    // Compute SIF for first pair
    nodePairLoc = analysisModel().domainManager().giveCoordinatesOf(_nodePairA[0]);
    dX = crackTipLoc - nodePairLoc;
    r = std::sqrt(dX.dot(dX));
    crackTangent = {dX(0)/r, dX(1)/r};
    crackNormal = {-dX(1)/r, dX(0)/r};
    
    uDof[0] = analysisModel().domainManager().giveNodalDof(_dispDofNum[0],_nodePairA[0]);
    uDof[1] = analysisModel().domainManager().giveNodalDof(_dispDofNum[1],_nodePairA[0]);
    posDisp = {analysisModel().dofManager().giveValueOfPrimaryVariableAt(uDof[0], converged_value),
               analysisModel().dofManager().giveValueOfPrimaryVariableAt(uDof[1], converged_value)};

    posDisp.print("Pair A pos disp", 4);

    uDof[0] = analysisModel().domainManager().giveNodalDof(_dispDofNum[0],_nodePairA[1]);
    uDof[1] = analysisModel().domainManager().giveNodalDof(_dispDofNum[1],_nodePairA[1]);
    negDisp = {analysisModel().dofManager().giveValueOfPrimaryVariableAt(uDof[0], converged_value),
               analysisModel().dofManager().giveValueOfPrimaryVariableAt(uDof[1], converged_value)};
    
    negDisp.print("Pair A neg disp", 4);

    double pairA_modeI_opening = 0.5*(posDisp - negDisp).dot(crackNormal);

    std::printf("Pair A data: r = %f, Du = %f\n", r, pairA_modeI_opening);
    double pairA_modeII_opening = 0.5*(posDisp - negDisp).dot(crackTangent);
    double KI_A = _E/((1. + _nu)*(1. + k))*std::sqrt(2.*PI/r)*pairA_modeI_opening;
    double KII_A = _E/((1. + _nu)*(1. + k))*std::sqrt(2.*PI/r)*pairA_modeII_opening;

    // Compute SIF for second pair
    nodePairLoc = analysisModel().domainManager().giveCoordinatesOf(_nodePairB[0]);
    dX = crackTipLoc - nodePairLoc;
    r = std::sqrt(dX.dot(dX));
    crackTangent = {dX(0)/r, dX(1)/r};
    crackNormal = {-dX(1)/r, dX(0)/r};
    
    uDof[0] = analysisModel().domainManager().giveNodalDof(_dispDofNum[0],_nodePairB[0]);
    uDof[1] = analysisModel().domainManager().giveNodalDof(_dispDofNum[1],_nodePairB[0]);
    posDisp = {analysisModel().dofManager().giveValueOfPrimaryVariableAt(uDof[0], converged_value),
               analysisModel().dofManager().giveValueOfPrimaryVariableAt(uDof[1], converged_value)};

    posDisp.print("Pair B pos disp", 4);

    uDof[0] = analysisModel().domainManager().giveNodalDof(_dispDofNum[0],_nodePairB[1]);
    uDof[1] = analysisModel().domainManager().giveNodalDof(_dispDofNum[1],_nodePairB[1]);
    negDisp = {analysisModel().dofManager().giveValueOfPrimaryVariableAt(uDof[0], converged_value),
               analysisModel().dofManager().giveValueOfPrimaryVariableAt(uDof[1], converged_value)};
    
    negDisp.print("Pair B neg disp", 4);

    double pairB_modeI_opening = 0.5*(posDisp - negDisp).dot(crackNormal);

    std::printf("Pair A data: r = %f, Du = %f\n", r, pairB_modeI_opening);

    double pairB_modeII_opening = 0.5*(posDisp - negDisp).dot(crackTangent);
    double KI_B = _E/((1. + _nu)*(1. + k))*std::sqrt(2.*PI/r)*pairB_modeI_opening;
    double KII_B = _E/((1. + _nu)*(1. + k))*std::sqrt(2.*PI/r)*pairB_modeII_opening;

    double KI = 4.*KI_A/3. - KI_B/3.;
    double KII = 4.*KII_A/3. - KII_B/3.;

    // Output to screen
    if ( _mode == 1 )
        std::printf("    Mode I SIF = %E\n", KI);
    else
        std::printf("    Mode II SIF = %E\n", KII);

    if ( _mode == 1 )
        return KI;
    else
        return KII;
}
// -----------------------------------------------------------------------------
void SIF_Elastic_QPE::initialize()
{
    // Get node associated with crack tip
    int crackTipPhysNum = analysisModel().domainManager().givePhysicalEntityNumberFor(_crackTipLabel);

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
                    _crackTipNode = bndCellNode[0];
            }
        }
    }
    if ( !crackTipFound )
        throw std::runtime_error("ERROR: Failed to find node corresponding to crack tip!\nSource: " + _name);

    // ---->
    std::printf("\nCrack tip found at node # %d\n", analysisModel().domainManager().giveIdOf(_crackTipNode));
    RealVector temp = analysisModel().domainManager().giveCoordinatesOf(_crackTipNode);
    temp.print("Crack tip coordinates: ", 4);
    // ---->>

    // Find relevant boundary cells on crack face
    int crackFacePhysNum[2];
    crackFacePhysNum[0] = analysisModel().domainManager().givePhysicalEntityNumberFor(_crackFaceLabel[0]);
    crackFacePhysNum[1] = analysisModel().domainManager().givePhysicalEntityNumberFor(_crackFaceLabel[1]);
    std::vector<Cell*> crackFaceCell[2];

    // A. "Positive" face
    for ( int i = 0; i < nBndCells; i++ )
    {
        Cell* curBndCell = analysisModel().domainManager().giveBoundaryCell(i);
        if ( analysisModel().domainManager().giveLabelOf(curBndCell) == crackFacePhysNum[0] )
        {
            std::vector<Node*> node = analysisModel().domainManager().giveNodesOf(curBndCell);
            for ( int j = 0; j < (int)node.size(); j++ )
                if ( node[j] == _crackTipNode )
                {
                    // Sanity check
                    if ( (int)node.size() != 3 )
                        throw std::runtime_error("ERROR: Boundary cell attached to crack tip must have 3 nodes!\nSource: " + _name);
                    
                    // ---->
                    std::printf("Positive boundary cell found at cell # %d.\n", analysisModel().domainManager().giveIdOf(curBndCell));
                    std::printf("Nodes: ");
                    auto bcellNode = analysisModel().domainManager().giveNodesOf(curBndCell);
                    for ( int k = 0; k < (int)bcellNode.size(); k++ )
                        std::printf("%d ", analysisModel().domainManager().giveIdOf(bcellNode[k]));
                    std::printf("\n");
                    for ( int k = 0; k < (int)bcellNode.size(); k++ )
                    {
                        RealVector bNodeCoor = analysisModel().domainManager().giveCoordinatesOf(bcellNode[k]);
                        bNodeCoor.print(("Coordinates of node " + std::to_string(k)).c_str(), 4);
                    }
                    // ---->>

                    if ( j == 0 )
                    {
                        _nodePairA[0] = node[2];
                        _nodePairB[0] = node[1];
                    }
                    else // (j == 1)
                    {
                        _nodePairA[0] = node[2];
                        _nodePairB[0] = node[0];
                    }
                }
        }
    }

    // B. "Negative" face
    for ( int i = 0; i < nBndCells; i++ )
    {
        Cell* curBndCell = analysisModel().domainManager().giveBoundaryCell(i);
        if ( analysisModel().domainManager().giveLabelOf(curBndCell) == crackFacePhysNum[1] )
        {
            std::vector<Node*> node = analysisModel().domainManager().giveNodesOf(curBndCell);
            for ( int j = 0; j < (int)node.size(); j++ )
                if ( node[j] == _crackTipNode )
                {
                    // Sanity check
                    if ( (int)node.size() != 3 )
                        throw std::runtime_error("ERROR: Boundary cell attached to crack tip must have 3 nodes!\nSource: " + _name);
                    
                    // ---->
                    std::printf("Negative boundary cell found at cell # %d.\n", analysisModel().domainManager().giveIdOf(curBndCell));
                    std::printf("Nodes: ");
                    auto bcellNode = analysisModel().domainManager().giveNodesOf(curBndCell);
                    for ( int k = 0; k < (int)bcellNode.size(); k++ )
                        std::printf("%d ", analysisModel().domainManager().giveIdOf(bcellNode[k]));
                    std::printf("\n");
                    for ( int k = 0; k < (int)bcellNode.size(); k++ )
                    {
                        RealVector bNodeCoor = analysisModel().domainManager().giveCoordinatesOf(bcellNode[k]);
                        bNodeCoor.print(("Coordinates of node " + std::to_string(k)).c_str(), 4);
                    }
                    // ---->>
                    
                    if ( j == 0 )
                    {
                        _nodePairA[1] = node[2];
                        _nodePairB[1] = node[1];
                    }
                    else // (j == 1)
                    {
                        _nodePairA[1] = node[2];
                        _nodePairB[1] = node[0];
                    }
                }
        }
    }
}
// -----------------------------------------------------------------------------
void SIF_Elastic_QPE::readDataFrom( FILE* fp )
{
    // Analysis mode
    _analysisMode = getStringInputFrom(fp, "Failed reading analysis mode from input file.", _name);
    if ( _analysisMode != "PlaneStress" && _analysisMode != "PlaneStrain" )
        throw std::runtime_error("ERROR: Invalid analysis mode '" + _analysisMode + "' detected!\nSource: " + _name);
    
    // Material parameters
    _E = getRealInputFrom(fp, "Failed to read Young's modulus from input file.", _name);
    _nu = getRealInputFrom(fp, "Failed to read Poisson ratio from input file.", _name);
    
    // Crack tip
    verifyKeyword(fp, "CrackTip", _name);
    _crackTipLabel = getStringInputFrom(fp, "Failed reading crack tip label from input file.", _name);
    
    // Crack faces
    verifyKeyword(fp, "CrackFaces", _name);
    _crackFaceLabel[0] = getStringInputFrom(fp, "Failed reading first crack face label from input file.", _name);
    _crackFaceLabel[1] = getStringInputFrom(fp, "Failed reading second crack face label from input file.", _name);

    // Displacement DOFs
    verifyKeyword(fp, "DisplacementDOFs", _name);
    std::string dofLabel;
    dofLabel = getStringInputFrom(fp, "Failed reading first displacement DOF label from input file.", _name);
    _dispDofNum[0] = analysisModel().dofManager().giveIndexForNodalDof(dofLabel);
    dofLabel = getStringInputFrom(fp, "Failed reading second displacement DOF label from input file.", _name);
    _dispDofNum[1] = analysisModel().dofManager().giveIndexForNodalDof(dofLabel);

    // SIF mode computation
    verifyKeyword(fp, "SIFMode", _name);
    _mode = getIntegerInputFrom(fp, "Failed reading SIF mode from input file.", _name);
    if ( _mode != 1 && _mode != 2 )
        throw std::runtime_error("ERROR: Invalid SIF mode '" + std::to_string(_mode) + "' detected!\nSource: " + _name);
}