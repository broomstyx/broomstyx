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
    // Get displacement components normal and tangent to the crack
    double disp_tn[4][2];
    RealVector disp;
        
    for ( int i = 0; i < 4; i++ )
    {
        disp = {analysisModel().dofManager().giveValueOfPrimaryVariableAt(_dof[i][0], converged_value),
                analysisModel().dofManager().giveValueOfPrimaryVariableAt(_dof[i][1], converged_value)};
        
        disp_tn[i][0] = disp.dot(_crackTangent);
        disp_tn[i][1] = disp.dot(_crackNormal);
    }

    // Compute stress intensity factors
    double K_I, K_II;
    double k;
    if ( _analysisMode == "PlaneStress" )
        k = (3. - _nu)/(1. + _nu);
    else
        k = 3. - 4.*_nu;
    
    double coef = _E/(6.*(1. + _nu)*(1. + k))*std::sqrt(2*PI/_h);
    K_I = coef*(8.*(disp_tn[0][0] - disp_tn[1][0]) - (disp_tn[2][0] - disp_tn[3][0]));
    K_II = coef*(8.*(disp_tn[0][1] - disp_tn[1][1]) - (disp_tn[2][1] - disp_tn[3][1]));

    // Output to screen
    if ( _mode == 1 )
        std::printf("\nMode I SIF = %f\n", K_I);
    else
        std::printf("\nMode II SIF = %f\n", K_II);

    if ( _mode == 1 )
        return K_I;
    else
        return K_II;
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
                    
                    if ( j == 0 )
                    {
                        // Compute crack orientation
                        RealVector coorA, coorB, dL;
                        coorA = analysisModel().domainManager().giveCoordinatesOf(node[0]);
                        coorB = analysisModel().domainManager().giveCoordinatesOf(node[1]);
                        dL = coorA - coorB;
                        double len = std::sqrt(dL.dot(dL));
                        _h = len/0.75;
                        _crackTangent = dL/len;
                        _crackNormal = {-_crackTangent(1), _crackTangent(0)};
                        

                        // Store DOFs
                        _dof[0][0] = analysisModel().domainManager().giveNodalDof(_dispDofNum[0], node[2]);
                        _dof[0][1] = analysisModel().domainManager().giveNodalDof(_dispDofNum[1], node[2]);
                        _dof[2][0] = analysisModel().domainManager().giveNodalDof(_dispDofNum[0], node[1]);
                        _dof[2][1] = analysisModel().domainManager().giveNodalDof(_dispDofNum[1], node[1]);
                    }
                    else // (j == 1)
                    {
                        // Compute crack orientation
                        RealVector coorA, coorB, dL;
                        coorA = analysisModel().domainManager().giveCoordinatesOf(node[1]);
                        coorB = analysisModel().domainManager().giveCoordinatesOf(node[0]);
                        dL = coorA - coorB;
                        double len = std::sqrt(dL.dot(dL));
                        _h = len/0.75;
                        _crackTangent = dL/len;
                        _crackNormal = {-_crackTangent(1), _crackTangent(0)};
                        

                        // Store DOFs
                        _dof[0][0] = analysisModel().domainManager().giveNodalDof(_dispDofNum[0], node[2]);
                        _dof[0][1] = analysisModel().domainManager().giveNodalDof(_dispDofNum[1], node[2]);
                        _dof[2][0] = analysisModel().domainManager().giveNodalDof(_dispDofNum[0], node[0]);
                        _dof[2][1] = analysisModel().domainManager().giveNodalDof(_dispDofNum[1], node[0]);
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
                    
                    if ( j == 0 )
                    {
                        // Store DOFs
                        _dof[1][0] = analysisModel().domainManager().giveNodalDof(_dispDofNum[0], node[2]);
                        _dof[1][1] = analysisModel().domainManager().giveNodalDof(_dispDofNum[1], node[2]);
                        _dof[3][0] = analysisModel().domainManager().giveNodalDof(_dispDofNum[0], node[1]);
                        _dof[3][1] = analysisModel().domainManager().giveNodalDof(_dispDofNum[1], node[1]);
                    }
                    else // (j == 1)
                    {
                        // Store DOFs
                        _dof[1][0] = analysisModel().domainManager().giveNodalDof(_dispDofNum[0], node[2]);
                        _dof[1][1] = analysisModel().domainManager().giveNodalDof(_dispDofNum[1], node[2]);
                        _dof[3][0] = analysisModel().domainManager().giveNodalDof(_dispDofNum[0], node[0]);
                        _dof[3][1] = analysisModel().domainManager().giveNodalDof(_dispDofNum[1], node[0]);
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