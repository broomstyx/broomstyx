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

#include "DomainManager.hpp"
#include <chrono>
#include <stdexcept>
#include <omp.h>

#include "AnalysisModel.hpp"
#include "Cell.hpp"
#include "DofManager.hpp"
#include "DomainManager.hpp"
#include "MaterialManager.hpp"
#include "Node.hpp"
#include "NumericsManager.hpp"
#include "SolutionManager.hpp"
#include <MeshReaders/MeshReader.hpp>
#include <Numerics/Numerics.hpp>
#include <Util/readOperations.hpp>

using namespace broomstyx;

// Constructor
DomainManager::DomainManager()
{
    _fieldsPerNode = -1;
    _constructFaces = false;
}

// Destructor
DomainManager::~DomainManager()
{
#ifdef VERBOSE_DESTRUCTION
    std::printf("\n  Destroying DomainManager... ");
    std::fflush(stdout);
#endif
    
    for ( auto curNode : _nodeList )
    {
        analysisModel().dofManager().destroyNodalDofsAt(curNode);
        delete curNode;
    }
    
    for ( auto curFace : _faceList )
    {
        analysisModel().dofManager().destroyFaceDofsAt(curFace);
        delete curFace;
    }
    
    for ( auto curDomCell : _domCellList )
    {
        analysisModel().dofManager().destroyCellDofsAt(curDomCell);
        this->giveNumericsFor(curDomCell)->deleteNumericsAt(curDomCell);
        
        delete curDomCell;
    }
    
    for ( auto curBndCell : _bndCellList )
        delete curBndCell;
    
#ifdef VERBOSE_DESTRUCTION
    std::printf("done.");
    std::fflush(stdout);
#endif
}

// Public methods
// ----------------------------------------------------------------------------
void DomainManager::createPhysicalEntity( int dim, int number, std::string label )
{
    PhysicalEntity newPhysEnt;
    
    newPhysEnt.dimension = dim;
    newPhysEnt.entityNumber = number;
    newPhysEnt.name = label;
    
    _physEnt.push_back(newPhysEnt);
}
// ----------------------------------------------------------------------------
DomainManager::PhysicalEntity DomainManager::giveDataForPhysicalEntity( int n )
{
    return _physEnt[n];
}
// ----------------------------------------------------------------------------
std::vector<Material*> DomainManager::giveMaterialSetForDomain( int label )
{
    std::string name = this->givePhysicalEntityNameFor(label);
    return _materialSet[name];
}
// ----------------------------------------------------------------------------
int DomainManager::giveNumberOfPhysicalNames()
{
    return _physEnt.size();
}
// ----------------------------------------------------------------------------
Numerics* DomainManager::giveNumericsForDomain( int label )
{
    std::string name = this->givePhysicalEntityNameFor(label);
    return _numerics[name];
}
// ----------------------------------------------------------------------------
std::string DomainManager::givePhysicalEntityNameFor( int physEntNum )
{
    std::string physEntName;
    bool success = false;
    
    for ( int i = 0; i < (int)_physEnt.size(); i++)
    {
        if ( physEntNum == _physEnt[i].entityNumber )
        {
            physEntName = _physEnt[i].name;
            success = true;
            break;
        }
    }
    
    if ( !success )
        throw std::runtime_error("Failed to find name corresponding to physical entity number '"
                + std::to_string(physEntNum) + "'!\nSource: DomainManager");
    
    return physEntName;
}
// ----------------------------------------------------------------------------
int DomainManager::givePhysicalEntityNumberFor( std::string name )
{
    int physEntNumber;
    bool success = false;
    
    for ( int i = 0; i < (int)_physEnt.size(); i++)
    {
        if ( name == _physEnt[i].name )
        {
            physEntNumber = _physEnt[i].entityNumber;
            success = true;
            break;
        }
    }
    
    if ( !success )
        throw std::runtime_error("Failed to find physical entity number corresponding to '" + name + "'!\nSource: DomainManager");
    
    return physEntNumber;
}
// ----------------------------------------------------------------------------
void DomainManager::readDomainAssignmentsFrom( FILE* fp )
{
    std::string key, src = "DomainManager";
    
    int nAssign = getIntegerInputFrom(fp, "Failed to read number of domain assignments from input file!", src);
    
    for ( int i = 0; i < nAssign; i++ )
    {
        // Read label
        std::string domainLabel = getStringInputFrom(fp, "Failed reading domain label from input file!", src);
                
        // Read numerics
        verifyKeyword(fp, key = "Numerics", src);
        int numericsLabel = getIntegerInputFrom(fp, "Failed reading numerics label from input file!", src);
        
        Numerics* numericsPtr = analysisModel().numericsManager().giveNumerics(numericsLabel);
        
        // Create map entry
        std::pair< std::map<std::string, Numerics*>::iterator, bool> entry;
        entry = _numerics.insert(std::pair<std::string, Numerics*>(domainLabel, numericsPtr));
        if ( !entry.second )
            throw std::runtime_error("Multiple declaration of numerics for label '" + domainLabel + "' detected in input file!\nSource: " + src);
        
        // Read material set
        int nMat = numericsPtr->requiredNumberOfMaterials();
        if ( nMat > 0 )
        {
            verifyKeyword(fp, key = "MaterialSet", src);
        
            std::vector<Material*> matSet;
            matSet.assign(nMat, nullptr);
            for ( int j = 0; j < nMat; j++ )
            {
                int matLabel = getIntegerInputFrom(fp, "Failed to read material label from input file.", src);
                matSet[j] = analysisModel().materialManager().giveMaterial(matLabel);
            }

            // Create map entry
            std::pair<std::map<std::string, std::vector<Material*> >::iterator,bool> tmp;
            tmp = _materialSet.insert(std::pair<std::string, std::vector<Material*> >(domainLabel, matSet));
            if ( !tmp.second )
                throw std::runtime_error("Multiple declaration of material sets for label '" + domainLabel + "' detected in input file!\nSource: " + src);
        }
    }
}
// ----------------------------------------------------------------------------
void DomainManager::countNodes()
{
    // Determine number of active nodes
    int nActiveNodes = 0;
    for ( auto curNode = _nodeList.begin(); curNode != _nodeList.end(); curNode++ )
        if ( (*curNode)->_isActive )
            (*curNode)->_id = nActiveNodes++;
    
    // Rewrite active node addresses into one contiguous array
    _node.assign(nActiveNodes, nullptr);
    int curCount = 0;
    for ( auto curNode = _nodeList.begin(); curNode != _nodeList.end(); curNode++ )
        if ( (*curNode)->_isActive )
            if ( (*curNode)->_isActive )
                _node[curCount++] = *curNode;
}
// ----------------------------------------------------------------------------
std::set<Cell*> DomainManager::giveAttachedDomainCellsOf( Node* node )
{
    return node->_attachedDomCell;
}
// ----------------------------------------------------------------------------
RealVector DomainManager::giveCoordinatesOf( Node* node )
{
    return node->_coordinates;
}
// ----------------------------------------------------------------------------
Dof* DomainManager::giveNodalDof( int dofNum, Node* node )
{
    return node->_dof[ dofNum ];
}
// ----------------------------------------------------------------------------
double DomainManager::giveFieldValueAt( Node* node, int fieldNum )
{
    if ( fieldNum == 0 )
        return 0.0;
    else
        return node->_fieldVal( fieldNum - 1 );
}
// ----------------------------------------------------------------------------
int DomainManager::giveIdOf( Node* node )
{
    return node->_id;
}
// ----------------------------------------------------------------------------
Node* DomainManager::giveNode( int nodeNum ) 
{
    return _node[nodeNum];
}
// ----------------------------------------------------------------------------
int DomainManager::giveNumberOfNodes()
{
    return _node.size();
}
// ----------------------------------------------------------------------------
void DomainManager::makeNewNodeAt( RealVector& location )
{
    std::string errmsg, src = "DomainManager";
    if ( _fieldsPerNode == -1 )
    {
        errmsg = "Cannot create new node due to undefined number of fields" + std::string(" per node!\nSource: ") + src;
        throw std::runtime_error(errmsg);
    }
        
    Node* newNode = new Node();
    
    if ( location.dim() == 2 )
    {
        newNode->_coordinates.init(3);
        newNode->_coordinates(0) = location(0);
        newNode->_coordinates(1) = location(1);
        newNode->_coordinates(2) = 0.0;
    }
    else if ( location.dim() == 3 )
        newNode->_coordinates = location;
    else
        throw std::runtime_error("Invalid size of vector input for nodal coordinates!\nSource: " + src);
    
    // Instantiate degrees of freedom for new node
    analysisModel().dofManager().createNodalDofsAt(newNode);
    
    // Initialize nodal fields
    if ( _fieldsPerNode > 0)
        newNode->_fieldVal.init( _fieldsPerNode );
    
    // Append node to list of nodes
    _nodeList.push_back(newNode);
}
// ----------------------------------------------------------------------------
#ifdef USING_DUNE_GRID_BACKEND
void DomainManager::makeNewNodeFrom( VertexSeedType seed )
{
    std::string errmsg, src = "DomainManager";
    if ( _fieldsPerNode == -1 )
    {
        errmsg = "Cannot create new node due to undefined number of fields" + std::string(" per node!\nSource: ") + src;
        throw std::runtime_error(errmsg);
    }
    
    Node* newNode = new Node();
    newNode->_vertexSeed = seed;

    // Instantiate degrees of freedom for new node
    analysisModel().dofManager().createNodalDofsAt(newNode);
    
    // Initialize nodal fields
    if ( _fieldsPerNode > 0)
        newNode->_fieldVal.init( _fieldsPerNode );
    
    // Append node to list of nodes
    _nodeList.push_back(newNode);
}
#endif
// ----------------------------------------------------------------------------
void DomainManager::performNodalPostProcessing()
{
    std::chrono::time_point<std::chrono::system_clock> tic, toc;
    std::chrono::duration<double> tictoc;
    
    std::printf("    %-40s", "Performing nodal post-processing ...");
    std::fflush(stdout);
    tic = std::chrono::high_resolution_clock::now();
    
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < (int)_node.size(); i++ )
    {
        // Get field values corresponding to DOF variables
        analysisModel().dofManager().writeConvergedDofValuesTo(_node[i]);
        
        // Post-process nodal field values from element data
//        analysisModel->solutionManager->recoverVariablesAt(_node[i]);
    }
    
    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    std::printf("done (time = %f sec.)\n", tictoc.count());
}
// ----------------------------------------------------------------------------
void DomainManager::readNumberOfFieldsPerNodeFrom( FILE* fp )
{
    std::string str, errmsg, src = "DomainManager";
    _fieldsPerNode = getIntegerInputFrom(fp, "\nFailed to read number of fields per node from input file", src);
}
// ----------------------------------------------------------------------------
void DomainManager::mustConstructFaces()
{
    _constructFaces = true;
}
// ----------------------------------------------------------------------------
void DomainManager::setFieldValueAt( Node* targetNode, int fieldNum, double val )
{
    targetNode->_fieldVal( fieldNum - 1 ) = val;
}

// Methods for cell access
// ----------------------------------------------------------------------------
void DomainManager::constructCellFaces()
{
    if ( _constructFaces )
    {
        std::chrono::time_point<std::chrono::system_clock> tic, toc;
        std::chrono::duration<double> tictoc;

        std::printf("  %-40s", "Constructing cell faces ...");
        std::fflush(stdout);
        tic = std::chrono::high_resolution_clock::now();
    
        // Initialize vectors to hold pointers to faces and their orientations
        // in each domain cell

        int nDomCells = _domCell.size();

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for ( int i = 0; i < nDomCells; i++ )
        {
            Cell* curCell = _domCell[i];
            int nNeighbors = curCell->_neighbor.size();
            curCell->_face.assign(nNeighbors, nullptr);
            curCell->_faceOrient.assign(nNeighbors, 0);
        }

        // Actual construction of faces

        // Note that this method requires domain cells to already have been counted,
        // and the cell-node/cell-cell connectivities to have been determined. 
        // It cannot be parallelized since the face cells are constructed on the
        // fly and appended to a linked list.

        if ( !_faceList.empty() )
            _faceList.clear();

        for ( int i = 0; i < nDomCells; i++ )
        {
            Cell* curCell = _domCell[i];
            int nFaces = analysisModel().meshReader().giveNumberOfFacesForElementType(curCell->_elType);

            // Sanity check
            if ( (int)curCell->_face.size() != nFaces )
                throw std::runtime_error("ERROR: Mismatch encountered in number of faces for element type!\nSource:DomainManager");

            for ( int j = 0; j < nFaces; j++ )
            {
                if ( curCell->_faceOrient[j] == 0 )
                {
                    // Create cell object representing face
                    this->makeNewFaceBetween(curCell, curCell->_neighbor[j], j);
                }
            }
        }
        this->countFaces();
        
        toc = std::chrono::high_resolution_clock::now();
        tictoc = toc - tic;
        std::printf("done (time = %f sec.)\n", tictoc.count());
        
        std::printf("\n    Domain cell faces = %d\n\n", (int)_faceList.size());
    }
}
// ----------------------------------------------------------------------------
void DomainManager::countBoundaryCells()
{
    _bndCell.assign(_bndCellList.size(), nullptr);
    int curCount = 0;
    for ( auto it = _bndCellList.begin(); it != _bndCellList.end(); it++ )
    {
        _bndCell[curCount] = *it;
        (*it)->_id = curCount++;
    }
}
// ----------------------------------------------------------------------------
void DomainManager::countDomainCells()
{
    _domCell.assign(_domCellList.size(), nullptr);
    int curCount = 0;
    for ( auto it = _domCellList.begin(); it != _domCellList.end(); it++ )
    {
        _domCell[curCount] = *it;
        (*it)->_id = curCount++;
    }
}
// ----------------------------------------------------------------------------
void DomainManager::countFaces()
{
    _face.assign(_faceList.size(), nullptr);
    int curCount = 0;
    for ( auto it = _faceList.begin(); it != _faceList.end(); it++ )
    {
        _face[curCount] = *it;
        (*it)->_id = curCount++;
    }
}
// ----------------------------------------------------------------------------
void DomainManager::finalizeCellData()
{
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < (int)_domCell.size(); i++)
    {
        Numerics* numerics = analysisModel().domainManager().giveNumericsForDomain(_domCell[i]->_label);
        numerics->finalizeDataAt(_domCell[i]);
    }
}
// ----------------------------------------------------------------------------
void DomainManager::findBoundaryAssociations()
{
    std::chrono::time_point<std::chrono::system_clock> tic, toc;
    std::chrono::duration<double> tictoc;
    
    std::printf("  %-40s", "Finding boundary associations ...");
    std::fflush(stdout);
    tic = std::chrono::high_resolution_clock::now();

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < (int)_bndCell.size(); i++ )
        this->findDomainCellsAssociatedWith(_bndCell[i]);

    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    std::printf("done (time = %f sec.)\n", tictoc.count());
}
// ----------------------------------------------------------------------------
void DomainManager::findDomainCellNeighbors()
{
    std::chrono::time_point<std::chrono::system_clock> tic, toc;
    std::chrono::duration<double> tictoc;
    
    std::printf("  %-40s", "Finding domain cell neighbors ...");
    std::fflush(stdout);
    tic = std::chrono::high_resolution_clock::now();
    
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < (int)_domCell.size(); i++ )
    {
        this->findNeighborsOf(_domCell[i]);
    }

    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    std::printf("done (time = %f sec.)\n", tictoc.count());
}
// ----------------------------------------------------------------------------
void DomainManager::findDomainCellsAssociatedWith( Cell *targetCell )
{    
    if ( targetCell->_isPartOfDomain )
        throw std::runtime_error("ERROR: Search request made for associated domain cell of a cell\nthat is already part of the domain!\nSource: DomainManager");
    
    // Form set of all elements attached to all the nodes of the 
    // current boundary element
    std::set<Cell*> candCell;
    int nNodes = targetCell->_node.size();
    
    for (  int i = 0; i < nNodes; i++ )
    {
        std::set<Cell*> temp = targetCell->_node[i]->_attachedDomCell;
        
        std::set<Cell*>::iterator it;
        for ( it = temp.begin(); it != temp.end(); it++ )
            candCell.insert(*it);
    }
    
    std::vector<bool> hasNode;
    std::set<Cell*>::iterator it;
    it = candCell.begin();
    
    for ( it = candCell.begin(); it != candCell.end(); ++it )
    {
        std::vector<Node*> ccNode = (*it)->_node;
        
        hasNode.assign(nNodes, false);
        bool found = true;
        for ( int i = 0; i < nNodes; i++ )
        {
            for ( int j = 0; j < (int)ccNode.size(); j++ )
                if ( targetCell->_node[i] == ccNode[j] )
                    hasNode[i] = true;
            
            if ( hasNode[i] == false )
                found = false;
        }
        
        if (found)
            targetCell->_assocDomCell.push_back(*it);
    }

    if ( targetCell->_assocDomCell.size() == 0 )
        throw std::runtime_error("ERROR: Failed finding domain cell association!\nSource: DomainManager");
}
// ----------------------------------------------------------------------------
void DomainManager::findNeighborsOf( Cell* targetCell )
{
    int nFaces = analysisModel().meshReader().giveNumberOfFacesForElementType(targetCell->_elType);
    targetCell->_neighbor.assign(nFaces, nullptr);

    for ( int j = 0; j < nFaces; j++ )
    {
        std::vector<int> faceNodes = analysisModel().meshReader().giveFaceNodeNumbersForElementType(targetCell->_elType, j);
        std::set<Cell*> candidate = analysisModel().domainManager().giveAttachedDomainCellsOf(targetCell->_node[faceNodes[0]]);
        for ( auto it = candidate.begin(); it != candidate.end(); ++it )
        {
            std::vector<Node*> adjNode = analysisModel().domainManager().giveNodesOf(*it);

            bool isNeighbor = true;
            for ( int k = 1; k < (int)faceNodes.size(); k++ )
            {
                bool hasFaceNode = false;
                for ( int m = 0; m < (int)adjNode.size(); m++ )
                    if ( adjNode[m] == targetCell->_node[faceNodes[k]] )
                        hasFaceNode = true;

                if ( !hasFaceNode )
                    isNeighbor = false;
            }

            if ( isNeighbor && (*it) != targetCell )
                targetCell->_neighbor[j] = *it;
        }
    }
}
// ----------------------------------------------------------------------------
void DomainManager::formDomainPartitions()
{
    // Determine number of partitions
    int nPartitions = 0;
    for ( int i = 0; i < (int)_domCell.size(); i++)
        if ( _domCell[i]->_partition > nPartitions )
            nPartitions = _domCell[i]->_partition;
    
    // Count number of cells in each partition
    std::vector<int> partitionCellCount;
    partitionCellCount.assign(nPartitions + 1, 0);
    for ( int i = 0; i < (int)_domCell.size(); i++)
        partitionCellCount[_domCell[i]->_partition] += 1;
    
    // Create cell address vector for each partition
    _partition.assign(nPartitions + 1, std::vector<Cell*>());
    for ( int i = 0; i <= nPartitions; i++)
    {
        _partition[i].assign(partitionCellCount[i], nullptr);
        partitionCellCount[i] = 0;
    }
    
    for ( int i = 0; i < (int)_domCell.size(); i++)
    {
        int cellPartition = _domCell[i]->_partition;
        _partition[cellPartition][partitionCellCount[cellPartition]++] = _domCell[i];
    }
}
// ----------------------------------------------------------------------------
Cell* DomainManager::giveBoundaryCell( int cellNum )
{
    return _bndCell[ cellNum ];
}
// ----------------------------------------------------------------------------
Dof* DomainManager::giveCellDof( int dofNum, Cell *targetCell )
{
    return targetCell->_dof[ dofNum ];
}
// ----------------------------------------------------------------------------
Cell* DomainManager::giveDomainCell( int cellNum )
{
    return _domCell[ cellNum ];
}
// ----------------------------------------------------------------------------
Cell* DomainManager::giveDomainCellInPartition( int partNum, int cellNum )
{
    return _partition[partNum][cellNum];
}
// ----------------------------------------------------------------------------
std::vector<Cell*> DomainManager::giveDomainCellsAssociatedWith( Cell *targetCell )
{
    return targetCell->_assocDomCell;
}
// ----------------------------------------------------------------------------
int DomainManager::giveElementTypeOf( Cell* targetCell )
{
    return targetCell->_elType;
}
// ----------------------------------------------------------------------------
int DomainManager::giveIdOf( Cell *targetCell )
{
    return targetCell->_id;
}
// ----------------------------------------------------------------------------
int DomainManager::giveLabelOf( Cell *targetCell )
{
    return targetCell->_label;
}
// ----------------------------------------------------------------------------
std::vector<Cell*> DomainManager::giveNeighborsOf(Cell *targetCell )
{
    return targetCell->_neighbor;
}
// ----------------------------------------------------------------------------
std::vector<Node*> DomainManager::giveNodesOf( Cell *targetCell )
{
    return targetCell->_node;
}
// ----------------------------------------------------------------------------
int DomainManager::giveNumberOfBoundaryCells()
{
    return _bndCell.size();
}
// ----------------------------------------------------------------------------
int DomainManager::giveNumberOfDomainCells()
{ 
    return _domCell.size();
}
// ----------------------------------------------------------------------------
int DomainManager::giveNumberOfDomainCellsInPartition( int partNum )
{
    return _partition[partNum].size();
}
// ----------------------------------------------------------------------------
int DomainManager::giveNumberOfNodesOf( Cell *targetCell )
{
    return targetCell->_node.size();
}
// ----------------------------------------------------------------------------
int DomainManager::giveNumberOfPartitions()
{
    return _partition.size();
}
// ----------------------------------------------------------------------------
Numerics* DomainManager::giveNumericsFor(Cell* targetCell)
{
    if ( targetCell->_isPartOfDomain )
        return this->giveNumericsForDomain(targetCell->_label);
    else
        return nullptr;
}
// ----------------------------------------------------------------------------
void DomainManager::initializeMaterialsAtCells()
{
    std::chrono::time_point<std::chrono::system_clock> tic, toc;
    std::chrono::duration<double> tictoc;
    
    std::printf("  %-40s", "Initializing material data at cells ...");
    std::fflush(stdout);
    tic = std::chrono::high_resolution_clock::now();

#ifdef _OPENMP    
#pragma omp parallel for
#endif
    for (int i = 0; i < (int)_domCell.size(); i++)
    {
        Numerics *numerics = analysisModel().domainManager().giveNumericsForDomain(_domCell[i]->_label);
        numerics->initializeMaterialsAt(_domCell[i]);
    }    

    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    std::printf("done (time = %f sec.)\n", tictoc.count());
}
// ----------------------------------------------------------------------------
void DomainManager::initializeNumericsAtCells()
{
    std::chrono::time_point<std::chrono::system_clock> tic, toc;
    std::chrono::duration<double> tictoc;
    
    std::printf("  %-40s", "Initializing numerics at cells ...");
    std::fflush(stdout);
    tic = std::chrono::high_resolution_clock::now();
    
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < (int)_domCell.size(); i++)
    {
        Numerics* numerics = this->giveNumericsFor(_domCell[i]);
        numerics->initializeNumericsAt( _domCell[i] );
    }
    
    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    std::printf("done (time = %f sec.)\n", tictoc.count());
}
// ----------------------------------------------------------------------------
#ifdef USING_DUNE_GRID_BACKEND
Cell* DomainManager::makeNewBoundaryCellFrom( DomainSeedType seed, int cellLabel )
{
    // Instantiate new cells
    Cell* newCell = new Cell();
    newCell->_domainSeed = seed;
    newCell->_label = cellLabel;

    // Retrieve numerics type based on cell label
    Numerics* numerics = this->giveNumericsForDomain(cellLabel);
    
    // Add new cell object to relevant list
    if ( numerics )
    {
        newCell->_isPartOfDomain = true;
        _domCellList.push_back(newCell);
        analysisModel().dofManager().createCellDofsAt(newCell);
        if ( _fieldsPerCell > 0 )
            newCell->cellData.init(_fieldsPerCell);
    }
    else
    {
        newCell->_isPartOfDomain = false;
        _bndCellList.push_back(newCell);
    }

    return newCell;
}
// ----------------------------------------------------------------------------
Cell* DomainManager::makeNewDomainCellFrom( DomainSeedType seed, int cellLabel )
{
    // Instantiate new cells
    Cell* newCell = new Cell();
    newCell->_domainSeed = seed;
    newCell->_label = cellLabel;

    // Retrieve numerics type based on cell label
    Numerics* numerics = this->giveNumericsForDomain(cellLabel);
    
    // Add new cell object to relevant list
    if ( numerics )
    {
        newCell->_isPartOfDomain = true;
        _domCellList.push_back(newCell);
        analysisModel().dofManager().createCellDofsAt(newCell);
        if ( _fieldsPerCell > 0 )
            newCell->cellData.init(_fieldsPerCell);
    }
    else
    {
        newCell->_isPartOfDomain = false;
        _bndCellList.push_back(newCell);
    }

    return newCell;
}
#endif
// ----------------------------------------------------------------------------
Cell* DomainManager::makeNewCellWithLabel( int cellLabel )
{
    // Instantiate new cells
    Cell* newCell = new Cell();
    newCell->_label = cellLabel;
    
    // Retrieve numerics type based on cell label
    Numerics* numerics = this->giveNumericsForDomain(cellLabel);
    
    // Add new cell object to relevant list
    if ( numerics )
    {
        newCell->_isPartOfDomain = true;
        _domCellList.push_back(newCell);
        analysisModel().dofManager().createCellDofsAt(newCell);
        if ( _fieldsPerCell > 0 )
            newCell->cellData.init(_fieldsPerCell);
    }
    else
    {
        newCell->_isPartOfDomain = false;
        _bndCellList.push_back(newCell);
    }

    return newCell;
}
// ----------------------------------------------------------------------------
void DomainManager::makeNewFaceBetween( Cell* posCell, Cell* negCell, int posFaceNum )
{
    Cell* newFace = new Cell();
    _faceList.push_back(newFace);
    
    // Store face nodes
    std::vector<int> faceNodeNum = analysisModel().meshReader().giveFaceNodeNumbersForElementType(posCell->_elType, posFaceNum);
    int nFaceNodes = faceNodeNum.size();
    newFace->_node.assign(nFaceNodes, nullptr);
    for ( int i = 0; i < nFaceNodes; i++ )
        newFace->_node[i] = posCell->_node[faceNodeNum[i]];
    
    // Set face orientations
    posCell->_faceOrient[posFaceNum] = 1;
    
    // Store adjoining cells (Note: negCell may be 'nullptr' if posCell is at the boundary)
    newFace->_neighbor.assign({posCell, negCell});
    
    // Register new face on neighbor cell with negative orientation
    if ( negCell )
    {
        int nFaces = analysisModel().meshReader().giveNumberOfFacesForElementType(negCell->_elType);
        bool assigned = false;
        for ( int i = 0; i < nFaces; i++ )
            if ( negCell->_neighbor[i] == posCell )
            {
                negCell->_face[i] = newFace;
                negCell->_faceOrient[i] = -1;
                assigned = true;
            }
        
        // Sanity check
        if ( !assigned )
            throw std::runtime_error("ERROR: Unable to assign face in neighbor cell with negative orientation!\nSource: DomainManager");
    }
    
    // Initialize face fields
    newFace->cellData.init(_fieldsPerFace);
}
// ----------------------------------------------------------------------------
void DomainManager::readNumberOfFieldsPerCellFrom( FILE* fp )
{
    _fieldsPerCell = getIntegerInputFrom(fp, "\nFailed to read number of fields per cell in input file!", "DomainManager");
}
// ----------------------------------------------------------------------------
void DomainManager::readNumberOfFieldsPerFaceFrom( FILE* fp )
{
    _fieldsPerFace = getIntegerInputFrom(fp, "\nFailed to read number of fields per face in input file!", "DomainManager");
    
    if ( _fieldsPerFace > 0 )
        this->_constructFaces = true;
}
// ----------------------------------------------------------------------------
void DomainManager::removeAllCellConstraints()
{
    for (int i = 0; i < (int)_domCell.size(); i++)
    {
        Numerics* numerics = analysisModel().domainManager().giveNumericsForDomain(_domCell[i]->_label);
        numerics->removeConstraintsOn(_domCell[i]);
    }
}
// ----------------------------------------------------------------------------
void DomainManager::reportDetailedStatus()
{
    this->reportStatus();
    
    std::printf("\n");
    for ( int i = 0; i < (int)this->_node.size(); i++ )
    {
        RealVector coor = _node[i]->_coordinates;
        std::printf("    node %d: x = %e, y = %e, z = %e, dofs = ", i, coor(0), coor(1), coor(2));
        for ( int j = 0; j < (int)_node[i]->_dof.size(); j++ )
            std::printf("%d ", analysisModel().dofManager().giveEquationNumberAt(_node[i]->_dof[j]));
        std::printf("\n");
    }
    
    for ( int i = 0; i < (int)this->_domCell.size(); i++ )
    {
        std::vector<Node*> node = _domCell[i]->_node;
        std::printf("    cell %d: nodes = ", i);
        for ( int j = 0; j < (int)node.size(); j++ )
            std::printf("%d ", node[j]->_id);
        std::printf(", dofs = ");
        for ( int j = 0; j < (int)_domCell[i]->_dof.size(); j++ )
            std::printf("%d ", analysisModel().dofManager().giveEquationNumberAt(_domCell[i]->_dof[j]));
        std::printf(", neighbors = ");
        for ( int j = 0; j < (int)_domCell[i]->_neighbor.size(); j++)
            if ( _domCell[i]->_neighbor[j] )
                std::printf("%d ", _domCell[i]->_neighbor[j]->_id);
        std::printf("\n");
    }
}
// ----------------------------------------------------------------------------
void DomainManager::reportStatus()
{
    std::printf("    Nodes          = %ld\n", _node.size() );
    std::printf("    Domain cells   = %ld\n", _domCell.size());
    std::printf("    Boundary cells = %ld\n\n", _bndCell.size());
    
    std::printf("    Number of domain partitions = %ld\n", _partition.size());
    for ( int i = 0; i < (int)_partition.size(); i++ )
        std::printf("      Cells in partition # %2d = %ld\n", i, _partition[i].size());
}
// ----------------------------------------------------------------------------
void DomainManager::setElementTypeOf( Cell* targetCell, int elemType )
{
    targetCell->_elType = elemType;
}
// ----------------------------------------------------------------------------
void DomainManager::setHaloOf(Cell *targetCell, std::vector<int>& halo )
{
    targetCell->_halo = halo;
}
// ----------------------------------------------------------------------------
void DomainManager::setNeighborsOf( Cell *targetCell, std::vector<Cell*>& neighbors )
{
    targetCell->_neighbor = neighbors;
}
// ----------------------------------------------------------------------------
void DomainManager::setNodesOf( Cell *targetCell, std::vector<int>& cellNodes )
{
    // Important: this method must not be called from within a loop that is
    // parallelized using OpenMP
    
    std::vector<Node*> node(cellNodes.size(), nullptr);
    
    for ( int i = 0; i < (int)cellNodes.size(); i++ ) 
    {
        node[i] = _node[cellNodes[i]];
        if ( targetCell->_isPartOfDomain )
            node[i]->_attachedDomCell.insert(targetCell);
        else
            node[i]->_attachedBndCell.insert(targetCell);
    }
    
    targetCell->_node = node;
}
// ----------------------------------------------------------------------------
#ifdef USING_DUNE_GRID_BACKEND
void DomainManager::setNodesOf( Cell* targetCell, std::vector<VertexSeedType>& vertices )
{
    // Important: this method must not be called from within a loop that is
    // parallelized using OpenMP

    std::vector<Node*> node(vertices.size(), nullptr);

    for ( int i = 0; i < (int)_node.size(); i++ )
    {
        VertexSeedType nodeSeed = _node[i]->_vertexSeed;
        int found = 0;
        for ( int j = 0; j < (int)vertices.size(); j++ )
        if ( nodeSeed == vertices[j] )
        {
            node[j] = _node[i];
            ++found;
        }
        
        if ( found == (int)vertices.size() )
            break;
    }

    for ( int i = 0; i < (int)node.size(); i++ ) 
    {
        if ( targetCell->_isPartOfDomain )
            node[i]->_attachedDomCell.insert(targetCell);
        else
            node[i]->_attachedBndCell.insert(targetCell);
    }
    
    targetCell->_node = node;
}
#endif
// ----------------------------------------------------------------------------
void DomainManager::setPartitionOf( Cell* targetCell, int partition )
{
    targetCell->_partition = partition;
}

// Methods helpful in debugging
// ----------------------------------------------------------------------------
void DomainManager::printInfoFor( Cell* targetCell )
{
    int cellId = targetCell->_id;
    int cellLabel = targetCell->_label;

    std::printf("\n  Cell ID = %d, label = %d, dim = %d\n", targetCell->_id, targetCell->_label, targetCell->_dim);
    std::printf("  Cell nodes: ");
    
    for ( int i = 0; i < (int)targetCell->_node.size(); ++i )
        std::printf("%d ", targetCell->_node[i]->_id);
    
    std::printf("\n  Neighbors: ");
    for ( int i = 0; i < (int)targetCell->_neighbor.size(); ++i )
        if ( targetCell->_neighbor[i] )
            std::printf("%d ", targetCell->_neighbor[i]->_id);
        else
            std::printf("null ");
    std::printf("\n");
    
    if ( !targetCell->_isPartOfDomain )
    {
        std::printf("  Number of associated domain cells = %d\n", (int)targetCell->_assocDomCell.size());
        std::printf("  Associated domain cell = ");
        for ( auto it = targetCell->_assocDomCell.begin(); it != targetCell->_assocDomCell.end(); ++it )
            std::printf("%d ", (*it)->_id);
        std::printf("\n");
    }
}