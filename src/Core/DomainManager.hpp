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

#ifndef DOMAINMANAGER_HPP
#define	DOMAINMANAGER_HPP

#include <cstdio>
#include <map>
#include <set>
#include <vector>
#include <list>

#include "TimeData.hpp"
#include "Util/RealVector.hpp"

namespace broomstyx
{
    class Cell;
    class Dof;
    class Material;
    class Node;
    class Numerics;

    class DomainManager final
    {
        friend class AnalysisModel;
        
    public:
        // Disable copy constructor and assignment operator
        DomainManager( const DomainManager& ) = delete;
        DomainManager& operator=( const DomainManager& ) = delete;
        
        struct PhysicalEntity
        {
            int dimension;
            int entityNumber;
            std::string name;
        };
        
        void                   createPhysicalEntity( int dim, int number, std::string label );
        PhysicalEntity         giveDataForPhysicalEntity( int n );
        std::vector<Material*> giveMaterialSetForDomain( int label );
        int                    giveNumberOfPhysicalNames();
        Numerics*              giveNumericsForDomain( int label );
        std::string            givePhysicalEntityNameFor( int physEntNum );
        int                    givePhysicalEntityNumberFor( std::string name );
        void                   readDomainAssignmentsFrom( FILE* fp );
        
        // Methods involving node access
        
        void   countNodes();
        std::set<Cell*> giveAttachedDomainCellsOf( Node* node );
        RealVector      giveCoordinatesOf( Node* node );
        Dof*   giveNodalDof( int dofNum, Node* node );
        double giveFieldValueAt( Node* node, int fieldNum );
        int    giveIdOf( Node* node );
        Node*  giveNode( int nodeNum );
        int    giveNumberOfNodes();
        void   makeNewNodeAt( RealVector& location );
        void   performNodalPostProcessing();
        void   readNumberOfFieldsPerNodeFrom( FILE* fp );
        void   setCoordinatesOf( Node* targetNode, const RealVector& coor );
        void   setFieldValueAt( Node* targetNode, int fieldNum, double val );
        
        // Methods involving cell access
        
        void   constructCellFaces();
        void   countBoundaryCells();
        void   countDomainCells();
        void   countFaces();
        void   finalizeCellDataAt( const TimeData& time );
        void   findBoundaryAssociations();
        void   findDomainCellNeighbors();
        void   findDomainCellsAssociatedWith( Cell* targetCell );
        void   findNeighborsOf( Cell *targetCell );
        void   formDomainPartitions();
        Cell*  giveBoundaryCell( int cellNum );
        Dof*   giveCellDof( int dofNum, Cell *targetCell );
        Cell*  giveDomainCell( int cellNum );
        std::vector<Cell*> 
               giveDomainCellsAssociatedWith( Cell *targetCell );
        Cell*  giveDomainCellInPartition( int partNum, int cellNum );
        int    giveElementTypeOf( Cell* targetCell );
        int    giveIdOf( Cell *targetCell );
        int    giveLabelOf( Cell *targetCell );
        
        std::vector<Cell*> giveNeighborsOf( Cell *targetCell );    
        std::vector<Node*> giveNodesOf( Cell *targetCell );
        
        int   giveNumberOfBoundaryCells();
        int   giveNumberOfDomainCells();
        int   giveNumberOfDomainCellsInPartition( int partNum );
        int   giveNumberOfNodesOf( Cell *targetCell );
        int   giveNumberOfPartitions();
        Numerics* giveNumericsFor( Cell* targetCell );
        void  initializeMaterialsAtCells();
        void  initializeNumericsAtCells();
        Cell* makeNewCellWithLabel( int cellLabel );
        void  makeNewFaceBetween( Cell* posCell, Cell* negCell, int posFaceNum );
        void  mustConstructFaces();
        void  readNumberOfFieldsPerCellFrom( FILE* fp );
        void  readNumberOfFieldsPerFaceFrom( FILE* fp );
        void  removeAllCellConstraints();
        void  reorderNodesOf( Cell* targetCell, std::vector<int>& reordering );
        void  reportDetailedStatus();
        void  reportStatus();
        void  setElementTypeOf( Cell* targetCell, int elemType );
        void  setHaloOf( Cell *targetCell, std::vector<int>& halo );
        void  setNeighborsOf( Cell *targetCell, std::vector<Cell*>& neighbors);
        void  setNodesOf( Cell *targetCell, std::vector<int>& cellNodes );
        void  setPartitionOf( Cell *targetCell, int partition );

    private:
        std::vector<PhysicalEntity> _physEnt;
        
        std::map<std::string, Numerics*> _numerics;
        std::map<std::string, std::vector<Material*> > _materialSet;
        
        int _fieldsPerNode;
        std::list<Node*> _nodeList;
        std::vector<Node*> _node;
                
        int _fieldsPerCell;
        std::list<Cell*> _bndCellList;
        std::vector<Cell*> _bndCell;
        
        std::list<Cell*> _domCellList;
        std::vector<Cell*> _domCell;
        
        // Cells representing faces of domain cells
        int _fieldsPerFace;
        bool _constructFaces;
        std::list<Cell*> _faceList;
        std::vector<Cell*> _face;
        
        std::vector<std::vector<Cell*> > _partition;
        
        DomainManager();
        virtual ~DomainManager();
    };
}

#endif	/* DOMAINMANAGER_HPP */
