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

#ifndef DOFMANAGER_HPP
#define	DOFMANAGER_HPP

#include <cstdio>
#include <map>
#include <vector>
#include <list>

#include "Util/RealVector.hpp"

#define UNASSIGNED -1

namespace broomstyx
{    
    enum ValueType
    {
        current_value = 1
        , incremental_value = 2
        , converged_value = 3
        , correction = 4
        , replacement_correction = 5
    };

    class Dof;
    class Cell;
    class Node;

    class DofManager final
    {
        friend class AnalysisModel;
        
    public:
        // Disable copy constructor and assignment operator
        DofManager( const DofManager& ) = delete;
        DofManager& operator=( const DofManager& ) = delete;
        
        void   addToSecondaryVariableAt( Dof* targetDof, double val );
        void   createCellDofsAt( Cell* targetCell );
        void   createFaceDofsAt( Cell* targetFace );
        void   createNodalDofsAt( Node* targetNode );
        Dof*   createNumericsDofWithGroup( int grp );
        void   destroyCellDofsAt( Cell* targetCell );
        void   destroyFaceDofsAt( Cell* targetFace );
        void   destroyNodalDofsAt( Node* targetNode );
        void   destroyNumericsDof( Dof*& targetDof );
        void   enslave( Dof* targetDof, Dof* masterDof );
        void   finalizeDofPrimaryValues();
        void   findActiveDofs();
        std::vector<Dof*> 
               giveActiveDofsAtStage( int stg );
        int    giveGroupNumberFor( Dof* targetDof );
        int    giveIndexForCellDof( const std::string& name );
        int    giveIndexForFaceDof( const std::string& name );
        int    giveIndexForNodalDof( const std::string& name );
        int    giveEquationNumberAt( Dof* targetDof );
        int    giveNumberOfActiveDofsAtStage( int stgNum );
        int    giveSubsystemNumberFor( Dof* targetDof );
        double giveValueOfConstraintAt( Dof* targetDof, ValueType valType );
        double giveValueOfPrimaryVariableAt( Dof* targetDof, ValueType valType );
        double giveValueOfSecondaryVariableAt( Dof* targetDof );
        double giveValueOfResidualAt( Dof* targetDof );
        void   imposeMultiFreedomConstraints();
        void   putDirichletConstraintOn( Dof* targetDof );
        void   readCellDofsFrom( FILE* fp );
        void   readFaceDofsFrom( FILE* fp );
        void   readMultiFreedomConstraintsFrom( FILE* fp );
        void   readNodalDofsFrom( FILE* fp );
        void   removeAllDofConstraints();
        void   reportNumberOfActiveDofs();
        void   resetDofCurrentPrimaryValues();
        void   resetSecondaryVariablesAtStage( int stage );
        void   setConstraintValueAt( Dof* targetDof, double val );
        void   setEquationNumberFor( Dof* targetDof, int eqNo );
        void   setStageFor( Dof* targetDof, int stgNum );
        void   setSubsystemFor( Dof* targetDof, int subsysNum );        
        void   updatePrimaryVariableAt( Dof* targetDof, double val, ValueType valType );
        void   updateResidualAt( Dof* targetDof, double val );
        void   writeConvergedDofValuesTo( Node* targetNode );

    private:
        struct DofInfo
        {
            std::string tag;
            int group;
            int primField;
            int secField;
        };
        
        struct MultiFreedomConstraint
        {
            std::string type;
            std::string masterTag;
            int         masterDofNum;
            std::string slaveTag;
            int         slaveDofNum;
        };

        std::vector<DofInfo> _cellDofInfo;
        std::vector<DofInfo> _faceDofInfo;
        std::vector<DofInfo> _nodalDofInfo;
        std::vector<Dof*> _numericsDof;

        std::map< int, int > _nActiveDof;
        std::map< int, int > _nInactiveDof;
        std::map< int, std::vector<Dof*> > _activeDof;
        std::map< int, std::vector<Dof*> > _inactiveDof;
        std::vector<MultiFreedomConstraint> _multiFreedomConstraint;
        
        DofManager();
        virtual ~DofManager();
        
        void imposeNodalDofSlaveConstraint( MultiFreedomConstraint& mfc );
        void readNodalDofSlaveConstraintDataFrom( FILE* fp, MultiFreedomConstraint& mfc );
    };
}

#endif	/* DOFMANAGER_HPP */