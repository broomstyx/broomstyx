/* --------------------------------------------------------------------------
 * Input file declaration:
 * 
 * ELEMENT_TYPES
 *    ...
 *    <tag> HeatTransfer_Fe_Tri3
 *            NodalDof <dof1>
 *            Stage <stg>
 *            CellFieldOutput <n2>
 *                   ...             ...
 *                   ...             ...
 * 
 * Field data tags:
 * 
 *  T_x : temperature gradient along x
 *  T_y : temperature gradient along y
 *
 * Material Set:
 *  material 1 --> density
 *  material 2 --> heat capacity
 *  material 3 --> thermal conductivity
 * 
 ****************************************************************************/

#ifndef HEATTRANSFER_FE_TRI3_HPP
#define	HEATTRANSFER_FE_TRI3_HPP

#include "Numerics/Numerics.hpp"
#include "Core/DofManager.hpp"
#include "Core/NumericsManager.hpp"
#include "BasisFunctions/Triangle_P1.hpp"

namespace broomstyx
{
    class MaterialStatus;
        
    class NumericsStatus_HeatTransfer_Fe_Tri3 final : public NumericsStatus
    {
        friend class HeatTransfer_Fe_Tri3;
        
    public:
        NumericsStatus_HeatTransfer_Fe_Tri3();
        virtual ~NumericsStatus_HeatTransfer_Fe_Tri3();
        
    private:
        RealVector _T; // This should be a double but is stored as a vector in order to avoid temporary
                       // instantiation of vectors during local calculations
        RealVector _gradT;
        RealMatrix _JmatInv;
        double     _Jdet;
        MaterialStatus* _materialStatus;
    };
    
    class HeatTransfer_Fe_Tri3 final : public Numerics
    {
    public:
        HeatTransfer_Fe_Tri3();
        virtual ~HeatTransfer_Fe_Tri3();

        void   deleteNumericsAt( Cell* targetCell ) override;
        void   finalizeDataAt( Cell* targetCell, const TimeData& time ) override;
        double giveCellFieldValueAt( Cell* targetCell, int fieldNum ) override;
        
        RealVector 
            giveCellNodeFieldValuesAt( Cell* targetCell, int fieldNum ) override;
        
        std::vector<RealVector> 
            giveEvaluationPointsFor( Cell* targetCell ) override;
        
        std::tuple< RealVector,RealVector > 
            giveFieldOutputAt( Cell* targetCell, const std::string& fieldTag ) override;
        
        std::tuple< std::vector<Dof*>
                  , std::vector<Dof*>
                  , RealVector >
            giveStaticCoefficientMatrixAt( Cell*           targetCell
                                         , int             stage
                                         , int             subsys
                                         , const TimeData& time ) override;
        
        std::tuple< std::vector<Dof*>, RealVector >
            giveStaticLeftHandSideAt( Cell*           targetCell
                                    , int             stage
                                    , int             subsys
                                    , const TimeData& time ) override;
        
        std::tuple< std::vector<Dof*>, RealVector >
            giveStaticRightHandSideAt( Cell*                    targetCell
                                     , int                      stage
                                     , int                      subsys
                                     , const BoundaryCondition& bndCond
                                     , const TimeData&          time ) override;

        std::tuple< std::vector<Dof*>, RealVector >
            giveStaticRightHandSideAt( Cell*                 targetCell
                                     , int                   stage
                                     , int                   subsys
                                     , const FieldCondition& fldCond
                                     , const TimeData&       time ) override;

        std::tuple< std::vector<Dof*>, std::vector<Dof*>, RealVector >
            giveTransientCoefficientMatrixAt( Cell*           targetCell
                                            , int             stage
                                            , int             subsys
                                            , const TimeData& time ) override;

        std::tuple< std::vector<Dof*>, RealVector >
            giveTransientLeftHandSideAt( Cell*           targetCell
                                       , int             stage
                                       , int             subsys
                                       , const TimeData& time
                                       , ValueType       valType ) override;
        
        void imposeConstraintAt( Cell*                    targetCell
                               , int                      stage
                               , const BoundaryCondition& bndCond
                               , const TimeData&          time ) override;
        
        void initializeNumericsAt( Cell* targetCell ) override;
        void setDofStagesAt( Cell* targetCell ) override;

    private:
        Triangle_P1 _basisFunction;
        RealVector  _basisFunctionValues;
        RealMatrix  _basisFunctionDerivatives;
        
        RealVector  _gpNatCoor;
        double      _wt;
        RealMatrix  _massMatrix;
        
        NumericsStatus_HeatTransfer_Fe_Tri3*
                   getNumericsStatusAt( Cell* targetCell );
        RealMatrix giveBmatAt( Cell* targetCell );
        RealMatrix giveJacobianMatrixAt( Cell* targetCell, const RealVector& natCoor );
        RealVector giveLocalValuesAt( std::vector<Dof*>& dof, ValueType valType );
        std::vector<Dof*> giveNodalDofsAt( Cell* targetCell );
    };
}

#endif	/* HEATTRANSFER_FE_TRI3_HPP */
