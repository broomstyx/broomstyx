// ---------------------------------------------------------------------------
//   Declaration in input file:
//
//   <n>  ConstantPermeability <AM> <k11, k22, ... >
//
//   where n  = material ID
//         AM = Analysis mode: '2D'/'3D'
//   
//         For 2D, 3 permeability components must be entered
//             <k_11, k_22, k_12>
//
//         For 3D, 6 permeability components must be entered
//             <k_11, k_22, k_33, k_12, k_23, k_12>
//
// ---------------------------------------------------------------------------
//   Requirements on function arguments:
//
//   For 2D:
//      conState = {dh_dx1, dh_dx2}^T
//
//   For 3D:
//      conState = {dh_dx1, dh_dx2, dh_dx3}^T
//
//   Argument 'matData' is not used.
// ---------------------------------------------------------------------------

#ifndef CONSTANTPERMEABILITY_HPP
#define	CONSTANTPERMEABILITY_HPP

#include "../../Material.hpp"

namespace broomstyx
{
    class ConstantPermeability final : public Material
    {
    public:
        ConstantPermeability();
        virtual ~ConstantPermeability();
        
        RealVector giveForceFrom( const RealVector& conState, const MaterialStatus* matStatus ) override;
        RealMatrix giveModulusFrom( const RealVector& conState, const MaterialStatus* matStatus ) override;
        double     giveParameter( const std::string& str ) override;
        void       readParamatersFrom( FILE* fp ) override;
        
    private:
        int        _analysisMode;
        RealMatrix _k;
    };
}
    
#endif	/* CONSTANTPERMEABILITY_HPP */