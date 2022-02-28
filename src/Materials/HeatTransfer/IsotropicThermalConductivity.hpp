// ---------------------------------------------------------------------------
//   Declaration in input file:
//
//   <n>  IsotropicThermalConductivity <AM> <k>
//
//   where n  = material ID
//         AM = Analysis mode: '2D'/'3D'
//   
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
// ---------------------------------------------------------------------------

#ifndef ISOTROPICTHERMALCONDUCTIVITY_HPP
#define	ISOTROPICTHERMALCONDUCTIVITY_HPP

#include "Materials/Material.hpp"

namespace broomstyx
{
    class IsotropicThermalConductivity final : public Material
    {
    public:
        IsotropicThermalConductivity();
        virtual ~IsotropicThermalConductivity();
        
        RealMatrix giveModulusFrom( const RealVector& conState, const MaterialStatus* matStatus ) override;
        double     giveParameter( const std::string& str ) override;
        void       readParamatersFrom( FILE* fp ) override;
        
    private:
        int        _analysisMode;
        RealMatrix _k;
    };
}
    
#endif	/* ISOTROPICTHERMALCONDUCTIVITY_HPP */
