// ---------------------------------------------------------------------------
//   Declaration in input file:
//
//   <n>  HeatCapacity <cp>
//
//   where   n = material ID
//          cp = density
//
// ---------------------------------------------------------------------------

#ifndef HEATCAPCITY_HPP
#define	HEATCAPCITY_HPP

#include "Materials/Material.hpp"

namespace broomstyx
{
    class HeatCapacity final : public Material
    {
    public:
        HeatCapacity();
        virtual ~HeatCapacity();

        double giveMaterialVariable( const std::string& str, const MaterialStatus* matStatus ) override;
        double giveParameter( const std::string& str ) override;
        void   readParamatersFrom( FILE* fp ) override;
        
    private:
        double _cp;
    };
}

#endif	/* HEATCAPCITY_HPP */

