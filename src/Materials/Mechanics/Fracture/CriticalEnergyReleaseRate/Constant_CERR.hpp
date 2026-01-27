#ifndef CONSTANT_CERR_HPP
#define CONSTANT_CERR_HPP

#include "Materials/Material.hpp"

namespace broomstyx
{
    class Constant_CERR final : public Material
    {
    public:
        Constant_CERR();
        virtual ~Constant_CERR();
        
        RealMatrix giveModulusFrom( const RealVector& conState, const MaterialStatus* matStatus) override;
        void       readParamatersFrom( FILE* fp ) override;

    private:
        double _Gc;
    };
}

#endif /* CONSTANT_CERR_HPP */
