#ifndef TWOCRACKSPOINTSOURCE_HPP
#define TWOCRACKSPOINTSOURCE_HPP

#include "../UserFunction.hpp"

namespace broomstyx
{
    class TwoCracksPointSource : public UserFunction
    {
    public:
        TwoCracksPointSource();
        virtual ~TwoCracksPointSource();
        
        double at( const RealVector& coor, const TimeData& time ) override;

    private:
    };
}

#endif /* TWOCRACKSPOINTSOURCE_HPP */

