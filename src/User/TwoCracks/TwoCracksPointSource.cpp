#include "TwoCracksPointSource.hpp"
#include "../../Core/ObjectFactory.hpp"

using namespace broomstyx;

registerBroomstyxObject(UserFunction, TwoCracksPointSource)


TwoCracksPointSource::TwoCracksPointSource() {}

TwoCracksPointSource::~TwoCracksPointSource() {}

double TwoCracksPointSource::at(const RealVector& coor, const TimeData& time)
{
    return 300. + 100*time.target;
}