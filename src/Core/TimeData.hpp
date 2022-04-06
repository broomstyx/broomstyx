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

#ifndef TIMEDATA_HPP
#define	TIMEDATA_HPP

#include <cmath>

namespace broomstyx
{
    class TimeData final
    {
    public:
        TimeData()
            : _start(0.)
            , _end(0.)
            , _current(0.)
            , _increment(0.)
            , _target(0.)
        {}

        void advanceTime()
        {
            _current = _target;
            _target += _increment;

            if ( _target > _end )
                _target = _end;
        }

        double giveCurrentTime() const { return _current; }
        double giveEndTime() const { return _end; }
        double giveStartTime() const { return _start; }
        double giveTargetTime() const { return _target; }
        double giveTimeIncrement() const { return _increment; }

        bool   hasReachedEnd() const
        {
            if ( std::fabs(_current - _end) < 1.0e-13 )
                return true;
            else
                return false;
        }

        // Note: modifying the current time or the time increment
        // causes a recalculation of the target time
        void setCurrentTimeTo( double val )
        {
            _current = val;
            _target = _current + _increment;
        }

        void setEndTimeTo( double val )
        {
            _end = val;
        }

        void setStartTimeTo( double val )
        {
            _start = val;
        }

        void setTimeIncrementTo( double val )
        {
            _increment = val;
            _target = _current + _increment;
        }

    private:
        double _start;
        double _end;
        double _current;
        double _increment;
        double _target;
    };
}

#endif	/* TIMEDATA_HPP */
