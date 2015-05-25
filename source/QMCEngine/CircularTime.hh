/******************************************************************\

Copyright 2015 Dimitrios Galanakis

This file is part of QMCSGF

QMCSGF is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 3.

QMCSGF is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with QMCSGF.  If not, see <http://www.gnu.org/licenses/>.

\*******************************************************************/

#ifndef __CIRCULARTIME__
#define __CIRCULARTIME__

#include <cmath>

namespace SGF {
/*
  Class Circular time

  represents a "special" kind of "double" which implements addition
  and subtraction modulo 1.0. For example if t1=0.3 and t2=0.8 and
  1.0=1.0 t1+t2=0.1, t1-t2=0.5 etc. The class has only one variable
  _time which is in the interval 0<_time<=1.0.

  To do this we define all possible additions/subtractions between
  CircilarTimes and doubles, such that the result is always of type
  Circular time.

  Note that we include 1.0 and exclude zero. This is because we
  don't want delta tau to ever be zero. A worse case scenario would
  occur when there are no operators in the string. In this case
  Delta Tau would be the time of the left minus the time of the right
  operator, which are all the same operator. Then delta tau would be
  zero, the time shift would be zero and the algorithm would get stuck.
	This is why the zero is excluded.

*/

typedef double circulartime_t;

class CircularTime {
  circulartime_t _time;
public:
  CircularTime() : _time(0) {}
  explicit CircularTime(circulartime_t t) : _time(t) {}

  inline circulartime_t  time() const {
    return _time+circulartime_t(1.0)-ceil(_time);
  }

  friend CircularTime operator+(const CircularTime &a,const CircularTime &b);
  friend CircularTime operator-(const CircularTime &a,const CircularTime &b);

};

inline CircularTime operator+(const CircularTime &a,const CircularTime &b) {return CircularTime(a._time+b._time);}
inline CircularTime operator-(const CircularTime &a,const CircularTime &b) {return CircularTime(a._time-b._time);}

std::ostream &operator<<(std::ostream &os,const CircularTime &o) { return os<<o.time(); }

}

#endif
