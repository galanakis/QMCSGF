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
public:
  static inline circulartime_t  Plus(const circulartime_t  t1,const circulartime_t  t2) {
    circulartime_t  sum=t1+t2;
    return (sum<=circulartime_t(1.0))?sum:sum-circulartime_t(1.0);
  }
  static inline circulartime_t  Minus(const circulartime_t  t1,const circulartime_t  t2) {
    circulartime_t  diff=t1-t2;
    return (diff<=circulartime_t(0))?diff+circulartime_t(1.0):diff;
  }

  /*  Shift is the key function that could replace Plus and Minus.
      It projects every double in the interval ]0,1.0]. To understand why
      we consider the function 1+x-ceiling(x), which projects every number
      in the interval ]0,1]. */
  static inline circulartime_t  Shift(const circulartime_t  t) { return t+circulartime_t(1.0)-ceil(t); }
  circulartime_t  _time;
public:
  CircularTime() : _time(1.0) {}

  CircularTime(const CircularTime &o) : _time(o._time) {}
  inline CircularTime &operator= (const CircularTime &o) { _time=o._time; return *this; }

  CircularTime(circulartime_t  t) { _time=Shift(t); }
  inline CircularTime &operator= (const circulartime_t  t) { _time=Shift(t); return *this; }

  inline circulartime_t  time() const {return _time;}

  inline CircularTime &operator+=(const CircularTime &o) { _time=Plus(_time,o._time);  return *this; }
  inline CircularTime &operator+=(const circulartime_t  t) { _time=Shift(t+_time);  return *this; }

  inline CircularTime &operator-=(const CircularTime &o) { _time=Minus(_time,o._time); return *this; }
  inline CircularTime &operator-=(const circulartime_t  t) { _time=Shift(_time-t); return *this; }

  friend CircularTime operator+(const CircularTime &,const CircularTime &);
  friend CircularTime operator+(const CircularTime &,const circulartime_t );
  friend CircularTime operator+(const circulartime_t ,const CircularTime &);

  friend CircularTime operator-(const CircularTime &,const CircularTime &);
  friend CircularTime operator-(const CircularTime &,const circulartime_t );
  friend CircularTime operator-(const circulartime_t ,const CircularTime &);

  friend bool operator==(const CircularTime &a,const CircularTime &b);
};

inline CircularTime operator+(const CircularTime &a,const CircularTime &b) {return CircularTime(CircularTime::Plus (a._time,b._time));}
inline CircularTime operator+(const CircularTime &a,const circulartime_t  t) {return CircularTime(a._time+t);}
inline CircularTime operator+(const circulartime_t  t,const CircularTime &a) {return CircularTime(t+a._time);}

inline CircularTime operator-(const CircularTime &a,const CircularTime &b) {return CircularTime(CircularTime::Minus(a._time,b._time));}
inline CircularTime operator-(const CircularTime &a,const circulartime_t  t) {return CircularTime(a._time-t);}
inline CircularTime operator-(const circulartime_t  t,const CircularTime &a) {return CircularTime(t-a._time);}

inline bool operator==(const CircularTime &a,const CircularTime &b) { return (a._time-b._time)*(a._time-b._time)<1e-14 && (a._time-b._time)*(a._time-b._time)<(a._time+b._time)*(a._time+b._time)*1e-14 ;}

}

#endif
