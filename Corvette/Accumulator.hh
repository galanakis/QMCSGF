#ifndef __ACCUMULATOR__
#define __ACCUMULATOR__

#include <cmath>
#include <vector>
#include <ostream>
#include "Conventions.hh"

namespace SGF {

/* 
  class Accumulator
  You push values to it, and it remembers the average or
  other statistical moments.
*/

template<const int MomentOrder,class T>
class Accumulator {
  T _sum[MomentOrder];
  unsigned long _count;
public:
  Accumulator() { reset(); }
  Accumulator(const Accumulator &o) {
    _count=o._count; 
    for(int i=0;i<MomentOrder;++i) _sum[i]=o._sum[i]; 
  }
  Accumulator & push(const T &data,const T &Weight) {
    T moment(Weight);
    for(int i=0;i<MomentOrder;++i) {
      _sum[i]=(_sum[i]*_count+moment)/(_count+1);
      moment*=data;
    }
    ++_count;
    return *this;
  }
  inline void reset() {
    _count=0; 
    for(int i=0;i<MomentOrder;++i) _sum[i]=0; 
  }
  inline long count() const {return _count;}
  inline T operator[](int i) const {return _sum[i]*_count;} 
  inline T operator()(int i) const {return _sum[i]/_sum[0];}
  inline int nmoments() const {return MomentOrder;}
};

template<class T>
class BinnedAccumulator : public Accumulator<3,T> {
  Accumulator<2,T> Buffer;
public:
  BinnedAccumulator() : Accumulator<3,T>(), Buffer() {}
  inline void flush() {
    Accumulator<3,T>::push(Buffer(1),1.0);
    Buffer.reset();
  }
  inline void flush(double Weight) {
    Accumulator<3,T>::push(Buffer[1]/Weight,1.0);
    Buffer.reset();
  }
  inline void push(T data,double Weight) { Buffer.push(data,Weight); }
  inline double average() const {return (*this)(1);}
  inline double sigma() const {return sqrt(fabs((*this)(2)-(*this)(1)*(*this)(1))/(this->count()-1));}
};

}

#endif
