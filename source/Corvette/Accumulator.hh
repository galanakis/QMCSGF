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
  unsigned long long _count;
public:
  Accumulator() { reset(); }
  Accumulator(const Accumulator &o) {
    _count=o._count; 
    for(int i=0;i<MomentOrder;++i) _sum[i]=o._sum[i]; 
  }
  Accumulator & push(const T &data,const T &Weight) {
    T moment(Weight);
    for(int i=0;i<MomentOrder;++i) {
      _sum[i]+=moment;
      moment*=data;
    }
    ++_count;
    return *this;
  }
  inline void reset() {
    _count=0; 
    for(int i=0;i<MomentOrder;++i) _sum[i]=T(0); 
  }
  inline unsigned long long count() const {return _count;}
  inline T operator[](int i) const {return _sum[i];} 
  inline T operator()(int i) const {return _sum[i]/_sum[0];}
  inline int nmoments() const {return MomentOrder;}
};

template<class T>
class BinnedAccumulator : public Accumulator<3,T> {
  Accumulator<2,T> Buffer;
	Accumulator<3,T> Bins;
public:
  BinnedAccumulator() : Accumulator<3,T>(), Buffer() {}
  inline void flush(T Weight) {
    Bins.push(Buffer[1]/Weight,1.0);
    Buffer.reset();
  }
  inline void push(T data,T Weight) { Buffer.push(data,Weight); }
  inline T average() const {return Bins(1);}
  inline T sigma() const {return sqrt(fabs(Bins(2)-Bins(1)*Bins(1))/(Bins.count()-1));}
};

}

#endif
