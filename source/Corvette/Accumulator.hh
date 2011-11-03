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
  Accumulator & push(const T &data,const T &Weight=T(1.0)) {
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
  inline T sum(int i) const {return _sum[i];} 
  inline T operator()(int i) const {return _sum[i]/_sum[0];}
  inline int nmoments() const {return MomentOrder;}
};

template<class T>
class KahanSum {
	T sum;
	T compensation;
public:
	KahanSum(const T _sum,const T _compensation) : sum(_sum), compensation(_compensation) {}
	KahanSum(const KahanSum &o) : sum(o._sum), compensation(o._compensation) {}
	inline KahanSum & operator+=(const T &input) {
		T y=input-compensation;
		T t=sum+y;
		compensation=(t-sum)-y;
		sum=t;
	}
};

template<class T>
class BinnedAccumulator {
public:
	T Base;                       // A constant value
	unsigned long long Bins_sum0;
	T Bins_sum1;
	T Bins_sum2;
	T Buffer;
public:
  BinnedAccumulator() : Base(0), Bins_sum0(0), Bins_sum1(0), Bins_sum2(0), Buffer(0) {}
  inline void flush(T Weight) {
    
		T data=Buffer/Weight;
		
		++Bins_sum0;
		Bins_sum1+=data;
		Bins_sum2+=data*data;

		Buffer=T(0);
  }

  inline void push(T data) { 
		Buffer+=data;
	}
	
	inline T &constant() {return Base;}
  inline T average() const {return Base+Bins_sum1/Bins_sum0;}
  inline T sigma() const {return sqrt(fabs(Bins_sum2/Bins_sum0-Bins_sum1*Bins_sum1/Bins_sum0/Bins_sum0)/(Bins_sum0-1));}
};

template<class T> inline std::ostream& operator<<(std::ostream& output, const BinnedAccumulator<T> &o) { 
return output<<o.average()<<" +/- "<<o.sigma(); 
} 

}

#endif
