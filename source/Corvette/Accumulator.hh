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
	inline T operator()(int i) const {return _sum[i]/_sum[0];}
	inline int nmoments() const {return MomentOrder;}
};


template<class T>
class BinnedAccumulator {
public:
	T Base;                       // A constant value
	T Buffer;
	Accumulator<3,T> Bins;
public:
	BinnedAccumulator() : Base(0), Buffer(0) {}
	inline void flush(T Weight) {
		Bins.push(Buffer/Weight);
		Buffer=T(0);
	}
  
	inline void push(T data) { Buffer+=data; }
	inline T &constant() {return Base;}
	inline T average() const {return Base+Bins(1);}
	inline T sigma() const {return sqrt(fabs(Bins(2)-Bins(1)*Bins(1))/(Bins.count()-1));}
};

template<class T> inline std::ostream& operator<<(std::ostream& output, const BinnedAccumulator<T> &o) { 
return output<<o.average()<<" +/- "<<o.sigma(); 
} 

}

#endif
