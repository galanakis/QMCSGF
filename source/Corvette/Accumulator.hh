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

	=========== Usage ==========
	Accumulator<3,double> acc;   // an accumulator of doubles that remembers the first 3 moments.
	acc.push(data,weigth);       // push a value and the corresponding weight.
	acc.push(data);              // the default weight is 1.0;
	double average=acc(1);       // the first moment is the average;
	double averagesquare=acc(2); // average of squares;
	int count=acc.count();       // how many numbers have been pushed.
	acc.reset();                 // clears everything.
	int nmoments=acc.nmoments(); // how many moments this accumulator remembers, in this example "3".

*/

template<const int MomentOrder,class T>
class Accumulator {
	T _sum[MomentOrder];
	_integer_counter _count;
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
	inline _integer_counter count() const {return _count;}
	inline T operator()(int i) const {return _sum[i]/_sum[0];}
	inline int nmoments() const {return MomentOrder;}
};


/* 

	class BinnedAccumulator
  An accumulator that has bins. When you push data to it, it stores them in
	a buffer until the bin is full in which case it pushes them in the bin.
	
	======== Usage ==========
	BinnedAccumulator<double> bacc; // a binned accumulator of doubles
	bacc.constant()=10;             // Constant offset.
	bacc.push(data);                // push a value in to the buffer (repeat many times).
	bacc.flush(Weight);             // push in to the bin the sum of the buffer divided by the Weight.
	double average=bacc.average();  // get the average of the bins plus the constant value.
	double sigma=bacc.sigma();      // get the standard deviation.
	std::cout<<bacc<<std::endl;     // the "<<" has bee properly overloaded.
	
*/

template<class T>
class BinnedAccumulator {
public:
	T Base;                       // A constant value
	Accumulator<3,T> Bins;
public:
	BinnedAccumulator() : Base(0) {}
	inline void push(const T &val) {
		Bins.push(val);
	}
  
	inline T &constant() {return Base;}
	inline T average() const {return Base+Bins(1);}
	inline T sigma() const {return sqrt(fabs(Bins(2)-Bins(1)*Bins(1))/(Bins.count()-1));}
};

template<class T> inline std::ostream& operator<<(std::ostream& output, const BinnedAccumulator<T> &o) { 
return output<<o.average()<<" +/- "<<o.sigma(); 
} 

}

#endif
