#ifndef __ACCUMULATOR__
#define __ACCUMULATOR__

#include "Conventions.hh"

namespace SGF {


/*

  class BinnedAccumulator
  An accumulator that has bins. When you push data to it, it stores them in
  a buffer until the bin is full in which case it pushes them in the bin.

  ======== Usage ==========
  BinnedAccumulator<double> bacc; // a binned accumulator of doubles
  bacc.push(data);                // push a value in to the buffer (is meaned to be called frequently).
  double average=bacc.average();  // get the average of the bins plus the constant value.
  double sigma=bacc.sigma();      // get the standard deviation.

*/


template<class T>
class BinnedAccumulator {
  _integer_counter n;
  T mean;
  T M2;
public:
  BinnedAccumulator() : n(0), mean(0), M2(0) {}
  inline void push(const T& x) {
    ++n;
    T delta = x - mean;
    T R = delta / n;
    mean = mean + R;
    M2 = M2 + delta * (x - mean);
  }
  inline T sigma() const {
    return sqrt(M2 / (n - 1));
  }
  inline T average() const {
    return mean;
  }

};

}

#endif
