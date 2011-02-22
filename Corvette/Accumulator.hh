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

/* 
class AutoAccumulator : public BinnedAccumulator {
  std::vector<double> VBins;
public:
  AutoAccumulator(int _size) : BinnedAccumulator(_size), VBins() {}
  AutoAccumulator &push(double data,double Weight) {
    BinnedAccumulator::push(data,Weight);
    if(full())
      VBins.push_back(bufferaverage());
    return *this;
  }
  
}; 
*/
/*
  class Alpha Accumulator
  This class is used to define a container for the Alpha coefficients
  together with the keep creating and keep destroying probabilities.
  It can return the values of Alpha[ADD/REMOVE] and in can accumulate
  the probabilities. Also it provides the method update which recalculates
  the Alphas from the Probabilities.
*/

class AlphaAccumulator {
  uint ndata[2];
  double SumProbability[2]; // Probability to keep creating and keep destroying
  double Alpha[2];
  double AlphaParameter;
public:
  AlphaAccumulator(double a) : AlphaParameter(a) {
    Alpha[ADD]=Alpha[REMOVE]=AlphaParameter;
    SumProbability[ADD]=SumProbability[REMOVE]=0.0;
    ndata[ADD]=ndata[REMOVE]=0;
  }
  AlphaAccumulator(AlphaAccumulator &a) : AlphaParameter(a.AlphaParameter) {
    SumProbability[ADD]=a.SumProbability[ADD];
    SumProbability[REMOVE]=a.SumProbability[REMOVE];
    Alpha[ADD]=a.Alpha[ADD];
    Alpha[REMOVE]=a.Alpha[REMOVE];
    ndata[ADD]=a.ndata[ADD];
    ndata[REMOVE]=a.ndata[REMOVE];
  }
  inline void push(int action,double Probability) { 
    SumProbability[action]+=Probability;
    ++ndata[action];
  }
  inline double &Parameter() {return AlphaParameter;}
  inline double operator[](int action) const { return Alpha[action]; }
  inline double Probability(int action) const {return SumProbability[action]/ndata[action];}
  inline double Equilibrium(int action) const {return AlphaParameter*Probability(!action)/Max(Probability(REMOVE),Probability(ADD));}
  inline void equilibrate(int action) {
    if(ndata[action]%1000==1) { 
      Alpha[action]=Equilibrium(action);
      ndata[action]=0;
      SumProbability[action]=0.0;
    }
  } 
};

}

#endif
