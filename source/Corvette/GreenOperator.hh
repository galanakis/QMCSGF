#ifndef __GREENOPERATOR__
#define __GREENOPERATOR__

#include <cmath>
#include <vector>
#include "Conventions.hh"

namespace SGF {


/*
class GreenOperator
This defines a functor which stores and returns the values
of the Green Operator. The values are saved in a cache.
If the requested argument is outside the cache, then
the cache is expanded.
The usage is very simple:

GreenOperator G(12120);
double a=G(13);
*/
template<class T>
class GreenOperator {
  std::vector<T> cache;
  unsigned long NSites;
	unsigned int LinesCutoff;
public:
  T _GreenOperator(int n) const { 
    double result=1.0;
    if(n <= LinesCutoff && n!=0) result=1.0/NSites;
    if(n > LinesCutoff ) result=pow(1.0/NSites,n);
    return result; 
  }
public:
  GreenOperator() : NSites(0), LinesCutoff(0) {}
  GreenOperator(int _nsites,int lines) : NSites(_nsites),LinesCutoff(lines) {}
  inline T operator()(int n) const { 
    return (n<cache.size() && n>=0)?cache[n]:0.0; 
  }

  void initialize(unsigned long _nsites,unsigned int _cutoff) {
    if(_nsites==0) {
      std::cout<<"Number of sites cannot be zero"<<std::endl;
      exit(123);
    }
    NSites=_nsites;
		LinesCutoff=_cutoff;
		std::cout<<"GreenOperator Initialize, nsites= "<<NSites<<", cutoff= "<<LinesCutoff<<std::endl;
    cache.clear();
    T value;
    int i=0;
    while((value=_GreenOperator(i++)))
      cache.push_back(value);
    
  }
};

}

#endif