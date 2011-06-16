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
  int NSites;
public:
  T _GreenOperator(int n) const { 
    double result=1.0;
    if(n<=4 && n!=0) result=1.0/NSites;
    if(n>4) result=pow(1.0/NSites,n);
    return result; 
  }
public:
  GreenOperator() : NSites(0) {}
  GreenOperator(int _nsites) : NSites(_nsites) {}
  inline T operator()(int n) const { 
    return (n<cache.size() && n>=0)?cache[n]:0.0; 
  }

  void initialize(uint _nsites) {
    if(_nsites==0) {
      std::cout<<"Number of sites cannot be zero"<<std::endl;
      exit(123);
    }
    NSites=_nsites;
    cache.clear();
    T value;
    int i=0;
    while((value=_GreenOperator(i++)))
      cache.push_back(value);
    
  }
};

}

#endif
