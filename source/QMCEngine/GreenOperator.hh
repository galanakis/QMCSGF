/******************************************************************\

Copyright 2015 Dimitrios Galanakis

This file is part of QMCSGF

QMCSGF is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 3.

QMCSGF is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with QMCSGF.  If not, see <http://www.gnu.org/licenses/>.

\*******************************************************************/

#ifndef __GREENOPERATOR__
#define __GREENOPERATOR__

#include <iostream>
#include <cmath>
#include <vector>

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
public:
  static T _GreenOperator(unsigned int n,unsigned long NSites,unsigned int LinesCutoff) {
    double result=1.0;
    if(n <= LinesCutoff && n!=0) result=1.0/NSites;
    if(n > LinesCutoff ) result=pow(1.0/NSites,1.0*n);
    return result;
  }
public:
  GreenOperator() {}
  GreenOperator(const GreenOperator &o) : cache(o.cache) {}
  GreenOperator(unsigned long _nsites,unsigned int _cutoff) { initialize(_nsites,_cutoff); }
  inline T operator()(unsigned int n) const {
    return (n<cache.size())?cache[n]:0.0;
  }

  void initialize(unsigned long _nsites,unsigned int _cutoff) {
    if(_nsites==0) {
      std::cerr<<"Number of sites cannot be zero"<<std::endl;
      exit(123);
    }
    cache.clear();
    T value;
    unsigned int i=0;
    while((value=_GreenOperator(i++,_nsites,_cutoff)))
      cache.push_back(value);

  }
};

}

#endif
