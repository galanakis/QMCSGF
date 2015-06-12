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

#ifndef __ATOMICTERM__
#define __ATOMICTERM__

#include <vector>
#include <cmath>
#include "Conventions.hh"

namespace SGF {


template<typename Float,typename OProd,int direction>
Float MultiplyMe(const typename std::vector<OProd> &p) {
  unsigned int result=1;
  for(typename std::vector<OProd>::const_iterator ip=p.begin(); ip<p.end(); ++ip)
    result*=ip->template amplitude<direction>();
  return sqrt(result);
}


template<typename Float,typename OProd>
class AtomicTermLight {
public:
  typedef typename std::vector<OProd> OProdVector;
  typedef typename OProdVector::const_iterator const_iterator;
private:
  Float _coefficient;
  OProdVector _product;
public:
  AtomicTermLight(Float c,const OProdVector &p) : _coefficient(c), _product(p) {}
  AtomicTermLight(const AtomicTermLight &o) : _coefficient(o._coefficient), _product(o._product) {}

  const_iterator begin() const {return _product.begin();}
  const_iterator end() const {return _product.end();}
  inline const OProdVector& product() const {return _product;}

  inline Float &coefficient() {return _coefficient;};
  inline const Float &coefficient() const {return _coefficient;};

  /* The matrix element after applying
     the operator to the left or the right */
  template<int direction>
  inline Float me() const {
    unsigned int result=1;
    for(const_iterator ip=_product.begin(); ip!=_product.end(); ++ip)
      result*=ip->template amplitude<direction>();
    return _coefficient*sqrt(result);
  }

  template<int direction,int action>
  inline void update_psi() const {
    for(const_iterator ip=_product.begin(); ip!=_product.end(); ++ip)
      ip->template update<direction,action>();
  }

};

/*

  class AtomicTerm

  This is the actual data structure that holds a kinetic/potential
  term, that is a product creation/annihilation operators for
  different indices.

  Note: the structure does not necessarily guarantee that there can't
  be two ProductElements for the same index. To do this one we
  should use a map instead of a vector. However the std::map costs a
  factor of 10 in performance. I am in serious need of an amazing
  implementation of a map or hash map for integers. Or I should just
  give up some performance.

  This class is defined as a template with two parameters
  the one is the type of the MatrixElement (usually double)
  and the other one a special class which contains the combination
  of a particle and an operator.

*/

template<typename Float,typename OProd>
class AtomicTerm : public AtomicTermLight<Float,OProd> {
public:

  using typename AtomicTermLight<Float,OProd>::OProdVector;
  using typename AtomicTermLight<Float,OProd>::const_iterator;
  using AtomicTermLight<Float,OProd>::begin;
  using AtomicTermLight<Float,OProd>::end;

  AtomicTerm(Float c,const OProdVector &p) : AtomicTermLight<Float,OProd>(c,p) {}

  AtomicTerm(const AtomicTermLight<Float,OProd> &o) : AtomicTermLight<Float,OProd>(o) {}

  AtomicTerm(const AtomicTerm &o) : AtomicTermLight<Float,OProd>(o) {}

  /* offset after adding/removing the present term */
  template<int action=ADD>
  inline int offset() const {
    int result=0;
    for(const_iterator ip=begin(); ip!=end(); ++ip)
      result+=ip->template offset<action>();
    return result;
  }

  inline int absdelta() const {
    int result=0;
    for(const_iterator ip=begin(); ip!=end(); ++ip)
      result+=Abs(ip->delta());
    return result;        
  }

};

}

#endif
