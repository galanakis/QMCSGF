#ifndef __ATOMICTERM__
#define __ATOMICTERM__

#include <vector>
#include <cmath>
#include "Conventions.hh"

namespace SGF {


template<typename IndexedProductElementClass,int direction>
MatrixElement MultiplyMe(const typename std::vector<IndexedProductElementClass> &p) {
  unsigned int result=1;
  for(typename std::vector<IndexedProductElementClass>::const_iterator ip=p.begin(); ip<p.end(); ++ip)
    result*=ip->template amplitude<direction>();
  return sqrt(result);
}

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

template<typename MatrixElementClass,typename IndexedProductElementClass>
class AtomicTerm {
  MatrixElementClass _coefficient;
  typedef typename std::vector<IndexedProductElementClass> IndexProductElementVector;
  IndexProductElementVector _product;
public:
  AtomicTerm() : _coefficient(1) {} // The constructor really depends on the interface

  AtomicTerm(const MatrixElementClass &c,const IndexedProductElementClass &p) : _coefficient(c), _product(1,p) {}

  AtomicTerm(const MatrixElementClass &c,const IndexedProductElementClass &p1,const IndexedProductElementClass &p2) : _coefficient(c) {
    _product.reserve(2);
    _product.push_back(p1);
    _product.push_back(p2);
  }

  AtomicTerm(MatrixElementClass c,const IndexProductElementVector &p) : _coefficient(c), _product(p) {}

  AtomicTerm(const IndexProductElementVector &p) : _coefficient(MatrixElementClass(1.0)), _product(p) {}

  AtomicTerm(const AtomicTerm &o) : _coefficient(o._coefficient), _product(o._product) {}

  typedef typename IndexProductElementVector::const_iterator const_iterator;

  const_iterator begin() const {return _product.begin();}
  const_iterator end() const {return _product.end();}

  inline MatrixElementClass &coefficient() {return _coefficient;};
  inline const MatrixElementClass &coefficient() const {return _coefficient;};
  inline const IndexProductElementVector& product() const {return _product;}

  /* offset after adding/removing the present term */
  template<int action>
  inline int offset() const {
    int result=0;
    for(const_iterator ip=_product.begin(); ip!=_product.end(); ++ip)
      result+=ip->template offset<action>();
    return result;
  }

  inline int offset() const {
    return offset<ADD>();
  }

  /* The matrix element after applying
     the operator to the left or the right */
  template<int direction>
  inline unsigned int amplitude() const {
    unsigned int result=1;
    for(const_iterator ip=_product.begin(); ip!=_product.end(); ++ip)
      result*=ip->template amplitude<direction>();
    return result;
  }

  /* The matrix element after applying
     the operator to the left or the right */
  template<int direction>
  inline MatrixElementClass me() const {
    return _coefficient*sqrt(amplitude<direction>());
  }

  template<int direction,int action>
  inline void update_psi() const {
    for(const_iterator ip=_product.begin(); ip!=_product.end(); ++ip)
      ip->template update<direction,action>();
  }

  inline int absdelta() const {
    int result=0;
    for(const_iterator ip=_product.begin(); ip!=_product.end(); ++ip)
      result+=Abs(ip->delta());
    return result;        
  }

};

}

#endif
