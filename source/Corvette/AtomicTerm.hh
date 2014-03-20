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
public:
  typedef typename std::vector<IndexedProductElementClass> IndexProductElementVector;
  typedef typename IndexProductElementVector::const_iterator const_iterator;
private:
  MatrixElementClass _coefficient;
  IndexProductElementVector _product;
public:
  AtomicTerm(MatrixElementClass c,const IndexProductElementVector &p) : _coefficient(c), _product(p) {}

  AtomicTerm(const AtomicTerm &o) : _coefficient(o._coefficient), _product(o._product) {}

  const_iterator begin() const {return _product.begin();}
  const_iterator end() const {return _product.end();}
  inline const IndexProductElementVector& product() const {return _product;}

  inline MatrixElementClass &coefficient() {return _coefficient;};
  inline const MatrixElementClass &coefficient() const {return _coefficient;};

  /* offset after adding/removing the present term */
  template<int action=ADD>
  inline int offset() const {
    int result=0;
    for(const_iterator ip=_product.begin(); ip!=_product.end(); ++ip)
      result+=ip->template offset<action>();
    return result;
  }

  /* The matrix element after applying
     the operator to the left or the right */
  template<int direction>
  inline MatrixElementClass me() const {
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

  inline int absdelta() const {
    int result=0;
    for(const_iterator ip=_product.begin(); ip!=_product.end(); ++ip)
      result+=Abs(ip->delta());
    return result;        
  }

};

}

#endif
