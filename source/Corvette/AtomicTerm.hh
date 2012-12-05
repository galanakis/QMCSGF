#ifndef __ATOMICTERM__
#define __ATOMICTERM__

#include <vector>
#include <cmath>
#include "Conventions.hh"

namespace SGF {

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
  std::vector<IndexedProductElementClass> _product;
public:
  AtomicTerm() : _coefficient(1) {} // The constructor really depends on the interface
  
  AtomicTerm(const MatrixElementClass &c,const IndexedProductElementClass &p) : _coefficient(c) { 
    _product.reserve(1);
    _product.push_back(p); 
  }

  AtomicTerm(const MatrixElementClass &c,const IndexedProductElementClass &p1,const IndexedProductElementClass &p2) : _coefficient(c) {
    _product.reserve(2); 
    _product.push_back(p1);
    _product.push_back(p2); 
  }  
  
	AtomicTerm(MatrixElementClass c,const std::vector<IndexedProductElementClass> &p) : _coefficient(c), _product(p) {}
 
  AtomicTerm(const std::vector<IndexedProductElementClass> &p) : _coefficient(MatrixElementClass(1.0)), _product(p) {}

  AtomicTerm(const AtomicTerm &o) : _coefficient(o._coefficient), _product(o._product) {}
  
  typedef typename std::vector<IndexedProductElementClass>::const_iterator iterator;
	typedef typename std::vector<IndexedProductElementClass>::size_type size_type;
  
	inline MatrixElementClass &coefficient() {return _coefficient;};
  inline const MatrixElementClass &coefficient() const {return _coefficient;}; 
  inline const std::vector<IndexedProductElementClass>& product() const {return _product;}
    
  inline bool diagonal() const { 
    int result=0;
    for(size_type i=0;i<_product.size();++i)
      result+=abs(_product[i].delta());
    return result==0;         
  }

  /* offset after adding/removing the present term */
  inline int offset(int action=ADD) const {
    int result=0;
    for(size_type i=0;i<_product.size();++i)
      result+=_product[i].offset(action);
    return result;      
  }

  /* The matrix element after applying 
     the operator to the left or the right */
  inline unsigned int amplitude(int direction) const {
    unsigned int result=1;
    for(size_type i=0;i<_product.size();++i)
      result*=_product[i].amplitude(direction);
    return result;
  }

  /* The matrix element after applying 
     the operator to the left or the right */
  inline MatrixElementClass me(int direction) const {
     return _coefficient*sqrt(amplitude(direction));
  }
                    
  inline void update_psi(int direction,int action) const {
    for(size_type i=0;i<_product.size();++i)
      _product[i].update(direction,action);
  }
  
  inline int maxoffset() const {
    int result=0;
    for(size_type i=0;i<_product.size();++i)
      result+=_product[i].maxoffset();
    return result;      
  }

  inline int minoffset() const {
    int result=0;
    for(size_type i=0;i<_product.size();++i)
      result+=_product[i].minoffset();
    return result;      
  }
    
};

}

#endif
