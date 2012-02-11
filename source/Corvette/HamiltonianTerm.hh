#ifndef __HAMILTONIANTERM__
#define __HAMILTONIANTERM__

#include "ProductElement.hh"
#include "Boson.hh"
#include <vector>
#include <cmath>
#include <map> 
#include <cstdlib>

namespace SGF {

/* 
  The kinetic or potential term represents a product of creation or annihilation
  operators in no particular order:
    c_1 c_2 c^+_1 c^+_3 c_2 c^+_3 c^+_2
  In the product there may be operators with the same indices appearing
  in different positions. We need to represent this product in such
  a way that the offset (number of broken lines it creates) and the 
  matrix element applied on the right or left occupancies is evaluated
  efficiently. We will represent the products of creation and 
  annihilation operators at the same index with the class ProductElement.
  
  In this case we can factor the term as a product of  ProductElements:
  c_1 c^+_1, c_2 c_2 c^+_2, c^+_3 c^+_3. We see that there are three
  ProductElement c c+, cc and c+ c+ applied to particle indices 1,2,3
  respectively.
      
*/


/*
  This class represents a product element with an associated
  index. The index is just a pointer to a particle, rather than
  an index to a vector of particles. This is a bit more efficient
  than keeping track of indices. The class provides a method
  to calculate the offset of an operator (given the occupancies
  of the boson) and the amplitude when the operator is applied
  to the right or left. To calculate the amplitude in particular
  we note the following: for an operator T (acting on a single site)
  <n+Delta | T |n> = f(n), where Delta is the number of creation
  minus the number of annihilation operators. If we act the operator
  on the right at occupancy n we will get f(n), but if we act
  on the left at occupancy n then we will get f(n-Delta). 
  Therefore when we act to the left (direction=0), we need to shift
  the occupancy.
  
*/

class IndexedProductElement : public ProductElement {
 	
	
protected:
  Boson* particle;
public:
	IndexedProductElement(const IndexedProductElement &o) : ProductElement(o), particle(o.particle) {}
  IndexedProductElement(const ProductElement peprod,Boson* const _p) : ProductElement(peprod),particle(_p) {}
  Boson* const &particle_id() const {return particle;}
  inline int n(int direction) const {return particle->n(direction)-delta()*!direction;}
  inline int offset(int action=ADD) const {return ProductElement::offset(Sign[action]*particle->delta());}
  inline uint amplitude(int direction) const {return ProductElement::amplitude(n(direction),particle->nmax());}
  inline void update(int direction,int action) const {

#ifdef DEBUG    
    if(particle->n(direction)==0 && delta()*Sign[action==direction]<0) {
      cout<<"Error: Negative occupancy encountered."<<std::endl;
      exit(2);
    }
#endif

    particle->n(direction) += delta()*Sign[action==direction];
     
  }
  
  inline int maxoffset() const { return Abs(delta()); }
  inline int minoffset() const {
    uint Nmax=particle->nmax();
    uint dn=Abs(delta()); 
    return (Nmax!=0)?(dn-2*Min(dn,Nmax)):(-dn);
  }
};


/* This is the actual data structure that holds a kinetic/potential 
  term, that is a product creation/annihilation operators for
  different indices.
  
  Note: the structure does not necessarily guarantee that there can't
  be two ProductElements for the same index. To do this one we
  should use a map instead of a vector. However the std::map costs a 
  factor of 10 in performance. I am in serious need of an amazing
  implementation of a map or hash map for integers. Or I should just
  give up some performance.
*/
class HamiltonianTerm {
  MatrixElement _coefficient;
  std::vector<IndexedProductElement> _product;
public:
  HamiltonianTerm() : _coefficient(1) {} // The constructor really depends on the interface
  
  HamiltonianTerm(MatrixElement c,const ProductElement &op,Boson *index) : _coefficient(c) {
    _product.reserve(1);
    _product.push_back(IndexedProductElement(op,index));
  }

  HamiltonianTerm(MatrixElement c,const ProductElement &op_i,Boson *index_i,const ProductElement &op_j,Boson *index_j) : _coefficient(c) {
    _product.reserve(2);
    _product.push_back(IndexedProductElement(op_i,index_i));
    _product.push_back(IndexedProductElement(op_j,index_j));
  }
   
  
  HamiltonianTerm(MatrixElement c,const std::map< Boson* , ProductElement > &m) : _coefficient(c) { 
    _product.reserve(m.size());
    std::map<Boson*,ProductElement>::const_iterator it;
    for(it=m.begin();it!=m.end();++it)
      _product.push_back(IndexedProductElement(it->second,it->first));
      
  }
  
	HamiltonianTerm(const std::vector<IndexedProductElement> &p) : _coefficient(MatrixElement(1.0)), _product(p) {}
 
  HamiltonianTerm(const HamiltonianTerm &o) : _coefficient(o._coefficient), _product(o._product) {}
  
  typedef std::vector<IndexedProductElement>::const_iterator iterator;
	typedef std::vector<IndexedProductElement>::size_type size_type;
  
	inline MatrixElement &coefficient() {return _coefficient;};
  inline const MatrixElement &coefficient() const {return _coefficient;}; 
  inline const std::vector<IndexedProductElement>& product() const {return _product;}
  
  /* The total number of creation minus the number of annihilation operators */
  inline int delta() const {
    int result=0;
    for(size_type i=0;i<_product.size();++i)
      result+=_product[i].delta();
    return result;         
  }
  
  inline bool constant() const { return product().size()==0; }
	inline bool atom() const {return length()==1;}
  
  inline bool diagonal() const { 
    int result=0;
    for(size_type i=0;i<_product.size();++i)
      result+=abs(_product[i].delta());
    return result==0;         
  }

  /* The total number of creation and annihilation operators */
  inline int length() const {
    int result=0;
    for(size_type i=0;i<_product.size();++i)
      result+=_product[i].length();
    return result;         
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
  inline MatrixElement me(int direction) const {
    int result=1;
    for(size_type i=0;i<_product.size();++i)
      result*=_product[i].amplitude(direction);
    return _coefficient*sqrt(result);
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

inline bool operator==(const IndexedProductElement &a,const IndexedProductElement &b) { return a.id() == b.id() && a.particle_id() == b.particle_id(); }

// Define the type of the Kinetic and Potential Part of the Hamiltonian
typedef std::vector<HamiltonianTerm> Hamiltonian;

}

#endif
