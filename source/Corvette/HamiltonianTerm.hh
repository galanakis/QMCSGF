#ifndef __HAMILTONIANTERM__
#define __HAMILTONIANTERM__

#include <vector>
#include <cmath>
#include "AtomicTerm.hh"
#include "Conventions.hh"

namespace SGF {


/*
  class ProductElement is a container of utilities
  that are used when dealing with products of creation (c)
  annihilation (a) operators. Such a product is represented
  by a binary number, caprod, with the 1's mapped to c and 
  the zeros to a. There is a placeholder "1" as a leading
  digit. For example 101101 represents accac.
  Such an operator product is applied to an occupancy n,
  either to the left or to the right. There are functions
  to calculate the amplitude, the change in occupancy, delta,
  and the maximum intermediate occupancy, maxdelta, reached 
  while applying one operator after another.

  For the implementation we define the static functions  
  Delta, MaxDelta and Amplitude
  only for the case of applying the operator to the right
  and we build all the other cases on top of those.
  
  We use the following naming convension
  uint caprod; // the binary number representing the operator product
  int n; // an occupancy state
  int direction; // Takes values 1 (right) or 0 (left)
  int nmax: // the maximum allowed occupancy. no limit --> nmax=0.
  
  Caching the values of delta and maxdelta offers some performance 
  gain, but not much (about 30%). However in order to keep this
  possibility open, I made the Delta, MaxDelta and Amplitude
  static, such that they can be tabulated more easily.
  
  CAVEAT: because the product element is represented by a single number, 
  the maximum number of operators that can be represented depend on the 
  hardware, but they are typically more than 31 and nowadays 63 
  (the first bit only marks the end of the chain). I don't think 
  this is a practical limitation.

  Also I added the static function Multiply (and the associated
  operator* ) which multiplies two products. For example it will 
  take 10011 and 110 and return 1001110.
  
*/ 


class ProductElement {
  inline static int Delta(uint caprod);
  inline static uint Amplitude(uint caprod,occupancy_t n,occupancy_t nmax);
  static uint Multiply(uint caprod1,uint caprod2);

  uint _caprod;  // This is the binary representation of a product of c/a operators and the only data in the class

public:
  
  ProductElement(int ca=1) : _caprod(ca) {}
  ProductElement(const ProductElement &o) : _caprod(o._caprod) {}
  
  inline const uint &id() const {return _caprod;};

  inline int delta() const {return Delta(_caprod);} 
  inline uint amplitude(int n,int nmax) const {return Amplitude(_caprod,n,nmax); }
  inline int offset(int dn) const {return Abs(dn+delta())-Abs(dn);} // The offset of the operator for occupancy difference dn=NR-NL
  
  inline ProductElement &operator*=(const ProductElement &o) {

     _caprod=Multiply(_caprod,o._caprod);
     return *this;

  }

  inline int maxoffset() const { return Abs(delta()); }
  inline int minoffset(uint nmax) const {
    uint dn=Abs(delta()); 
    return (nmax!=0)?(dn-2*Min(dn,nmax)):(-dn);
  }


};

inline ProductElement operator*(const ProductElement &a,const ProductElement &b) { return ProductElement(a)*=b;}


uint ProductElement::Amplitude(uint caprod,occupancy_t n,const occupancy_t nmax) {

  if(nmax!=0 && n>nmax)
    return 0;
  uint amplitude=1;
  while(caprod>1) {

    amplitude*=n+(caprod&1);

    n += Sign[caprod&1];
    if(nmax!=0 && n>nmax)
      return 0;
    caprod>>=1;
  }
  return amplitude;

}

 

/* The population of 1's minus the population of 0's after the leading digit. */
int ProductElement::Delta(uint caprod) {
  int result=0;
  while(caprod>1) {
    result += Sign[caprod&1];
    caprod>>=1;
  }
  return result;
}
 
 

/* The product of two products. For example 10*11=101. */
uint ProductElement::Multiply(uint caprod1,uint caprod2) {
  uint log=1;
  while(2*log<=caprod2) log<<=1;
  return caprod1*log+(caprod2&(log-1));
}


/* Shortcuts from which we can create any ProductElement by using the "*" operator. */
const ProductElement A(2);
const ProductElement C(3);



/*
  Class Boson
  It represents a single degree of freedom of a boson.
  In the context of SGF it is described by two occupancies,
  one for the left and one for the right wave function
  and a common maximum occupancy.
  
*/

class Boson {
protected:
  occupancy_t _maximum_occupancy; // 1 for hardcore bosons and infinity for bosons
  occupancy_t _occupancy[2];      // The left and right occupancy 
public:
  Boson() : _maximum_occupancy(0) {
    _occupancy[0]=_occupancy[1]=0;
  };
  Boson(const Boson &o) { *this=o; }
  Boson &operator=(const Boson &o) {
    _maximum_occupancy=o._maximum_occupancy;
    _occupancy[0]=o._occupancy[0];
    _occupancy[1]=o._occupancy[1];
    return *this;    
  }
  inline occupancy_t &n(int direction) {return _occupancy[direction];};
  inline const occupancy_t &n(int direction) const {return _occupancy[direction];};
  inline occupancy_t &nmax() {return _maximum_occupancy;};
  inline const occupancy_t &nmax() const {return _maximum_occupancy;};
  inline int delta() const {return _occupancy[RIGHT]-_occupancy[LEFT];};
  
};
 





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
  

  Concerning the minimum and maximum offsets:

  For an operator product at the same index, the offset is defined as
  Abs(dn+NR-NL)-Abs(NR-NL) and for a given dn it takes values -dn,-dn+2,...,dn-2,dn.
  If we have products at different indices, the offset is additive. This reminds a
  bit addition of spins. The maximum offset is the Abs(delta())
  (# of creation -# of annihilation operators) of the product.
  The step of the offset is always 2. Now if we have a lot of terms, some of them
  will have an even and some of them will have an odd and some others and even
  offset. If all of them are even then we can set the offset step=2, otherwise
  it will be 1.

  All this does not assume the case of hard core bosons, in which the number
  of offsets is limitted.
  In this case we need to consider Abs(dn+DN)-Abs(DN), -Nmax<=(DN=NR-NL)<=Nmax.
  It is easy tos show that given the constraint:
  |dn|-2*Min(|dn|,Nmax) <= Abs(|dn|+DN)-Abs(DN) <= |dn|

  This is handled by IndexedProductElement, which can return the minimum and
  maximum offset.
  Given that we can scan through all kinetic terms and generate a list of
  offsets (between min and max with step 2). Then we can map them to
  sequential integers.




*/

class IndexedProductElement : public ProductElement {
 	
	
protected:
  Boson* particle;
public:
	IndexedProductElement(const IndexedProductElement &o) : ProductElement(o), particle(o.particle) {}
  IndexedProductElement(const ProductElement &peprod,Boson* const _p) : ProductElement(peprod),particle(_p) {}
  Boson* const &particle_id() const {return particle;}
  inline int offset(int action=ADD) const {return ProductElement::offset(Sign[action]*particle->delta());}
  inline uint amplitude(int direction) const {return ProductElement::amplitude(particle->n(direction)-delta()*!direction,particle->nmax());}
  inline void update(int direction,int action) const {

#ifdef DEBUG    
    if(particle->n(direction)==0 && delta()*Sign[action==direction]<0) {
      std::cerr<<"Error: Negative occupancy encountered."<<std::endl;
      exit(2);
    }

    if(particle->nmax()!=0 && particle->n(direction)==particle->nmax() && delta()*Sign[action==direction]>0) {
      std::cerr<<"Error: Overpopulation encountered."<<std::endl;
      exit(2);
    }

#endif

    particle->n(direction) += delta()*Sign[action==direction];
     
  }

  inline int minoffset() const { return ProductElement::minoffset(particle->nmax()); }

};

inline bool operator==(const IndexedProductElement &a,const IndexedProductElement &b) { return a.id() == b.id() && a.particle_id() == b.particle_id(); }

typedef AtomicTerm<MatrixElement,IndexedProductElement> HamiltonianTerm;
// Define the type of the Kinetic and Potential Part of the Hamiltonian
typedef std::vector<HamiltonianTerm> Hamiltonian;

}

#endif
