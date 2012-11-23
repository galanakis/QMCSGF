#ifndef __PRODUCT_ELEMENT__
#define __PRODUCT_ELEMENT__

#include "Conventions.hh"

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
  
  Update: I added the static function Length which returns the 
  number of creation and annihilation operators in the product.
  Also I added the static function Multiply (and the associated
  operator* ) which multiplies two products. For example it will 
  take 10011 and 110 and return 1001110.
  
*/ 

#include <iostream>
 
namespace SGF {

class ProductElement {
  static uint Length(uint caprod);
  static int Delta(uint caprod);
  static int MaxDelta(uint caprod);
  static uint Amplitude(uint caprod,occupancy_t n);
public:
  static uint Multiply(uint caprod1,uint caprod2);

  uint _caprod;  // This is the binary representation of a product of c/a operators and the only data in the class
  
public:
  ProductElement(int ca=1) : _caprod(ca) {}
  
  inline uint &id() {return _caprod;};
  inline const uint &id() const {return _caprod;};

  inline uint length() const {return Length(_caprod);}
  inline int delta() const {return Delta(_caprod);} 
  inline int maxdelta() const {return MaxDelta(_caprod);}
  inline uint count(int optype) {return (length()+delta()*Sign[optype])/2;}
  inline uint amplitude(int n,int nmax) const {return (nmax!=0 && maxdelta()+n>nmax)?0:ProductElement::Amplitude(_caprod,n);}
  inline int offset(int dn) const {return Abs(dn+delta())-Abs(dn);} // The offset of the operator for occupancy difference dn=NR-NL
  
  friend inline ProductElement operator*(const ProductElement &a,const ProductElement &b);
    
};

inline ProductElement operator*(const ProductElement &a,const ProductElement &b) { return ProductElement::Multiply(a._caprod,b._caprod);}


uint ProductElement::Length(uint caprod) {
  uint result=0;
  while(caprod>1) {
    ++result;
    caprod>>=1;
  }
  return result;
}


uint ProductElement::Amplitude(uint caprod,occupancy_t n) {
  uint amplitude=1;
  while(caprod>1) {
    amplitude*=(n+(caprod&1)); 
    n += Sign[caprod&1];
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
 


/* The maximum intermediate occupancy occurs at the same operator
   regardless of direction. */
int ProductElement::MaxDelta(uint caprod) {
  int max=0;
  int n=0;
  while(caprod>1) {
    n += Sign[caprod&1];
    if(n>max) max=n;
    caprod>>=1;
  } 
  return max;
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

}

#endif
