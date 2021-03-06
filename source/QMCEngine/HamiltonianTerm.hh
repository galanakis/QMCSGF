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

#ifndef __HAMILTONIANTERM__
#define __HAMILTONIANTERM__

#include <vector>
#include <cmath>
#include "AtomicTerm.hh"
#include "Conventions.hh"
#include <string>

namespace SGF {



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

  template<int direction>
  inline occupancy_t &n() {return _occupancy[direction];};
  inline occupancy_t &nL() {return n<LEFT>();};
  inline occupancy_t &nR() {return n<RIGHT>();};
  template<int direction>
  inline const occupancy_t &n() const {return _occupancy[direction];};
  inline const occupancy_t &nL() const {return n<LEFT>();};
  inline const occupancy_t &nR() const {return n<RIGHT>();};
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
  offset=Abs(dn+DN)-Abs(DN),
  where DN=NR-NL, is the difference in occupancy and
  dn is the number of creation minus the number of annihilation operators.
  We can show that if Abs[DN]<=Nmax and Abd[dn]<=Nmax, it takes values
  -dn,-dn+2,...,dn-2,dn,
  so it ranges between -dn and dn with a step of 2.
  If we have products at different indices, the offset is additive.

  This is handled by IndexedProductElement, which can return the minimum and
  maximum offset.
  Given that we can scan through all kinetic terms and generate a list of
  offsets (between min and max with step 2). Then we can map them to
  sequential integers.


*/


/*
  A product of creation (c) and annihilation (a) operators is
  represented by a binary number, caprod, with the 1's mapped
  to c and the zeros to a. There is a placeholder "1" as a leading
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


class ProductElement {

// The population of 1's minus the population of 0's after the leading digit.
  inline static int Delta(uint caprod) {
    int result=0;
    while(caprod>1) {
      result += Sign[caprod&1];
      caprod>>=1;
    }
    return result;
  }

  inline static uint Amplitude(uint caprod,occupancy_t n,occupancy_t nmax) {

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


};


*/



// Deltas: The population of 1's minus the population of 0's after the leading digit.
// MaxDel: Maximum occupancy
// MinDel: Minimum occupancy
//                       0,  1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
typedef enum      {INVALID,  I,   A,   C,  AA,  AC,  CA,  CC, AAA, AAC, ACA, ACC, CAA, CAC, CCA, CCC,AAAA,AAAC,AACA,AACC,ACAA,ACAC,ACCA,ACCC,CAAA,CAAC,CACA,CACC,CCAA,CCAC,CCCA,CCCC} ProductCode;
const int Deltas[32]=   {0,  0,  -1,   1,  -2,   0,   0,   2,  -3,  -1,  -1,   1,  -1,   1,   1,   3,  -4,  -2,  -2,   0,  -2,   0,   0,   2,  -2,   0,   0,   2,   0,   2,   2,   4};
const int MaxDel[32]=   {0,  0,   0,   1,   0,   1,   0,   2,   0,   1,   0,   2,   0,   1,   1,   3,   0,   1,   0,   2,   0,   1,   1,   3,   0,   1,   0,   2,   0,   2,   2,   4};
const int MinDel[32]=   {0,  0,  -1,   0,  -2,   0,  -1,   0,  -3,  -1,  -1,   0,  -2,   0,  -1,   0,  -4,  -2,  -2,   0,  -2,   0,  -1,   0,  -3,  -1,  -1,   0,  -2,   0,  -1,   0};
std::string CAnames[32]= {"", "I","A", "C","AA","AC","CA","CC","AAA","AAC","ACA","ACC","CAA","CAC","CCA","CCC","AAAA","AAAC","AACA","AACC","ACAA","ACAC","ACCA","ACCC","CAAA","CAAC","CACA","CACC","CCAA","CCAC","CCCA","CCCC"};

class IndexedProductElement {


  inline static int Amplitude(uint caprod,occupancy_t n,const occupancy_t nmax) {

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

  ProductCode code;  // This is the binary representation of a product of c/a operators and the only data in the class
  Boson * particle;

public:
  IndexedProductElement(const IndexedProductElement &o) : code(o.code), particle(o.particle) {}
  IndexedProductElement(const ProductCode _code,Boson* const _p) : code(_code),particle(_p) {}
  Boson* const &particle_id() const {return particle;}
  inline const ProductCode &code_id() const {return code;};
  std::string code_name() const {return CAnames[code];}

  template<int action>
  inline int offset() const {
    int dn=particle->n<action>()-particle->n<!action>();
    return Abs(dn+delta())-Abs(dn);
  }

  inline const int &delta() const {
    return Deltas[code];
  }

  template<int direction>
  inline uint amplitude(occupancy_t n,occupancy_t nmax) const {
    return Amplitude(code,n-delta()*!direction,nmax);
  }

  template<int direction>
  inline uint amplitude() const {
    return amplitude<direction>(particle->n<direction>(),particle->nmax());
  }

  template<int direction,int action>
  inline void update() const {

#ifdef DEBUG
    if(particle->n<direction>()==0 && delta()*Sign[action==direction]<0) {
      std::cerr<<std::endl<<"Error in IndexedProductElement::update()"<<std::endl<<"Negative occupancy encountered."<<std::endl;
      exit(20);
    }

    if(particle->nmax()!=0 && particle->n<direction>()==particle->nmax() && delta()*Sign[action==direction]>0) {
      std::cerr<<std::endl<<"Error in IndexedProductElement::update()"<<std::endl<<"Error: Overpopulation encountered."<<std::endl;
      exit(21);
    }

#endif

    particle->n<direction>() += delta()*Sign[action==direction];

  }

};



inline bool operator==(const IndexedProductElement &a,const IndexedProductElement &b) { return a.code_id() == b.code_id() && a.particle_id() == b.particle_id(); }


IndexedProductElement C_(Boson* const _p) {
  return IndexedProductElement(C,_p);
}

IndexedProductElement A_(Boson* const _p) {
  return IndexedProductElement(A,_p);
}

struct IndexedProduct {
  std::vector<IndexedProductElement> vector;
  IndexedProduct(const IndexedProductElement &e) {
    vector.push_back(e);
    vector.shrink_to_fit();
  }
private:
  IndexedProduct() {}
  IndexedProduct(const IndexedProduct &o) : vector(o.vector) {}
  friend struct HTerm;
  friend IndexedProduct operator*(const IndexedProductElement &e,const IndexedProduct &p);
};

struct HTerm : public IndexedProduct {
  MatrixElement coefficient;
  HTerm(const IndexedProduct &p) : IndexedProduct(p), coefficient(1.0) {}
private:
  HTerm() : coefficient(1.0) {}
  friend HTerm operator*(const IndexedProduct &p,MatrixElement c);
};


ProductCode CodeMultiply(ProductCode codeL,ProductCode codeR) {
  int i=0;
  while(codeR>>i > 1)
    ++i;
  return static_cast<ProductCode> ((codeL<<i) | (codeR ^ (1<<i)));
}

IndexedProduct operator*(const IndexedProductElement &e,const IndexedProduct &p) {
  IndexedProduct result(p);
  std::vector<IndexedProductElement>::const_iterator it=result.vector.begin();
  while( it!=result.vector.end() && e.particle_id() > it->particle_id() )
    it++;

  if( it!=result.vector.end() && e.particle_id() == it->particle_id() )
    result.vector.insert(it,IndexedProductElement(CodeMultiply(e.code_id(),it->code_id()),it->particle_id()));
  else 
    result.vector.insert(it,e);

  result.vector.shrink_to_fit();

  return result;

}

HTerm operator*(const IndexedProduct &p,MatrixElement c) {
  HTerm result(p);
  result.coefficient=c;
  return result;
}

}

#endif
