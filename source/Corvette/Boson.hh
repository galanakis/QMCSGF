#ifndef __PARTICLE__
#define __PARTICLE__
 
#include "Conventions.hh"
#include <vector>

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
  inline occupancy_t &n(int direction) {return _occupancy[direction];};
  inline occupancy_t &nmax() {return _maximum_occupancy;};
  inline int delta() const {return _occupancy[RIGHT]-_occupancy[LEFT];};
  
};
 

}

#endif
