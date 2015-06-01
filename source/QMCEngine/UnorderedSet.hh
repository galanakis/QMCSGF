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

#ifndef __UNORDEREDSET__
#define __UNORDEREDSET__

#include <algorithm>

namespace SGF {

class UnorderedSet {

protected:
  typedef unsigned long size_type;

  size_type _capacity;
  size_type *data_begin, *data_end;
  size_type **map;

public:

  typedef size_type* iterator;

  inline iterator begin() const {return data_begin;}
  inline iterator end() const {return data_end;}

  inline size_type capacity() const {return _capacity;}
  inline size_type size() const {return data_end-data_begin;}

  UnorderedSet() : _capacity(0), data_begin(0), data_end(0), map(0) {}

  UnorderedSet(size_type size) {initialize(size);}

  void initialize(size_type size) {

    _capacity=size;
    map=new size_type*[_capacity];
    data_end=data_begin=new size_type[_capacity];

    for(size_type i=0; i<_capacity; ++i) {
      data_begin[i]=0;
      map[i]=0;
    }

  }


  void clear() {
    for(iterator it=begin(); it!=end(); ++it) {
      map[*it]=0;
    }
    data_end=data_begin;

  }

  ~UnorderedSet() {
    delete [] data_begin;
    delete [] map;
  }


  inline size_type operator[](size_type i) const {
    return *(data_begin+i);
  }

  inline void insert(size_type i) {

    if(!map[i]) {

      *data_end=i;
      map[i]=data_end;
      ++data_end;

    }
  }

  inline void erase(size_type i) {

    if(map[i]) {

      --data_end;
      *map[i]=*data_end;
      map[*data_end]=map[i];
      map[i]=0;

    }

  }

  inline void sort() {
    std::sort(begin(),end());
    for(iterator it=begin(); it!=end(); ++it)
      map[*it]=it;
  }


};

}

#endif