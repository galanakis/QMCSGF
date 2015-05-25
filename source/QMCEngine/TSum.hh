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

#ifndef __TSUM__
#define __TSUM__

#include "RandomNumberGenerator.hh"
#include <vector>
#include "UnorderedSet.hh"


namespace SGF {

/*
class TSum
This is a helper class which consists of a vector of matrix elements.
The class interprets this vector as a balanced binary tree. To get the
path of the i^th element we just look at the binary representation of i+1,
where the most significant digit corresponds to ancestor nodes.
The number of elements of the tree are meant to be fixed. Each
element contains a value which is interpreted as a relative probability
At each node we don't need to store that relative probability, but
only the total sum of the probability of the node and all its children.
The probability can be deduced by subtracting the sum of a node
minus the sums of its children. When a probability is changed this change
is propagated to the node's parents with logarithmic complexity.
Selecting a term according to its relative probability is also done
with logarithmic complexity algorithm. We start from the root of the
tree and select either a left or right branch or the root node itself
(three way selection). We continue in the branch we selected or return
the root.

In this class N terms are stored in a tree with 2*P-1 where P=2^d is the
smallest power of 2 which is greater than N. The partial sums are stored
in the first P-1 terms and matrix elements of the terms are in the last P
terms. The size of the tree is thus size=2*P-1 and the terms start
at size/2=(2*(P-1)+1)/2=P-1.

The key functions are:
update(index,me): change the probability of the index by "me".
choose(): randomly chose an index using its relative probability.
norm(): returns the sum of all relative probabilities

========== Usage =========
TSum t;
t.resize(100);          // Keep the relative probabilities of 100 indices.
t.update(10,25);        // set the relative probability of index 10 to be 25.
double p=t.element(10); // return the relative probability of index 10.
double norm=t.norm();   // the sum of all relative probabilities
int index=t.choose();   // pick an index at random using their relative probabilities.
t.reset();              // clear everything.

*/


template<typename _Float, typename _DoubleFloat>
class TSumBase {
  typedef std::vector<_Float> FloatArray;
  typedef typename FloatArray::iterator Iterator;
public:
  typedef _Float Float;
  typedef _DoubleFloat DoubleFloat;
  typedef typename FloatArray::size_type IndexType;
private:
  FloatArray _sums;
  Float* _elements;
  IndexType _nsums, _base;
  DoubleFloat _norm;

  IndexType _nterms;
  bool* _guard;
  IndexType* _buffer, *_buffer_entry;

  /* Makes all elements zero */
  inline void reset() {
    _buffer_entry = _buffer;
    _norm = 0;
    for (IndexType i = 0; i < _nterms; ++i)
      _guard[i] = true;
    for (Iterator it = _sums.begin(); it != _sums.end(); ++it)
      *it = Float(0);
  }

  inline void flush() {

    while (_buffer_entry > _buffer) {
      --_buffer_entry;

      _guard[*_buffer_entry] = true;
      IndexType index = *_buffer_entry / 2 + _base;
      Float newme = _sums[2 * index - 1] + _sums[2 * index];
      Float oldme = _sums[index - 1];

      if (newme != oldme) {
        while (index > 1) {
          _sums[index - 1] =  newme;
          newme += _sums[(index ^ 1) - 1];
          index /= 2;
        }
        _sums[0] = newme;
      }



    }

    _norm = _sums[0];
  }

public:

  TSumBase() : _sums(), _elements(0), _norm(0), _guard(0) {}
  ~TSumBase() {
    delete [] _guard;
    delete [] _buffer;
  }

  /* resizes the vector and sets all elements to zero */
  inline void resize(IndexType NTerms) {
    IndexType _size = 1;
    while (_size < NTerms) _size <<= 1;
    _sums.resize(2 * _size - 1, Float(0));
    _nsums = _size - 1;
    _base = (1 + _nsums) / 2;
    _elements = &_sums[_nsums];
    _nterms = NTerms;
    _guard = new bool[_nterms];
    _buffer = new IndexType[_nterms];
    _buffer_entry = _buffer;
    _norm = 0;

    reset();
  }

  inline void update(IndexType index, Float me) {
    if (me != _elements[index]) {
      if (_guard[index]) {
        *_buffer_entry = index;
        ++_buffer_entry;
        _guard[index] = false;
      }
      _norm += me - _elements[index];
      _elements[index] = me;
    }
  }

  inline const Float &element(IndexType index) const {
    return _elements[index];
  }

  inline IndexType choose() {
    flush();

    Float w, wl;
    IndexType index = 0;
    while (index < _nsums) {
      w = _sums[index];
      index = 2 * index + 1;
      wl = _sums[index];
      index += !(w * RNG::Uniform() < wl);
    }

    return index - _nsums;

  }

  const DoubleFloat &norm() const {
    return _norm;
  }

};

template<typename Float, typename DoubleFloat>
class TSumHC {

  UnorderedSet list;
  DoubleFloat _norm;
  std::vector<Float>  _elements;
  typedef unsigned long IndexType;

  /* Makes all elements zero */
  inline void reset() {
    _norm = 0;
    for (IndexType i = 0; i < _elements.size(); ++i)
      _elements[i] = 0;
    list.clear();
  }

public:
  TSumHC() : _norm(0), list(10000) {}
  inline void resize(IndexType NTerms) {
    _elements.resize(NTerms);
    reset();
  }
  const DoubleFloat &norm() const {
    return _norm;
  }
  inline IndexType choose() {
    return list.element(RNG::Uniform() * list.size());
  }

  inline const Float &element(IndexType index) const {
    return _elements[index];
  }

  inline void update(IndexType index, Float me) {
    if (me != _elements[index]) {
      _norm += me - _elements[index];
      _elements[index] = me;
      if (me == 0) {
        list.erase(index);
      } else {
        list.insert(index);
      }
    }
  }

};

}

#endif
