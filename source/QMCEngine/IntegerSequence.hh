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

#ifndef __INTEGERSEQUENCE__
#define __INTEGERSEQUENCE__

#include <vector>
#include <algorithm>
#include <climits>


/*
  class IntegerSequence

  It holds N pairs of integers and their corresponding index,
  where the index is a sequencial number from 0 to N-1.

  It is designed to hold integers which are not very sparce.

  The index(integer) function will return the location of integer
  or an invalid number if the integer is not in the list.

  The operator [i] with return the i^th smaller integer.

  Both functions don't do any error or bounds checking.

*/

class IntegerSequence {
public:
  typedef std::vector<int> IntegerArray;
  typedef IntegerArray::size_type size_type;
  typedef std::vector<size_type> IntegerMap;
private:
  const IntegerArray sequence;    // Keeps a sorted list of the offsets
  const IntegerMap   indices;     // Map from an offset to consecutive integers

  // Helper function which generates the indices.
  static IntegerMap getIndices(const IntegerArray& o) {

    int Min = *std::min_element(o.begin(), o.end());
    int Max = *std::max_element(o.begin(), o.end());

    size_type nindex = Max - Min + 1;

    IntegerMap result;
    result.reserve(nindex);

    for (size_type i = 0; i < nindex; ++i)
      result.push_back(UINT_MAX);

    int count = 0;
    for (IntegerArray::const_iterator it = o.begin(); it != o.end(); ++it)
      result[(*it) - Min] = count++;

    return result;

  }

public:

  IntegerSequence(const IntegerArray& o) : sequence(o), indices(getIndices(o)) {}

  // Number of different offsets
  inline size_type size() const {
    return sequence.size();
  }
  // IntegerSequence.index(int offset) will return a consecutive integer
  inline const size_type& index(int offset) const {
    return indices[offset - sequence[0]];
  }
  // IntegerSequence[int i] will return the i^th smallest offset
  inline const int& operator[](int i) const {
    return sequence[i];
  }

};

#endif // __INTEGERSEQUENCE__

