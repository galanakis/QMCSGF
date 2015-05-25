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

#ifndef __CIRCDLIST__
#define __CIRCDLIST__

#include <deque>
#include "Conventions.hh"

namespace SGF {

/*
  class CircDList

  It represents the operator string holding N terms as follows:

  (RIGHT) : 0, 1, 2, ... , N (LEFT)

  where RIGHT and LEFT is relative to the Green Operator,
  that is element 0 of the string is the first element in
  the immediate right of the Green Operator.

  The class is implemented as a wrapper to deque, which replaces
  push_back, push_front, etc with push(direction,data).

  The class members are
  pop(direction): remove an element from the right or left.
  push(direction): place a new element to the right or left
  top(direction): read the top element in the right or left
  top(direction,depth): read an element a distance depth from the right or left
  length(): how many elements there are in the string
  empty(): returns true if there are no elements
  operator[][i]: access directly the elements of the deque.

*/

template<class T>
class CircDList {
protected:
  std::deque<T> que;
public:
  typedef typename std::deque<T>::size_type string_size_type;
  CircDList() : que() {};

  template<int direction>
  inline void pop() { direction==RIGHT ? que.pop_front() : que.pop_back(); }

  template<int direction>
  inline void push(T&& args) { direction==RIGHT ?  que.emplace_front(args) : que.emplace_back(args); }

  template<int direction>
  inline const T& top() const { return direction==RIGHT ? que.front() : que.back(); }

  template<int direction>
  inline const T& top(int depth) const {return direction==RIGHT ? que[depth] : que[que.size()-depth-1];}

  inline string_size_type length() const {return que.size();}
  inline bool empty() const {return que.empty();}
  inline const T &operator[](string_size_type i) const {return que[i];}
  inline T &operator[](string_size_type i) {return que[i];}

  friend struct SGFBase;

};



}


#endif