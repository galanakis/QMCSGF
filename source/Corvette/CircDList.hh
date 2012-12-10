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
  inline void pop(int direction) { direction==RIGHT ? que.pop_front() : que.pop_back(); }
  inline void push(int direction,const T& data) { direction==RIGHT ?  que.push_front(data) : que.push_back(data); }
  inline const T& top(int direction) const { return direction==RIGHT ? que.front() : que.back(); }
	inline const T& top(int direction,int depth) const {return direction==RIGHT ? que[depth] : que[que.size()-depth-1];}
  inline string_size_type length() const {return que.size();}
  inline bool empty() const {return que.empty();}
	inline const T &operator[](string_size_type i) const {return que[i];}
}; 



}


#endif