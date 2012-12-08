#ifndef __CIRCDLIST__
#define __CIRCDLIST__

#include <deque>

#include "CircularTime.hh"
#include "HamiltonianTerm.hh"

namespace SGF {

/*
  class CircDList
  This is a wrapper to deque, which replaces push_back, push_front, etc
  with push(direction,data). 
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


  
struct Operator {
  CircularTime Time;
  const HamiltonianTerm* Term;
	_float_accumulator Energy;
  
  Operator(const Operator &o) : Time(o.Time), Term(o.Term), Energy(o.Energy) {}
	Operator &operator=(const Operator &o) {Time=o.Time; Term=o.Term; Energy=o.Energy; return *this;}
  Operator(const CircularTime &_time,const HamiltonianTerm *_term,const _float_accumulator &_energy) : Time(_time), Term(_term), Energy(_energy) {}
};                       


typedef CircDList<Operator> OperatorCircDlist;



}


#endif