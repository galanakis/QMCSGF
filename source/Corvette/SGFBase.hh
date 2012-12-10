#ifndef __SGFBASE__
#define __SGFBASE__

#include "CircularTime.hh"
#include "HamiltonianTerm.hh"
#include "CircDList.hh"
#include "GreenOperator.hh"
#include "Conventions.hh"

namespace SGF {

/*
	This struct represents a an entry in the operator
	string. It contains
	- a pointer to the Hamiltonian Term
	- The time index of the Hamiltonian Term
	- the energy of the time slice4
*/

struct Operator {
  CircularTime Time;
  const HamiltonianTerm* Term;
	_float_accumulator Energy;
  
  Operator(const Operator &o) : Time(o.Time), Term(o.Term), Energy(o.Energy) {}
	Operator &operator=(const Operator &o) {Time=o.Time; Term=o.Term; Energy=o.Energy; return *this;}
  Operator(const CircularTime &_time,const HamiltonianTerm *_term,const _float_accumulator &_energy) : Time(_time), Term(_term), Energy(_energy) {}
};                       


typedef CircDList<Operator> OperatorCircDlist;

//
// This class contains all the independent data
//
struct SGFBase {

   std::vector<Boson> Psi;
   OperatorCircDlist OperatorCDL;
   Hamiltonian T,V;
   GreenOperator<long double> g;
   double Alpha;
   double Beta;
   ensemble_t Ensemble;

};

}

#endif