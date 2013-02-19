#ifndef __SGFBASE__
#define __SGFBASE__

#include "CircularTime.hh"
#include "HamiltonianTerm.hh"
#include "CircDList.hh"
#include "GreenOperator.hh"
#include "Conventions.hh"
#include <iostream>
#include <fstream>
#include <iomanip>

namespace SGF {

/*
	This struct represents an entry in the operator
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

  void write(const std::string &fname) {

    const unsigned int precision=17;

    std::ofstream output;
    output.open(fname.c_str());

    output<<"{"<<std::endl;
    output<<"  \"nsites\": "<<Psi.size()<<","<<std::endl;
    output<<"  \"configuration\": ["<<std::endl;
    std::vector<Boson>::size_type i;
    for(i=0; i<Psi.size()-1; ++i) {
      output<<"    "<<Psi[i].nR()<<","<<std::endl;
    }
    output<<"    "<<Psi[i].nR()<<std::endl;

    output<<"  ],"<<std::endl;

    output<<"  \"nterms\": "<<T.size()<<","<<std::endl;

    output<<"  \"operators\": ["<<std::endl;

    Operator *o;

    for(i=0; i+1<OperatorCDL.que.size(); ++i) {

      o=&OperatorCDL.que[i];

      output<<std::right<<"    ["<<std::setw(6)<<o->Term-&T[0]<<", "<<std::fixed<<std::setprecision(precision)<<std::setw(precision+5)<<std::left<<o->Time.time()<<", "<<std::setw(21)<<o->Energy<<"],"<<std::endl;

    }

    o=&OperatorCDL.que[i];

    output<<std::right<<"    ["<<std::setw(6)<<o->Term-&T[0]<<", "<<std::fixed<<std::setprecision(precision)<<std::setw(precision+5)<<std::left<<o->Time.time()<<", "<<std::setw(21)<<o->Energy<<"]"<<std::endl;

    output<<"  ]"<<std::endl;

    output<<"}"<<std::endl;


  }

};

}

#endif