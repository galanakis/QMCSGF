#ifndef __SGFBASE__
#define __SGFBASE__

#include "HamiltonianTerm.hh"
#include "CircDList.hh"
#include "GreenOperator.hh"
#include "Conventions.hh"

namespace SGF {

//
// This class contains all the independent data
//
struct SGFBase {

   std::vector<Boson> Psi;
   CircDList<Operator> OperatorCDList;
   Hamiltonian T,V;
   GreenOperator<long double> g;
   double Alpha;
   double Beta;
   ensemble_t Ensemble;

};


}

#endif