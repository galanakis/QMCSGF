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
#include <set>
#include <map>


namespace SGF {

/*
	This struct represents an entry in the operator
	string. It contains
	- a pointer to the Hamiltonian Term
	- The time index of the Hamiltonian Term
	- the energy of the time slice4
*/


typedef AtomicTerm<MatrixElement,IndexedProductElement> HamiltonianTerm;
// Define the type of the Kinetic and Potential Part of the Hamiltonian
typedef std::vector<HamiltonianTerm> Hamiltonian;



typedef std::vector<int> int_vector_t;
typedef std::vector<Boson*> boson_vector_t;
typedef std::set<Boson*> boson_set_t;
typedef std::vector<const HamiltonianTerm *> term_vector_t;
typedef term_vector_t::const_iterator adjacency_iterator;
typedef std::pair<adjacency_iterator,adjacency_iterator> range_type;
typedef std::vector<term_vector_t > AdjacencyList;

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

  //
  // Input: a list of kinetic terms
  // Output: a softed list containing all the possible offsets
  //
  static int_vector_t GetOffsets(const Hamiltonian & T);
  //
  // Input: a list of kinetic terms
  // Output: a softed list containing all the possible bosons
  //
  static boson_vector_t GetConfiguration(const Hamiltonian & T);
  //
  // Input: a list of all the kinetic terms
  // Output: a list of the "extra" kinetic terms needed for the Grand canonical ensemble
  //
  static Hamiltonian GetExtraTerms(const Hamiltonian &T);
  //
  // Input: a softed list containing all the possible bosons
  // Output: a list of the "extra" kinetic terms needed for the Grand canonical ensemble
  //
  static Hamiltonian GetExtraTerms(const boson_vector_t &psi);
  //
  // Input: a vector of boson pointers
  // Output: the number of broken lines in this vector
  //
  static unsigned int CountBrokenLines(const boson_vector_t &psi);
  //
  // It generates the number operator
  //
  static Hamiltonian GenerateNumberOperator(std::vector<Boson> &psi); 
  //
  // It generates the density matrix rho
  //
  static Hamiltonian GenerateDensityMatrix(std::vector<Boson> &psi);
  //
  // Input: two lists of terms, one called row and another called column
  // Output: a vector of size Trow, which contains the indices of the operators
  //         in Tcol which share a common boson
  //
  static AdjacencyList GetAdjacencyList(const Hamiltonian &Trow,const Hamiltonian &Tcol);
  //
  // Input: reference to a Hamiltonian and reference to a Hamiltonian term
  // Output: None.
  // It appends the term to the Hamiltonian. 
  //
  static void AppendHamiltonianTerm(Hamiltonian &H,const HamiltonianTerm &term) {
  if(term.coefficient()!=0)
    H.push_back(term);
  }
  //
  // Creates a Hamiltonian term with a single bosonic component
  //
  static const HamiltonianTerm CreateHamiltonianTerm(MatrixElement c,const ProductCode code,Boson* const p) {
    std::vector<IndexedProductElement> v;
    v.reserve(1);
    v.push_back(IndexedProductElement(code,p));
    return HamiltonianTerm(c,v);
  }
  //
  // Creates a Hamiltonian term with two bosonic components
  // Note: There will be undefined behavior if p1=p2.
  //
  static const HamiltonianTerm CreateHamiltonianTerm(MatrixElement c,const ProductCode code1,Boson* const p1,const ProductCode code2,Boson* const p2) {
    std::vector<IndexedProductElement> v;
    v.reserve(2);
    v.push_back(IndexedProductElement(code1,p1));
    v.push_back(IndexedProductElement(code2,p2));
    return HamiltonianTerm(c,v);
  }

  //
  // Member functions
  //


  //
  // Input: a filename
  // Output: It stores the all the data in the container
  //
  void write(const std::string &fname);





};





//
// Input: a list of kinetic terms
// Output: a softed list containing all the possible offsets
//
int_vector_t SGFBase::GetOffsets(const Hamiltonian & T) {
  int_vector_t offsets;
  // Scan all kinetic terms to find all the lengths
  std::set<int> set_offsets;
  for (Hamiltonian::size_type i = 0; i < T.size(); ++i)
    for (int offset = -T[i].absdelta(); offset <= T[i].absdelta(); offset += 2)
      set_offsets.insert(offset);

  offsets.clear();
  offsets.insert(offsets.begin(), set_offsets.begin(), set_offsets.end());

  return offsets;
}

//
// Input: a list of kinetic terms
// Output: a softed list containing all the possible bosons
//
boson_vector_t SGFBase::GetConfiguration(const Hamiltonian & T) {
  boson_vector_t indices;
  boson_set_t indexset;
  for (Hamiltonian::size_type i = 0; i < T.size(); ++i)
    for (HamiltonianTerm::const_iterator jt = T[i].begin(); jt != T[i].end(); ++jt)
      indexset.insert( jt->particle_id());

  indices.insert(indices.begin(), indexset.begin(), indexset.end());

  return indices;

}

//
// Input: a softed list containing all the possible bosons
// Output: a list of the "extra" kinetic terms needed for the Grand canonical ensemble
//
Hamiltonian SGFBase::GetExtraTerms(const boson_vector_t &psi) {
  Hamiltonian result;
  for(boson_vector_t::size_type i=0; i<psi.size(); ++i) {
    result.push_back(CreateHamiltonianTerm(1.0,A,psi[i]));
    result.push_back(CreateHamiltonianTerm(1.0,C,psi[i]));
  }

  return result;

}

//
// Input: a list of all the kinetic terms
// Output: a list of the "extra" kinetic terms needed for the Grand canonical ensemble
//
Hamiltonian SGFBase::GetExtraTerms(const Hamiltonian &T) {

  return GetExtraTerms(GetConfiguration(T));

}


//
// Input: a vector of boson pointers
// Output: the number of broken lines in this vector
//
unsigned int SGFBase::CountBrokenLines(const boson_vector_t &o) {
  int result=0;
  for(boson_vector_t::size_type i=0; i<o.size(); ++i)
    result+=Abs(o[i]->delta());
  return result;

}


//
// Input: two lists of terms, one called row and another called column
// Output: a vector of size Trow, which contains the indices of the operators
//         in Tcol which share a common boson
//
AdjacencyList SGFBase::GetAdjacencyList(const Hamiltonian &Trow,const Hamiltonian &Tcol)  {

  AdjacencyList adjacency; // The adjacency list is stored here

  // First, categorize the Tcol terms by index.
  std::map<Boson*,std::set<Hamiltonian::size_type> > map_to_set;
  for(Hamiltonian::size_type i=0; i<Tcol.size(); ++i)
    for(HamiltonianTerm::const_iterator jt=Tcol[i].begin(); jt!=Tcol[i].end(); ++jt)
      map_to_set[jt->particle_id()].insert(i);

  adjacency.clear();
  adjacency.resize(Trow.size());

  /* For each term in Trow, merge the sets corresponding to it's indices
     Then, copy the set elements to a vector. */
  for(Hamiltonian::size_type i=0; i<Trow.size(); ++i) {
    std::set<Hamiltonian::size_type> merged;
    for(HamiltonianTerm::const_iterator jt=Trow[i].begin(); jt!=Trow[i].end(); ++jt) {
      if(jt->delta()!=0) {
        std::set<Hamiltonian::size_type> &s=map_to_set[jt->particle_id()];
        merged.insert(s.begin(),s.end());
      }
    }

    adjacency[i].reserve(merged.size());
    for(std::set<Hamiltonian::size_type>::const_iterator sit=merged.begin(); sit!=merged.end(); ++sit)
      adjacency[i].push_back(&Tcol[*sit]);

  }

  return adjacency;

}

//
// It generates the number operator
//
Hamiltonian SGFBase::GenerateNumberOperator(std::vector<Boson> &psi) {
  Hamiltonian result;
  for(boson_vector_t::size_type i=0; i<psi.size(); ++i)
    result.push_back(CreateHamiltonianTerm(1.0,CA,&psi[i]));
  return result;
}

//
// It generates the density matrix rho
//
Hamiltonian SGFBase::GenerateDensityMatrix(std::vector<Boson> &psi) {
  Hamiltonian result;
  for(boson_vector_t::size_type i=0; i<psi.size(); ++i) {
    for(boson_vector_t::size_type j=0; j<psi.size(); ++j) {
      result.push_back(CreateHamiltonianTerm(1.0,C,&psi[i],A,&psi[j]));
    }
  }
  return result;
}


void SGFBase::write(const std::string &fname) {

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



}

#endif