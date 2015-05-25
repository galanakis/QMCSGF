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
#include <limits>
#include <cassert>
#include <sstream>


namespace SGF {

/*
  This struct represents an entry in the operator
  string. It contains
  - a pointer to the Hamiltonian Term
  - The time index of the Hamiltonian Term
  - the energy of the time slice4
*/


typedef AtomicTerm<MatrixElement, IndexedProductElement> HamiltonianTerm;
// Define the type of the Kinetic and Potential Part of the Hamiltonian
typedef std::vector<HamiltonianTerm> Hamiltonian;



typedef std::vector<int> int_vector_t;
typedef std::vector<Boson*> boson_vector_t;
typedef std::set<Boson*> boson_set_t;
typedef std::vector<const HamiltonianTerm*> term_vector_t;
typedef term_vector_t::const_iterator adjacency_iterator;
typedef std::pair<adjacency_iterator, adjacency_iterator> range_type;
typedef std::vector<term_vector_t > AdjacencyList;

struct Operator {
  CircularTime Time;
  const HamiltonianTerm* Term;
  _float_accumulator Energy;

  Operator(const CircularTime& _time, const HamiltonianTerm* _term, const _float_accumulator& _energy) : Time(_time), Term(_term), Energy(_energy) {}
};


typedef CircDList<Operator> OperatorCircDlist;



//
// This class contains all the independent data
//
struct SGFBase {


  std::vector<Boson> Psi;
  OperatorCircDlist OperatorCDL;
  Hamiltonian T, V;
  GreenOperator<long double> g;
  double Alpha;
  double Beta;
  ensemble_t Ensemble;
  double ExtraTermProbability;

  //
  // Input: a list of kinetic terms
  // Output: a softed list containing all the possible offsets
  //
  static int_vector_t GetOffsets(const Hamiltonian& T);
  //
  // Input: a list of kinetic terms
  // Output: a softed list containing all the possible bosons
  //
  static boson_vector_t GetConfiguration(const Hamiltonian& T);
  //
  // Input: a list of all the kinetic terms
  // Output: a list of the "extra" kinetic terms needed for the Grand canonical ensemble
  //
  static Hamiltonian GetExtraTerms(const Hamiltonian& T);
  //
  // Input: a softed list containing all the possible bosons
  // Output: a list of the "extra" kinetic terms needed for the Grand canonical ensemble
  //
  static Hamiltonian GetExtraTerms(const boson_vector_t& psi);
  //
  // Input: a vector of boson pointers
  // Output: the number of broken lines in this vector
  //
  static unsigned int CountBrokenLines(const boson_vector_t& psi);
  //
  // It generates the number operator
  //
  static Hamiltonian GenerateNumberOperator(std::vector<Boson>& psi);
  //
  // It generates the density matrix rho
  //
  static Hamiltonian GenerateDensityMatrix(std::vector<Boson>& psi);
  //
  // It generates the density-density matrix <n_i n_j>
  //
  static Hamiltonian GenerateDensityDensityMatrix(std::vector<Boson>& psi);
  //
  // From the list of bosons it generates the extra terms
  // C^{\dagger}, C
  //
  static Hamiltonian GenerateExtraTerms(std::vector<Boson>& psi);
  //
  // Input: a vector of boson pointers
  // Output: the number of broken lines in this vector
  //
  static unsigned int CountBrokenLines2(const std::vector<Boson>& o);
  //
  // Input: two lists of terms, one called row and another called column
  // Output: a vector of size Trow, which contains the indices of the operators
  //         in Tcol which share a common boson
  //
  static AdjacencyList GetAdjacencyList(const Hamiltonian& Trow, const Hamiltonian& Tcol);
  //
  // Input: reference to a Hamiltonian and reference to a Hamiltonian term
  // Output: None.
  // It appends the term to the Hamiltonian.
  //
  static void AppendHamiltonianTerm(Hamiltonian& H, const HamiltonianTerm& term) {
    if (term.coefficient() != 0)
      H.push_back(term);
  }
  //
  // Creates a Hamiltonian term with a single bosonic component
  //
  static const HamiltonianTerm CreateHamiltonianTerm(MatrixElement c, const ProductCode code, Boson* const p) {
    std::vector<IndexedProductElement> v;
    v.reserve(1);
    v.push_back(IndexedProductElement(code, p));
    return HamiltonianTerm(c, v);
  }
  //
  // Creates a Hamiltonian term with two bosonic components
  // Note: There will be undefined behavior if p1=p2.
  //
  static const HamiltonianTerm CreateHamiltonianTerm(MatrixElement c, const ProductCode code1, Boson* const p1, const ProductCode code2, Boson* const p2) {
    std::vector<IndexedProductElement> v;
    v.reserve(2);
    assert(p1 != p2); // Error if both operators act on the same boson
    v.push_back(IndexedProductElement(code1, p1));
    v.push_back(IndexedProductElement(code2, p2));
    return HamiltonianTerm(c, v);
  }

  //
  // Member functions
  //


  //
  // Input: a filename
  // Output: It stores the all the data in the container
  //

  void write(const std::string& fname) {

    std::ofstream o;
    o.open(fname.c_str());

    o << "{" << std::endl;
    o << "  \"nsites\": " << Psi.size() << "," << std::endl;
    o << "  \"nterms\": " << T.size() << "," << std::endl;
    json_print(o, "configuration", Psi, 1);
    o << "," << std::endl << std::endl;
    json_print(o, "String", OperatorCDL, 1);
    o << std::endl;
    o << "}" << std::endl;

  }

  //
  // Input: a filename
  // Output: It stores the all the data in the container including the Hamiltonian
  //
  void json_print_hamiltonian(std::ostream& o) {
    o << "{" << std::endl;
    json_print(o, "Kinetic", T, 1);
    o << "," << std::endl << std::endl;
    json_print(o, "Potential", V, 1);
    o << "}" << std::endl;
  }


  std::string json_print(const HamiltonianTerm& term) {

    std::stringstream ss;
    ss << "[" << std::setw(7) << term.coefficient() << ",";
    HamiltonianTerm::const_iterator ip = term.begin();
    while (ip != term.end()) {
      ss << std::setw(7) << '\"' + ip->code_name() + '\"' << "," << std::setw(6) << ip->particle_id() - &Psi[0];
      ++ip;
      if (ip != term.end())
        ss << ",";
    }
    ss << "]";
    return ss.str();
  }

  std::string json_print(const Operator* op) {
    std::stringstream o;
    const unsigned int timeprecision = std::numeric_limits<circulartime_t>::digits10 + 1;
    unsigned int timewidth = timeprecision + 4;

    const unsigned int energyprecision = std::numeric_limits<_float_accumulator>::digits10 + 1;
    unsigned int energywidth = energyprecision + 4;

    unsigned int intwidth = std::numeric_limits<Hamiltonian::size_type>::digits10;

    o  << std::right << "  ["
       << std::setw(intwidth) << op->Term - &T[0] << ", "
       << std::setprecision(timeprecision) << std::setw(timewidth) << std::left << op->Time << ", "
       << std::setprecision(energyprecision) << std::setw(energywidth) << op->Energy
       << "]";

    return o.str();

  }


  std::ostream& json_print(std::ostream& o, std::string tag, const Hamiltonian& H, unsigned int depth) {

    std::string indent(2 * depth, ' ');
    o << indent << '\"' << tag << '\"' << ": [" << std::endl;
    Hamiltonian::const_iterator term = H.begin();
    while (term != H.end()) {
      o << indent << "  " << json_print(*term);
      ++term;
      if (term != H.end())
        o << ",";
      o << std::endl;
    }
    o << indent << "]";
    return o;

  }

  std::ostream& json_print(std::ostream& o, std::string tag, std::vector<Boson>& psi, unsigned int depth) {

    std::string indent(2 * depth, ' ');
    o << indent << '\"' + tag + '\"' << ": [" << std::endl;
    std::vector<Boson>::const_iterator it = psi.begin();
    while (it != psi.end()) {
      o << indent << "  [ " << it->nR() << ", " << it->nmax() << " ]";
      ++it;
      if (it != psi.end())
        o << ',';
      o << std::endl;
    }
    o << indent << "]";
    return o;
  }

  std::ostream& json_print(std::ostream& o, std::string tag, OperatorCircDlist& ostring, unsigned int depth) {

    std::string indent(2 * depth, ' ');
    o << indent << '\"' + tag + '\"' << ": [" << std::endl;
    for (unsigned long i = 0; i < OperatorCDL.que.size(); ++i) {
      o << indent << json_print(&ostring.que[i]);
      if (i + 1 != ostring.que.size() )
        o << ",";
      o << indent << std::endl;
    }

    o << indent << "]";

    return o;
  }

  void write_verbose(const std::string& fname, unsigned int depth = 0)  {
    std::ofstream output;
    output.open(fname.c_str());
    write_verbose(output, depth);
  }

  void write_verbose(std::ofstream& o, unsigned int depth = 0) {

    std::string indent(2 * depth, ' ');
    const unsigned int prec = std::numeric_limits<double>::digits10 + 1;
    o << indent << "{" << std::endl;
    o << indent << std::setprecision(prec) << "  \"Temperature\": " << 1.0 / Beta << "," << std::endl;
    o << std::endl;
    json_print(o, "configuration", Psi, 1 + depth);
    o << "," << std::endl;
    o << std::endl;
    json_print(o, "Kinetic", T, 1 + depth);
    o << "," << std::endl;
    o << std::endl;
    json_print(o, "Potential", V, 1 + depth);
    o << "," << std::endl;
    o << std::endl;
    json_print(o, "String", OperatorCDL, 1 + depth);
    o << std::endl;
    o << indent << "}";

  }

  void write_as_input(const std::string& fname)  {
    std::ofstream output;
    output.open(fname.c_str());
    write_as_input(output);
  }

  void write_as_input(std::ofstream& o) {
    o << "{" << std::endl;
    o << "  \"Model\": ";
    write_verbose(o, 1);
    o << "," << std::endl;
    o << "  \"SGF\": {" << std::endl;

    const unsigned int prec = std::numeric_limits<double>::digits10;
    o << "    \"Alpha\": " << std::setprecision(prec) << Alpha;
    o << "," << std::endl;
    o << "    \"Enseble\": ";
    if (Ensemble == Canonical)
      o << "\"Canonical\"";
    else
      o << "\"GrandCanonical\"";
    o << "," << std::endl;
    o << "    \"ExtraTermProbability\": " << std::setprecision(prec) << ExtraTermProbability;
    o << std::endl;
    o << "  }" << std::endl;

    o << "}" << std::endl;


  }


};





//
// Input: a list of kinetic terms
// Output: a softed list containing all the possible offsets
//
int_vector_t SGFBase::GetOffsets(const Hamiltonian& T) {
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
boson_vector_t SGFBase::GetConfiguration(const Hamiltonian& T) {
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
Hamiltonian SGFBase::GetExtraTerms(const boson_vector_t& psi) {
  Hamiltonian result;
  for (boson_vector_t::size_type i = 0; i < psi.size(); ++i) {
    result.push_back(CreateHamiltonianTerm(1.0, A, psi[i]));
    result.push_back(CreateHamiltonianTerm(1.0, C, psi[i]));
  }

  return result;

}

//
// Input: a list of all the kinetic terms
// Output: a list of the "extra" kinetic terms needed for the Grand canonical ensemble
//
Hamiltonian SGFBase::GetExtraTerms(const Hamiltonian& T) {

  return GetExtraTerms(GetConfiguration(T));

}


//
// Input: a vector of boson pointers
// Output: the number of broken lines in this vector
//
unsigned int SGFBase::CountBrokenLines(const boson_vector_t& o) {
  int result = 0;
  for (boson_vector_t::size_type i = 0; i < o.size(); ++i)
    result += Abs(o[i]->delta());
  return result;

}


//
// Input: two lists of terms, one called row and another called column
// Output: a vector of size Trow, which contains the indices of the operators
//         in Tcol which share a common boson
//
AdjacencyList SGFBase::GetAdjacencyList(const Hamiltonian& Trow, const Hamiltonian& Tcol)  {

  AdjacencyList adjacency; // The adjacency list is stored here

  // First, categorize the Tcol terms by index.
  std::map<Boson*, std::set<Hamiltonian::size_type> > map_to_set;
  for (Hamiltonian::size_type i = 0; i < Tcol.size(); ++i)
    for (HamiltonianTerm::const_iterator jt = Tcol[i].begin(); jt != Tcol[i].end(); ++jt)
      map_to_set[jt->particle_id()].insert(i);

  adjacency.clear();
  adjacency.resize(Trow.size());

  /* For each term in Trow, merge the sets corresponding to it's indices
     Then, copy the set elements to a vector. */
  for (Hamiltonian::size_type i = 0; i < Trow.size(); ++i) {
    std::set<Hamiltonian::size_type> merged;
    for (HamiltonianTerm::const_iterator jt = Trow[i].begin(); jt != Trow[i].end(); ++jt) {
      if (jt->delta() != 0) {
        std::set<Hamiltonian::size_type>& s = map_to_set[jt->particle_id()];
        merged.insert(s.begin(), s.end());
      }
    }

    adjacency[i].reserve(merged.size());
    for (std::set<Hamiltonian::size_type>::const_iterator sit = merged.begin(); sit != merged.end(); ++sit)
      adjacency[i].push_back(&Tcol[*sit]);

  }

  return adjacency;

}

//
// It generates the number operator
//
Hamiltonian SGFBase::GenerateNumberOperator(std::vector<Boson>& psi) {
  Hamiltonian result;
  for (boson_vector_t::size_type i = 0; i < psi.size(); ++i)
    result.push_back(CreateHamiltonianTerm(1.0, CA, &psi[i]));
  return result;
}

//
// It generates the density matrix rho
//
Hamiltonian SGFBase::GenerateDensityMatrix(std::vector<Boson>& psi) {
  Hamiltonian result;
  for (boson_vector_t::size_type i = 0; i < psi.size(); ++i) {
    for (boson_vector_t::size_type j = 0; j < psi.size(); ++j) {
      if (i != j)
        result.push_back(CreateHamiltonianTerm(1.0, C, &psi[i], A, &psi[j]));
      else
        result.push_back(CreateHamiltonianTerm(1.0, CA, &psi[i]));
    }
  }
  return result;
}

//
// It generates the density-density matrix n_i n_j
//
Hamiltonian SGFBase::GenerateDensityDensityMatrix(std::vector<Boson>& psi) {
  Hamiltonian result;
  for (boson_vector_t::size_type i = 0; i < psi.size(); ++i) {
    for (boson_vector_t::size_type j = 0; j < psi.size(); ++j) {
      if (i != j)
        result.push_back(CreateHamiltonianTerm(1.0, CA, &psi[i], CA, &psi[j]));
      else
        result.push_back(CreateHamiltonianTerm(1.0, CACA, &psi[i]));
    }
  }
  return result;
}


//
// It generates the extra terms for the canonical ensemble
//
Hamiltonian SGFBase::GenerateExtraTerms(std::vector<Boson>& psi) {
  Hamiltonian Extra;
  for (std::vector<Boson>::size_type i = 0; i < psi.size(); ++i) {
    Extra.push_back(SGFBase::CreateHamiltonianTerm(1.0, C, &psi[i]));
    Extra.push_back(SGFBase::CreateHamiltonianTerm(1.0, A, &psi[i]));
  }
  return Extra;
}

//
// Input: a vector of boson pointers
// Output: the number of broken lines in this vector
//
unsigned int SGFBase::CountBrokenLines2(const std::vector<Boson>& o) {
  int result = 0;
  for (boson_vector_t::size_type i = 0; i < o.size(); ++i)
    result += Abs(o[i].delta());
  return result;

}




}

#endif