#ifndef __PROBABILITIES__
#define __PROBABILITIES__

#include "RandomNumberGenerator.hh"
#include "GreenOperator.hh"
#include "SGFBase.hh"
#include "UnorderedSet.hh"
#include "Forest.hh"
#include "IntegerSequence.hh"
#include "ListSum.hh"

namespace SGF {



// It represents a type for the Hamiltonian terms
//typedef HamitonianTerm TermType;
// It represents a collection of terms
//typedef Hamiltonian TermList;
// It represents the type of mantrix element container
//typedef TSum TermValueHash;


class CanonicalProbabilities {
public:
  class UpdatableObject {
  public:
    UpdatableObject(CanonicalProbabilities& o) {
      o.insert_update(this);
    }
    virtual inline void update(const HamiltonianTerm*, int, int) = 0;
  };

  //
  // Input: a vector of boson pointers
  // Output: the number of broken lines in this vector
  //
  static unsigned int CountBrokenLines(const std::vector<Boson>& o) {
    int result = 0;
    for (boson_vector_t::size_type i = 0; i < o.size(); ++i)
      result += Abs(o[i].delta());
    return result;

  }

protected:
  typedef TSumBase<MatrixElement, _float_accumulator> TSum;
  typedef ForestType<TSum> ForestContainer;
  typedef ListSum<MatrixElement, _float_accumulator> PotentialEnergyContainer;
  std::vector<UpdatableObject*> UpdatableObjects;

  std::vector<Boson>& _indices;           // A list of all the bosons
  int _NBWL;                              // The number of all broken world lines
  GreenOperator<long double>& GF;         // Defines the Green operator function

  const Hamiltonian& Potential;           // Local reference of the potential operators
  const AdjacencyList pot_adjacency;      // For each term it holds a list of other kinetic terms with one common site

  PotentialEnergyContainer PotentialEnergy[2];

  const Hamiltonian& Kinetic;             // Local copy of the kinetic operators
  const AdjacencyList kin_adjacency;      // For each term it holds a list of other kinetic terms with one common site

  const IntegerSequence offsets;          // Contains the offsets
  ForestContainer forest;                 // Contains two trees one for each direction

  inline boson_vector_t::size_type BosonHash(const Boson* p) const {
    return p - &_indices[0];
  }

  inline Hamiltonian::size_type KinHash(const HamiltonianTerm* term) const {
    return term - &Kinetic[0];
  }
  inline const HamiltonianTerm* KinTerm(Hamiltonian::size_type iterm) const {
    return &Kinetic[iterm];
  }

  inline const HamiltonianTerm* PotTerm(Hamiltonian::size_type iterm) const {
    return &Potential[iterm];
  }
  inline Hamiltonian::size_type PotHash(const HamiltonianTerm* term) const {
    return term - &Potential[0];
  }

  inline Hamiltonian::size_type PotSize() const {
    return Potential.size();
  }

  inline Hamiltonian::size_type KinSize() const {
    return Kinetic.size();
  }


  template<int rl>
  inline void update_trees(const term_vector_t& kin) {
    for (term_vector_t::const_iterator nbr = kin.begin(); nbr != kin.end(); ++nbr)
      forest.template update<rl>(offsets.index((*nbr)->offset()), KinHash(*nbr), (*nbr)->me<rl>());
  }

  // It updates the Energies[2], EnergyME[2][].
  template<int rl>
  inline void update_energies(const term_vector_t& pot) {
    for (term_vector_t::const_iterator nbr = pot.begin(); nbr != pot.end(); ++nbr)
      PotentialEnergy[rl].update(PotHash(*nbr), (*nbr)->me<rl>());
  }

  void insert_update(UpdatableObject* p) {
    UpdatableObjects.push_back(p);
  }


public:

  CanonicalProbabilities(SGFBase& base) :

    UpdatableObjects(),
    _indices(base.Psi),
    _NBWL(CountBrokenLines(_indices)),
    GF(base.g),
    Potential(base.V),
    pot_adjacency(SGFBase::GetAdjacencyList(base.T, base.V)),
    Kinetic(base.T),
    kin_adjacency(SGFBase::GetAdjacencyList(base.T, base.T)),
    offsets(SGFBase::GetOffsets(base.T)),
    forest(base.T.size(), offsets.size())

  {

    PotentialEnergy[0].resize(PotSize());
    PotentialEnergy[1].resize(PotSize());

    for (Hamiltonian::size_type i = 0; i < PotSize(); ++i) {
      PotentialEnergy[0].update(i, PotTerm(i)->me<0>());
      PotentialEnergy[1].update(i, PotTerm(i)->me<1>());
    }

    for (Hamiltonian::size_type i = 0; i < KinSize(); ++i) {
      forest.set<0>(offsets.index(KinTerm(i)->offset()), i, KinTerm(i)->me<0>());
      forest.set<1>(offsets.index(KinTerm(i)->offset()), i, KinTerm(i)->me<1>());
    }

  }

  ~CanonicalProbabilities() {}

  template<int rl>
  inline const _float_accumulator& Energy() const {
    return PotentialEnergy[rl].value();
  }

  inline MatrixElement G(int offset = 0) const {
    return GF(NBrokenLines() + offset); // The value of the Green Operator given the total broken lines and the offset.
  }

  template<int rl>
  inline double weight(IntegerSequence::size_type i) const {
    return G(offsets[i]) * forest.template norm<rl>(i);
  }

  template<int rl>
  inline double weight() const {
    _float_accumulator s = 0.0;
    for ( IntegerSequence::size_type i = 0; i < forest.size(); ++i)
      s += weight<rl>(i);
    return s;
  }

  inline const int& NBrokenLines() const {
    return _NBWL;
  }

  /*

     If the ensemble is GrandCanonical (ensemble!=0)
     AND if there is no extra term (ensemble->WormInit!=0)
     AND if there are no broken lines (NBrokenLines()==0)
     then pick a regular term with a fixed probability.

  */

  template<int rl>
  const HamiltonianTerm* choose() const {

    double R = RNG::Uniform() * weight<rl>();
    IntegerSequence::size_type i = 0;
    while ((R -= weight<rl>(i)) >= 0)
      ++i;
    return KinTerm(forest.template choose<rl>(i));


  }


  template<int rl, int arflag>
  inline void update_config(const HamiltonianTerm* term) {

    term->update_psi<rl, arflag>();

    for (std::vector<UpdatableObject*>::size_type i = 0; i < UpdatableObjects.size(); ++i)
      UpdatableObjects[i]->update(term, rl, arflag);

    /* Updating the number of broken lines.
      if the update is before the configuration change then
      _NBWL+=term->offset(arflag);
      otherwise it is
      _NBWL-=term->offset(!arflag); */

    _NBWL -= term->offset < !arflag > ();

  }

  template<int rl, int arflag>
  inline void update(const HamiltonianTerm* term) {

    update_config<rl, arflag>(term);

    unsigned long iterm = KinHash(term);
    update_trees<rl>(kin_adjacency[iterm]);
    update_energies<rl>(pot_adjacency[iterm]);

  }

  friend class BrokenLines;

};



class GrandProbabilities : public CanonicalProbabilities {

  const double ExtraTermProbability;

  Hamiltonian Extra;                            // The list of extra terms
  AdjacencyList extra_kin_adjacency;            // Adjacency list between the extra terms and the kinetic terms
  AdjacencyList extra_pot_adjacency;            // Adjacency list between the extra terms and the potential terms
  UnorderedSet available[2];                    // contains a list of the permitted terms (non zero matrix element) in each direction
  bool LockedTerms;                             // false if there is already an extra term.
  const HamiltonianTerm* WormInit;              // It contains a pointer to the extra term in the opertor string or zero if there is none

public:


  GrandProbabilities(SGFBase& base) : CanonicalProbabilities(base), WormInit(0), LockedTerms(false), ExtraTermProbability(0.05) {
    //
    // The extra terms appear in pairs with the same index.
    //
    for (boson_vector_t::size_type i = 0; i < _indices.size(); ++i) {
      Extra.push_back(SGFBase::CreateHamiltonianTerm(1.0, C, &_indices[i]));
      Extra.push_back(SGFBase::CreateHamiltonianTerm(1.0, A, &_indices[i]));
    }

    extra_kin_adjacency = SGFBase::GetAdjacencyList(Extra, base.T);
    extra_pot_adjacency = SGFBase::GetAdjacencyList(Extra, base.V);

    available[LEFT].initialize(Extra.size());
    available[RIGHT].initialize(Extra.size());
    for (Hamiltonian::size_type i = 0; i < Extra.size(); ++i) {
      update_extra<LEFT>(i);
      update_extra<RIGHT>(i);
    }

  }

  template<int rl>
  void update_extra(unsigned int i) {
    if (Extra[i].me<rl>() == 0)
      available[rl].erase(i);
    else
      available[rl].insert(i);
  }

  template<int rl>
  const HamiltonianTerm* choose() {

    double Weight = weight<rl>();

    if ( Weight == 0 || (!LockedTerms && NBrokenLines() == 0 && RNG::Uniform() < ExtraTermProbability) )
      return WormInit = &Extra[available[rl].element(RNG::UniformInteger(available[rl].size()))];
    else {
      double R = RNG::Uniform() * Weight;
      IntegerSequence::size_type i = 0;
      while ((R -= weight<rl>(i)) >= 0)
        ++i;
      return KinTerm(forest.template choose<rl>(i));
    }


  }

  template<int rl, int arflag>
  inline void update(const HamiltonianTerm* term) {

    update_config<rl, arflag>(term);

    // If there is no reference to the GrandCanonicalContainer, we work in the canonical ensemble
    if ( WormInit != term ) {
      unsigned long iterm = KinHash(term);
      update_trees<rl>(kin_adjacency[iterm]);
      update_energies<rl>(pot_adjacency[iterm]);
    } else {
      unsigned long iterm = term - &Extra[0];
      update_trees<rl>(extra_kin_adjacency[iterm]);
      update_energies<rl>(extra_pot_adjacency[iterm]);
      LockedTerms = !LockedTerms;
    }

    for (HamiltonianTerm::const_iterator it = term->begin(); it != term->end(); ++it) {
      unsigned int ind = BosonHash(it->particle_id());
      for (unsigned int u = 0; u < 2; ++u)
        update_extra<rl>(2 * ind + u);
    }

  }


};



/*
  class BrokenLines

  It holds a list of the indices with broken lines. It is only used for measurements
  and this is why it is UpdatableObject.

*/


typedef std::vector<std::pair<Boson*, int> > BosonDeltaMapType;

// This one is ok to be slow since it is called only in the initializer
inline BosonDeltaMapType GetTermHash(const std::vector<IndexedProductElement>& p) {
  std::map<Boson*, int> indices;
  for (std::vector<IndexedProductElement>::const_iterator it = p.begin(); it != p.end(); ++it) {
    int delta = -it->delta();
    if (delta != 0) indices[it->particle_id()] = delta;
  }

  BosonDeltaMapType vmap;
  vmap.reserve(indices.size());
  std::map<Boson*, int>::const_iterator it;
  for (it = indices.begin(); it != indices.end(); ++it)
    vmap.push_back(*it);

  return vmap;
}


class BrokenLines : public CanonicalProbabilities::UpdatableObject {
  UnorderedSet _broken_lines;
  Boson* Psi0;
public:


  BrokenLines(CanonicalProbabilities& c) :  CanonicalProbabilities::UpdatableObject(c), _broken_lines(c._indices.size()) {

    Psi0 = &c._indices[0];
    for (std::vector<Boson>::size_type i = 0; i < c._indices.size(); ++i) {
      if (c._indices[i].delta() != 0)
        _broken_lines.insert(i);
    }
  }

  inline void update(const HamiltonianTerm* term, int, int) {
    for (HamiltonianTerm::const_iterator it = term->begin(); it != term->end(); ++it) {
      Boson* pind = it->particle_id();
      if (pind->delta() != 0)
        _broken_lines.insert(pind - Psi0);
      else
        _broken_lines.erase(pind - Psi0);
    }
  }

  // This one must be very fast
  inline const BosonDeltaMapType map() {
    BosonDeltaMapType vmap;
    _broken_lines.sort();
    vmap.reserve(_broken_lines.size());
    for (UnorderedSet::iterator it = _broken_lines.begin(); it != _broken_lines.end(); ++it) {
      int delta = Psi0[*it].delta();
      if (delta != 0) vmap.push_back(std::pair<Boson*, int>(Psi0 + *it, delta));
    }
    return vmap;
  }

  inline const BosonDeltaMapType operator()() {
    return map();
  }


};


}

#endif
