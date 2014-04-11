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

template<class Vector>
class VectorPointerHash {
public:
  typedef typename Vector::size_type Index;
  typedef typename Vector::pointer Pointer;
  typedef typename Vector::const_pointer ConstPointer;
  typedef typename Vector::const_reference Reference;
private:
  Pointer first;
  Index n;
public:
  VectorPointerHash(Vector& v) : first(&v[0]), n(v.size()) {}
  VectorPointerHash(const VectorPointerHash& o) : first(o.first), n(o.n) {}

  inline Pointer pointer(const Index& i) const {
    return first + i;
  }
  inline Index index(const Pointer& p) const {
    return p - first;
  }
  inline Index index(const ConstPointer& p) const {
    return p - first;
  }
  inline Reference reference(const Index& i) const {
    return &*(first + i);
  }
  inline Index size() const {
    return n;
  }
  inline bool contains(const Pointer& p) const {
    return p >= first && p < first + n;
  }
  inline bool contains(const ConstPointer& p) const {
    return p >= first && p < first + n;
  }

};


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
    UpdatableObject(CanonicalProbabilities& o) {}
    virtual inline void update(const HamiltonianTerm*, int, int) = 0;
  };


protected:
  typedef TSumBase<MatrixElement, _float_accumulator> TSum;
  typedef ForestType<TSum> ForestContainer;
  typedef ListSum<MatrixElement, _float_accumulator> PotentialEnergyContainer;
  typedef IntegerSequence OffsetContainer;
  std::vector<UpdatableObject*> UpdatableObjects;

  std::vector<Boson>& _indices;                // A local list of all the bosons

  const VectorPointerHash<std::vector<Boson> > &boson;
  const VectorPointerHash<Hamiltonian> kinetic;
  const VectorPointerHash<Hamiltonian> potential;

  const AdjacencyList kin_pot_adjacency;       // For each term it holds a list of other kinetic terms with one common boson
  const AdjacencyList kin_kin_adjacency;       // For each term it holds a list of other potential terms with one common boson
  const OffsetContainer offsets;               // Contains the offsets
  ForestContainer forest;                      // Contains two trees one for each direction
  PotentialEnergyContainer PotentialEnergy[2]; // Contains the potential energy matrix elements and their sum

  int _NBWL;                                   // The number of all broken world lines
  GreenOperator<long double>& GF;              // Defines the Green operator function

  template<int rl>
  void init() {
    PotentialEnergy[rl].resize(potential.size());
    for (Hamiltonian::size_type i = 0; i < potential.size(); ++i)
      PotentialEnergy[rl].update(i, potential.pointer(i)->me<rl>());
    for (Hamiltonian::size_type i = 0; i < kinetic.size(); ++i)
      forest.set<rl>(offsets.index(kinetic.pointer(i)->offset()), i, kinetic.pointer(i)->me<rl>());
  }

  template<int rl>
  inline double weight(OffsetContainer::size_type i) const {
    return G(offsets[i]) * forest.template norm<rl>(i);
  }

  template<int rl>
  inline void update_trees(const term_vector_t& kin) {
    for (term_vector_t::const_iterator nbr = kin.begin(); nbr != kin.end(); ++nbr)
      forest.template update<rl>(offsets.index((*nbr)->offset()), kinetic.index(*nbr), (*nbr)->me<rl>());
  }

  template<int rl>
  inline void update_energies(const term_vector_t& pot) {
    for (term_vector_t::const_iterator nbr = pot.begin(); nbr != pot.end(); ++nbr)
      PotentialEnergy[rl].update(potential.index(*nbr), (*nbr)->me<rl>());
  }

public:

  CanonicalProbabilities(SGFBase& base) :

    UpdatableObjects(),
    _indices(base.Psi),
    boson(base.Psi),
    kinetic(base.T),
    potential(base.V),
    kin_pot_adjacency(SGFBase::GetAdjacencyList(base.T, base.V)),
    kin_kin_adjacency(SGFBase::GetAdjacencyList(base.T, base.T)),
    offsets(SGFBase::GetOffsets(base.T)),
    forest(base.T.size(), offsets.size()),
    _NBWL(SGFBase::CountBrokenLines2(base.Psi)),
    GF(base.g)

  {

    init<0>();
    init<1>();

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
  inline double weight() const {
    _float_accumulator s = 0.0;
    for ( OffsetContainer::size_type i = 0; i < forest.size(); ++i)
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
    OffsetContainer::size_type i = 0;
    while ((R -= weight<rl>(i)) >= 0)
      ++i;
    return kinetic.pointer(forest.template choose<rl>(i));


  }

  template<int rl>
  inline void update_structures(const HamiltonianTerm* term) {
    unsigned long iterm = kinetic.index(term);
    update_trees<rl>(kin_kin_adjacency[iterm]);
    update_energies<rl>(kin_pot_adjacency[iterm]);
  }

  void insert_update(UpdatableObject* p) {
    UpdatableObjects.push_back(p);
  }

  template<int rl, int arflag>
  inline void update_others(const HamiltonianTerm* term) const {
    for (std::vector<UpdatableObject*>::size_type i = 0; i < UpdatableObjects.size(); ++i)
      UpdatableObjects[i]->update(term, rl, arflag);
  }

  template<int rl, int arflag>
  inline void update_lines(const HamiltonianTerm* term) {
    /* Updating the number of broken lines.
      if the update is before the configuration change then
      _NBWL+=term->offset(arflag);
      otherwise it is
      _NBWL-=term->offset(!arflag); */

    _NBWL -= term->offset < !arflag > ();

  }

  template<int rl, int arflag>
  inline void update_config(const HamiltonianTerm* term) {

    term->update_psi<rl, arflag>();

    update_others<rl, arflag>(term);

    update_lines<rl, arflag>(term);
  }

  template<int rl, int arflag>
  inline void update(const HamiltonianTerm* term) {

    update_config<rl, arflag>(term);

    update_structures<rl>(term);

  }

  friend class BrokenLines;

};



class GrandProbabilities : public CanonicalProbabilities {

  Hamiltonian Extra;                            // The list of extra terms (C dagger, C) is stored locally
  VectorPointerHash<Hamiltonian> extra;         // a term hash for the extra terms
  const AdjacencyList extra_kin_adjacency;      // Adjacency list between the extra terms and the kinetic terms
  const AdjacencyList extra_pot_adjacency;      // Adjacency list between the extra terms and the potential terms
  const AdjacencyList kin_extra_adjacency;      // Adjacency list between the kinetic terms and the extra terms
  const AdjacencyList extra_extra_adjacency;    // Adjacency list between the kinetic terms and the extra terms
  UnorderedSet available[2];                    // contains a list of the permitted terms (non zero matrix element) in each direction
  bool LockedTerms;                             // false if there is already an extra term.
  const double ExtraTermProbability;            // The probability to chose an extra term

  template<int rl>
  void init() {
    available[rl].initialize(Extra.size());
    for (Hamiltonian::size_type i = 0; i < extra.size(); ++i)
      update_extra<rl>(i);
  }

  template<int rl>
  void update_extra(unsigned int i) {
    if (Extra[i].me<rl>() == 0)
      available[rl].erase(i);
    else
      available[rl].insert(i);
  }

  template<int rl>
  inline void update_extra(const term_vector_t& kin) {
    for (term_vector_t::const_iterator nbr = kin.begin(); nbr != kin.end(); ++nbr)
      update_extra<rl>(extra.index(*nbr));
  }


public:


  GrandProbabilities(SGFBase& base) :
    CanonicalProbabilities(base),
    Extra(SGFBase::GenerateExtraTerms(base.Psi)),
    extra(Extra),
    extra_kin_adjacency(SGFBase::GetAdjacencyList(Extra, base.T)),
    extra_pot_adjacency(SGFBase::GetAdjacencyList(Extra, base.V)),
    kin_extra_adjacency(SGFBase::GetAdjacencyList(base.T, Extra)),
    extra_extra_adjacency(SGFBase::GetAdjacencyList(Extra, Extra)),
    LockedTerms(false),
    ExtraTermProbability(base.ExtraTermProbability)
  {

    init<0>();
    init<1>();

  }

  template<int rl>
  const HamiltonianTerm* choose() {

    double Weight = weight<rl>();

    if ( Weight == 0 || (!LockedTerms && NBrokenLines() == 0 && ExtraTermProbability>0 && RNG::Uniform() < ExtraTermProbability) )
      return extra.pointer(available[rl].element(RNG::UniformInteger(available[rl].size())));
    else
      return CanonicalProbabilities::choose<rl>();

  }

  template<int rl, int arflag>
  inline void update(const HamiltonianTerm* term) {

    update_config<rl, arflag>(term);

    // If there is no reference to the GrandCanonicalContainer, we work in the canonical ensemble
    if ( extra.contains(term) ) {
      unsigned long iterm = extra.index(term);
      update_trees<rl> (extra_kin_adjacency[iterm]);
      update_extra<rl> (extra_extra_adjacency[iterm]);
      update_energies<rl>(extra_pot_adjacency[iterm]);
      LockedTerms = !LockedTerms;
    } else {
      unsigned long iterm = kinetic.index(term);
      update_trees<rl> (kin_kin_adjacency[iterm]);
      update_extra<rl> (kin_extra_adjacency[iterm]);
      update_energies<rl>(kin_pot_adjacency[iterm]);
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
