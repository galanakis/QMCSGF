#ifndef __PROBABILITIES__
#define __PROBABILITIES__

#include "RandomNumberGenerator.hh"
#include "GreenOperator.hh"
#include "TSum.hh"
#include "SGFBase.hh"
#include "UnorderedSet.hh"


namespace SGF {


/*
  class OffsetMap
  Helper class. It defines a map that can convert offsets into consecutive
  integers and vice versa.
*/

class OffsetMap {
public:
  typedef int_vector_t::size_type size_type;
private:
  int_vector_t _offsets;   // Keeps a sorted list of the offsets
  std::vector<size_type> _index;    // Map from an offset to consecutive integers
public:

  OffsetMap(const int_vector_t &o) : _offsets(o), _index() {

    int MinOffset=_offsets[0];
    int MaxOffset=_offsets[_offsets.size()-1];
    size_type nindex=MaxOffset-MinOffset+1;
    _index.reserve(nindex);

    for(size_type i=0; i<nindex; ++i)
      _index.push_back(Infinity);

    int_vector_t::const_iterator sit;
    int count=0;
    for(sit=_offsets.begin(); sit!=_offsets.end(); ++sit)
      _index[(*sit)-MinOffset]=count++;

  }

  // Number of different offsets
  inline size_type size() const {
    return _offsets.size();
  }
  // OffsetMap(int offset) will return a consecutive integer
  inline const size_type &operator()(int offset) const {
    return _index[offset-_offsets[0]];
  }
  // OffsetMap[int i] will return the i^th smallest offset
  inline const int &operator[](int i) const {
    return _offsets[i];
  }

};


    // This structure stores one tree in every direction
  struct TreeType {
    TSum tsum[2];
    inline void resize(Hamiltonian::size_type nterms) {
      tsum[0].resize(nterms);
      tsum[1].resize(nterms);
    }
    inline void reset() {
      tsum[0].reset();
      tsum[1].reset();
    }


  };

//
// This Functor overloads the operator()() and it returns
// to return a pointer to the tree in which a particular
// term currently belongs (say based on it's offset).
//

class ForestType {
  const OffsetMap offsets;            // This will map the offsets to consecutive integers
  TreeType *trees;                    // Holds the probability trees
public:
  TreeType *operator()(const HamiltonianTerm *term) const {
    return &trees[offsets(term->offset())];
  }
  TreeType *operator()(unsigned int i) const {
    return &trees[i];
  }
  inline OffsetMap::size_type size() const {
    return offsets.size();
  }
  inline const int & offset(int i) const {
    return offsets[i];
  }

  ForestType(unsigned int nterms,int_vector_t _offsets) : offsets(_offsets) {
    trees=new TreeType[size()];
    for(OffsetMap::size_type i=0; i<size(); ++i)
      trees[i].resize(nterms);
    for(OffsetMap::size_type i=0; i<size(); ++i)
      trees[i].reset();

  }
  ~ForestType() {
    delete [] trees;
  }
};


class CanonicalProbabilities {
public:
  class UpdatableObject {
  public:
    UpdatableObject(CanonicalProbabilities &o) {
      o.insert_update(this);
    }
    virtual inline void update(const HamiltonianTerm*,int,int) = 0;
  };

  //
  // Input: a vector of boson pointers
  // Output: the number of broken lines in this vector
  //
  static unsigned int CountBrokenLines(const std::vector<Boson> &o) {
    int result=0;
    for(boson_vector_t::size_type i=0; i<o.size(); ++i)
      result+=Abs(o[i].delta());
    return result;

  }

protected:
  std::vector<UpdatableObject*> UpdatableObjects;

  std::vector<Boson> &_indices;           // A list of all the bosons
  int _NBWL;                              // The number of all broken world lines
  GreenOperator<long double> &GF;         // Defines the Green operator function

  const Hamiltonian &Potential;           // Local reference of the potential operators
  const AdjacencyList pot_adjacency;      // For each term it holds a list of other kinetic terms with one common site

  const unsigned int RebuildPeriod;       // Number of updates before a rebuild.
  _integer_counter NUpdates;              // Number of updates since last rebuild. This is only used to fix accumulating floating point errors for the energies
  _float_accumulator Energies[2];         // energy of right and left state. It is an accumulator
  std::vector<MatrixElement> EnergyME[2]; // All the potential energies

  const Hamiltonian &Kinetic;             // Local copy of the kinetic operators
  const AdjacencyList kin_adjacency;      // For each term it holds a list of other kinetic terms with one common site
  std::vector<TreeType*> tree_cache;      // remembers the offset of each operator
  ForestType forest;                      // Contains two trees one for each direction


  void print_probabilities() {

    std::cout<<"NumIndices= "<<_indices.size()<<std::endl;
    for(unsigned long i=0; i<_indices.size(); ++i)
      std::cout<<"  "<<_indices[i].nL()<<"\t"<<_indices[i].nR()<<"\t"<<std::endl;

    std::cout<<"EnergyL= "<<Energies[LEFT]<<"\tEnergyR= "<<Energies[RIGHT]<<std::endl;
    std::cout<<"All Energies "<<EnergyME[LEFT].size()<<"\t"<<EnergyME[RIGHT].size()<<std::endl;
    for(unsigned long i=0; i<EnergyME[LEFT].size(); ++i) {
      std::cout<<"   "<<EnergyME[LEFT][i]<<"\t"<<EnergyME[RIGHT][i]<<std::endl;
    }


  }

  inline boson_vector_t::size_type BosonHash(const Boson *p) const {
    return p-&_indices[0];
  }

  inline Hamiltonian::size_type KinHash(const HamiltonianTerm* term) const {
    return term-&Kinetic[0];
  }
  inline const HamiltonianTerm* KinTerm(Hamiltonian::size_type iterm) const {
    return &Kinetic[iterm];
  }

  inline const HamiltonianTerm* PotTerm(Hamiltonian::size_type iterm) const {
    return &Potential[iterm];
  }
  inline Hamiltonian::size_type PotHash(const HamiltonianTerm* term) const {
    return term-&Potential[0];
  }

  // It changes _tsum[2][] and tree_cache.
  // Implementing this as a template results in a small but noticeable performance improvement
  template<int rl>
  inline void update_trees(const term_vector_t &kin) {
    for(term_vector_t::const_iterator nbr=kin.begin(); nbr!=kin.end(); ++nbr) {
      Hamiltonian::size_type fndex= KinHash(*nbr);
      TreeType *i_tree = tree_cache[fndex];
      TreeType *f_tree = forest(*nbr);
      MatrixElement fme = (*nbr)->me<rl>();
      f_tree->tsum[ rl].update(fndex,fme);
      if(i_tree!=f_tree) {
        MatrixElement jme = i_tree->tsum[!rl].element(fndex);
        i_tree->tsum[ rl].update(fndex,0);
        i_tree->tsum[!rl].update(fndex,0);
        f_tree->tsum[!rl].update(fndex,jme);
        tree_cache[fndex] = f_tree;
      }
    }
  }

  // It updates the Energies[2], EnergyME[2][].
  template<int rl>
  inline void update_energies(const term_vector_t &pot) {
    for(term_vector_t::const_iterator nbr=pot.begin(); nbr!=pot.end(); ++nbr) {
      MatrixElement me=(*nbr)->me<rl>();
      std::vector<MatrixElement>::size_type i=PotHash(*nbr);
      Energies[rl]+=me-EnergyME[rl][i];
      EnergyME[rl][i]=me;
    }
    if((++NUpdates % RebuildPeriod)==0)
      rebuild_energies();

  }

  inline void rebuild_energies() {
    NUpdates=0;
    Energies[LEFT]=0;
    Energies[RIGHT]=0;
    for(Hamiltonian::size_type i=0; i<Potential.size(); ++i) {
      Energies[LEFT ]+=EnergyME[LEFT ][i]=PotTerm(i)->me<LEFT>();
      Energies[RIGHT]+=EnergyME[RIGHT][i]=PotTerm(i)->me<RIGHT>();
    }

  }

  void insert_update(UpdatableObject *p) {
    UpdatableObjects.push_back(p);
  }


public:

  CanonicalProbabilities(SGFBase &base) :

    UpdatableObjects(),
    _indices(base.Psi),
    _NBWL(CountBrokenLines(_indices)),
    GF(base.g),
    Potential(base.V),
    pot_adjacency(SGFBase::GetAdjacencyList(base.T,base.V)),
    RebuildPeriod(1000000),
    NUpdates(0),
    Kinetic(base.T),
    kin_adjacency(SGFBase::GetAdjacencyList(base.T,base.T)),
    forest(base.T.size(),SGFBase::GetOffsets(base.T))

  {

    EnergyME[0].resize(base.V.size());
    EnergyME[1].resize(base.V.size());
    rebuild_energies();

    tree_cache.resize(base.T.size());

    for(Hamiltonian::size_type i=0; i<base.T.size(); ++i) {
      TreeType *tree = forest(&base.T[i]);
      tree_cache[i] = tree;
      tree->tsum[RIGHT].update(i,base.T[i].me<RIGHT>());
      tree->tsum[LEFT].update(i,base.T[i].me<LEFT>());
    }

  }

  ~CanonicalProbabilities() {}

  template<int rl>
  inline const _float_accumulator &Energy() const {
    return Energies[rl];
  }

  inline MatrixElement G(int offset=0) const {
    return GF(NBrokenLines()+offset);   // The value of the Green Operator given the total broken lines and the offset.
  }

  template<int rl>
  inline double weight(OffsetMap::size_type i) const {
    return G(forest.offset(i))*forest(i)->tsum[rl].norm();
  }

  template<int rl>
  inline double weight() const {
    _float_accumulator s=0.0;
    for( OffsetMap::size_type i=0; i<forest.size(); ++i)
      s+=weight<rl>(i);
    return s;
  }

  inline const int &NBrokenLines() const {
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

    double R=RNG::Uniform()*weight<rl>();
    OffsetMap::size_type i=0;
    while((R-=weight<rl>(i))>=0)
      ++i;
    return KinTerm(forest(i)->tsum[rl].choose());


  }


  template<int rl,int arflag>
  inline void update_config(const HamiltonianTerm* term) {

    term->update_psi<rl,arflag>();

    for(std::vector<UpdatableObject*>::size_type i=0; i<UpdatableObjects.size(); ++i)
      UpdatableObjects[i]->update(term,rl,arflag);

    /* Updating the number of broken lines.
      if the update is before the configuration change then
      _NBWL+=term->offset(arflag);
      otherwise it is
      _NBWL-=term->offset(!arflag); */

    _NBWL-=term->offset<!arflag>();

  }

  template<int rl,int arflag>
  inline void update(const HamiltonianTerm* term) {

    update_config<rl,arflag>(term);

    unsigned long iterm=KinHash(term);
    update_trees<rl>(kin_adjacency[iterm]);
    update_energies<rl>(pot_adjacency[iterm]);

  }

  friend class BrokenLines;

};


const double ExtraTermProbability=0.05;

class GrandProbabilities : public CanonicalProbabilities {

  Hamiltonian Extra;                            // The list of extra terms
  AdjacencyList extra_kin_adjacency;            // Adjacency list between the extra terms and the kinetic terms
  AdjacencyList extra_pot_adjacency;            // Adjacency list between the extra terms and the potential terms
  UnorderedSet available[2];                    // contains a list of the permitted terms (non zero matrix element) in each direction
  bool LockedTerms;                             // false if there is already an extra term.
  const HamiltonianTerm* WormInit;              // It contains a pointer to the extra term in the opertor string or zero if there is none

public:


  GrandProbabilities(SGFBase &base) : CanonicalProbabilities(base), WormInit(0), LockedTerms(false) {
    //
    // The extra terms appear in pairs with the same index.
    //
    for(boson_vector_t::size_type i=0; i<_indices.size(); ++i) {
      Extra.push_back(SGFBase::CreateHamiltonianTerm(1.0,C,&_indices[i]));
      Extra.push_back(SGFBase::CreateHamiltonianTerm(1.0,A,&_indices[i]));
    }

    extra_kin_adjacency=SGFBase::GetAdjacencyList(Extra,base.T);
    extra_pot_adjacency=SGFBase::GetAdjacencyList(Extra,base.V);

    available[LEFT].initialize(Extra.size());
    available[RIGHT].initialize(Extra.size());
    for(Hamiltonian::size_type i=0; i<Extra.size(); ++i) {
      update_extra<LEFT>(i);
      update_extra<RIGHT>(i);
    }

  }

  template<int rl>
  void update_extra(unsigned int i) {
    if(Extra[i].me<rl>()==0)
      available[rl].erase(i);
    else
      available[rl].insert(i);
  }

  template<int rl>
  const HamiltonianTerm* choose() {

    double Weight=weight<rl>();

    if( Weight==0 || (!LockedTerms && NBrokenLines()==0 && RNG::Uniform()<ExtraTermProbability) )
      return WormInit=&Extra[available[rl].element(RNG::UniformInteger(available[rl].size()))];
    else {
      double R=RNG::Uniform()*Weight;
      OffsetMap::size_type i=0;
      while((R-=weight<rl>(i))>=0)
        ++i;
      return KinTerm(forest(i)->tsum[rl].choose());
    }


  }

  template<int rl,int arflag>
  inline void update(const HamiltonianTerm* term) {

    update_config<rl, arflag>(term);

    // If there is no reference to the GrandCanonicalContainer, we work in the canonical ensemble
    if( WormInit!=term ) {
      unsigned long iterm=KinHash(term);
      update_trees<rl>(kin_adjacency[iterm]);
      update_energies<rl>(pot_adjacency[iterm]);
    } else {
      update_trees<rl>(extra_kin_adjacency[term-&Extra[0]]);
      update_energies<rl>(extra_pot_adjacency[term-&Extra[0]]);
      LockedTerms=!LockedTerms;
    }

    for(HamiltonianTerm::const_iterator it=term->begin(); it!=term->end(); ++it) {
      unsigned int ind=BosonHash(it->particle_id());
      for(unsigned int u=0; u<2; ++u)
        update_extra<rl>(2*ind+u);
    }

  }


};



/*
  class BrokenLines

  It holds a list of the indices with broken lines. It is only used for measurements
  and this is why it is UpdatableObject.

*/


typedef std::vector<std::pair<Boson*,int> > BosonDeltaMapType;

// This one is ok to be slow since it is called only in the initializer
inline BosonDeltaMapType GetTermHash(const std::vector<IndexedProductElement> &p) {
  std::map<Boson*,int> indices;
  for(std::vector<IndexedProductElement>::const_iterator it=p.begin(); it!=p.end(); ++it) {
    int delta=-it->delta();
    if(delta!=0) indices[it->particle_id()]=delta;
  }

  BosonDeltaMapType vmap;
  vmap.reserve(indices.size());
  std::map<Boson*,int>::const_iterator it;
  for(it=indices.begin(); it!=indices.end(); ++it)
    vmap.push_back(*it);

  return vmap;
}


class BrokenLines : public CanonicalProbabilities::UpdatableObject {
  UnorderedSet _broken_lines;
  Boson* Psi0;
public:


  BrokenLines(CanonicalProbabilities &c) :  CanonicalProbabilities::UpdatableObject(c),_broken_lines(c._indices.size()) {

    Psi0=&c._indices[0];
    for(std::vector<Boson>::size_type i=0; i<c._indices.size(); ++i) {
      if(c._indices[i].delta()!=0)
        _broken_lines.insert(i);
    }
  }

  inline void update(const HamiltonianTerm* term,int,int) {
    for(HamiltonianTerm::const_iterator it=term->begin(); it!=term->end(); ++it) {
      Boson *pind=it->particle_id();
      if(pind->delta()!=0)
        _broken_lines.insert(pind-Psi0);
      else
        _broken_lines.erase(pind-Psi0);
    }
  }

  // This one must be very fast
  inline BosonDeltaMapType map() {
    BosonDeltaMapType vmap;
    _broken_lines.sort();
    vmap.reserve(_broken_lines.size());
    for(UnorderedSet::iterator it=_broken_lines.begin(); it!=_broken_lines.end(); ++it) {
      int delta=Psi0[*it].delta();
      if(delta!=0) vmap.push_back(std::pair<Boson*,int>(Psi0+*it,delta));
    }
    return vmap;
  }

  inline const BosonDeltaMapType operator()() {
    return map();
  }


};


}

#endif
