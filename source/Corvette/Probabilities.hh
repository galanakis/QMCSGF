#ifndef __PROBABILITIES__
#define __PROBABILITIES__

#include "RandomNumberGenerator.hh"
#include "HamiltonianTerm.hh"
#include "GreenOperator.hh"
#include "TSum.hh"
#include <vector>
#include <list>
#include <map>
#include <set>
#include <queue>

namespace SGF {


/*

  This class contains some utilities that I don't know where to put.
  They are all static functions with a specific input
  and a specific output

*/

typedef std::vector<int> int_vector_t;
typedef std::vector<Boson*> boson_vector_t;
typedef std::set<Boson*> boson_set_t;
typedef std::vector<const HamiltonianTerm *> term_vector_t;
typedef term_vector_t::const_iterator adjacency_iterator;
typedef std::pair<adjacency_iterator,adjacency_iterator> range_type;
typedef std::vector<term_vector_t > AdjacencyList;


class Orphans {
public:

   //
   // Input: a list of kinetic terms
   // Output: a softed list containing all the possible offsets
   //
   static int_vector_t GetOffsets(const Hamiltonian & T) {
      int_vector_t offsets;
      // Scan all kinetic terms to find all the lengths
      std::set<int> set_offsets;
      for (Hamiltonian::size_type i = 0; i < T.size(); ++i)
         for (int offset = T[i].minoffset(); offset <= T[i].maxoffset(); offset += 2)
            set_offsets.insert(offset);

      offsets.clear();
      offsets.insert(offsets.begin(), set_offsets.begin(), set_offsets.end());

      return offsets;
   }

   //
   // Input: a list of kinetic terms
   // Output: a softed list containing all the possible bosons
   //
   static boson_vector_t GetConfiguration(const Hamiltonian & T) {
      boson_vector_t indices;
      boson_set_t indexset;
      for (Hamiltonian::size_type i = 0; i < T.size(); ++i)
         for (HamiltonianTerm::size_type j = 0; j < T[i].product().size(); ++j)
            indexset.insert(T[i].product()[j].particle_id());

      indices.insert(indices.begin(), indexset.begin(), indexset.end());

      return indices;

   }

   //
   // Input: a softed list containing all the possible bosons
   // Output: a list of the "extra" kinetic terms needed for the Grand canonical ensemble
   //
   static Hamiltonian GetExtraTerms(const boson_vector_t &psi) {
      Hamiltonian result;
      for(boson_vector_t::size_type i=0;i<psi.size();++i) {
         result.push_back(HamiltonianTerm(1.0,IndexedProductElement(A,psi[i])));
         result.push_back(HamiltonianTerm(1.0,IndexedProductElement(C,psi[i])));
      }

      return result;

   }

   //
   // Input: a list of all the kinetic terms
   // Output: a list of the "extra" kinetic terms needed for the Grand canonical ensemble
   //
   static Hamiltonian GetExtraTerms(const Hamiltonian &T) {

      return GetExtraTerms(GetConfiguration(T));

   }


   //
   // Input: a vector of boson pointers
   // Output: the number of broken lines in this vector
   //
   static unsigned int CountBrokenLines(const boson_vector_t &o) {
      int result=0;
      for(boson_vector_t::size_type i=0; i<o.size(); ++i)
         result+=Abs(o[i]->delta());
      return result;

   }

   //
   // Input: a vector of boson pointers
   // Output: a set of the boson pointers where the lines are broken
   //
   static  boson_set_t GetListBrokenLines(const boson_vector_t &o) {
      boson_set_t broken;
      for(boson_vector_t::const_iterator it=o.begin(); it!=o.end(); ++it) {
         if((*it)->delta() != 0)
            broken.insert(*it);
      }
      return broken;
   }


   //
   // Input: two lists of terms, one called row and another called column
   // Output: a vector of size Trow, which contains the indices of the operators
   //         in Tcol which share a common boson
   //
   static AdjacencyList GetAdjacencyList(const Hamiltonian &Trow,const Hamiltonian &Tcol)  {

      AdjacencyList adjacency; // The adjacency list is stored here

      // First, categorize the Tcol terms by index.
      std::map<Boson*,std::set<Hamiltonian::size_type> > map_to_set;
      for(Hamiltonian::size_type i=0; i<Tcol.size(); ++i)
         for(Hamiltonian::size_type j=0; j<Tcol[i].product().size(); ++j)
            map_to_set[Tcol[i].product()[j].particle_id()].insert(i);

      adjacency.clear();
      adjacency.resize(Trow.size());

      /* For each term in Trow, merge the sets corresponding to it's indices
         Then, copy the set elements to a vector. */
      for(Hamiltonian::size_type i=0; i<Trow.size(); ++i) {
         std::set<Hamiltonian::size_type> merged;
         for(Hamiltonian::size_type j=0; j<Trow[i].product().size(); ++j) {
            if(Trow[i].product()[j].delta()!=0) {
               Boson* pid=Trow[i].product()[j].particle_id();
               std::set<Hamiltonian::size_type> &s=map_to_set[pid];
               merged.insert(s.begin(),s.end());
            }
         }

         adjacency[i].reserve(merged.size());
         for(std::set<Hamiltonian::size_type>::const_iterator sit=merged.begin(); sit!=merged.end(); ++sit)
            adjacency[i].push_back(&Tcol[*sit]);

      }

      return adjacency;

   }

};



/*
  class Configuration {
   boson_vector_t _indices;
   std::vector<UpdatableObject*> UpdatableObjects;
  }

  holds the indices of the Bosons, together the number of broken world lines.
  When updated it also updates all UpdatableObjects.

  sub-class UpdatableObject

  This is class to be used as the pure base class of objects that update themselves
  after each addition/insertion of an operator.
  The constructor of the class will push the new object in a static
  list of pointers. Then the Configuration object which is a
  friend, will read this list and call the update function of
  the objects after each update. This way, all objects that
  are declared UpdatableObject will be automatically updated.
  This relies on the virtual function mechanism which is slightly
  slow, so I only use it for measurable quantities.


*/

class Configuration {
public:
   class UpdatableObject {
   public:
      UpdatableObject(Configuration &o) { o.insert_update(this); }
      virtual inline void update(const HamiltonianTerm*,int,int) = 0;
   };

protected:
   boson_vector_t _indices;
   std::vector<UpdatableObject*> UpdatableObjects;
public:
   Configuration(const Hamiltonian &T) : _indices(Orphans::GetConfiguration(T)), UpdatableObjects() {}
   Configuration(const Configuration &o) : _indices(o._indices), UpdatableObjects() {}
   int CountBrokenLines() const {
      return Orphans::CountBrokenLines(_indices);
   }
   std::set<Boson*> GetListBrokenLines() const {
      return Orphans::GetListBrokenLines(_indices);
   }
   void insert_update(UpdatableObject *p) {
      UpdatableObjects.push_back(p);
   }

   //
   // This function updates the configuration and calls whatever needs
   // updating.
   //
   inline void update(const HamiltonianTerm* term,int rl,int arflag) {

      term->update_psi(rl,arflag);

      for(std::vector<UpdatableObject*>::size_type i=0; i<UpdatableObjects.size(); ++i)
         UpdatableObjects[i]->update(term,rl,arflag);

   }


};




/*
  class BrokenLines

  It holds a list of the indices with broken lines. It is only used for measurements
  and this is why it is UpdatableObject.

*/


class BrokenLines : public Configuration::UpdatableObject {
   std::set<Boson*> _broken_lines;
public:

   typedef std::vector<std::pair<Boson*,int> > BosonDeltaMapType;
   BrokenLines(const std::set<Boson*> &o,Configuration &c) :  Configuration::UpdatableObject(c), _broken_lines(o) {}

   inline void update(const std::vector<IndexedProductElement> &p) {
      for(std::vector<IndexedProductElement>::const_iterator it=p.begin(); it!=p.end(); ++it) {
         Boson *pind=it->particle_id();
         if(pind->delta()!=0)
            _broken_lines.insert(pind);
         else
            _broken_lines.erase(pind);
      }
   }

   inline void update(const HamiltonianTerm* term,int,int) {
      update(term->product());
   }

   // This one must be very fast
   inline BosonDeltaMapType map() const {
      BosonDeltaMapType vmap;
      vmap.reserve(_broken_lines.size());
      std::set<Boson*>::const_iterator it;
      for(it=_broken_lines.begin(); it!=_broken_lines.end(); ++it) {
         int delta=(*it)->delta();
         if(delta!=0) vmap.push_back(std::pair<Boson*,int>(*it,delta));
      }
      return vmap;
   }

   inline const BosonDeltaMapType operator()() const {
      return map();
   }

   // This one is ok to be slow since it is called only in the initializer
   static inline BosonDeltaMapType map(const std::vector<IndexedProductElement> &p) {
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

};

/*
  class TermCount

  measures how many times each operator gets inserted minus the times
  it is removed from the right of the GreenOperator.

*/

class TermCount : public Configuration::UpdatableObject {

   const HamiltonianTerm * Kinetic0;
   typedef std::vector<long> count_type;
   typedef count_type::iterator iterator;
   count_type _count;

public:
   TermCount(const Hamiltonian &T,Configuration &o) : Configuration::UpdatableObject(o), Kinetic0(&T[0]), _count(T.size()) {
      reset();
   }

   inline void update(const HamiltonianTerm* term,int rl,int ar) {
      if(rl==RIGHT) _count[term-Kinetic0]+=Sign[ar];
   }

   inline void reset() {
      for(iterator ptr=_count.begin(); ptr!=_count.end(); ++ptr)
         *ptr=0;
   }

};

/*

  class PotentialEnergies {
   static unsigned int RebuildPeriod;
   const Hamiltonian &Potential; 
   const AdjacencyList pot_adjacency;
   const HamiltonianTerm* Kinetic0;
   _integer_counter NUpdates;
   _float_accumulator Energies[2]; 
   std::vector<MatrixElement> EnergyME[2];
  }

  it holds and updates then right and left potential energy.

*/

class PotentialEnergies {
   const Hamiltonian &Potential;           // Local reference of the potential operators
   const AdjacencyList pot_adjacency;      // For each term it holds a list of other kinetic terms with one common site
   const HamiltonianTerm* Kinetic0;

   const unsigned int RebuildPeriod;       // Number of updates before a rebuild.
   _integer_counter NUpdates;              // Number of updates since last rebuild. This is only used to fix accumulating floating point errors for the energies
   _float_accumulator Energies[2];         // energy of right and left state. It is an accumulator
   std::vector<MatrixElement> EnergyME[2];

   inline const HamiltonianTerm * Term(Hamiltonian::size_type iterm) const {
      return &Potential[iterm];
   }
   inline Hamiltonian::size_type Hash(const HamiltonianTerm* term) const {
      return term-&Potential[0];
   }


   inline void rebuild() {
      NUpdates=0;
      for(int direction=0; direction<2; ++direction) {
         Energies[direction]=0;
         for(HamiltonianTerm::size_type i=0; i<Potential.size(); ++i) {
            MatrixElement me=Term(i)->me(direction);
            EnergyME[direction][i]=me;
            Energies[direction]+=me;
         }
      }
   }


public:

   PotentialEnergies(const Hamiltonian &T,const Hamiltonian &P) : Potential(P), pot_adjacency(Orphans::GetAdjacencyList(T,P)), Kinetic0(&T[0]), NUpdates(0), RebuildPeriod(1000000) {
      EnergyME[0].resize(Potential.size());
      EnergyME[1].resize(Potential.size());
      rebuild();
   }

   PotentialEnergies(const PotentialEnergies &o) : Potential(o.Potential), pot_adjacency(o.pot_adjacency), Kinetic0(o.Kinetic0), NUpdates(o.NUpdates), RebuildPeriod(o.RebuildPeriod) {
      EnergyME[0]=o.EnergyME[0];
      EnergyME[1]=o.EnergyME[1];
      Energies[0]=o.Energies[0];
      Energies[1]=o.Energies[1];
   }

   inline const _float_accumulator &operator()(int rl) const {
      return Energies[rl];
   }

   // It updates the Energies[2], EnergyME[2][].
   inline void update(const HamiltonianTerm* const &term,int rl) {
      const term_vector_t &pot=pot_adjacency[term-Kinetic0];
      for(term_vector_t::const_iterator nbr=pot.begin(); nbr!=pot.end(); ++nbr) {
         MatrixElement me=(*nbr)->me(rl);
         std::vector<MatrixElement>::size_type i=Hash(*nbr);
         Energies[rl]+=me-EnergyME[rl][i];
         EnergyME[rl][i]=me;
      }
      if((++NUpdates % RebuildPeriod)==0)
         rebuild();

   }

};


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

   ForestType(const Hamiltonian &K) : offsets(Orphans::GetOffsets(K)) {
      trees=new TreeType[size()];
      for(OffsetMap::size_type i=0; i<size(); ++i)
         trees[i].resize(K.size());
      for(OffsetMap::size_type i=0; i<size(); ++i)
         trees[i].reset();

   }
   ~ForestType() {
      delete [] trees;
   }
};



/*

  class KineticProbabilities {

   const Hamiltonian &Kinetic;
   const AdjacencyList kin_adjacency;
   std::vector<TreeType*> tree_cache;
   ForestType forest;
  
  }

  holds and updates the Trees for each direction.

*/



class KineticProbabilities {
   const Hamiltonian &Kinetic;         // Local copy of the kinetic operators
   const AdjacencyList kin_adjacency;  // For each term it holds a list of other kinetic terms with one common site
   inline Hamiltonian::size_type Hash(const HamiltonianTerm* term) const {
      return term-&Kinetic[0];
   }
   inline const HamiltonianTerm * Term(Hamiltonian::size_type iterm) const {
      return &Kinetic[iterm];
   }
   std::vector<TreeType*> tree_cache;      // remembers the offset of each operator

   ForestType forest;

public:

   inline OffsetMap::size_type size() const {
      return forest.size();
   }
   inline const int &offset(int i) const {
      return forest.offset(i);
   }

   KineticProbabilities(const Hamiltonian &T) : Kinetic(T), kin_adjacency(Orphans::GetAdjacencyList(Kinetic,Kinetic)), forest(Kinetic) {

      tree_cache.resize(Kinetic.size());

      for(int direction=0; direction<2; ++direction) {
         for(Hamiltonian::size_type i=0; i<Kinetic.size(); ++i) {
            TreeType *tree = forest(&Kinetic[i]);
            tree_cache[i] = tree;
            tree->tsum[direction].update(i,Kinetic[i].me(direction));
         }
      }


   }

   ~KineticProbabilities() {
   }


   inline const _float_accumulator &norm(int rl,int i) const {
      return forest(i)->tsum[rl].norm();
   }
   inline const HamiltonianTerm * choose(int rl,int i) const {
      return Term(forest(i)->tsum[rl].choose());
   }

   // It changes _tsum[2][] and tree_cache.
   inline void update(const HamiltonianTerm* const &term,int rl) {

      const term_vector_t &kin=kin_adjacency[Hash(term)];
      for(term_vector_t::const_iterator nbr=kin.begin(); nbr!=kin.end(); ++nbr) {
         Hamiltonian::size_type fndex= Hash(*nbr);
         TreeType *i_tree = tree_cache[fndex];
         TreeType *f_tree = forest(*nbr);
         MatrixElement fme = (*nbr)->me(rl);
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

};


class Probabilities : public Configuration {
protected:
   KineticProbabilities Trees;
   int _NBWL;
   GreenOperator<long double> GF;                // Defines the Green operator function
   PotentialEnergies Energy;


public:
   Probabilities(const Hamiltonian &T,const Hamiltonian &P,const GreenOperator<long double> &g) : Configuration(T), Trees(T), GF(g), Energy(T,P) {
      _NBWL=CountBrokenLines();
   }

   inline MatrixElement G(int offset=0) const {
      return GF(NBrokenLines()+offset);   // The value of the Green Operator given the total broken lines and the offset.
   }

   inline double weight(int rl) const {
      _float_accumulator s=0.0;
      for( OffsetMap::size_type i=0; i<Trees.size(); ++i)
         s+=G(Trees.offset(i))*Trees.norm(rl,i);
      return s;
   }

   /* Choose the offset first. */
   const HamiltonianTerm* choose(int rl) const {

      double R=RNG::Uniform()*weight(rl);
      OffsetMap::size_type i=0;
      while((R-=G(Trees.offset(i))*Trees.norm(rl,i))>=0)
         ++i;

      return Trees.choose(rl,i);
   }

   inline const int &NBrokenLines() const {
      return _NBWL;
   }


   inline void update(const HamiltonianTerm* term,int rl,int arflag) {

      Configuration::update(term,rl,arflag);

      /* Updating the number of broken lines.
        if the update is before the configuration change then
        _NBWL+=term->offset(arflag);
        otherwise it is
        _NBWL-=term->offset(!arflag); */

      _NBWL-=term->offset(!arflag);
      Trees.update(term,rl);
      Energy.update(term,rl);

   }

};

}

#endif
