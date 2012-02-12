#ifndef __PROBABILITIES__
#define __PROBABILITIES__

#include "RandomNumberGenerator.hh"
#include "HamiltonianTerm.hh"
#include "GreenOperator.hh"
#include "AdjacencyList.hh"
#include "TSum.hh"
#include <vector>
#include <list>
#include <set>
#include <queue>

namespace SGF {


/* Crop small doubles. Use them when you expect a finite double
   in which case you can disregard the small ones as numerical
   errors */


/* 
  class OffsetMap
  Helper class. It defines a map that can convert offsets into consecutive
  integers and vice versa.
  
  For an operator product at the same index, the offset is defined as
  Abs(dn+NR-NL)-Abs(NR-NL) and for a given dn it takes values -dn,-dn+2,...,dn-2,dn.
  If we have products at different indices, the offset is additive. This reminds a 
  bit addition of spins. The maximum offset is the Abs(delta()) 
  (# of creation -# of annihilation operators) of the product.
  The step of the offset is always 2. Now if we have a lot of terms, some of them
  will have an even and some of them will have an odd and some others and even
  offset. If all of them are even then we can set the offset step=2, otherwise
  it will be 1.
  
  All this does not assume the case of hard core bosons, in which the number
  of offsets is limitted.
  In this case we need to consider Abs(dn+DN)-Abs(DN), -Nmax<=(DN=NR-NL)<=Nmax.
  It is easy tos show that given the constraint:
  |dn|-2*Min(|dn|,Nmax) <= Abs(|dn|+DN)-Abs(DN) <= |dn|
  
  This is handled by IndexedProductElement, which can return the minimum and
  maximum offset.
  Given that we can scan through all kinetic terms and generate a list of 
  offsets (between min and max with step 2). Then we can map them to
  sequential integers.

*/

class OffsetMap {
public:
  typedef std::vector<int>::size_type size_type;
private:
  std::vector<size_type> _index;    // Map from an offset to consecutive integers
  std::vector<int> _offsets;   // Keeps a sorted list of the offsets
	size_type _size;
  
    
  void initialize(const Hamiltonian &T) {
    // Scan all kinetic terms to find all the lengths
    std::set<int> set_offsets;
    for(Hamiltonian::size_type i=0;i<T.size();++i)
      for(int offset=T[i].minoffset();offset<=T[i].maxoffset();offset+=2)
        set_offsets.insert(offset);
    
		_size=set_offsets.size();
		
    _offsets.clear();
    _offsets.insert(_offsets.begin(),set_offsets.begin(),set_offsets.end());

		int MinOffset=_offsets[0];
	  int MaxOffset=_offsets[_size-1];
	  size_type nindex=MaxOffset-MinOffset+1;
	  
    _index.clear();
    _index.reserve(nindex);
    for(size_type i=0;i<nindex;++i)
      _index.push_back(Infinity); // The default value is set to infinity

    std::vector<int>::const_iterator sit;
    int count=0;
    for(sit=_offsets.begin();sit!=_offsets.end();++sit)
      _index[(*sit)-MinOffset]=count++;    

  }
  
public:

  OffsetMap(const Hamiltonian &T) {initialize(T);}
    
  // Number of different offsets
  inline const size_type &size() const {return _size;}
  // OffsetMap(int offset) will return a consecutive integer
  inline const size_type &operator()(int offset) const {return _index[offset-_offsets[0]];}
  // OffsetMap[int i] will return the i^th smallest offset
  inline const int &operator[](int i) const {return _offsets[i];}
  
};

/*
	class UpdatableObject
	
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



class UpdatableObject {
	static std::vector<UpdatableObject*> Objects;   // This lists contains pointers to all UpdatableObjects.
public:
	typedef std::vector<UpdatableObject*>::size_type size_type;
	UpdatableObject() {Objects.push_back(this);}
	virtual inline void update(const HamiltonianTerm*,int,int) = 0;

	friend class Configuration;
};

std::vector<UpdatableObject*> UpdatableObject::Objects;
 


/*
	class Configuration
	
	holds the indices of the Bosons, together the number of broken world lines.
	When updated it also updates all UpdatableObjects.
*/

class Configuration {
protected:	
	typedef std::vector<Boson*> boson_vector;
	boson_vector _indices;

public:
	Configuration(const Hamiltonian &T) {
		std::set<Boson*> indexset;
		for(Hamiltonian::size_type i=0;i<T.size();++i)
			for(HamiltonianTerm::size_type j=0;j<T[i].product().size();++j)
			indexset.insert(T[i].product()[j].particle_id());

		_indices.insert(_indices.begin(),indexset.begin(),indexset.end());		
 		
	}
	
	Configuration(const Configuration &o) : _indices(o._indices) {}

	
	int CountBrokenLines() const {
		int result=0;
		for(std::vector<Boson*>::size_type i=0;i<_indices.size();++i) 
     result+=Abs(_indices[i]->delta()); 
	  return result;
	}

	std::set<Boson*> GetListBrokenLines() const {
		std::set<Boson*> _broken_lines;
		for(std::vector<Boson*>::const_iterator it=_indices.begin();it!=_indices.end();++it) {
			if((*it)->delta() != 0) 
				_broken_lines.insert(*it);
		}
		return _broken_lines;		
	}

  
	// call this update you update the occupancies
	inline void update(const HamiltonianTerm* term,int rl,int arflag) {
		/*********************************\  
			The configuration changes here
		\*********************************/  
    term->update_psi(rl,arflag); 

		/* Updating the number of broken lines. 
			if the update is before the configuration change then
			_NBWL+=term->offset(arflag);
			otherwise it is
			_NBWL-=term->offset(!arflag); 
		*/

		for(UpdatableObject::size_type i=0;i<UpdatableObject::Objects.size();++i)
			UpdatableObject::Objects[i]->update(term,rl,arflag);

	}

   
};


/*
	class BrokenLines

	It holds a list of the indices with broken lines. It is only used for measurements
	and this is why it is UpdatableObject.

*/


class BrokenLines : public UpdatableObject {
	typedef std::set<Boson*> boson_set;
	boson_set _broken_lines;
 
public:
  
	typedef std::vector<std::pair<Boson*,int> > BosonDeltaMapType;

	BrokenLines(const std::set<Boson*> &o) :  _broken_lines(o) {}
 
	inline void update(const HamiltonianTerm* term,int,int) {
		for(HamiltonianTerm::size_type i=0;i<term->product().size();++i) {
			Boson *pind=term->product()[i].particle_id();
			if(pind->delta()!=0)
				_broken_lines.insert(pind);
			else
				_broken_lines.erase(pind);
		}		
	}

		// This one must be very fast
	inline BosonDeltaMapType map() const {
		BosonDeltaMapType vmap;
		vmap.reserve(_broken_lines.size());
		std::set<Boson*>::const_iterator it;
		for(it=_broken_lines.begin();it!=_broken_lines.end();++it) {
			int delta=(*it)->delta();
			if(delta!=0) vmap.push_back(std::pair<Boson*,int>(*it,delta));
		}
		return vmap;
	}	
  
	inline const BosonDeltaMapType operator()() const { return map(); } 

		// This one is ok to be slow since it is called only in the initializer
	static inline BosonDeltaMapType map(const HamiltonianTerm &term) {
		std::map<Boson*,int> indices;
		for(unsigned int i=0;i<term.product().size();++i) {
			int delta=-term.product()[i].delta();
			if(delta!=0) indices[term.product()[i].particle_id()]=delta;
		}

		BosonDeltaMapType vmap;
		vmap.reserve(indices.size());
		std::map<Boson*,int>::const_iterator it;
		for(it=indices.begin();it!=indices.end();++it)
			vmap.push_back(*it);

		return vmap;
	}

};

/*
	class TermCount
	
	measures how many times each operator gets inserted minus the times 
	it is removed from the right of the GreenOperator.
	
*/
 
class TermCount : public UpdatableObject {

	const HamiltonianTerm * Kinetic0;
	typedef std::vector<long> count_type;
	typedef count_type::iterator iterator;
	count_type _count;

public:
	TermCount(const Hamiltonian &T) : _count(T.size()), Kinetic0(&T[0]) { reset(); }
	
	inline void update(const HamiltonianTerm* term,int rl,int ar) { if(rl==RIGHT) _count[term-Kinetic0]+=Sign[ar]; }

	inline void reset() {
		for(iterator ptr=_count.begin(); ptr!=_count.end(); ++ptr)
			*ptr=0;
	}
	
};


/*

	class PotentialEnergies
	
	it holds and updates then right and left potential energy.
	
*/

class PotentialEnergies {
	static unsigned int RebuildPeriod;
  const Hamiltonian &Potential;       		// Local reference of the potential operators    
  const AdjacencyList pot_adjacency;  		// For each term it holds a list of other kinetic terms with one common site
	const HamiltonianTerm* Kinetic0;

  _integer_counter NUpdates;     					// Number of updates since last rebuild. This is only used to fix accumulating floating point errors for the energies
  _float_accumulator Energies[2];       	// energy of right and left state. It is an accumulator
	std::vector<MatrixElement> EnergyME[2];

	inline const HamiltonianTerm * Term(Hamiltonian::size_type iterm) const {return &Potential[iterm];}
	inline Hamiltonian::size_type Hash(const HamiltonianTerm* term) const {return term-&Potential[0];} 

	typedef AdjacencyList::const_iterator adjacency_iterator;
	typedef std::pair<adjacency_iterator,adjacency_iterator> adjancency_range;
	inline  const adjancency_range &potential_adjancency_range(const HamiltonianTerm *term) const { return pot_adjacency.range(term-Kinetic0); }

	inline void rebuild() {
    NUpdates=0;
    for(int direction=0;direction<2;++direction) {
			Energies[direction]=0;
			for(HamiltonianTerm::size_type i=0;i<Potential.size();++i) {
				MatrixElement me=Term(i)->me(direction);
				EnergyME[direction][i]=me;
        Energies[direction]+=me;
			}
		}
	}


public:
 
 	PotentialEnergies(const Hamiltonian &T,const Hamiltonian &P) : Potential(P), pot_adjacency(T,P), Kinetic0(&T[0]), NUpdates(0) {
		EnergyME[0].resize(Potential.size());
		EnergyME[1].resize(Potential.size());	
		rebuild();
	}

	PotentialEnergies(const PotentialEnergies &o) : Potential(o.Potential), pot_adjacency(o.pot_adjacency), Kinetic0(o.Kinetic0), NUpdates(o.NUpdates) {
		EnergyME[0]=o.EnergyME[0];
		EnergyME[1]=o.EnergyME[1];
		Energies[0]=o.Energies[0];
		Energies[1]=o.Energies[1];
	}

	inline const _float_accumulator &operator()(int rl) const {return Energies[rl];}
	
  // It updates the Energies[2], EnergyME[2][]. It needs a potential_adjancency_range and PotentialHash
	inline void update(const HamiltonianTerm* const &term,int rl) {
		const adjancency_range &pot=potential_adjancency_range(term);
    for(adjacency_iterator nbr=pot.first;nbr!=pot.second;++nbr) {
			MatrixElement me=(*nbr)->me(rl);
			std::vector<MatrixElement>::size_type i=Hash(*nbr);
			Energies[rl]+=me-EnergyME[rl][i];
			EnergyME[rl][i]=me;
    }
    if((++NUpdates % RebuildPeriod)==0)
			rebuild();

	}	

};

uint PotentialEnergies::RebuildPeriod=1000000;

 
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


class ForestType {
  const OffsetMap offsets;            // This will map the offsets to consecutive integers
  TreeType *trees; 										// Holds the probability trees 
public:
	TreeType *operator()(const HamiltonianTerm *term) const {return &trees[offsets(term->offset())]; }
	TreeType *operator()(unsigned int i) const {return &trees[i];}
	inline const OffsetMap::size_type &size() const {return offsets.size();}
	inline const int & offset(int i) const {return offsets[i];}

	ForestType(const Hamiltonian &K) : offsets(K) {
		trees=new TreeType[size()];
		for(OffsetMap::size_type i=0;i<size();++i)
			trees[i].resize(K.size());
		for(OffsetMap::size_type i=0;i<size();++i) 
			trees[i].reset();
		
	}
	~ForestType() {
		delete [] trees;
	}
};



/*
  
	class KineticProbabilities
	
	holds and updates the Trees for each direction.

*/



class KineticProbabilities {
  const Hamiltonian &Kinetic;         // Local copy of the kinetic operators
  const AdjacencyList kin_adjacency;  // For each term it holds a list of other kinetic terms with one common site
  inline Hamiltonian::size_type Hash(const HamiltonianTerm* term) const {return term-&Kinetic[0];}
	inline const HamiltonianTerm * Term(Hamiltonian::size_type iterm) const {return &Kinetic[iterm];}
	std::vector<TreeType*> tree_cache; 			// remembers the offset of each operator

	ForestType forest;

public:
	
	inline const OffsetMap::size_type &size() const {return forest.size();}
	inline const int &offset(int i) const {return forest.offset(i);}
	
	typedef AdjacencyList::const_iterator adjacency_iterator;
	typedef std::pair<adjacency_iterator,adjacency_iterator> adjancency_range;
	inline  const adjancency_range &kinetic_adjancency_range(const HamiltonianTerm *term) const { return kin_adjacency.range(Hash(term)); } 
	KineticProbabilities(const Hamiltonian &T) : Kinetic(T), kin_adjacency(Kinetic), forest(Kinetic) {
    
		tree_cache.resize(Kinetic.size());

    for(int direction=0;direction<2;++direction) {
      for(Hamiltonian::size_type i=0;i<Kinetic.size();++i) {
				TreeType *tree = forest(&Kinetic[i]);
				tree_cache[i] = tree;
				tree->tsum[direction].update(i,Kinetic[i].me(direction));
			}
		}		

		
	}

	~KineticProbabilities() {
	}

 
	inline const _float_accumulator &norm(int rl,int i) const {return forest(i)->tsum[rl].norm();}
	inline const HamiltonianTerm * choose(int rl,int i) const {return Term(forest(i)->tsum[rl].choose());}

  // It changes _tsum[2][] and tree_cache. It needs kinetic_adjancency_range and KineticHash
  inline void update(const HamiltonianTerm* const &term,int rl) { 
    
		const adjancency_range &kin=kinetic_adjancency_range(term);
    for(adjacency_iterator nbr=kin.first;nbr!=kin.second;++nbr) {
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
protected:
	PotentialEnergies Energy;

		
public:
	Probabilities(const Hamiltonian &T,const Hamiltonian &P) : Configuration(T), Trees(T), Energy(T,P) {
		_NBWL=CountBrokenLines();
	}
 
	inline void GreenInit(int nsites,int cutoff) {GF.initialize(nsites,cutoff);}
  inline MatrixElement G(int offset=0) const {return GF(NBrokenLines()+offset);}  // The value of the Green Operator given the total broken lines and the offset.


 
  inline double weight(int rl) const {
    _float_accumulator s=0.0;
    for( OffsetMap::size_type i=0;i<Trees.size();++i) 
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
  
	inline const int &NBrokenLines() const {return _NBWL;}
	
	
  inline void update(const HamiltonianTerm* term,int rl,int arflag) {
      
		Configuration::update(term,rl,arflag);
		_NBWL-=term->offset(!arflag);
    Trees.update(term,rl);
    Energy.update(term,rl);

	}

  
 
};

}

#endif
