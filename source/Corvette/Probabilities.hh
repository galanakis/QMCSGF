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
  std::vector<uint> _index;    // Map from an offset to consecutive integers
  std::vector<int> _offsets;   // Keeps a sorted list of the offsets

  
  inline int MinOffset() const {return _offsets[0];}    // Minimum possible offset
  inline int MaxOffset() const {return _offsets[_offsets.size()-1];}    // Maximum possible offset

  inline uint nindex() const {return MaxOffset()-MinOffset()+1;}

    
  void initialize(const Hamiltonian &T) {
    // Scan all kinetic terms to find all the lengths
    std::set<int> set_offsets;
    for(Hamiltonian::size_type i=0;i<T.size();++i)
      for(int offset=T[i].minoffset();offset<=T[i].maxoffset();offset+=2)
        set_offsets.insert(offset);

    _index.clear();
    _offsets.clear();
    _offsets.insert(_offsets.begin(),set_offsets.begin(),set_offsets.end());
     
    _index.reserve(nindex());
    for(uint i=0;i<nindex();++i)
      _index.push_back(Infinity); // The default value is set to infinity

    std::vector<int>::const_iterator sit;
    int count=0;
    for(sit=_offsets.begin();sit!=_offsets.end();++sit)
      _index[(*sit)-MinOffset()]=count++;    

  }
  
public:

  typedef std::vector<int>::size_type size_type;

  OffsetMap(const Hamiltonian &T) {initialize(T);}
    
  // Number of different offsets
  inline size_type size() const {return _offsets.size();}
  // OffsetMap(int offset) will return a consecutive integer
  inline uint operator()(int offset) const {return _index[offset-MinOffset()];}
  // OffsetMap[int i] will return the i^th smallest offset
  inline int operator[](int i) const {return _offsets[i];}
  
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
	int _NBWL;

public:
	Configuration(const Hamiltonian &T) {
		std::set<Boson*> indexset;
		for(Hamiltonian::size_type i=0;i<T.size();++i)
			for(HamiltonianTerm::size_type j=0;j<T[i].product().size();++j)
			indexset.insert(T[i].product()[j].particle_id());

		_indices.insert(_indices.begin(),indexset.begin(),indexset.end());		
 
		_NBWL=CountBrokenLines();
		
	}
	
	Configuration(const Configuration &o) : _indices(o._indices), _NBWL(o._NBWL) {}


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
		_NBWL-=term->offset(!arflag);

		for(int i=0;i<UpdatableObject::Objects.size();++i)
			UpdatableObject::Objects[i]->update(term,rl,arflag);

	}

	inline int NBrokenLines() const {return _NBWL;}


};


/*
	class BrokenLines

	It holds a list of the indices with broken lines. It is only used for measurements
	and this is why it is UpdatableObject.

*/

class BrokenLines : public UpdatableObject {
 
public:

	typedef std::map< Boson*,int > BosonDeltaMapType;
	BosonDeltaMapType _broken_lines;
	
	BrokenLines(const std::set<Boson*> &o) : _broken_lines(map(o)) {}
 
	inline void update(const HamiltonianTerm* term) {
		for(HamiltonianTerm::size_type i=0;i<term->product().size();++i) {
			Boson *pind=term->product()[i].particle_id();
			int delta=pind->delta();
			if(delta!=0)
				_broken_lines.insert(std::pair<Boson*,int>(pind,delta));
			else
				_broken_lines.erase(pind);
		}		
	}
  
	inline void update(const HamiltonianTerm* term,int,int) {update(term);}

		// This one must be very fast
	inline const BosonDeltaMapType &operator()() const { return _broken_lines; }
 
	static inline BosonDeltaMapType map(const std::set<Boson*> &o) {
		BosonDeltaMapType vmap;
		std::set<Boson*>::const_iterator it;
		for(it=o.begin();it!=o.end();++it) {
			int delta=(*it)->delta();
			if(delta!=0) vmap.insert(std::pair<Boson*,int>(*it,delta));
		}
		return vmap;		
	}

		// This one is ok to be slow since it is called only in the initializer
	static inline BosonDeltaMapType map(const HamiltonianTerm &term) {
		std::map<Boson*,int> indices;
		for(unsigned int i=0;i<term.product().size();++i) {
			int delta=-term.product()[i].delta();
			if(delta!=0) indices[term.product()[i].particle_id()]=delta;
		}
		return indices;
 	}

};

/*

	class PotentialEnergies
	
	it holds and updates then right and left potential energy.
	
*/

class PotentialEnergies {
	static uint RebuildPeriod;
  const AdjacencyList pot_adjacency;  		// For each term it holds a list of other kinetic terms with one common site
  const Hamiltonian &Potential;       		// Local reference of the potential operators    
	const HamiltonianTerm* Kinetic0;

  _float_accumulator Energies[2];       	// energy of right and left state. It is an accumulator
	std::vector<MatrixElement> EnergyME[2];
  _integer_counter NUpdates;     					// Number of updates since last rebuild. This is only used to fix accumulating floating point errors for the energies

	inline const HamiltonianTerm * Term(Hamiltonian::size_type iterm) const {return &Potential[iterm];}
	inline Hamiltonian::size_type Hash(const HamiltonianTerm* term) const {return term-&Potential[0];} 
	inline Hamiltonian::size_type size() const {return Potential.size();}

	typedef AdjacencyList::const_iterator adjacency_iterator;
	typedef std::pair<adjacency_iterator,adjacency_iterator> adjancency_range;
	inline  adjancency_range potential_adjancency_range(const HamiltonianTerm *term) { return pot_adjacency.range(term-Kinetic0); }

	inline void rebuild() {
    NUpdates=0;
    for(int direction=0;direction<2;++direction) {
			Energies[direction]=0;
			for(uint i=0;i<size();++i) {
				MatrixElement me=Term(i)->me(direction);
				EnergyME[direction][i]=me;
        Energies[direction]+=me;
			}
		}
	}


public:
 
 	PotentialEnergies(const Hamiltonian &T,const Hamiltonian &P) : Potential(P), pot_adjacency(T,P), NUpdates(0), Kinetic0(&T[0]) {
		EnergyME[0].resize(Potential.size());
		EnergyME[1].resize(Potential.size());	
		rebuild();
	}

	inline _float_accumulator operator()(int rl) const {return Energies[rl];}
	
  // It updates the Energies[2], EnergyME[2][]. It needs a potential_adjancency_range and PotentialHash
	inline void update(const HamiltonianTerm* term,int rl) {
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

	inline void update(const HamiltonianTerm* term,int rl,int) {update(term,rl);}

};

uint PotentialEnergies::RebuildPeriod=1000000;


/*
  
	class KineticProbabilities
	
	holds and updates the Trees for each direction.

*/

class KineticProbabilities {
  const Hamiltonian &Kinetic;         // Local copy of the kinetic operators
  const OffsetMap offsets;            // This will map the offsets to consecutive integers
  const AdjacencyList kin_adjacency;  // For each term it holds a list of other kinetic terms with one common site
  inline OffsetMap::size_type noffsets() const {return offsets.size();}
  inline Hamiltonian::size_type Hash(const HamiltonianTerm* term) const {return term-&Kinetic[0];}
	inline const HamiltonianTerm * Term(Hamiltonian::size_type iterm) const {return &Kinetic[iterm];}

  TSum *_tsum[2]; 										// Holds the probability trees
	std::vector<int> offset_cache; 			// remembers the offset of each operator

public:
	
	inline OffsetMap::size_type size() const {return noffsets();}
	inline int offset(int i) const {return offsets[i];}
	
	typedef AdjacencyList::const_iterator adjacency_iterator;
	typedef std::pair<adjacency_iterator,adjacency_iterator> adjancency_range;
	inline  adjancency_range kinetic_adjancency_range(const HamiltonianTerm *term) { return kin_adjacency.range(Hash(term)); } 
	KineticProbabilities(const Hamiltonian &T) : Kinetic(T), offsets(Kinetic), kin_adjacency(Kinetic) {

		for(int direction=0;direction<2;++direction) {
			_tsum[direction]=new TSum[noffsets()];
			for(int count=0;count<noffsets();++count) 
				_tsum[direction][count].resize(Kinetic.size());
		}
		offset_cache.resize(Kinetic.size());

		for(int direction=0;direction<2;++direction)
			for(int count=0;count<noffsets();++count) 
				_tsum[direction][count].reset();
				
    for(int direction=0;direction<2;++direction) {
      for(uint i=0;i<Kinetic.size();++i) {
				int ioffset = offsets(Kinetic[i].offset());
				offset_cache[i] = ioffset;
				_tsum[direction][ ioffset ].update(i,Kinetic[i].me(direction));
			}
		}		

		
	}

	~KineticProbabilities() {
		delete [] _tsum[0];
		delete [] _tsum[1];
	}

 
	inline double norm(int rl,int i) const {return _tsum[rl][i].norm();}
	inline const HamiltonianTerm * choose(int rl,int i) const {return Term(_tsum[rl][i].choose());}

  // It changes _tsum[2][] and offset_cache. It needs kinetic_adjancency_range and KineticHash
  inline void update(const HamiltonianTerm* term,int rl) { 
    
		const adjancency_range &kin=kinetic_adjancency_range(term);
    for(adjacency_iterator nbr=kin.first;nbr!=kin.second;++nbr) {
      Hamiltonian::size_type fndex= Hash(*nbr);
      int ioffset = offset_cache[fndex];
      int foffset = offsets((*nbr)->offset());
      MatrixElement fme = (*nbr)->me(rl);
			_tsum[ rl][foffset].update(fndex,fme);
			if(ioffset!=foffset) {
				offset_cache[fndex] = foffset;
				MatrixElement jme = _tsum[!rl][ioffset].element(fndex); 
				_tsum[ rl][ioffset].update(fndex,0);
				_tsum[!rl][ioffset].update(fndex,0);
				_tsum[!rl][foffset].update(fndex,jme); 
			}
    }    
  }

	inline void update(const HamiltonianTerm* term,int rl,int) { update(term,rl); }
 
};

 
class Probabilities : public Configuration {
	
	KineticProbabilities Trees;
	PotentialEnergies Energies;

 	GreenOperator<long double> GF;                // Defines the Green operator function

		
public:
	Probabilities(const Hamiltonian &T,const Hamiltonian &P) : Configuration(T), Trees(T), Energies(T,P) {}
 
	inline void GreenInit(int nsites,int cutoff) {GF.initialize(nsites,cutoff);}
  inline MatrixElement G(int offset=0) const {return GF(NBrokenLines()+offset);}  // The value of the Green Operator given the total broken lines and the offset.

	inline _float_accumulator Energy(int rl) const {return Energies(rl);}

 
  inline double weight(int rl) const {
    _float_accumulator s=0.0;
    for(uint i=0;i<Trees.size();++i) 
      s+=G(Trees.offset(i))*Trees.norm(rl,i);
    return s;
  }

	/* Choose the offset first. */
	const HamiltonianTerm* choose(int rl) const {

		double R=RNG::Uniform()*weight(rl);
		int i=0;
		while((R-=G(Trees.offset(i))*Trees.norm(rl,i))>=0)
			++i;

		return Trees.choose(rl,i);
	}
  

  inline void update(const HamiltonianTerm* term,int rl,int arflag) {
      
		Configuration::update(term,rl,arflag);
    Trees.update(term,rl);
		Energies.update(term,rl);

	}

  
 
};

}

#endif
