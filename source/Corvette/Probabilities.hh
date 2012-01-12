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


class Configuration {
protected:	
	typedef std::vector<Boson*> boson_vector;
	boson_vector _indices;
	int _NBWL;
 	GreenOperator<long double> GF;                // Defines the Green operator function

public:
	Configuration(const Hamiltonian &T) {
		std::set<Boson*> indexset;
		for(Hamiltonian::size_type i=0;i<T.size();++i)
			for(HamiltonianTerm::size_type j=0;j<T[i].product().size();++j)
			indexset.insert(T[i].product()[j].particle_id());

		_indices.insert(_indices.begin(),indexset.begin(),indexset.end());		
 
		rebuild();
		
	}
	
	Configuration(const Configuration &o) : _indices(o._indices) {}

	inline void GreenInit(int nsites,int cutoff) {GF.initialize(nsites,cutoff);}
  inline MatrixElement G() const {return GF(NBrokenLines());};                  // The value of the Green Operator for the given broken lines
  inline MatrixElement G(int offset) const {return GF(NBrokenLines()+offset);}  // The value of the Green Operator given the total broken lines and the offset.

	void rebuild() {
		_NBWL=0;
		for(std::vector<Boson*>::size_type i=0;i<_indices.size();++i) 
     _NBWL+=Abs(_indices[i]->delta()); 
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
	}

	inline int NBrokenLines() const {return _NBWL;}

	inline const boson_vector &configuration() const {return _indices;}

};


class BrokenLines {
	typedef std::set<Boson*> boson_set;
	boson_set _broken_lines;
 
	void initialize(const std::vector<Boson*> &_indices) {
		for(std::vector<Boson*>::const_iterator it=_indices.begin();it!=_indices.end();++it) {
			if((*it)->delta() != 0) 
				_broken_lines.insert(*it);
		}		
	}
public:

	typedef std::vector< std::pair<Boson*,int> > BosonDeltaMapType;

	BrokenLines(const std::vector<Boson*> &_indices) { initialize(_indices); }
	BrokenLines(const Configuration &o) { initialize(o.configuration()); }
 
	inline void update(const HamiltonianTerm* term) {
		for(HamiltonianTerm::size_type i=0;i<term->product().size();++i) {
			Boson *pind=term->product()[i].particle_id();
			if(pind->delta()!=0)
				_broken_lines.insert(pind);
			else
				_broken_lines.erase(pind);
		}		
	}

	inline const boson_set &broken_lines() const { return _broken_lines;}

		// This one must be very fast
	inline BosonDeltaMapType operator()() const {
		BosonDeltaMapType vmap;
		vmap.reserve(_broken_lines.size());
		std::set<Boson*>::const_iterator it;
		for(it=_broken_lines.begin();it!=_broken_lines.end();++it) {
			int delta=(*it)->delta();
			if(delta!=0) vmap.push_back(std::pair<Boson*,int>(*it,delta));
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

		BosonDeltaMapType vmap;
		vmap.reserve(indices.size());
		std::map<Boson*,int>::const_iterator it;
		for(it=indices.begin();it!=indices.end();++it)
			vmap.push_back(*it);

		return vmap;
	}

};


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
};

uint PotentialEnergies::RebuildPeriod=1000000;

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

		rebuild();
		
	}

	~KineticProbabilities() {
		delete [] _tsum[0];
		delete [] _tsum[1];
	}

	// Note: it changes offset_cache and _tsum[2][]. It needs the Kinetic[i].me() and Kinetic[i].offset()
	inline void rebuild() {

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
 
};

 
class Probabilities {
	
	KineticProbabilities Trees;
	PotentialEnergies Energies;
	Configuration Psi; 
	BrokenLines BLines;
		
public:
	Probabilities(const Hamiltonian &T,const Hamiltonian &P) : Psi(T), Trees(T), Energies(T,P), BLines(Psi) {}
 
	inline int NBrokenLines() const {return Psi.NBrokenLines();}

  inline MatrixElement G() const {return Psi.G();};                  // The value of the Green Operator for the given broken lines
  const std::set<Boson*> &ListBrokenLines() const {return BLines.broken_lines();};  // A set of the boson indices that are broken.
  inline void GreenInit(int nsites,int cutoff) {Psi.GreenInit(nsites,cutoff);}
	inline _float_accumulator Energy(int rl) const {return Energies(rl);}

 
  inline double weight(int rl) const {
    _float_accumulator s=0.0;
    for(uint i=0;i<Trees.size();++i) 
      s+=Psi.G(Trees.offset(i))*Trees.norm(rl,i);
    return s;
  }

	/* Choose the offset first. */
	const HamiltonianTerm* choose(int rl) const {

		double R=RNG::Uniform()*weight(rl);
		int i=0;
		while((R-=Psi.G(Trees.offset(i))*Trees.norm(rl,i))>=0)
			++i;

		return Trees.choose(rl,i);
	}
  

  inline void update(const HamiltonianTerm* term,int rl,int arflag) {
      
		Psi.update(term,rl,arflag);
    Trees.update(term,rl);
		Energies.update(term,rl);
		BLines.update(term);

	}

  
 
};

}

#endif
