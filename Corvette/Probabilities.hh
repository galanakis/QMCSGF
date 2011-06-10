#ifndef __PROBABILITIES__
#define __PROBABILITIES__

#include "RandomNumberGenerator.hh"
#include "HamiltonianTerm.hh"
#include "GreenOperator.hh"
#include "AdjacencyList.hh"
#include <vector>
#include <list>
#include <set>

namespace SGF {

uint RebuildFrequency=1000000;

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
    for(int i=0;i<T.size();++i)
      for(int offset=T[i].minoffset();offset<=T[i].maxoffset();offset+=2)
        set_offsets.insert(offset);

    _index.clear();
    _offsets.clear();
    _offsets.insert(_offsets.begin(),set_offsets.begin(),set_offsets.end());
     
    _index.reserve(nindex());
    for(int i=0;i<nindex();++i)
      _index.push_back(Infinity); // The default value is set to infinity

    std::vector<int>::const_iterator sit;
    int count=0;
    for(sit=_offsets.begin();sit!=_offsets.end();++sit)
      _index[(*sit)-MinOffset()]=count++;    

  }
  
public:

  OffsetMap(const Hamiltonian &T) {initialize(T);}
    
  // Number of different offsets
  inline uint size() const {return _offsets.size();}
  // OffsetMap(int offset) will return a consecutive integer
  inline uint operator()(int offset) const {return _index[offset-MinOffset()];}
  // OffsetMap[int i] will return the i^th smallest offset
  inline int operator[](int i) const {return _offsets[i];}
  
};

/*
	class TSum
	This is a helper class which consists of a vector of matrix elements.
	The class interprets this vector as a balanced binary tree. To get the 
	path of the i^th element we just look at the binary representation of i+1,
	where the most significant digit corresponds to ancestor nodes.
	The number of elements of the tree are meant to be fixed. Each
	element contains a value which is interpreted as a relative probability
	At each node we don't need to store that relative probability, but
	only the total sum of the probability of the node and all its children.
	The probability can be deduced by subtracting the sum of a node
	minus the sums of its children.	When a probability is changed this change 
	is propagated to the node's parents with logarithmic complexity. 
	Selecting a term according to its relative probability is also done
	with logarithmic complexity algorithm. We start from the root of the
	tree and select either a left or right branch or the root node itself
	(three way selection). We continue in the branch we selected or return
	the root.
	
	The key functions are:
	update(index,me): change the probability of the index by "me".
	choose(): randomly chose an index using its relative probability (sum-sum_Left-sum_Right)
	norm(): returns the sum of all relative probabilities
	
*/

class TSum {
	std::vector<MatrixElement> _sums;
	long double _head; // a copy of the _sums[0] is stored here. The redundancy is used for error tracking.
	
	// Chop off small numerical values. The thresshold is Tolerance which is set as the minimum coefficient of any kinetic operator.
	long double Crop(long double x) const { return (fabs(x)>Tolerance)?x:0; }

public:
  static MatrixElement Tolerance;

	TSum() {}
	TSum(uint N) : _sums(N,MatrixElement(0)),_head(0) {}
	TSum(const TSum &o) : _sums(o._sums),_head(0) {}
	~TSum() {}
	
	/* resizes the vector */
	inline void resize(uint N) {_sums.resize(N,MatrixElement(0));}
	
	inline void update(uint index,MatrixElement me) {
		if(me!=0) {
			index++;
			while(index>0) {
			  _sums[index-1]+=me;
        index>>=1;
	    }
			_head+=me;
    }
	}

	inline uint choose() const {
		int _nterms=_sums.size();
		int index=0;
		while(index<_nterms) {
			int indr=(index+1)<<1;
			int indl=indr-1;

			double w =Crop(_sums[index]);
			double wr=(indr<_nterms) ? Crop(_sums[indr]) : MatrixElement(0);
			double wl=(indl<_nterms) ? Crop(_sums[indl]) : MatrixElement(0);

			if(w*RNG::Uniform()>=wl+wr)
				return index; 
			else
				index=((wl+wr)*RNG::Uniform()>=wr)?indl:indr;
		}
		std::cout<<"TSum::choose has reached a dead end. Cannot choose term"<<std::endl;
		exit(13);
		return _nterms;

	}

	/* Makes all elements zero */
	inline void reset() {
		_head=0;
		for(uint i=0;i<_sums.size();++i)
			_sums[i]=MatrixElement(0);	
	}
	
	/* Returns the sum of all the partial probabilities which is stored at the root */
	inline MatrixElement norm() const { return Crop(_head); }
	inline MatrixElement error() const {return fabs(_head-_sums[0]);}

};

MatrixElement TSum::Tolerance;

class ProbabilitiesBase {
  static MatrixElement GetMinCoefficient(const Hamiltonian &T) {
      /* Finding the minimum of all Kinetic operator coefficients */
      MatrixElement result=T[0].coefficient();
      for(int term=0;term<T.size();++term)
        result=Min(result,T[term].coefficient());
      return result;	
  }
	
	
  // Scan all kinetic terms to find all the indices
  static std::vector<Boson*> GetIndices(const Hamiltonian &T) {
    std::set<Boson*> indexset;
    for(int i=0;i<T.size();++i)
      for(int j=0;j<T[i].product().size();++j)
        indexset.insert(T[i].product()[j].particle_id());

	std::vector<Boson*> result;
    result.insert(result.begin(),indexset.begin(),indexset.end());
	return result;
	
}
	
protected:
  const Hamiltonian &Kinetic;         // Local copy of the kinetic operators
  const Hamiltonian &Potential;       // Local copy of the potential operators    
  const AdjacencyList kin_adjacency;  // For each term it holds a list of other kinetic terms with one common site
  const AdjacencyList pot_adjacency;  // For each term it holds a list of other kinetic terms with one common site

  const std::vector<Boson*> _indices; // Vector containing all indices appearing in the row terms  
  const OffsetMap offsets;            // This will map the offsets to consecutive integers

  GreenOperator GF;                // Defines the Green operator function


  inline uint noffsets() const {return offsets.size();}
  inline uint term_index(const HamiltonianTerm* term) const {return term-&Kinetic[0];}
	inline const HamiltonianTerm * index_term(uint iterm) const {return &Kinetic[iterm];}

	uint _ensemble;
	
	uint nextra() const {return 2*_indices.size()*_ensemble;}
	uint nregular() const {return Kinetic.size()-nextra();}
public:
  inline void GreenInit(int nsites) {GF.initialize(nsites);}

  ProbabilitiesBase(const Hamiltonian &T,const Hamiltonian &V) :  _indices(GetIndices(T)), Kinetic(T), offsets(Kinetic), Potential(V), kin_adjacency(Kinetic), pot_adjacency(Kinetic,Potential) {
		/* Initialize the Green operator function.
		The number of sites is just the number of different
		indices appearing in the Kinetic operators */
		GF.initialize(_indices.size());
	  TSum::Tolerance=0.5*GetMinCoefficient(T);
	
		// The SGFContainer class has put the extra kinetic terms at the end. So we only count to indirectly extract the ensemble.
		uint _nextra=0;
		for(int i=0;i<Kinetic.size();++i) {
			if(Kinetic[i].atom()) _nextra++;
		}
		
		_ensemble=(_nextra==2*_indices.size())?GrandCanonical:Canonical;

}

};

/*
  class Probabilities
  This is the core class of the algorithm. It is able to evalute
  the probabilities that are necessary for the SGF algorithm
  with directed updates.
  
  To use the class one needs to declare the variable Kinetic
  which is of type Hamiltonian and stores all Kinetic energy terms.
  Then it can be called to provide all sorts of probabilities.
  
  The class provices
  int choose(): randomly chooses a term
  update(int,int): which adds or removes a term
  weight(): returns the total probability of the tree (the weight of the root)

  The class is itinialized only by an array of type Hamiltonian.

  Typical Usage of the class:
  
  static Probabilites::Kinetic=Kinetic_Terms;   // Gives the path to the Kinetic energy Terms;
  Probabilities p();                            // Initializes the internal structures
  unsigned long p.choose(RIGHT);                // pick a term for RIGHT addition. Similarly for LEFT addition.
  p.update(term,RIGHT or LEFT,ADD or REMOVE);   // update the probabilities, BEFORE updating the Psi's!
  p.weight(RIGHT or LEFT);                      // gives the probability normalization constant
  p.rebuild();                                  // Fixes accumulated errors, by starting from scratch. This is slow
  
*/


class Probabilities : public ProbabilitiesBase {
  double Energies[2];              // energy of right and left state
  std::set<Boson*> _broken_lines;  // A set of the boson indices that are broken.
  long long NUpdates;              // Number of updates since last rebuild

	bool ExtraLock;                  // whether to ignore the extra terms or not. 0 means ignore, 1 means accept.

  TSum *_tsum;
	TSum & tsum(int direction,int indoffset) const {return _tsum[direction*noffsets()+indoffset];}
public:

  inline double Energy(int direction) const {return Energies[direction];}

  inline double G() const {return GF(NBrokenLines());};                  // The value of the Green Operator for the given broken lines
  inline double G(int offset) const {return GF(NBrokenLines()+offset);}  // The value of the Green Operator given the total broken lines and the offset.
  inline int NBrokenLines() const {return _broken_lines.size();}
  const std::set<Boson*> &ListBrokenLines() const {return _broken_lines;};  // A set of the boson indices that are broken.


  inline double weight(int rl) const {
    double s=0.0;
    for(int i=0;i<noffsets();++i) 
      s+=G(offsets[i])*tsum(rl,i).norm();
    return s;
  }


	Probabilities(const Hamiltonian &T,const Hamiltonian &V) : ProbabilitiesBase(T,V),ExtraLock(1) {
		_tsum=new TSum[2*noffsets()];
		for(int rl=0;rl<2;++rl) 
			for(int i=0;i<noffsets();++i) {
				tsum(rl,i).resize(Kinetic.size());
			}
			
		rebuild();
	}

	~Probabilities() {delete [] _tsum;}


	/* Choose the offset first. This version
	calls TSum rather than weight and so it's faster */
	const HamiltonianTerm* choose(int rl) const {

		double R=RNG::Uniform()*weight(rl);
		int i=0;
		while((R-=G(offsets[i])*tsum(rl,i).norm())>=0)
			++i;

		return index_term(tsum(rl,i).choose());
	}

	/* makes everything zero */
	inline void reset() {
		_broken_lines.clear();
		for(int rl=0;rl<2;++rl) 
    	for(int i=0;i<noffsets();++i)
				tsum(rl,i).reset();
    Energies[LEFT]=Energies[RIGHT]=0;  
    NUpdates=0;
		
	}

  // Evaluates the matrix elements and populates the trees
  void rebuild() {
    
		reset();
		
    for(int i=0;i<_indices.size();++i) 
    	if(_indices[i]->delta()!=0) 
    		_broken_lines.insert(_indices[i]);

    for(int rl=0;rl<2;++rl) {
      for(uint i=0;i<Kinetic.size();++i)
				tsum(rl,offsets(Kinetic[i].offset())).update(i,Kinetic[i].me(rl));
			for(uint i=0;i<Potential.size();++i)
        Energies[rl]+=Potential[i].me(rl);
		}
        

  }


  inline void update(const HamiltonianTerm* term,int rl,int arflag) {
		// If an extra term appears, the toggle the lock.
		if(term->atom()) ExtraLock!=ExtraLock;
	
    /* This needs to be confirmed: the tree needs to be rebuilt from time
      to time to fix accumulated floating points errors */
    if((++NUpdates % RebuildFrequency)==0) {
		// Slow update: update the state and rebuild
			term->update_psi(rl,arflag);
			rebuild();
		}
		else
		 	fast_update(term,rl,arflag);
		
	}
  
  inline void fast_update(const HamiltonianTerm* term,int rl,int arflag) { 
		uint index=term_index(term);

    adjacency_list_t::const_iterator nbr;
    for(nbr=kin_adjacency[index].begin();nbr!=kin_adjacency[index].end();++nbr) {
      int fndex= term_index(nbr->term());
      int ioffset = nbr->offset();
      int foffset = nbr->offset(arflag);
      MatrixElement ime = nbr->me(rl);
      MatrixElement fme = nbr->me(rl,arflag);
      // The if statements decrease the time by 30%
      if(ioffset!=foffset) {
        MatrixElement jme =nbr->me(!rl);
        // Note that those statements are parallelizable
				tsum( rl,offsets(ioffset)).update(fndex,-ime);
				tsum( rl,offsets(foffset)).update(fndex,+fme);
				tsum(!rl,offsets(ioffset)).update(fndex,-jme);
				tsum(!rl,offsets(foffset)).update(fndex,+jme);
      }
      else 
				tsum( rl,offsets(ioffset)).update(fndex,fme-ime);

    }
    
    for(nbr=pot_adjacency[index].begin();nbr!=pot_adjacency[index].end();++nbr) 
      Energies[rl]+=nbr->me(rl,arflag)-nbr->me(rl);
      
    term->update_psi(rl,arflag);   
    
    for(int i=0;i<term->product().size();++i) {
        Boson *pind=term->product()[i].particle_id();
        if(pind->delta()!=0)
           _broken_lines.insert(pind);
        else
           _broken_lines.erase(pind);
    }
    
  }

};

}

#endif
