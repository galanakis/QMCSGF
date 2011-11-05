#ifndef __TSUM__
#define __TSUM__

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
	
	In this class N terms are stored in a tree with 2*P-1 where P=2^d is the 
	smallest power of 2 which is greater than N. The partial sums are stored 
	in the first P-1 terms and matrix elements of the terms are in the last P 
	terms. The size of the tree is thus size=2*P-1 and the terms start
	at size/2=(2*(P-1)+1)/2=P-1.

	The key functions are:
	update(index,me): change the probability of the index by "me".
	choose(): randomly chose an index using its relative probability (sum-sum_Left-sum_Right)
	norm(): returns the sum of all relative probabilities
	
*/


  
class TSum {
	std::vector<MatrixElement> _sums;
 
  typedef std::vector<MatrixElement>::size_type index_type;
	inline index_type nsums() const {return _sums.size()/2;} // This is where the terms start.
 
	// Chop off small numerical values. The thresshold is Tolerance which is set as the minimum coefficient of any kinetic operator.
	long double Crop(long double x) const { return (fabs(x)>Tolerance)?x:0; }

public:
  static MatrixElement Tolerance;

	TSum() : _sums() {}
	TSum(const TSum &o) : _sums(o._sums) {}
	~TSum() {}
	
	/* resizes the vector */
	inline void resize(index_type NTerms) {
		index_type _size=1;
		while(_size<2*NTerms) _size<<=1;
		_sums.resize(_size-1,MatrixElement(0));
	}

	
	inline void update(index_type index,MatrixElement me) {
		if(me!=0) {
			index+=1+nsums();
			while(index>0) {
			  _sums[index-1]+=me;
        index>>=1;
	    }
    }
	}

	inline index_type choose() const {
		index_type _nterms=nsums();
		index_type index=0;
		while(index<_nterms) {
			
			index_type indr=(index+1)<<1;
			index_type indl=indr-1;

			MatrixElement wr=Crop(_sums[indr]);
			MatrixElement wl=Crop(_sums[indl]);
			
			index=((wr+wl)*RNG::Uniform()<wl) ? indl : indr;
		}
		
		return index-_nterms;
		
	}


	/* Makes all elements zero */
	inline void reset() {
		for(std::vector<MatrixElement>::iterator it=_sums.begin();it!=_sums.end();++it)
			*it=MatrixElement(0);
	}
	
	/* Returns the sum of all the partial probabilities which is stored at the root */
	inline MatrixElement norm() const { return Crop(_sums[0]); }

};

MatrixElement TSum::Tolerance;

}

#endif