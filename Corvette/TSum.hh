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
	
	The key functions are:
	update(index,me): change the probability of the index by "me".
	choose(): randomly chose an index using its relative probability (sum-sum_Left-sum_Right)
	norm(): returns the sum of all relative probabilities
	
*/

class TSum {
	std::vector<MatrixElement> _sums;
	long double _head; // a copy of the _sums[0] is stored here. The redundancy is used for error tracking.
	
	struct _buffer_data {
		uint index;
		MatrixElement me;
		_buffer_data(uint i,MatrixElement m) : index(i),me(m) {}
	};
	
	struct _buffer_comparison { inline bool operator()(_buffer_data a1,_buffer_data a2) {return a1.index<a2.index;} };
	
	std::priority_queue<_buffer_data,std::vector<_buffer_data>,_buffer_comparison> _buffer;
	
	// Chop off small numerical values. The thresshold is Tolerance which is set as the minimum coefficient of any kinetic operator.
	long double Crop(long double x) const { return (fabs(x)>Tolerance)?x:0; }

public:
  static MatrixElement Tolerance;

	TSum() : _buffer(), _sums(), _head(0) {}
	//TSum(uint N) : _sums(N,MatrixElement(0)),_head(0),_buffer() {std::cout<<"Tsum init"<<std::endl;}
	//TSum(const TSum &o) : _sums(o._sums),_head(0),_buffer(o._buffer) {std::cout<<"Tsum copy"<<std::endl;}
	~TSum() {}
	
	/* resizes the vector */
	inline void resize(uint N) {_sums.resize(N,MatrixElement(0));}

	
	inline void push(uint index,MatrixElement me) { 
		if(me!=0) {
			_head+=me;
			_buffer.push(_buffer_data(index,me));
		}
	}
	
	
	
	
	inline void flush() {
		
		//std::cout<<"Begin flush"<<_buffer.size()<<std::endl;
		
		while(!_buffer.empty()) {
			_buffer_data data=_buffer.top();
			_buffer.pop();

			while(!_buffer.empty() && _buffer.top().index==data.index) {
				data.me+=_buffer.top().me;
				_buffer.pop();
			}
			
			_sums[data.index]+=data.me;

			if(data.index!=0) {
			 data.index=((data.index+1)>>1)-1;
			 _buffer.push(data);
		  }

		}
		
	}
	
	
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


	inline uint choose() {
		
		flush();
		
		int _nterms=_sums.size();
		int index=0;
		while(index<_nterms) {
			int indr=(index+1)<<1;
			int indl=indr-1;

			MatrixElement w =Crop(_sums[index]);
			MatrixElement wr=(indr<_nterms) ? Crop(_sums[indr]) : MatrixElement(0);
			MatrixElement wl=(indl<_nterms) ? Crop(_sums[indl]) : MatrixElement(0);

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
		flush();
		_head=0;
		for(uint i=0;i<_sums.size();++i)
			_sums[i]=MatrixElement(0);	
	}
	
	/* Returns the sum of all the partial probabilities which is stored at the root */
	inline MatrixElement norm() const { return Crop(_head); }

};

MatrixElement TSum::Tolerance;

}

#endif