#ifndef __TSUM__
#define __TSUM__

#include "RandomNumberGenerator.hh"
#include "Conventions.hh"
#include <vector>

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
choose(): randomly chose an index using its relative probability.
norm(): returns the sum of all relative probabilities

========== Usage =========
TSum t;
t.resize(100);          // Keep the relative probabilities of 100 indices.
t.update(10,25);        // set the relative probability of index 10 to be 25.
double p=t.element(10); // return the relative probability of index 10.
double norm=t.norm();   // the sum of all relative probabilities
int index=t.choose();   // pick an index at random using their relative probabilities.
t.reset();              // clear everything.

*/



class TSum {
	typedef std::vector<MatrixElement>::size_type index_type;

	std::vector<MatrixElement> _sums;
	MatrixElement *_elements;
	index_type _nsums,_base;
	_float_accumulator _norm;

	index_type _nterms;
	bool *_guard;
	index_type *_buffer,*_buffer_entry;


	inline void flush() {

		while(_buffer_entry>_buffer) {
			--_buffer_entry;

			_guard[*_buffer_entry]=true;
			index_type index=*_buffer_entry/2+_base;
			MatrixElement newme=_sums[2*index-1]+_sums[2*index];
			MatrixElement oldme=_sums[index-1];

			if(newme!=oldme) {
				while(index) {
					_sums[index-1] =  newme;
					newme += _sums[index-2*(index&1)];
					index/=2;
				}
			}



		}

		_norm=_sums[0];
	}

public:

	TSum() : _sums(), _elements(0), _guard(0) {}
	~TSum() {
		delete [] _guard;
		delete [] _buffer;
	}

	/* resizes the vector */
	inline void resize(index_type NTerms) {
		index_type _size=1;
		while(_size<NTerms) _size<<=1;
		_sums.resize(2*_size-1,MatrixElement(0));
		_nsums=_size-1;
		_base=(1+_nsums)/2;
		_elements=&_sums[_nsums];
		_nterms=NTerms;
		_guard=new bool[_nterms];
		_buffer=new index_type[_nterms];
		_buffer_entry=_buffer;
		_norm=0;
	}

	inline void update(index_type index,MatrixElement me) {
		if(me!=_elements[index]) {
			if(_guard[index]) {
				*_buffer_entry=index;
				++_buffer_entry;
				_guard[index]=false;
			}
			_norm+=me-_elements[index];
			_elements[index]=me;
		}
	}

	inline const MatrixElement &element(index_type index) const {return _elements[index];}

	inline index_type choose() {
		flush();

		MatrixElement w,wl;
		index_type index=0;
		while(index<_nsums) {
			w =_sums[index];
			index=2*index+1;
			wl=_sums[index];
			index+=!(w*RNG::Uniform()<wl);
		}

		return index-_nsums;

	}

	const _float_accumulator &norm() const {return _norm;}
	/* Makes all elements zero */
	inline void reset() {
		_buffer_entry=_buffer;
		_norm=0;
		for(index_type i=0;i<_nterms;++i)
			_guard[i]=true;
		for(std::vector<MatrixElement>::iterator it=_sums.begin();it!=_sums.end();++it)
			*it=MatrixElement(0);
	}

};


}

#endif
