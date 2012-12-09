#ifndef __MEASURABLE__
#define __MEASURABLE__

#include <map>
#include <vector>
#include <set>
#include <string>
#include <cmath>
#include <algorithm>

#include "HamiltonianTerm.hh" 
#include "Accumulator.hh"
#include "OperatorString.hh"

namespace SGF {


typedef BinnedAccumulator<_float_accumulator> BinnedAccumulatorME;


class MeasurableTerms  {

	struct TermBuffer {
		_float_accumulator buffer;
		HamiltonianTerm term;
		TermBuffer(const HamiltonianTerm &_term) : buffer(0), term(_term) {}
		TermBuffer(const TermBuffer &o) : buffer(o.buffer), term(o.term) {}
	};

	typedef BrokenLines::BosonDeltaMapType KeyType;

	typedef std::multimap<KeyType,TermBuffer> multimap_type;
	typedef std::pair<multimap_type::iterator,multimap_type::iterator> equal_range_type;
	multimap_type _multimap;
public:
	MeasurableTerms() : _multimap() {}

	_float_accumulator* insert(const HamiltonianTerm &term) {

		const KeyType key=BrokenLines::map(term.product());

		equal_range_type equal_range=_multimap.equal_range(key);
		multimap_type::iterator it=equal_range.first;
		while(it!=equal_range.second && it->second.term.product() != term.product() ) 
			++it;

		if(it==equal_range.second) {
			it=_multimap.insert(std::pair<KeyType,TermBuffer>( key, TermBuffer(term.product()) ) );
		}

		return &(it->second.buffer);

	}

	inline void measure(const KeyType &key,double Weight) {
		equal_range_type equal_range=_multimap.equal_range(key);
		for(multimap_type::iterator it=equal_range.first; it!=equal_range.second; ++it) {
			it->second.buffer += it->second.term.me(RIGHT)*Weight; 
		}
	}

	inline void reset() {
		for(multimap_type::iterator it=_multimap.begin(); it!=_multimap.end(); ++it) 
			it->second.buffer=0;
	}

};


class MeasurableFunction {

public:
	MeasurableFunction() {}
	virtual inline _float_accumulator evaluate() const = 0;
	virtual ~MeasurableFunction() {}
};

class MeasurableSum : public MeasurableFunction {
	typedef std::vector<std::pair<_float_accumulator*,MatrixElement> > vector_pair_t;
	vector_pair_t data;
public:
	MeasurableSum() : MeasurableFunction() {}
	inline _float_accumulator evaluate() const {
		_float_accumulator result=0;
		for(vector_pair_t::const_iterator it=data.begin(); it!=data.end(); ++it)
			result += *it->first * it->second;
		return result;
	}
	void push_back(_float_accumulator *a,MatrixElement m) {
		data.push_back(std::pair<_float_accumulator*,MatrixElement>(a,m));
	}

};

class MeasurableNumber : public MeasurableFunction {
	_float_accumulator* data;
	MatrixElement coefficient;
public:
	MeasurableNumber(_float_accumulator* a,MatrixElement c) : data(a), coefficient(c) {}
	inline _float_accumulator evaluate() const {return *data * coefficient; }	
};

/*
	Notes for MPI parallelization.
	
	The BinnedAccumulatorME variables should only be defined in the Root node.
 	Also the print and flush methods should only be called by the root.
	During the flush operation the mpi reduce should be called for every bin.
	
*/

class Measurable  {
	const OperatorStringType &OperatorString;

	BrokenLines BrokenLineTracer;

	_float_accumulator buffer_BoltzmannWeight;
	_float_accumulator buffer_kinetic;
	_float_accumulator buffer_potential;
	MeasurableTerms _TermBuffers;

	std::vector<MeasurableFunction*> _Meas_Ptr;

	BinnedAccumulatorME _Kinetic;
	BinnedAccumulatorME _Potential;
	BinnedAccumulatorME _TotalEnergy;
	std::vector<std::string> _Tags;
	std::vector<BinnedAccumulatorME*> _Bins;

public:
	Measurable(OperatorStringType &OS) : OperatorString(OS), BrokenLineTracer(OS.GetListBrokenLines(),OS) { reset(); }

	~Measurable() {
		for(std::vector<MeasurableFunction*>::iterator it=_Meas_Ptr.begin(); it!=_Meas_Ptr.end(); ++it)
			delete *it;
		for(std::vector<BinnedAccumulatorME*>::iterator it=_Bins.begin(); it!=_Bins.end(); ++it)
			delete *it;
	}

	void insert_bin(const std::string &tag) {
		BinnedAccumulatorME *_bin_ptr=new BinnedAccumulatorME;
		_Bins.push_back(_bin_ptr);
		_Tags.push_back(tag);		
	}

	void insert(const std::string &tag,const Hamiltonian &Operator) {
		MeasurableSum *_meas_ptr=new MeasurableSum();
		for(Hamiltonian::const_iterator term_ptr=Operator.begin(); term_ptr!=Operator.end(); ++term_ptr) {
			_meas_ptr->push_back( buffer_pointer(*term_ptr), term_ptr->coefficient() );
		}
		insert_measurable(tag,_meas_ptr);
	}

	void insert(const std::string &tag,const HamiltonianTerm &term) {
		MeasurableNumber *_meas_ptr=new MeasurableNumber(buffer_pointer(term), term.coefficient());
		insert_measurable(tag,_meas_ptr);
	}

	_float_accumulator *buffer_pointer(const HamiltonianTerm &term) {
		return _TermBuffers.insert(term);
	}

	void insert_measurable(const std::string &tag,MeasurableFunction *_meas_ptr) {
		insert_bin(tag);
		_Meas_Ptr.push_back(_meas_ptr);
	}


	inline void measure() {

		_TermBuffers.measure(BrokenLineTracer(),OperatorString.BoltzmannWeight());

		if(OperatorString.NBrokenLines()==0) {
			const _float_accumulator Weight=OperatorString.BoltzmannWeight();
			buffer_BoltzmannWeight+=Weight;
			_float_accumulator Kinetic=-static_cast<_float_accumulator>(OperatorString.length())/OperatorString.Beta();
			_float_accumulator Potential=OperatorString.DiagonalEnergy();
			buffer_kinetic+=Kinetic*Weight;
			buffer_potential+=Potential*Weight;       
		}

	}

	std::vector<_float_accumulator> buffers() {
		std::vector<_float_accumulator> result;
		for(std::vector<MeasurableFunction*>::size_type i=0; i<_Meas_Ptr.size(); ++i) 
			result.push_back(_Meas_Ptr[i]->evaluate());
		return result;		
	}
  
	inline void flush() {

		std::vector<_float_accumulator> buffer_functions;
		for(std::vector<MeasurableFunction*>::size_type i=0; i<_Meas_Ptr.size(); ++i) 
			buffer_functions.push_back(_Meas_Ptr[i]->evaluate());

#ifdef USEMPI

		_float_accumulator send_buffer_kinetic=buffer_kinetic;
		_float_accumulator send_buffer_potential=buffer_potential;
		_float_accumulator send_buffer_BoltzmannWeight=buffer_BoltzmannWeight;
		std::vector<_float_accumulator> send_buffer_functions=buffer_functions;

		MPI_Reduce(&send_buffer_kinetic,&buffer_kinetic,1,MPI_LONG_DOUBLE,MPI_SUM,Master,MPI_COMM_WORLD);
		MPI_Reduce(&send_buffer_potential,&buffer_potential,1,MPI_LONG_DOUBLE,MPI_SUM,Master,MPI_COMM_WORLD);
		MPI_Reduce(&send_buffer_BoltzmannWeight,&buffer_BoltzmannWeight,1,MPI_LONG_DOUBLE,MPI_SUM,Master,MPI_COMM_WORLD);
		MPI_Reduce(&send_buffer_functions[0],&buffer_functions[0],buffer_functions.size(),MPI_LONG_DOUBLE,MPI_SUM,Master,MPI_COMM_WORLD);
		
	 	if(Rank==Master) {
			for(std::vector<MeasurableFunction*>::size_type i=0; i<_Meas_Ptr.size(); ++i) 
				_Bins[i]->push(buffer_functions[i]/buffer_BoltzmannWeight);

			_Kinetic.push( buffer_kinetic/buffer_BoltzmannWeight );
			_Potential.push( buffer_potential/buffer_BoltzmannWeight );
			_TotalEnergy.push( (buffer_kinetic+buffer_potential)/buffer_BoltzmannWeight );		
		
		} 
		 

	 
#else
		
		for(std::vector<MeasurableFunction*>::size_type i=0; i<_Meas_Ptr.size(); ++i) 
			_Bins[i]->push(buffer_functions[i]/buffer_BoltzmannWeight);

		_Kinetic.push( buffer_kinetic/buffer_BoltzmannWeight );
		_Potential.push( buffer_potential/buffer_BoltzmannWeight );
		_TotalEnergy.push( (buffer_kinetic+buffer_potential)/buffer_BoltzmannWeight );		
		

#endif
   

		reset();
	}
  
	inline void reset() {
		_TermBuffers.reset();
		buffer_BoltzmannWeight=0;
		buffer_kinetic=0;
		buffer_potential=0;		
	}
  

	std::ostream& print(std::ostream &o) const {


		o << "  ***********************************************************************************\n";
		o << "  * Energies (obtained from operator string length and Green operator state energy) *\n";
		o << "  ***********************************************************************************\n\n";
		o << "    Total energy: " << _TotalEnergy << "\n";
		o << "    Diagonal energy: " << _Potential << "\n";           
		o << "    Non-diagonal energy: " << _Kinetic << "\n\n";
		o << "  ******************************\n";
		o << "  * User's defined measurables *\n";
		o << "  ******************************\n\n";

		for(std::vector<MeasurableFunction*>::size_type i=0;i<_Meas_Ptr.size();++i) 
			o<<"    "<<_Tags[i]<<": "<<*_Bins[i]<<std::endl;

		return o;
	}

};  

inline std::ostream& operator<<(std::ostream& o,const Measurable &measurables) {
	return measurables.print(o);
}

}


#endif
