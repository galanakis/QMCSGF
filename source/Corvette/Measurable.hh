#ifndef __MEASURABLE__
#define __MEASURABLE__

#include <map>
#include <vector>
#include <set>
#include <string>
#include <cmath>
#include <algorithm>

#include "HamiltonianTerm.hh" 
#include "Boson.hh"
#include "Accumulator.hh"
#include "OperatorString.hh"

namespace SGF {

/*
* class MeasureDefaults

* Measures some basic quantities, such as the Kinetic, Potential and TotalEnergy and the
* count of occurences for each broken line number. This is a very light measurement
* class which can be used during thermalization.

* ======= Usage ======
* MeasureDefaults md;
* md.measure(OperatorString); // Do a measurement.
* md.flush();                 // push into the bins.
* cout<<md.KineticEnergy()<<std::endl;
* cout<<md.PotentialEnergy()<<std::endl;
* cout<<md.TotalEnergy()<<std::endl;
* int count=md.count(i);      // how many times i lines were broken. 

*/

typedef BinnedAccumulator<_float_accumulator> BinnedAccumulatorME;


class MeasurableTerms  {

	struct TermBuffer {
		_float_accumulator buffer;
		HamiltonianTerm term;
		TermBuffer(const HamiltonianTerm &_term) : buffer(0), term(_term) {}
		TermBuffer(const TermBuffer &o) : buffer(o.buffer), term(o.term) {}
	};

	typedef BrokenLines::BosonDeltaMapType KeyType;

	inline KeyType get_key(const HamiltonianTerm &term) const {return BrokenLines::map(term.product()); }

	typedef std::multimap<KeyType,TermBuffer> multimap_type;
	typedef std::pair<multimap_type::iterator,multimap_type::iterator> equal_range_type;
	multimap_type _multimap;
public:
	MeasurableTerms() : _multimap() {}

	_float_accumulator* insert(const HamiltonianTerm &term) {

		const KeyType key=get_key(term);

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
#ifdef DEBUG
      if(! it->second.term.match()) {
				std::cerr<<"Operator does not match with the broken lines"<<std::endl;
				exit(333);
			}
#endif			
			
			it->second.buffer += it->second.term.me(RIGHT)*Weight; 
		}
	}

	inline void reset() {
		for(multimap_type::iterator it=_multimap.begin(); it!=_multimap.end(); ++it) 
			it->second.buffer=0;
	}

};







struct FunctionFunctor {
	typedef std::vector<_float_accumulator*> data_type;
	typedef std::vector<MatrixElement> coeff_type;
	virtual _float_accumulator operator()(const data_type &_data,coeff_type &_coefficients) =0;
};

struct PlusFunctor : public FunctionFunctor {
	_float_accumulator operator()(const data_type &_data,coeff_type &_coefficients) {
		_float_accumulator result=0;
		for(data_type::size_type i=0;i<_data.size();++i)
			result+= *_data[i] * _coefficients[i];
		return result;
	}
};

PlusFunctor PlusData;

class MeasurableFunction {

	std::vector<_float_accumulator*> _data;
	std::vector<MatrixElement> _coefficients;
	FunctionFunctor* _function;
public:
	MeasurableFunction() :  _data(), _function(&PlusData) {}
	MeasurableFunction(const MeasurableFunction &o) : _data(o._data), _function(o._function) {
		for(std::vector<_float_accumulator*>::size_type i=0;i<o._data.size();++i)
			_data.push_back(o._data[i]);
	}


	inline _float_accumulator evaluate() { return (*_function)(_data,_coefficients); }
	inline void push_back(_float_accumulator* acc,MatrixElement me) { 
		_data.push_back(acc);
		_coefficients.push_back(me);
	} 

	FunctionFunctor* &function() {return _function;}

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

	std::vector<std::string> _Tags;
	BinnedAccumulatorME _Kinetic;
	BinnedAccumulatorME _Potential;
	BinnedAccumulatorME _TotalEnergy;
	std::vector<BinnedAccumulatorME*> _Bins;

public:
	Measurable(const OperatorStringType &OS) : OperatorString(OS), BrokenLineTracer(OS.GetListBrokenLines()) { reset(); }

	~Measurable() {
		for(std::vector<MeasurableFunction*>::iterator it=_Meas_Ptr.begin(); it!=_Meas_Ptr.end(); ++it)
			delete *it;
		for(std::vector<BinnedAccumulatorME*>::iterator it=_Bins.begin(); it!=_Bins.end(); ++it)
			delete *it;
	}

	void insert(const std::vector<std::string> &taglist,const std::vector<Hamiltonian> &OperatorList) {
		for(std::vector<Hamiltonian>::size_type i=0;i<OperatorList.size();++i) 
			insert_operator(taglist[i],OperatorList[i]);
	}

	void insert_bin(const std::string &tag) {
		BinnedAccumulatorME *_bin_ptr=new BinnedAccumulatorME;
		_Bins.push_back(_bin_ptr);
		_Tags.push_back(tag);		
	}

	void insert_operator(const std::string &tag,const Hamiltonian &Operator) {
    
		insert_bin(tag);

		MeasurableFunction *_meas_ptr=new MeasurableFunction;
		_Meas_Ptr.push_back(_meas_ptr);
   	
		for(Hamiltonian::const_iterator term_ptr=Operator.begin(); term_ptr!=Operator.end(); ++term_ptr) {
			_float_accumulator *buffer=_TermBuffers.insert(*term_ptr);
			_meas_ptr->push_back( buffer, term_ptr->coefficient() );
		}

	}

	void insert_function(const std::string &tag,FunctionFunctor * functor,const Hamiltonian &Operator) {
   
		insert_bin(tag);
		
		MeasurableFunction *_meas_ptr=new MeasurableFunction;
		_Meas_Ptr.push_back(_meas_ptr);		
		_meas_ptr->function()=functor;

		for(Hamiltonian::const_iterator term_ptr=Operator.begin(); term_ptr!=Operator.end(); ++term_ptr) {
			_float_accumulator *buffer=_TermBuffers.insert(*term_ptr);
			_meas_ptr->push_back( buffer, term_ptr->coefficient() );
		}		
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
