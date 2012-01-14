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
* std::cout<<md.KineticEnergy()<<std::endl;
* std::cout<<md.PotentialEnergy()<<std::endl;
* std::cout<<md.TotalEnergy()<<std::endl;
* int count=md.count(i);      // how many times i lines were broken. 

*/

typedef BinnedAccumulator<_float_accumulator> BinnedAccumulatorME;


class MeasurableTerms  {

	class TermBuffer {
		_float_accumulator _buffer;
		HamiltonianTerm _term;
	public:

		TermBuffer(const HamiltonianTerm &term) : _buffer(0), _term(term) {}
		TermBuffer(const TermBuffer &o) : _buffer(o._buffer), _term(o._term) {}
		inline _float_accumulator &buffer() {return _buffer; }
		inline const HamiltonianTerm &term() {return _term; }

	};

	typedef BrokenLines::BosonDeltaMapType KeyType;
	
	inline KeyType get_key(const HamiltonianTerm &term) const {return BrokenLines::map(term.product()); }
	
	typedef std::multimap<KeyType,TermBuffer> multimap_type;
	typedef pair<multimap_type::iterator,multimap_type::iterator> equal_range_type;
	multimap_type _multimap;
public:
	MeasurableTerms() : _multimap() {}

	_float_accumulator* insert(const HamiltonianTerm &term) {
    
		const KeyType key=get_key(term);

		equal_range_type equal_range=_multimap.equal_range(key);
		multimap_type::iterator it=equal_range.first;
		while(it!=equal_range.second && it->second.term().product() != term.product() ) 
			++it;

		if(it==equal_range.second) {
			it=_multimap.insert(pair<KeyType,TermBuffer>( key, TermBuffer(term.product()) ) );
		}

		return &(it->second.buffer());

	}

	inline void measure(const KeyType &key,double Weight) {
		equal_range_type equal_range=_multimap.equal_range(key);
		for(multimap_type::iterator it=equal_range.first; it!=equal_range.second; ++it) 
			it->second.buffer() += it->second.term().me(RIGHT)*Weight;
	}

	inline void reset() {
		for(multimap_type::iterator it=_multimap.begin(); it!=_multimap.end(); ++it) 
			it->second.buffer()=0;
	}

};



class MeasureDefaults {
public:
	typedef std::map<unsigned int,_integer_counter> BrokenHistogramType;
private:
	BinnedAccumulatorME _Kinetic;
	BinnedAccumulatorME _Potential;
	BinnedAccumulatorME _TotalEnergy;

	BrokenHistogramType _BrokenHistogram;

protected:
	const OperatorStringType &OperatorString;
	MatrixElement BoltzmannWeight;
public:
	MeasureDefaults(const OperatorStringType &OS) : _Kinetic(), _Potential(), _TotalEnergy(), BoltzmannWeight(0), OperatorString(OS) {}

	inline void flush() {
		_Kinetic.flush(BoltzmannWeight);
		_Potential.flush(BoltzmannWeight);
		_TotalEnergy.flush(BoltzmannWeight);
		BoltzmannWeight=0;
	}

	void measure() {


		_BrokenHistogram[OperatorString.NBrokenLines()]+=1;

		if(OperatorString.NBrokenLines()==0) {
			const _float_accumulator Weight=OperatorString.BoltzmannWeight();
			BoltzmannWeight+=Weight;
			_float_accumulator Kinetic=-static_cast<_float_accumulator>(OperatorString.length())/OperatorString.Beta();
			_float_accumulator Potential=OperatorString.DiagonalEnergy();        
			_Kinetic.push( Kinetic*Weight );
			_Potential.push( Potential*Weight );
			_TotalEnergy.push( (Kinetic+Potential)*Weight );
		}

	}		

	inline const BrokenHistogramType &BrokenHistogram() const {return _BrokenHistogram;}

	inline _float_accumulator BrokenNormalization() const {
		_float_accumulator Normalization(0);
		for(BrokenHistogramType::const_iterator it=_BrokenHistogram.begin();it!=_BrokenHistogram.end();++it)
			Normalization+=it->second;
		return Normalization;
	}

	inline _integer_counter count(unsigned int i=0) { return _BrokenHistogram[i]; }
	inline const BinnedAccumulatorME &KineticEnergy() const {return _Kinetic;}
	inline const BinnedAccumulatorME &PotentialEnergy() const {return _Potential;}
	inline const BinnedAccumulatorME &TotalEnergy() const {return _TotalEnergy;}

	void print_defaults() {

		std::cout<<::std::endl;
		std::cout<<"  *******************************\n";
		std::cout<<"  * Broken worldlines histogram *\n";
		std::cout<<"  *******************************\n\n";
		std::cout<<"    N lines\tCount\tProbability\n\n";


		double Normalization=BrokenNormalization();
		for(BrokenHistogramType::const_iterator it=_BrokenHistogram.begin();it!=_BrokenHistogram.end();++it)
			std::cout<<"    "<<it->first<<"\t\t"<<it->second<<"\t"<<it->second/Normalization<<std::endl;
		std::cout << endl;

		std::cout << "  ***********************************************************************************\n";
		std::cout << "  * Energies (obtained from operator string length and Green operator state energy) *\n";
		std::cout << "  ***********************************************************************************\n\n";
		std::cout << "    Total energy: " << TotalEnergy() << "\n";
		std::cout << "    Diagonal energy: " << PotentialEnergy() << "\n";           
		std::cout << "    Non-diagonal energy: " << KineticEnergy() << "\n\n";

	}

};





struct FunctionFunctor {
	typedef std::vector<_float_accumulator*> data_type;
	typedef std::vector<MatrixElement> coeff_type;
	virtual _float_accumulator operator()(const data_type &_data,coeff_type &_coefficients) =0;
};

struct PlusFunctor : public FunctionFunctor {
	_float_accumulator operator()(const data_type &_data,coeff_type &_coefficients) {
		_float_accumulator _buffer=0;
		for(int i=0;i<_data.size();++i)
			_buffer+= *_data[i] * _coefficients[i];
		return _buffer;
	}
};

PlusFunctor PlusData;

class MeasurableFunction {
	
	BinnedAccumulatorME bin_acc;
	std::vector<_float_accumulator*> _data;
	std::vector<MatrixElement> _coefficients;
	FunctionFunctor* _function;
public:
	MeasurableFunction() : bin_acc(), _function(&PlusData), _data() {}
	MeasurableFunction(const MeasurableFunction &o) : bin_acc(o.bin_acc), _function(o._function), _data(o._data) {
		for(std::vector<_float_accumulator*>::size_type i=0;i<o._data.size();++i)
			_data.push_back(o._data[i]);
	}
 
	inline void flush(_float_accumulator BoltzmannWeight) {
		bin_acc.push( evaluate() );
		bin_acc.flush( BoltzmannWeight );
	}

	const BinnedAccumulatorME & bin() const { return bin_acc; } 

	inline _float_accumulator evaluate() { return (*_function)(_data,_coefficients); }
	inline void push_back(_float_accumulator* acc,MatrixElement me) { 
		_data.push_back(acc);
		_coefficients.push_back(me);
	} 

	FunctionFunctor* &function() {return _function;}

};



class Measurable : public MeasureDefaults {

	std::vector<std::string> _Tags;
	std::vector<MeasurableFunction*> _Meas_Ptr;
	
	MeasurableTerms _TermBuffers;

	BrokenLines BrokenLineTracer;

public:
	void insert(const std::vector<std::string> &taglist,const std::vector<Hamiltonian> &OperatorList) {
		for(std::vector<Hamiltonian>::size_type i=0;i<OperatorList.size();++i) 
			insert_operator(taglist[i],OperatorList[i]);
	}

	void insert_operator(const std::string &tag,const Hamiltonian &Operator) {

		MeasurableFunction *_meas_ptr=new MeasurableFunction;
		_Meas_Ptr.push_back(_meas_ptr);
		_Tags.push_back(tag);

		for(Hamiltonian::const_iterator term_ptr=Operator.begin(); term_ptr!=Operator.end(); ++term_ptr) {
			_float_accumulator *buffer=_TermBuffers.insert(*term_ptr);
			_meas_ptr->push_back( buffer, term_ptr->coefficient() );
		}

	}

	void insert_function(const std::string &tag,FunctionFunctor * functor,const Hamiltonian &Operator) {

		MeasurableFunction *_meas_ptr=new MeasurableFunction;
		_Meas_Ptr.push_back(_meas_ptr);
		_Tags.push_back(tag);
		_meas_ptr->function()=functor;

		for(Hamiltonian::const_iterator term_ptr=Operator.begin(); term_ptr!=Operator.end(); ++term_ptr) {
			_float_accumulator *buffer=_TermBuffers.insert(*term_ptr);
			_meas_ptr->push_back( buffer, term_ptr->coefficient() );
		}		
	}


	inline void flush(MatrixElement BoltzmannW) {
		for(std::vector<MeasurableFunction*>::iterator it=_Meas_Ptr.begin(); it!=_Meas_Ptr.end(); ++it)
			(*it)->flush(BoltzmannW);
		_TermBuffers.reset();
	}


	inline void measure() {

		_TermBuffers.measure(BrokenLineTracer(),OperatorString.BoltzmannWeight());
		MeasureDefaults::measure();

	}

	inline void flush() {

		flush(BoltzmannWeight);
		MeasureDefaults::flush();

	}



	Measurable(const OperatorStringType &OS) : MeasureDefaults(OS), BrokenLineTracer(OS.GetListBrokenLines()) {}
	~Measurable() {
		for(std::vector<MeasurableFunction*>::iterator it=_Meas_Ptr.begin(); it!=_Meas_Ptr.end(); ++it)
			delete *it;
	}

	void print() {

		print_defaults();

		std::cout << "  ******************************\n";
		std::cout << "  * User's defined measurables *\n";
		std::cout << "  ******************************\n\n";

		for(std::vector<MeasurableFunction*>::size_type i=0;i<_Meas_Ptr.size();++i) 
			std::cout<<"    "<<_Tags[i]<<": "<<(_Meas_Ptr[i]->bin())<<std::endl;

	}

};  

}
