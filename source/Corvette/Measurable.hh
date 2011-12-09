#include <map>
#include <vector>
#include <set>
#include <string>
#include <cmath>

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

	class MeasureDefaults {
	public:
		typedef BinnedAccumulator<_float_accumulator> BinnedAccumulatorME;
		typedef std::map<unsigned int,_integer_counter> BrokenHistogramType;

	private:
    BinnedAccumulatorME _Kinetic;
    BinnedAccumulatorME _Potential;
    BinnedAccumulatorME _TotalEnergy;

		BrokenHistogramType _BrokenHistogram;
  
 	protected:
		MatrixElement BoltzmannWeight;
	public:
		MeasureDefaults() : _Kinetic(), _Potential(), _TotalEnergy(), BoltzmannWeight(0) {}

    inline void flush() {
      _Kinetic.flush(BoltzmannWeight);
      _Potential.flush(BoltzmannWeight);
      _TotalEnergy.flush(BoltzmannWeight);
			BoltzmannWeight=0;
		}
		
    void measure(const OperatorStringType &OperatorString) {
      
      const _float_accumulator Weight=OperatorString.BoltzmannWeight();

      _BrokenHistogram[OperatorString.NBrokenLines()]+=1;
      
      if(OperatorString.NBrokenLines()==0) {
				BoltzmannWeight+=Weight;
        _float_accumulator Kinetic=-static_cast<_float_accumulator>(OperatorString.length())/OperatorString.Beta();
        _float_accumulator Potential=OperatorString.Energy(RIGHT);        
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

	};

 


  class Measurable : public MeasureDefaults {


		typedef std::vector< std::pair<Boson*,int> > BosonDeltaMapType;

			// This one is ok to be slow since it is called only in the initializer
			inline BosonDeltaMapType map(const HamiltonianTerm &term) {
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


			// This one must be very fast
			inline BosonDeltaMapType map(const std::set<Boson*> &s) {
				BosonDeltaMapType vmap;
				vmap.reserve(s.size());
				std::set<Boson*>::const_iterator it;
				for(it=s.begin();it!=s.end();++it) {
					int delta=(*it)->delta();
					if(delta!=0) vmap.push_back(std::pair<Boson*,int>(*it,delta));
				}
				return vmap;
			}	

		  struct MeasurableData {
				MatrixElement coefficient;
		    BinnedAccumulatorME *sum;
				MeasurableData(const MeasurableData &o) : coefficient(o.coefficient), sum(o.sum) {}
		    MeasurableData(MatrixElement c,BinnedAccumulatorME * const s) : coefficient(c), sum(s) {}
		  };

		  class TermBuffer {
				typedef std::vector<MeasurableData> MeasurableDataVectorType;
				MeasurableDataVectorType _measurable_data_vector;
				_float_accumulator buffer;
			public:

		    TermBuffer() : buffer(0), _measurable_data_vector() {}
				TermBuffer(const TermBuffer &o) : buffer(o.buffer), _measurable_data_vector(o._measurable_data_vector) {}
				inline void push_me(_float_accumulator value) { buffer += value; }

				inline void push_buffers() {
					for(MeasurableDataVectorType::iterator it=_measurable_data_vector.begin(); it!=_measurable_data_vector.end(); ++it) 
						it->sum->push(buffer * it->coefficient);
					buffer=0;
				}

				void push_back(MatrixElement c,BinnedAccumulatorME * const s) {
					_measurable_data_vector.push_back(MeasurableData(c,s));
				}

		  };


    std::vector<BinnedAccumulatorME> Sums;

		typedef std::map<const HamiltonianTerm,TermBuffer> MapTermBufferType;
		MapTermBufferType _MapTermBuffer;

		typedef std::vector<MapTermBufferType::iterator> TermBufferVectorType; 
		typedef std::map<BosonDeltaMapType,TermBufferVectorType> AccType; 
    AccType _Acc;

	
  public:
    Measurable(const std::vector<Hamiltonian> &HList) : MeasureDefaults(), Sums(HList.size()) {

			for(unsigned int i=0;i<HList.size();++i)
      for(unsigned int j=0;j<HList[i].size();++j) {
				const HamiltonianTerm &term=HList[i][j];
				if(term.constant()) 
					Sums[i].constant()+=term.coefficient();
				else {
					_MapTermBuffer[term.product()].push_back(term.coefficient(),&Sums[i]);
					_Acc[map(term.product())].push_back(_MapTermBuffer.find(term.product()));
				}
      }
      
    }

		inline void flush(MatrixElement BoltzmannW) {
			for(MapTermBufferType::iterator it=_MapTermBuffer.begin(); it!=_MapTermBuffer.end(); ++it) 
				it->second.push_buffers();

			for(std::vector<BinnedAccumulatorME>::iterator it=Sums.begin();it!=Sums.end();++it) 
				it->flush(BoltzmannW);
			
		}

   	inline void measure(const std::set<Boson*> &ListBrokenLines,double Weight) {

			AccType::iterator v_it=_Acc.find(map(ListBrokenLines));
			if(v_it!=_Acc.end())
				for(TermBufferVectorType::iterator it=v_it->second.begin();it!=v_it->second.end();++it) {
					MatrixElement me=(*it)->first.me(RIGHT);
					(*it)->second.push_me( Weight * me );
				}
	
		}


    inline void measure(const OperatorStringType &OperatorString) {

			measure(OperatorString.ListBrokenLines(),OperatorString.BoltzmannWeight());
			MeasureDefaults::measure(OperatorString);

		}

    inline void flush() {
      
			flush(BoltzmannWeight);
			MeasureDefaults::flush();

    }

		 
		typedef std::vector<BinnedAccumulatorME>::size_type size_type;
    inline size_type size() const {return Sums.size();}
    inline const BinnedAccumulatorME &operator[](std::vector<BinnedAccumulatorME>::size_type i) const {return Sums[i];}
                                                  
  };  

}
