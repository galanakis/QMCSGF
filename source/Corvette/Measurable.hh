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

	class MeasureDefaults {
	public:
		typedef long double _accumulator_float;
		typedef unsigned long long _counter_int;
		typedef BinnedAccumulator<_accumulator_float> BinnedAccumulatorME;
		typedef std::map<unsigned int,_counter_int> BrokenHistogramType;

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
      
      const _accumulator_float Weight=OperatorString.BoltzmannWeight();

      _BrokenHistogram[OperatorString.NBrokenLines()]+=1;
      
      if(OperatorString.NBrokenLines()==0) {
				BoltzmannWeight+=Weight;
        _accumulator_float Kinetic=-static_cast<_accumulator_float>(OperatorString.length())/OperatorString.Beta();
        _accumulator_float Potential=OperatorString.Energy(RIGHT);        
        _Kinetic.push( Kinetic*Weight );
        _Potential.push( Potential*Weight );
        _TotalEnergy.push( (Kinetic+Potential)*Weight );
      }

		}		
				
		inline const BrokenHistogramType &BrokenHistogram() const {return _BrokenHistogram;}

		inline _accumulator_float BrokenNormalization() const {
			_accumulator_float Normalization=0;
      for(BrokenHistogramType::const_iterator it=_BrokenHistogram.begin();it!=_BrokenHistogram.end();++it)
        Normalization+=it->second;
			return Normalization;
		}

		inline _counter_int count(unsigned int i=0) { return _BrokenHistogram[i]; }
    inline const BinnedAccumulatorME &KineticEnergy() const {return _Kinetic;}
    inline const BinnedAccumulatorME &PotentialEnergy() const {return _Potential;}
    inline const BinnedAccumulatorME &TotalEnergy() const {return _TotalEnergy;}

	};


  class Measurable : public MeasureDefaults {
 
    std::vector<BinnedAccumulatorME> Sums;

	  struct MeasAccumulators {
	    const HamiltonianTerm *term;
	    BinnedAccumulatorME *sum;
	    MeasAccumulators(const HamiltonianTerm * const t,BinnedAccumulatorME * const s) : term(t), sum(s) {} 
	  };


		typedef std::vector< std::pair<Boson*,int> > BosonDeltaMapType; 
    std::map<BosonDeltaMapType,std::vector<MeasAccumulators> > _Acc;

		// This one is ok to be slow since it is called only in the initializer
		inline BosonDeltaMapType map(const HamiltonianTerm * const term) {
			std::map<Boson*,int> indices;
	    for(unsigned int i=0;i<term->product().size();++i) {
				int delta=-term->product()[i].delta();
				if(delta!=0) indices[term->product()[i].particle_id()]=delta;
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

	
  public:
    Measurable(const std::vector<Hamiltonian> &HList) : MeasureDefaults(), Sums(HList.size()) {

      for(unsigned int i=0;i<HList.size();++i)
      for(unsigned int j=0;j<HList[i].size();++j) {
        const HamiltonianTerm *term=&HList[i][j];
				_Acc[map(term)].push_back(MeasAccumulators(term,&Sums[i]));
      }

    }
    
    inline void flush() {
      
			for(std::vector<BinnedAccumulatorME>::iterator it=Sums.begin();it!=Sums.end();++it) 
				it->flush(BoltzmannWeight);

			MeasureDefaults::flush();

    }

    void measure(const OperatorStringType &OperatorString) {
			
      const _accumulator_float Weight=OperatorString.BoltzmannWeight();
			const std::vector<MeasAccumulators> &v=_Acc[map(OperatorString.ListBrokenLines())];
			for(std::vector<MeasAccumulators>::const_iterator it=v.begin();it!=v.end();++it) 
				it->sum->push( it->term->me(RIGHT) * Weight );

			MeasureDefaults::measure(OperatorString);

		}
		 
		typedef std::vector<BinnedAccumulatorME>::size_type size_type;
    inline size_type size() const {return Sums.size();}
    inline const BinnedAccumulatorME &operator[](unsigned long i) const {return Sums[i];}
                                                  
  };  

}
