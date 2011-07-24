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

	typedef long double _accumulator_float;
	typedef unsigned long long _counter_int;

	typedef BinnedAccumulator<_accumulator_float> BinnedAccumulatorME;

  struct Measurement {
    _accumulator_float value;
    _accumulator_float error;
    Measurement(_accumulator_float v,_accumulator_float e) : value(v), error(e) {}
  };

  inline std::ostream& operator<<(std::ostream& output, const Measurement &o) { return output<<o.value<<" +/- "<<o.error; }

   
  class BosonDeltaMap {
    std::map<Boson*,int> indices;
  public:
    BosonDeltaMap(const HamiltonianTerm * const term) {
      for(int i=0;i<term->product().size();++i)
        indices[term->product()[i].particle_id()]=-term->product()[i].delta();
    }
    BosonDeltaMap(const std::set<Boson*> &s) {
      std::set<Boson*>::const_iterator it;
      for(it=s.begin();it!=s.end();++it)
        indices[*it]=(*it)->delta();
    }
    const std::map<Boson*,int> &map() const {return indices;}
  };


  class MeasAccumulators {
    const HamiltonianTerm *term;
    BinnedAccumulatorME *sum;
  public:
    MeasAccumulators() {};
    MeasAccumulators(const HamiltonianTerm * const t,BinnedAccumulatorME * const s) : term(t), sum(s) {} 
    MeasAccumulators(const MeasAccumulators &o) : term(o.term), sum(o.sum) {}
    MeasAccumulators &operator=(MeasAccumulators &o) {
      term=o.term;
      sum=o.sum; 
      return *this;
    }
    inline void measure(_accumulator_float Weight) const { sum->push( term->me(RIGHT) ,Weight); }
  };


  class Measurable {
    _counter_int _ndiagonal;
    _counter_int _nflush;
    std::vector<BinnedAccumulatorME> Sums;
    std::vector<MatrixElement> Constants;

    const unsigned long _size;
    std::vector<MeasAccumulators> Diag_Acc;
    std::map<std::map<Boson*,int>,std::vector<MeasAccumulators> > OD_Acc;
    BinnedAccumulatorME _Kinetic;
    BinnedAccumulatorME _Potential;
    BinnedAccumulatorME _TotalEnergy;
    _accumulator_float _BoltzmannWeight;
    
    std::map<unsigned int,_counter_int> BrokenHistogram;

  public:
    Measurable(const std::vector<Hamiltonian> &HList) : _ndiagonal(0), _nflush(0), Sums(HList.size()), Constants(HList.size()), _size(HList.size()), _Kinetic(), _TotalEnergy(), _BoltzmannWeight(0) {

      for(unsigned int i=0;i<HList.size();++i)
      for(unsigned int j=0;j<HList[i].size();++j) {
        const HamiltonianTerm *term=&HList[i][j];

        if(term->constant())
          Constants[i]+=term->coefficient(); 
        else if(term->diagonal())
          Diag_Acc.push_back(MeasAccumulators(term,&Sums[i]));
        else
          OD_Acc[BosonDeltaMap(term).map()].push_back(MeasAccumulators(term,&Sums[i]));
      }
    }
                              
    inline void reset() {
      _ndiagonal=0;
      _nflush=0;
      _Kinetic.reset();
      _Potential.reset(); 
      _TotalEnergy.reset();

			std::vector<BinnedAccumulatorME>::iterator it;
			for(it=Sums.begin();it!=Sums.end();++it) it->reset();
			
      _BoltzmannWeight=0;
      BrokenHistogram.clear();
    }
    
    inline void flush() {
      ++_nflush;
      _Kinetic.flush(_BoltzmannWeight);
      _Potential.flush(_BoltzmannWeight);
      _TotalEnergy.flush(_BoltzmannWeight);

			std::vector<BinnedAccumulatorME>::iterator it;
			for(it=Sums.begin();it!=Sums.end();++it) it->flush(_BoltzmannWeight);
      _BoltzmannWeight=0;
    }

    void measure(const OperatorStringType &OperatorString) {
      
      const std::set<Boson*> &list=OperatorString.ListBrokenLines();
      const _accumulator_float Weight=OperatorString.BoltzmannWeight();

      BrokenHistogram[OperatorString.NBrokenLines()]+=1;
      
      if(OperatorString.NBrokenLines()!=0) {  // Off diagonal configurations
        std::map<std::map<Boson*,int>,std::vector<MeasAccumulators> >::const_iterator it(OD_Acc.find(BosonDeltaMap(list).map()));
        if(it!=OD_Acc.end()) {
          const std::vector<MeasAccumulators> &v=it->second;
					std::vector<MeasAccumulators>::const_iterator v_it;
					for(v_it=v.begin();v_it!=v.end();++v_it) v_it->measure(Weight);
        }
      }
      else {  // diagonal configuration
        ++_ndiagonal;
        _BoltzmannWeight+=Weight;
        _accumulator_float MeasKinetic=-static_cast<_accumulator_float>(OperatorString.length())/OperatorString.Beta();
        _accumulator_float MeasPotential=OperatorString.Energy(RIGHT);        
        _Kinetic.push(MeasKinetic,Weight);
        _Potential.push(MeasPotential,Weight);
        _TotalEnergy.push(MeasKinetic+MeasPotential,Weight);
				
				std::vector<MeasAccumulators>::const_iterator it;
				for(it=Diag_Acc.begin();it!=Diag_Acc.end();++it) it->measure(Weight);
      }

    }

    inline unsigned long size() const {return _size;}
    inline _counter_int count() const {return _ndiagonal;} 
    inline const Measurement operator[](unsigned long i) const {return Measurement(Constants[i]+Sums[i].average(),Sums[i].sigma());}
    inline const Measurement KineticEnergy() const {return Measurement(_Kinetic.average(),_Kinetic.sigma());}
    inline const Measurement PotentialEnergy() const {return Measurement(_Potential.average(),_Potential.sigma());}
    inline const Measurement TotalEnergy() const {return Measurement(_TotalEnergy.average(),_TotalEnergy.sigma());}
    
    void print_histogram() const {
      std::cout<<"  *******************************\n";
      std::cout<<"  * Broken worldlines histogram *\n";
      std::cout<<"  *******************************\n\n";
      std::cout<<"    N lines\tCount\tProbability\n\n";

      std::map<unsigned int,_counter_int>::const_iterator it;
      _accumulator_float Normalization=0;
      for(it=BrokenHistogram.begin();it!=BrokenHistogram.end();++it)
        Normalization+=it->second;

      for(it=BrokenHistogram.begin();it!=BrokenHistogram.end();++it)
        std::cout<<"    "<<it->first<<"\t\t"<<it->second<<"\t"<<it->second/Normalization<<std::endl;
    }
  };  

}  
