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

  struct Measurement {
    double value;
    double error;
    Measurement(double v,double e) : value(v), error(e) {}
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

	/* Adding a lot of numbers results in floating point errors. Therefore the accumulator 
	   should have higher precition than the data that are added */
	typedef SGF::BinnedAccumulator<long double> BinnedAccumulatorME;
  //typedef SGF::BinnedAccumulator<MatrixElement> BinnedAccumulatorME;

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
    inline void measure(double Weight) const { sum->push(me(),Weight); }
    inline MatrixElement me() const {return term->me(RIGHT);}
  };


  class Measurable {
    long _ndiagonal;
    long _nflush;
    std::vector<BinnedAccumulatorME> Sums;
    std::vector<double> Constants;
    const long _size;
    std::vector<MeasAccumulators> Diag_Acc;
    std::map<std::map<Boson*,int>,std::vector<MeasAccumulators> > OD_Acc;
    BinnedAccumulatorME _Kinetic;
    BinnedAccumulatorME _Potential;
    BinnedAccumulatorME _TotalEnergy;
    double _BoltzmannWeight;
    
    std::map<unsigned int,unsigned int> BrokenHistogram;

  public:
    Measurable(const std::vector<Hamiltonian> &HList) : _ndiagonal(0), _nflush(0), Sums(HList.size()), Constants(HList.size()), _size(HList.size()), _Kinetic(), _TotalEnergy(), _BoltzmannWeight(0) {

      for(int i=0;i<HList.size();++i)
      for(int j=0;j<HList[i].size();++j) {
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
      for(int i=0;i<Sums.size();++i) 
        Sums[i].reset();
      _BoltzmannWeight=0;
      BrokenHistogram.clear();
    }
    
    inline void flush() {
      ++_nflush;
      _Kinetic.flush(_BoltzmannWeight);
      _Potential.flush(_BoltzmannWeight);
      _TotalEnergy.flush(_BoltzmannWeight);
      for(int i=0;i<Sums.size();++i) 
        Sums[i].flush(_BoltzmannWeight);
      _BoltzmannWeight=0;
    }

    void measure(const OperatorStringType &OperatorString) {
      
      const std::set<Boson*> &list=OperatorString.ListBrokenLines();
      const double Weight=OperatorString.BoltzmannWeight();
      
      if(list.size()==0 && OperatorString.NBrokenLines()!=0 || list.size()!=0 && OperatorString.NBrokenLines()==0 ) {
        std::cout<<"The list of broken lines is not consistent"<<std::endl;
        exit(234);
      }

      BrokenHistogram[OperatorString.NBrokenLines()]+=1;
      
      if(OperatorString.NBrokenLines()!=0) {  // Off diagonal configurations
        std::map<std::map<Boson*,int>,std::vector<MeasAccumulators> >::const_iterator it(OD_Acc.find(BosonDeltaMap(list).map()));
        if(it!=OD_Acc.end()) {
          const std::vector<MeasAccumulators> &v=it->second;
          for(int i=0;i<v.size();++i) 
            v[i].measure(Weight);
        }
      }
      else {  // diagonal configuration
        ++_ndiagonal;
        _BoltzmannWeight+=Weight;
        double MeasKinetic=-OperatorString.length()/OperatorString.Beta();
        double MeasPotential=OperatorString.Energy(SGF::RIGHT);        
        _Kinetic.push(MeasKinetic,Weight);
        _Potential.push(MeasPotential,Weight);
        _TotalEnergy.push(MeasKinetic+MeasPotential,Weight);
        for(int i=0;i<Diag_Acc.size();++i)
          Diag_Acc[i].measure(Weight);
      }

    }

    inline long size() const {return _size;}
    inline long count() const {return _ndiagonal;} 
    inline long nflush() const {return _nflush;}
    inline const Measurement operator[](long i) const {return Measurement(Constants[i]+Sums[i](1),Sums[i].sigma());}
    inline const Measurement KineticEnergy() const {return Measurement(_Kinetic.average(),_Kinetic.sigma());}
    inline const Measurement PotentialEnergy() const {return Measurement(_Potential.average(),_Potential.sigma());}
    inline const Measurement TotalEnergy() const {return Measurement(_TotalEnergy.average(),_TotalEnergy.sigma());}
    
    void print_histogram() const {
      std::cout<<"  *******************************\n";
      std::cout<<"  * Broken worldlines histogram *\n";
      std::cout<<"  *******************************\n\n";
      std::cout<<"    N lines\tCount\tProbability\n\n";
      std::map<unsigned int,unsigned int>::const_iterator it;
      double Normalization=0;
      for(it=BrokenHistogram.begin();it!=BrokenHistogram.end();++it)
        Normalization+=it->second;

      for(it=BrokenHistogram.begin();it!=BrokenHistogram.end();++it)
        std::cout<<"    "<<it->first<<"\t\t"<<it->second<<"\t"<<it->second/Normalization<<std::endl;
    }
  };  

}  
