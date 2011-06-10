#include "EasyMathExpression.h" // For access to the interface
#include <set>

class SGFContainer {
  std::vector<SGF::Boson> _Psi;
  SGF::Hamiltonian _Kinetic,_Potential;
  double _ConstantEnergy;
  double _Beta;
  int _Ensemble;
  int _Seed;

  std::vector<SGF::Hamiltonian> _MeasurableOperators;

public: 

  inline const SGF::Hamiltonian &Kinetic() const {return _Kinetic;}
  inline const SGF::Hamiltonian &Potential() const {return _Potential;}
  inline double ConstantEnergy() const {return _ConstantEnergy;}
  inline double Beta() const {return _Beta;}
  inline int Ensemble() const {return _Ensemble;}
  inline int Seed() const {return _Seed;}
  inline int Nsites() const {return _Psi.size();} 
  
  inline unsigned int boson2index(SGF::Boson *a) const {return a-&_Psi[0];}
  
  inline void print_term(const SGF::HamiltonianTerm &term) const {
    char name[100][10]={"","","A","C","AA","AC","CA","CC","AAA","AAC","ACA","ACC","CAA","CAC","CCA","CCC","AAAA","AAAC","AACA","AACC","ACAA","ACAC","ACCA","ACCC","CAAA","CAAC","CACA","CACC","CCAA","CCAC","CCCA","CCCC",};
    std::cout<<(term.coefficient())<<"*";
    for(int i=0;i<term.product().size();++i) {
      int index=term.product()[i].particle_id()-&_Psi[0];
      int prodelem=term.product()[i].id();
      std::cout<<name[prodelem]<<"["<<index<<"]";
    }
    char dname[2][30]={"Non diagonal","Diagonal"};
    std::cout<<"\t"<<dname[term.diagonal()]<<std::endl;
  }
 
  inline void explicit_print_term(const SGF::HamiltonianTerm &term)  {
    char name[100][10]={"","","A","C","AA","AC","CA","CC","AAA","AAC","ACA","ACC","CAA","CAC","CCA","CCC","AAAA","AAAC","AACA","AACC","ACAA","ACAC","ACCA","ACCC","CAAA","CAAC","CACA","CACC","CCAA","CCAC","CCCA","CCCC",};
    std::cout<<(term.coefficient())<<"*";
    for(int i=0;i<term.product().size();++i) {
      int index=term.product()[i].particle_id()-&_Psi[0];
      int prodelem=term.product()[i].id();
      std::cout<<name[prodelem]<<"["<<index<<"->("<<_Psi[index].n(SGF::LEFT)<<","<<_Psi[index].n(SGF::RIGHT)<<","<<_Psi[index].nmax()<<")] ";
    }
    char dname[2][30]={"Non diagonal","Diagonal"};
    std::cout<<"\t(meL,meR): "<<term.me(SGF::LEFT)<<","<<term.me(SGF::RIGHT)<<")\t"<<dname[term.diagonal()]<<std::endl;
  }



  inline const std::vector<SGF::Hamiltonian> &MeasurableOperators() const {return _MeasurableOperators;}

  SGFContainer() {
    _Psi=EasyMathExpression::GeneratePsi();
    _Beta=EasyMathExpression::Beta();
    _Ensemble=EasyMathExpression::Ensemble();
    _Seed=EasyMathExpression::Seed();

    EasyMathExpression::OperatorIterator it(_Psi,"#Hamiltonian");
    
    while(it.increment()) {
      SGF::HamiltonianTerm Term=it.Term();
      if(Term.product().size()==0)
        _ConstantEnergy+=Term.coefficient();
      else if(Term.diagonal())
        _Potential.push_back(it.Term());
      else {
        Term.set_coefficient(-Term.coefficient());
        _Kinetic.push_back(Term);
      }

    }

		// Append the extra terms for the grand canonical ensemble
		if(_Ensemble==SGF::GrandCanonical) {
			// Scan all kinetic terms to find all the indices
		    std::set<SGF::Boson*> indexset;
		    for(int i=0;i<_Kinetic.size();++i)
		      for(int j=0;j<_Kinetic[i].product().size();++j)
		        indexset.insert(_Kinetic[i].product()[j].particle_id());

				std::set<SGF::Boson*>::const_iterator it;
				for(it=indexset.begin();it!=indexset.end();++it) {
					_Kinetic.push_back(SGF::HamiltonianTerm(1.0/NSites(),SGF::A,*it));
					_Kinetic.push_back(SGF::HamiltonianTerm(1.0/NSites(),SGF::C,*it));
				}
		}


    std::vector<std::string> opnames=EasyMathExpression::MeasurableNameList();
    _MeasurableOperators.resize(opnames.size());

    for(int i=0;i<opnames.size();++i) {
      EasyMathExpression::OperatorIterator oit(_Psi,opnames[i]);

      while(oit.increment())
        _MeasurableOperators[i].push_back(oit.Term());

    } 
    
  }
    
  inline int NSites() const {return EasyMathExpression::NSites();}
  void print_matrix_elements() {
    for(int i=0;i<Kinetic().size();++i) {
      std::cout<<"i="<<i<<" LEFT "<<Kinetic()[i].me(SGF::LEFT)<<" RIGHT "<<Kinetic()[i].me(SGF::RIGHT)<<std::endl;
    }
  }

  void print() {
    std::cout<<"Beta:\t"<<_Beta<<std::endl;
    std::cout<<"Ensemble:\t"<<_Ensemble<<std::endl;
    std::cout<<"Seed:\t"<<_Seed<<std::endl; 

    std::cout<<"Number of _Psi indices:\t"<<_Psi.size()<<std::endl;
    std::cout<<"Number of _Kinetic Terms:\t"<<_Kinetic.size()<<std::endl;
    std::cout<<"Number of _Potential Terms:\t"<<_Potential.size()<<std::endl;
    std::cout<<"Number of measurable Terms:\t"<<_MeasurableOperators.size()<<std::endl;
    for(int i=0;i<_MeasurableOperators.size();++i)
      std::cout<<"Size of Operator ["<<i<<"] "<<_MeasurableOperators[i].size()<<std::endl;


  }


};
