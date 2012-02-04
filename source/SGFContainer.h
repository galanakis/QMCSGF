#include "EasyMathExpression.h" // For access to the interface
#include <sstream>
#include <set>

class SGFContainer {
  std::vector<SGF::Boson> _Psi;
  SGF::Hamiltonian _Kinetic,_Potential;
  double _ConstantEnergy;
  double _Beta;
  int _Ensemble;
  int _Seed;
	
	int _GreenOperatorLines;
	double _AlphaParameter;
	unsigned int _RebuildPeriod;

  std::vector<SGF::Hamiltonian> _MeasurableOperators;

public: 

  inline const SGF::Hamiltonian &Kinetic() const {return _Kinetic;}
  inline const SGF::Hamiltonian &Potential() const {return _Potential;}
  inline double ConstantEnergy() const {return _ConstantEnergy;}
  inline double Beta() const {return _Beta;}
  inline int Ensemble() const {return _Ensemble;}
  inline int Seed() const {return _Seed;}
	inline int GreenOperatorLines() const {return _GreenOperatorLines;}
	inline double AlphaParameter() const {return _AlphaParameter;}

  inline std::vector<SGF::Boson>::size_type Nsites() const {return _Psi.size();} 
  
  
  inline void print_term(const SGF::HamiltonianTerm &term) const {
    char name[100][10]={"","","A","C","AA","AC","CA","CC","AAA","AAC","ACA","ACC","CAA","CAC","CCA","CCC","AAAA","AAAC","AACA","AACC","ACAA","ACAC","ACCA","ACCC","CAAA","CAAC","CACA","CACC","CCAA","CCAC","CCCA","CCCC",};
    std::cout<<(term.coefficient())<<"*";
    for(SGF::HamiltonianTerm::size_type i=0;i<term.product().size();++i) {
      long index=term.product()[i].particle_id()-&_Psi[0];
      int prodelem=term.product()[i].id();
      std::cout<<name[prodelem]<<"["<<index<<"]";
    }
    char dname[2][30]={"Non diagonal","Diagonal"};
    std::cout<<"\t"<<dname[term.diagonal()]<<std::endl;
  }
 
  inline void explicit_print_term(const SGF::HamiltonianTerm &term)  {
    char name[100][10]={"","","A","C","AA","AC","CA","CC","AAA","AAC","ACA","ACC","CAA","CAC","CCA","CCC","AAAA","AAAC","AACA","AACC","ACAA","ACAC","ACCA","ACCC","CAAA","CAAC","CACA","CACC","CCAA","CCAC","CCCA","CCCC",};
    std::cout<<(term.coefficient())<<"*";
    for(SGF::HamiltonianTerm::size_type i=0;i<term.product().size();++i) {
      long index=term.product()[i].particle_id()-&_Psi[0];
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
		_GreenOperatorLines=EasyMathExpression::GreenOperatorLines();
		_AlphaParameter=EasyMathExpression::AlphaParameter();

    EasyMathExpression::OperatorIterator it(_Psi,"#Hamiltonian");
    
    while(it.increment()) {
      SGF::HamiltonianTerm Term=it.Term();
      if(Term.product().size()==0)
        _ConstantEnergy+=Term.coefficient();
      else if(Term.diagonal())
        _Potential.push_back(it.Term());
      else {
        Term.coefficient()*=-1.0;
        _Kinetic.push_back(Term);
      }

    }

    std::set<SGF::Boson*> indexset;
    for(SGF::Hamiltonian::size_type i=0;i<_Kinetic.size();++i)
      for(SGF::HamiltonianTerm::size_type j=0;j<_Kinetic[i].product().size();++j)
        indexset.insert(_Kinetic[i].product()[j].particle_id());


		// Append the extra terms for the grand canonical ensemble
		if(_Ensemble==SGF::GrandCanonical) {
			// Scan all kinetic terms to find all the indices
				std::set<SGF::Boson*>::const_iterator it;
				for(it=indexset.begin();it!=indexset.end();++it) {
					_Kinetic.push_back(SGF::HamiltonianTerm(1.0/NSites(),SGF::A,*it));
					_Kinetic.push_back(SGF::HamiltonianTerm(1.0/NSites(),SGF::C,*it));
				}
		}


    std::vector<std::string> &opnames=EasyMathExpression::MeasurableNameList();
    _MeasurableOperators.resize(opnames.size());

    for(std::vector<std::string>::size_type i=0;i<opnames.size();++i) {
      EasyMathExpression::OperatorIterator oit(_Psi,opnames[i]);

      while(oit.increment())
        _MeasurableOperators[i].push_back(oit.Term());

    }

		/* Adding a default measurable set */
		 
		/*
		int nsites=_Psi.size();
		_MeasurableOperators.resize(opnames.size()+nsites*nsites);
		int i=opnames.size();



		// Additing the density matrix
		for(int a=0;a<nsites;++a)
			for(int b=0;b<nsites;++b) {
				std::stringstream ss;
				ss<<"<C[ "<<a<<" ]A[ "<<b<<" ]>";
				opnames.push_back(ss.str());
				if(a!=b) {
					SGF::HamiltonianTerm Term(1.0,SGF::C,&_Psi[a],SGF::A,&_Psi[b]);
					_MeasurableOperators[i].push_back(Term);
				}
				else {
					SGF::HamiltonianTerm Term(1.0,SGF::C*SGF::A,&_Psi[a]);
					_MeasurableOperators[i].push_back(Term);
				}
					
				++i;
			}

     */

		/* Determining conserved charges */
		
		std::vector<SGF::Boson*> indices;
		indices.insert(indices.begin(),indexset.begin(),indexset.end());
    std::map<SGF::Boson*,unsigned int> indexmap;
		for(unsigned int i=0;i<indices.size();++i)
			indexmap[indices[i]]=i;
		
		
		std::vector< std::map<unsigned int,int> > SparceRows;
		
		std::vector< std::vector<int> > SparceMatrix;
		

		for(SGF::Hamiltonian::size_type i=0;i<_Kinetic.size();++i) {
			std::map<unsigned int,int> temp_map;
			for(SGF::HamiltonianTerm::size_type n=0;n<_Kinetic[i].product().size();++n) {
				unsigned int index=indexmap[_Kinetic[i].product()[n].particle_id()];
				temp_map[index]=_Kinetic[i].product()[n].delta();
				//SparceMatrix[index].push_back(i);
			}
			//SparceRows.push_back(temp_map);
		}


  }
    
  inline int NSites() const {return EasyMathExpression::NSites();}
  void print_matrix_elements() {
    for(SGF::Hamiltonian::size_type i=0;i<Kinetic().size();++i) {
      std::cout<<"i="<<i<<" LEFT "<<Kinetic()[i].me(SGF::LEFT)<<" RIGHT "<<Kinetic()[i].me(SGF::RIGHT)<<std::endl;
    }
  }


};
