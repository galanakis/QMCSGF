#include <map>
#include <vector> 
#include <stack>
#include <iostream>
#include <string>

#include "HamiltonianTerm.hh"

class EasyMathExpression {
public:
  static inline double Beta() {return MathExpression::Find("#InverseTemperature")->Expression->Root->Evaluate().Re(); }
  static inline int Ensemble() {return MathExpression::Find("#Ensemble")->Expression->Root->Evaluate().Re(); }
  static inline int Seed() {return MathExpression::Find("#Seed")->Expression->Root->Evaluate().Re(); }
  static inline std::vector<string> &MeasurableNameList() {return MathExpression::GetMeasurableList(); }
  static inline int NSites() {return MathExpression::GetNumIndices()/MathExpression::GetNumSpecies();}
	static inline int GreenOperatorLines() {return MathExpression::GetValue("GreenOperatorLines").Re();}
	static inline unsigned int RebuildPeriod() {return MathExpression::GetValue("RebuildPeriod").Re();}
	static inline double AlphaParameter() {return MathExpression::GetValue("AlphaParameter").Re();}

  class Iterator {
  public:
    stack<MathExpression::Node*> pointers;
    Iterator(std::string token) { push(MathExpression::Find(token.c_str())->Expression->Root); }     
    MathExpression::Node* value() const {return pointers.top(); } 
    inline bool end() { return pointers.empty(); }
    inline void push(MathExpression::Node *Node) {
      while(Node!=NULL) {
        pointers.push(Node);
        Node=Node->Son1;
      }   
    }
    inline void increment() {
      if(!end()) {
        MathExpression::Node *Node=value()->Son2;
        pointers.pop();
        push(Node);
      }
    }
  };

  
  class OperatorIterator {
    std::map<SGF::Boson*,SGF::ProductElement> map;
    SGF::MatrixElement Coefficient;
    Iterator it;
    std::vector<SGF::Boson> &Psi;
  public:
    OperatorIterator(std::vector<SGF::Boson> &_Psi,std::string token) : map(), Coefficient(0), it(token), Psi(_Psi) {}
    
    inline bool end() {return it.end();}
    bool increment() {
      
      bool flag=!it.end();
      
      map.clear();
      Coefficient=0;
      

      if(!it.end()) { 

        Coefficient=it.value()->Value.Re();
        it.increment();

        while(!it.end() && it.value()->Type!=MathExpression::ArithmeticOperator1) {
          it.increment();
          SGF::Boson* boson=&Psi[it.value()->Value.Re()];
          it.increment();
          SGF::ProductElement F=(it.value()->Type==MathExpression::CreationOperator) ? SGF::C : SGF::A;
          map[boson]=(map.find(boson)==map.end())?F:map[boson]*F;
          it.increment();
        }
        it.increment();
      }
      
      return flag;
            
    }
    
    SGF::HamiltonianTerm Term() const {return SGF::HamiltonianTerm(Coefficient,map);}
    
  };

  static void build_terms() {

    std::vector<SGF::Boson> Psi=GeneratePsi();
    
    SGF::Hamiltonian Kinetic,Potential;
    
    for(std::vector<SGF::Boson>::size_type i=0;i<Psi.size();++i) {
      std::cout<<"Psi("<<i<<")="<<"("<<Psi[i].n(SGF::LEFT)<<","<<Psi[i].n(SGF::RIGHT)<<","<<Psi[i].nmax()<<")"<<std::endl;
    }


    OperatorIterator it(Psi,"#Hamiltonian");
    while(it.increment()) {

      SGF::HamiltonianTerm Term(it.Term());
      if(Term.diagonal()) 
        Potential.push_back(Term);
      else
        Kinetic.push_back(Term);

      std::cout<<Term.delta()<<", "<<Term.diagonal()<<", "<<Term.length()<<std::endl;

    }
 
    std::cout<<"__ Number of Kinetic Terms "<<Kinetic.size()<<std::endl;
    std::cout<<"__ Number of Potential Terms "<<Potential.size()<<std::endl;


  }
  
  static std::vector<SGF::Boson> GeneratePsi() {
    
    std::vector<SGF::Boson> _Psi;
    
    _Psi.resize(MathExpression::GetNumIndices());
    for(int i=0;i<MathExpression::GetNumIndices();++i)
      _Psi[i]=SGF::Boson(0,0,MathExpression::GetNmax(i));

    int NSites=MathExpression::GetNumIndices()/MathExpression::GetNumSpecies();

    for(int species=0;species<MathExpression::GetNumSpecies();++species) {
      for(int particle=0;particle<MathExpression::GetPopulation(species);++particle) {
        int i=NSites*species+particle%NSites;
        _Psi[i].n(SGF::LEFT)++;
        _Psi[i].n(SGF::RIGHT)++;
      }
    }  
      
    return _Psi;
    
  }

  static void PrintTerms(MathExpression::Node *Father) {

    switch (Father->Type) {
      case MathExpression::ComplexNumber:

      cout << endl << Father->Value.Re();
      break;

      // * We found a + or a * operator, so we treat the subtrees by induction *
      case MathExpression::ArithmeticOperator1: // Plus
      case MathExpression::ArithmeticOperator3: // Times

      PrintTerms(Father->Son1);
      PrintTerms(Father->Son2);

      break;

        // *We found a creation operator *
      case MathExpression::CreationOperator:

      cout << "C[" << (int) Father->Son1->Value.Re() << "]"; 
      break;

      // * We found an annihilation operator *
      case MathExpression::AnnihilationOperator:

      cout << "A[" << (int) Father->Son1->Value.Re() << "]";
      break;

      default:
      cout << "Whoops... This is probably a bug: Found a type of node that is not supposed to occur in this function!\n";
      cout << "Type=" << Father->Type << endl;
    }
  }

public:

      // *********************************************************************
      // * This function looks for the Hamiltonian tree in the symbol table, *
      // * and displays the list of terms.                                   *
      // *********************************************************************

  static void ListTerms(void) {


    build_terms();
    PrintTerms(MathExpression::Find("#Hamiltonian")->Expression->Root); 

  }
};



