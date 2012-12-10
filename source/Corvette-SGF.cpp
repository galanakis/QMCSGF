#include <set>
#include <iostream>
#include <cstring>
#include <limits>
#include <map>
#include <vector> 
#include <stack>
#include <iostream>
#include <string>

#include <Parser.h>
#include <MathExpression.h>
#include <ParserSGF.h>
#include <ScrEx.h>
#include <CheckCommandLine.h>
#include <OperatorString.hh>

#include <Simulation.hh> 
#include "HamiltonianTerm.hh"
#include "SGFBase.hh"





std::ostream cout(std::cout.rdbuf());
using std::endl;
using std::string;
using std::fstream;
using std::ios;
using std::streamoff;
using std::ostream;
using std::numeric_limits;

class Iterator {
public:
  std::stack<MathExpression::Node*> pointers;
  Iterator(const std::string &token) { push(MathExpression::Find(token.c_str())->Expression->Root); }     
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


  SGF::HamiltonianTerm Term() const {
    std::vector<SGF::IndexedProductElement> temp;
    std::map<SGF::Boson*,SGF::ProductElement>::const_iterator it;
    for(it=map.begin();it!=map.end();++it)
      temp.push_back(SGF::IndexedProductElement(it->second,it->first));

    return SGF::HamiltonianTerm(Coefficient,temp);
  }

};


void PrintTokens(std::ostream &o) {

    o << "*******************************************************************************************\n";
    o << "* This is a quantum Monte Carlo simulation performed by the \"Corvette SGF engine\".        *\n";
    o << "*                                                                                         *\n";
    o << "* For informations on the SGF algorithm:                                                  *\n";
    o << "*   Physical Review E 77, 056705 (2008)                                                   *\n";
    o << "*   Physical Review E 78, 056707 (2008)                                                   *\n";
    o << "*                                                                                         *\n";
    o << "* Dr Valy G. Rousseau and Dr Dimitris Galanakis";
    for (unsigned int i=0;i<32-strlen("");i++) o << " ";
    o << "*\n";
    o << "*******************************************************************************************\n\n";
    o << "***************************\n";
    o << "* User's input parameters *\n";
    o << "***************************\n\n  ";
    Parser::TokenHandle Input=Parser::First();

    while (Input)
    {
      if (Input->Type()==Parser::Number)
        o << *(double *) Input->Value();

      else
      {
        char *C=(char *) Input->Value();

        if (Input->Type()==Parser::String)
          o << "\"" << C << "\"";
        else
          o << C;

        if (MathExpression::IsKeyword(C))
          o << " ";

        else if (*C==';')
        {
          o << "\n";

          if (Input->NextToken())
            o << "  ";
        }
      }

      Input=Input->NextToken();
    }

    o << std::endl;


}

void Simulator() {



  unsigned int NumIndices=MathExpression::GetNumIndices();
  unsigned int NumSpecies=MathExpression::GetNumSpecies();
  unsigned int NumSites=NumIndices/NumSpecies;
  unsigned int NSites=MathExpression::GetNumIndices()/MathExpression::GetNumSpecies();
  unsigned int GreenOperatorLines=MathExpression::GetValue("GreenOperatorLines").Re();
  unsigned long WarmTime=MathExpression::GetValue("#WarmTime").Re();
  unsigned long WarmIterations=(MathExpression::Find("WarmIterations")!=NULL) ? static_cast<unsigned long>(MathExpression::GetValue("WarmIterations").Re()) : std::numeric_limits<unsigned long>::max();
  unsigned long MeasTime=MathExpression::GetValue("#MeasTime").Re();
  unsigned long MeasIterations=(MathExpression::Find("MeasIterations")!=NULL) ? static_cast<unsigned long>(MathExpression::GetValue("MeasIterations").Re()) : std::numeric_limits<unsigned long>::max();      
  unsigned long NBins=MathExpression::GetValue("#Bins").Re();

  SGF::SGFBase Container;

  Container.Beta=MathExpression::GetValue("#InverseTemperature").Re();
  Container.Alpha=MathExpression::GetValue("AlphaParameter").Re();
  Container.Ensemble=MathExpression::GetValue("#Ensemble").Re() == 0 ? SGF::Canonical : SGF::GrandCanonical;


  Container.g.initialize(NSites,GreenOperatorLines);

  // Initializing the configuration
  Container.Psi.resize(NumIndices);
  for(std::vector<SGF::Boson>::size_type i=0;i<Container.Psi.size();++i)
    Container.Psi[i].nmax()=MathExpression::GetNmax(i);

  for(unsigned int species=0;species<NumSpecies;++species) {
    for(int particle=0;particle<MathExpression::GetPopulation(species);++particle) {
      unsigned int i=NumSites*species+particle%NumSites;
      Container.Psi[i].nL()++;
      Container.Psi[i].nR()++;
    }
  }  


  // Initializing the kinetic and potential operators
  OperatorIterator it(Container.Psi,"#Hamiltonian");

   while(it.increment()) {
      SGF::HamiltonianTerm Term=it.Term();
      if(Term.product().size()!=0) {
         if(Term.diagonal())
            Container.V.push_back(it.Term());
         else {
            Term.coefficient()*=-1.0;
            Container.T.push_back(Term);
         }
      }

   }

  // Building the list of measurable operators
  std::vector<SGF::Hamiltonian> _MeasurableOperators;

  std::vector<std::string> &opnames=MathExpression::GetMeasurableList();
  _MeasurableOperators.resize(opnames.size());

  for(std::vector<std::string>::size_type i=0;i<opnames.size();++i) {
    OperatorIterator oit(Container.Psi,opnames[i]);

    while(oit.increment())
      _MeasurableOperators[i].push_back(oit.Term());

  }


  SGF::OperatorStringType OperatorString(Container);

	/* Initializing the simulation. Thermalize, Measure and pring the results */
  Simulation simul(MathExpression::GetSimulName(),cout);


	// We start warm up iterations
  simul.Thermalize(OperatorString,WarmIterations,WarmTime);

	// This defines the measurable objects some of which delay updates even if not measured.
	// This is why I declare the measurable operators after the thermalization.
  SGF::Measurable MeasuredOperators(OperatorString);

  std::vector<std::string> _MeasurableNameList = MathExpression::GetMeasurableList();
  if(_MeasurableNameList.size()!=_MeasurableOperators.size()) {
    std::cout<<"The number of of measurable operators does not match the number of the labels"<<std::endl;
    exit(3);
  }

  for(unsigned int i=0;i<_MeasurableOperators.size();++i) {
    MeasuredOperators.insert(_MeasurableNameList[i],_MeasurableOperators[i]);
  }


	//We start measurement iterations
  simul.Measure(OperatorString,MeasuredOperators,NBins,MeasIterations,MeasTime);
  
  // We display the tokens and parameters used in the simulation
  PrintTokens(cout);

  // We diplay the results of the simulation
  simul.Results(MeasuredOperators);  

}

// ***************************
// * Here starts the program *
// ***************************

int main(int NumArg,char **Arg)
  {
    Simulation::InitializeEnvironment();

    // *********************************************************************************
    // * We check the command line. If it is not correct, a help message is displayed. *
    // *********************************************************************************
    
    char *File=NULL;
    
    switch (CheckCommandLine(NumArg,Arg))
      {
        case 0: return 0;               // Simulation is aborded.
        case 1: File=Arg[1]; break;     // Command line is correct and requests to start the simulation.
        case 2: File=Arg[2]; break;     // Command line is correct and requests to display Hamiltonian terms.
        case 3: File=Arg[3]; break;     // Command line is correct and requests to display a specified quantity.
      }
    
    // ******************************************************************
    // * We read the input file and transform it into a list of tokens. *
    // ******************************************************************
    
    char *Error;

    if ((Error=Parser::ReadFile(File)))
      {
        cout << Error << endl;
        return 0;
      }
      
    // ********************************************************
    // * We transform the tokens into a list of SGF commands. *
    // ********************************************************

    ParserSGF ScriptSGF;
    
    if ((Error=ScriptSGF.ReadTokens()))
      {
        cout << Error << endl;
        return 0;
      }

    MathExpression::Initialize();	// We initialize the table of symbols with keywords and constants.
    
    // *******************************
    // * We execute the SGF commands *
    // *******************************
    
    if ((Error=ScriptSGF.ExecuteCommands()))
      {
        cout << Error << endl;
        return 0;
      }
    
    // ************************************************
    // * We build a suitable form of the Hamiltonian. *
    // ************************************************
    
    if ((Error=MathExpression::BuildHamiltonian()))
      {
        cout << Error << endl;
        return 0;
      }
      
    if (File==Arg[1])
      {
        // *************************
        // * Start the simulation. *
        // *************************

        if ((Error=MathExpression::SignOrPhaseProblem()))
          {
            cout << Error << endl;      // A sign or a phase problem has been detected. A warning message
            return 0;                   // is displayed and the simulation is aborded.
          }

					Simulator();

      }
      
    else if (File==Arg[2])
      MathExpression::ListHamiltonianTerms();                   // Display hamiltonian terms and abord the simulation.
      
    else
      {
        // ************************************************************************
        // * Display the quantity specified by the user and abord the simulation. *
        // ************************************************************************

        string StrError=MathExpression::DisplayQuantity(Arg[2]);
        
        if (StrError.length()!=0)
          cout << StrError << endl;
      }

      Simulation::FinalizeEnvironment();
    return 0;
  }
