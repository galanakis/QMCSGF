#include <iostream>
#include <cstring>
#include <limits>
#include <set>
#include <map>
#include <vector> 
#include <stack>
#include <string>

#include "HamiltonianTerm.hh"
#include <OperatorString.hh>
#include <Simulation.hh>
#include "Conventions.hh"


#ifdef USEMPI
// *******************************************
// * Declaration of global variables for MPI *
// *******************************************
#include <mpi.h>
#define Master 0
int NumProcessors,Rank,NameLength;
char ProcessorName[MPI_MAX_PROCESSOR_NAME];

#endif



std::ostream cout(std::cout.rdbuf());


void BoseHubbardPeriodic1D() {


  int Seed=34715;
  SGF::MatrixElement U=8.0;
  SGF::MatrixElement t=1.0;
  double Beta=1.0/1.02;
  double AlphaParameter=0.95;
  unsigned int GreenOperatorLines=4;
  unsigned long WarmTime=40000;
  unsigned long WarmIterations=100000000;
  unsigned long MeasTime=40000;
  unsigned long MeasIterations=100000000;
  unsigned long NBins=20;
  unsigned int Ensemble=SGF::Canonical;


  RNG::Initialize(Seed);

  std::vector<SGF::Boson> _Psi;

    unsigned int NSites=20;
    unsigned int Population=20;
    
    _Psi.resize(NSites);

    for(unsigned int i=0;i<NSites;++i)
      _Psi[i].nmax()=0;

    for(int p=0;p<Population;++p) {
      unsigned int i=p%NSites;
      _Psi[i].n(0)++;
      _Psi[i].n(1)++;
    }


  SGF::Hamiltonian T,V;
  
  for(unsigned int i=0;i<NSites;++i) {
    const SGF::IndexedProductElement atom(SGF::C*SGF::C*SGF::A*SGF::A,&_Psi[i]);
    V.push_back(SGF::HamiltonianTerm(U/2.0,atom));
  }

  for(unsigned int i=0;i<NSites;++i) {
    unsigned int j=(1+i)%NSites;
    const SGF::IndexedProductElement ci(SGF::C,&_Psi[i]);
    const SGF::IndexedProductElement cj(SGF::C,&_Psi[j]);
    const SGF::IndexedProductElement ai(SGF::A,&_Psi[i]);
    const SGF::IndexedProductElement aj(SGF::A,&_Psi[j]);

    T.push_back(SGF::HamiltonianTerm(t,ci,aj));
    T.push_back(SGF::HamiltonianTerm(t,cj,ai));

  }


  // Building the list of measurable operators
  std::vector<SGF::Hamiltonian> _MeasurableOperators;
  std::vector<std::string> _MeasurableNameList;

  _MeasurableOperators.push_back(V);
  _MeasurableNameList.push_back("Potential Energy");

  _MeasurableOperators.push_back(T);
  _MeasurableNameList.push_back("Atom Kinetic energy");




	
  SGF::GreenOperator<long double> g(NSites,GreenOperatorLines);

	SGF::OperatorStringType OperatorString(T,V,Beta,g,AlphaParameter);

	/* Initializing the simulation. Thermalize, Measure and pring the results */
	Simulation simul("BoseHubbard",cout);


  

	// We start warm up iterations
	simul.Thermalize(OperatorString,WarmIterations,WarmTime);

	// This defines the measurable objects some of which delay updates even if not measured.
	// This is why I declare the measurable operators after the thermalization.
	SGF::Measurable MeasuredOperators(OperatorString);
	MeasuredOperators.insert(_MeasurableNameList,_MeasurableOperators);
  
	//We start measurement iterations
	simul.Measure(OperatorString,MeasuredOperators,NBins,MeasIterations,MeasTime);
  
  // We diplay the results of the simulation
	simul.Results(MeasuredOperators);  
	
}

// ***************************
// * Here starts the program *
// ***************************

int main() {

#ifdef USEMPI
    // ***********************************
    // * Initialization of MPI functions *
    // ***********************************
    
    MPI_Init(&NumArg,&Arg);
    MPI_Comm_size(MPI_COMM_WORLD,&NumProcessors);
    MPI_Comm_rank(MPI_COMM_WORLD,&Rank);
    MPI_Get_processor_name(ProcessorName,&NameLength);

		// Silence the other nodes. Only the root node can print
		cout.rdbuf( Rank==Master ? std::cout.rdbuf() : 0 );
 
#endif

    BoseHubbardPeriodic1D();

    
#ifdef USEMPI
  return MPI_Finalize();
#endif
  return 0;
  }
