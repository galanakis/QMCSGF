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
#include "SGFBase.hh"


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
   unsigned int GreenOperatorLines=4;
   unsigned long WarmTime=40000;
   unsigned long WarmIterations=100000000;
   unsigned long MeasTime=40000;
   unsigned long MeasIterations=100000000;
   unsigned long NBins=20;
   unsigned int NSites=20;
   unsigned int Population=20;

   SGF::SGFBase Container;

   Container.Beta=1.0/1.02;
   Container.Alpha=0.95;
   Container.Ensemble=SGF::Canonical;


   Container.g.initialize(NSites,GreenOperatorLines);


   RNG::Initialize(Seed);


   Container.Psi.resize(NSites);

   for(unsigned int i=0; i<NSites; ++i)
      Container.Psi[i].nmax()=0;

   for(int p=0; p<Population; ++p) {
      unsigned int i=p%NSites;
      Container.Psi[i].n(0)++;
      Container.Psi[i].n(1)++;
   }




   for(unsigned int i=0; i<NSites; ++i) {
      const SGF::IndexedProductElement atom(SGF::C*SGF::C*SGF::A*SGF::A,&Container.Psi[i]);
      Container.V.push_back(SGF::HamiltonianTerm(U/2.0,atom));
   }

   for(unsigned int i=0; i<NSites; ++i) {
      unsigned int j=(1+i)%NSites;
      const SGF::IndexedProductElement ci(SGF::C,&Container.Psi[i]);
      const SGF::IndexedProductElement cj(SGF::C,&Container.Psi[j]);
      const SGF::IndexedProductElement ai(SGF::A,&Container.Psi[i]);
      const SGF::IndexedProductElement aj(SGF::A,&Container.Psi[j]);

      Container.T.push_back(SGF::HamiltonianTerm(t,ci,aj));
      Container.T.push_back(SGF::HamiltonianTerm(t,cj,ai));

   }


   // Building the list of measurable operators
   std::vector<SGF::Hamiltonian> _MeasurableOperators;
   std::vector<std::string> _MeasurableNameList;

   _MeasurableOperators.push_back(Container.V);
   _MeasurableNameList.push_back("Potential Energy");

   _MeasurableOperators.push_back(Container.T);
   _MeasurableNameList.push_back("Atom Kinetic energy");





   

   SGF::OperatorStringType OperatorString(Container);

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
