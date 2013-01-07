#include <iostream>
#include <cstring>
#include <limits>
#include <set>
#include <map>
#include <vector>
#include <stack>
#include <string>
#include <sstream>

#include "SGF.hh"
#include "Simulation.hh"
#include "ExtraMeasurables.hh"


namespace SGF {


void BoseHubbardPeriodic1DMakeContainer(SGFBase &Container) {



  // Model specific parameters
  MatrixElement U=8.0;
  MatrixElement t=1.0;

  unsigned int Lx=20;

  unsigned int NSites=Lx;
  unsigned int Population=NSites;
  unsigned int Nmax=0;


  Container.Beta=1.0/1.02;
  Container.Alpha=0.95;
  Container.Ensemble=Canonical;


  unsigned int GreenOperatorLines=4;
  Container.g.initialize(NSites,GreenOperatorLines);


  Container.Psi.resize(NSites);

  for(unsigned int i=0; i<NSites; ++i)
    Container.Psi[i].nmax() = Nmax;

  for(int p=0; p<Population; ++p) {
    unsigned int i=p%NSites;
    Container.Psi[i].nL()++;
    Container.Psi[i].nR()++;
  }




  for(unsigned int i=0; i<NSites; ++i) {
    const IndexedProductElement atom(C*C*A*A,&Container.Psi[i]);
    Container.V.push_back(HamiltonianTerm(U/2.0,atom));
  }


  list_links_t links=links_square(Lx,periodic);

  for(list_links_t::const_iterator it=links.begin(); it!=links.end(); ++it) {
    unsigned int i=it->first;
    unsigned int j=it->second;
    const IndexedProductElement ci(C,&Container.Psi[i]);
    const IndexedProductElement cj(C,&Container.Psi[j]);
    const IndexedProductElement ai(A,&Container.Psi[i]);
    const IndexedProductElement aj(A,&Container.Psi[j]);

    Container.T.push_back(HamiltonianTerm(t,ci,aj));
    Container.T.push_back(HamiltonianTerm(t,cj,ai));


  }



}

void BoseHubbardPeriodic1D() {


  int Seed=34715;

  unsigned long WarmTime=40000;
  unsigned long WarmIterations=100000000;
  unsigned long MeasTime=40000;
  unsigned long MeasIterations=100000000;
  unsigned long NBins=20;

  SGFBase Container;

  BoseHubbardPeriodic1DMakeContainer(Container);


  RNG::Initialize(Seed);

  OperatorStringType OperatorString(Container);

  /* Initializing the simulation. Thermalize, Measure and pring the results */
  Simulation simul("BoseHubbard");

  // We start warm up iterations
  simul.Thermalize(OperatorString,WarmIterations,WarmTime);

  // This defines the measurable objects some of which delay updates even if not measured.
  // This is why I declare the measurable operators after the thermalization.
  Measurable MeasuredOperators(OperatorString);

  InsertOperator("Potential Energy",Container.V,MeasuredOperators);
  InsertOperator("Atom Kinetic energy",Container.T,MeasuredOperators);
  InsertOperator("Number of Particles",Orphans::GenerateNumberOperator(Container.Psi),MeasuredOperators);
  InsertLocalDensity("Local Density",Container.Psi, MeasuredOperators);


  //We start measurement iterations
  simul.Measure(OperatorString,MeasuredOperators,NBins,MeasIterations,MeasTime);

  std::cout<<std::endl;
  // We diplay the results of the simulation
  simul.Results(MeasuredOperators);

}

}

// ***************************
// * Here starts the program *
// ***************************

int main(int argc,char *argv[]) {

  SGF::InitializeEnvironment(argc,argv);
  SGF::BoseHubbardPeriodic1D();
  SGF::FinalizeEnvironment();

}
