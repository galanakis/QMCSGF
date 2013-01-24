#include "Parameters.hh"
#include "ExtraMeasurables.hh"
#include "Simulation.hh"

namespace SGF {


void Simulator(const SGF::Parameters &p) {

  SGFBase Container;

  p.MakeContainer(Container);


  int Seed=p.Seed;

  unsigned long WarmTime=p.WarmTime;
  unsigned long WarmIterations=p.WarmIterations;
  unsigned long MeasTime=p.MeasTime;
  unsigned long MeasIterations=p.MeasIterations;
  unsigned long NBins=p.NBins;

  std::string Label=p.label();


  RNG::Initialize(Seed);

  OperatorStringType OperatorString(Container);

  /* Initializing the simulation. Thermalize, Measure and pring the results */
  Simulation simul(Label);

  // We start warm up iterations
  simul.Thermalize(OperatorString,WarmIterations,WarmTime);

  // This defines the measurable objects some of which delay updates even if not measured.
  // This is why I declare the measurable operators after the thermalization.
  Measurable MeasuredOperators;

  p.MakeMeasurables(Container,MeasuredOperators);

  //We start measurement iterations
  simul.Measure(OperatorString,MeasuredOperators,NBins,MeasIterations,MeasTime);


  p.print(std::cout);



  // We diplay the results of the simulation
  simul.Results(MeasuredOperators);


}


}

int main(int argc, char *argv[]) {

  SGF::InitializeEnvironment(argc,argv);

  if(argc!=2) {
    std::cout<<"No input file"<<std::endl;
    exit(1);
  }

  std::string json=readInputFile(argv[1]);

  SGF::Parameters p(json);

  Simulator(p);
  SGF::FinalizeEnvironment();

}
