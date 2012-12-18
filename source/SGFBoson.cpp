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
  unsigned long NBins=20;

  std::string Label=p.label();


  RNG::Initialize(Seed);

  OperatorStringType OperatorString(Container);

  /* Initializing the simulation. Thermalize, Measure and pring the results */
  Simulation simul(Label);

  // We start warm up iterations
  simul.Thermalize(OperatorString,WarmIterations,WarmTime);

  // This defines the measurable objects some of which delay updates even if not measured.
  // This is why I declare the measurable operators after the thermalization.
  Measurable MeasuredOperators(OperatorString);

  MeasuredOperators.insert("Potential Energy",Container.V);
  MeasuredOperators.insert("Kinetic Energy",Container.T);


  if(p.HasMeasurable("Number")) {
    MeasuredOperators.insert("Particle Number",Orphans::GenerateNumberOperator(Container.Psi));
  }

  if(p.HasMeasurable("LocalDensity")) {
    InsertLocalDensity(Container.Psi, MeasuredOperators);
  }

  if(p.HasMeasurable("DensityMatrixEigenSystem")) {
    InsertDensityMatrixEigenSystem(Container.Psi, MeasuredOperators);
  }

  if(p.HasMeasurable("DensityMatrix")) {
    InsertDensityMatrix(Container.Psi, MeasuredOperators);
  }


  //We start measurement iterations
  simul.Measure(OperatorString,MeasuredOperators,NBins,MeasIterations,MeasTime);


  std::cout<< "******************************\n";
  std::cout<< "* Model and parameters        *\n";
  std::cout<< "******************************\n\n";

  p.print(std::cout);

  std::cout<<std::endl;

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
