#include "Parameters.hh"
#include "ExtraMeasurables.hh"
#include <iostream>

#include <fstream>
#include <ctime>

#include "OperatorString.hh"
#include "Measurable.hh"
#include "Utilities.hh"

namespace SGF {



// Maximum number of broken lines is only used for the BrokenLineHistogram.
#define MAXNUMBROKENLINES 100

class Simulation {
public:

  typedef  std::vector<unsigned long> BrokenHistogramType;
  BrokenHistogramType BrokenHistorgram;

  unsigned long WarmTime;
  double ActualWarmTime;
  unsigned long WarmIterations;
  unsigned long NumWarmUpdates;
  unsigned long NumDirectedWarmUpdates;


  unsigned long MeasTime;
  double ActualMeasTime;
  unsigned long MeasIterations;
  unsigned long NumMeasUpdates;
  unsigned long NumDirectedMeasUpdates;

  unsigned long NumBins;
  int Seed;

  std::string SimulName;

  SGFBase Container;
  Measurable MeasuredOperators;

  std::string outputConfiguration;


  Simulation(const SGF::Parameters &p) {


    int Seed=p.Seed;

    RNG::Initialize(Seed);

    WarmTime=p.WarmTime;
    ActualWarmTime=0;
    WarmIterations=p.WarmIterations;
    NumWarmUpdates=0;
    NumDirectedWarmUpdates=0;

    MeasTime=p.MeasTime;
    ActualMeasTime=0;
    MeasIterations=p.MeasIterations;
    NumMeasUpdates=0;
    NumDirectedMeasUpdates=0;
    NumBins=p.NBins;

    SimulName=p.modelID;

    BrokenHistorgram.resize(MAXNUMBROKENLINES);

    p.MakeContainer(Container);
    p.MakeMeasurables(Container,MeasuredOperators);

    outputConfiguration=p.outputConfiguration;

  }

  template<class OperatorStringType>
  void run() {
    OperatorStringType OpString(Container);

    ProgressBar warmProgress(SimulName+std::string(": Thermalizing "),&NumWarmUpdates,WarmIterations,WarmTime);
    Timer time(&NumWarmUpdates,WarmIterations,WarmTime);

    do {

      do {
        NumWarmUpdates+=OpString.directed_update();
        ++NumDirectedWarmUpdates;
      } while(OpString.NBrokenLines()!=0);

      warmProgress.Update();

    } while ( time.Continue() );

    ActualWarmTime=time.ElapsedSeconds();


    warmProgress.reset_cout();

    ProgressBar measProgress(SimulName+std::string(": Measuring    "),&NumMeasUpdates,MeasIterations,MeasTime);

    NumMeasUpdates=0;
    NumDirectedMeasUpdates=0;

    SGF::BrokenLines BrokenLineTracer(OpString);  // It traces the list of broken lines



    for (unsigned int i=0; i<NumBins; ++i) {

      Timer time(&NumMeasUpdates,MeasIterations/NumBins,MeasTime/NumBins);

      do {

        do {
          NumMeasUpdates+=OpString.directed_update();   // Perform an update.
          ++NumDirectedMeasUpdates;
          MeasuredOperators.measure(OpString,BrokenLineTracer());          // Perform measurements.
          BrokenHistorgram[OpString.NBrokenLines()]+=1;

        } while(OpString.NBrokenLines()!=0);

        measProgress.Update();

      } while ( time.Continue() );

      MeasuredOperators.flush();                               // Bin the data.

    }

    ActualMeasTime=measProgress.ElapsedSeconds();

    measProgress.reset_cout();

  }


  void run() {

    if(Container.Ensemble==Canonical)
      run<CanonicalOperatorString>();
    else
      run<GrandOperatorString>();


  }

  void Results(std::ostream &o) {


#ifdef USEMPI

    unsigned long send_NumWarmUpdates(NumWarmUpdates);
    unsigned long send_NumDirectedWarmUpdates(NumDirectedWarmUpdates);
    double send_ActualWarmTime(ActualWarmTime);

    MPI_Reduce(&send_NumWarmUpdates,        &NumWarmUpdates,        1,MPI_UNSIGNED_LONG,MPI_SUM,Master,MPI_COMM_WORLD);
    MPI_Reduce(&send_NumDirectedWarmUpdates,&NumDirectedWarmUpdates,1,MPI_UNSIGNED_LONG,MPI_SUM,Master,MPI_COMM_WORLD);
    MPI_Reduce(&send_ActualWarmTime,        &ActualWarmTime,        1,MPI_DOUBLE,       MPI_SUM,Master,MPI_COMM_WORLD);

    unsigned long send_NumMeasUpdates(NumMeasUpdates);
    unsigned long send_NumDirectedMeasUpdates(NumDirectedMeasUpdates);
    double send_ActualMeasTime(ActualMeasTime);

    MPI_Reduce(&send_NumMeasUpdates,        &NumMeasUpdates,        1,MPI_UNSIGNED_LONG,MPI_SUM,Master,MPI_COMM_WORLD);
    MPI_Reduce(&send_NumDirectedMeasUpdates,&NumDirectedMeasUpdates,1,MPI_UNSIGNED_LONG,MPI_SUM,Master,MPI_COMM_WORLD);
    MPI_Reduce(&send_ActualMeasTime,        &ActualMeasTime,        1,MPI_DOUBLE,       MPI_SUM,Master,MPI_COMM_WORLD);


#endif



    o<<"\nResults:\n";

    unsigned int w=30;
    unsigned int w2=13;
    o<<"\n";
    o<<"  Thermalization:\n";
    o<<"    Iterations:                     "<<NumWarmUpdates<<"\n";
    o<<"    Time:                           "<<ActualWarmTime<<std::endl;
    o<<"    Iterations per sec:             "<<std::setprecision(2)<<std::fixed<<NumWarmUpdates/ActualWarmTime<<"\n";
    o<<"    Directed Updates:               "<<NumDirectedWarmUpdates<<"\n";
    o<<"    Update Length:                  "<<double(NumWarmUpdates)/NumDirectedWarmUpdates<<"\n";
    o<<"\n";
    o<<"  Measurements:\n";
    o<<"    Iterations:                     "<<NumMeasUpdates<<"\n";
    o<<"    Time:                           "<<ActualMeasTime<<"\n";
    o<<"    Iterations per sec:             "<<std::setprecision(2)<<std::fixed<<NumMeasUpdates/ActualMeasTime<<"\n";
    o<<"    Directed Updates:               "<<NumDirectedMeasUpdates<<"\n";
    o<<"    Update Length:                  "<<double(NumMeasUpdates)/NumDirectedMeasUpdates<<"\n";
    o<<"    Measurement Count:              "<< BrokenHistorgram[0]<<"\n";
    o<<"    Iterations per measurement:     "<<std::fixed<< double(NumMeasUpdates)/BrokenHistorgram[0]<<"\n";

#ifdef USEMPI
    o<<"    Number of Processors:           "<<NumProcessors<<"\n";
#endif

    o<<"\n";
    o<<"  Broken worldlines:\n\n";

    double Normalization=0;
    for(int i=0; i<BrokenHistorgram.size(); ++i)
      Normalization+=BrokenHistorgram[i];

    for(int i=0; i<BrokenHistorgram.size(); ++i)
      if(BrokenHistorgram[i]!=0)
        o<<"    - [ "<<std::right<<std::fixed<<std::setw(3)<<i<<","<<std::setw(12)<<BrokenHistorgram[i]<<std::setprecision(9)<<std::left<<","<<std::setw(13)<<std::right<<BrokenHistorgram[i]/Normalization<<" ]\n";

    o<< std::endl;

    o<< MeasuredOperators.TotalEnergy() << "\n";
    o<< MeasuredOperators.Potential() << "\n";
    o<< MeasuredOperators.Kinetic() << "\n\n";

    for(int i=0; i<MeasuredOperators.size(); ++i)
      o<<MeasuredOperators.Quantity(i)<<std::endl;


    if(outputConfiguration!="")
      Container.write(outputConfiguration);


  }

};



}

int main(int argc, char *argv[]) {

  SGF::InitializeEnvironment(argc,argv);

  if(argc!=2) {
    std::cerr<<"No input file"<<std::endl;
    exit(1);
  }

  std::string json=readInputFile(argv[1]);

  SGF::Parameters p(json);

  SGF::Simulation S(p);
  S.run();
  p.print(std::cout);
  S.Results(std::cout);
  SGF::FinalizeEnvironment();

}
