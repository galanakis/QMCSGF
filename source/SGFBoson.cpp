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

  struct SimCounters {
    unsigned long Time;
    double ActualTime;
    unsigned long Iterations;
    unsigned long Updates;
    unsigned long DirectedUpdates;
    SimCounters() : Time(0), ActualTime(0), Iterations(0), Updates(0), DirectedUpdates(0) {}

    std::ostream& print(std::ostream& o, unsigned int depth) {
      std::string indent(2 * depth, ' ');
      reduce();
      o << indent << "Iterations:                     " << Updates << "\n";
      o << indent << "Time:                           " << ActualTime << std::endl;
      o << indent << "Iterations per sec:             " << std::setprecision(2) << std::fixed << Updates / ActualTime << "\n";
      o << indent << "Directed Updates:               " << DirectedUpdates << "\n";
      o << indent << "Update Length:                  " << double(Updates) / DirectedUpdates << "\n";

      return o;
    }

    void reduce() {
#ifdef USEMPI
      unsigned long send_Updates(Updates);
      unsigned long send_DirectedUpdates(DirectedUpdates);
      double send_ActualTime(ActualTime);

      MPI_Reduce(&send_Updates,          &Updates,         1, MPI_UNSIGNED_LONG, MPI_SUM, Master, MPI_COMM_WORLD);
      MPI_Reduce(&send_DirectedUpdates,  &DirectedUpdates, 1, MPI_UNSIGNED_LONG, MPI_SUM, Master, MPI_COMM_WORLD);
      MPI_Reduce(&send_ActualTime,       &ActualTime,      1, MPI_DOUBLE,        MPI_SUM, Master, MPI_COMM_WORLD);
#endif
    }


    inline void increment(unsigned int length) {
      Updates += length;
      ++DirectedUpdates;
    }

  };

  typedef  std::vector<unsigned long> BrokenHistogramType;
  BrokenHistogramType BrokenHistorgram;

  SimCounters Warm, Meas;

  unsigned long NumBins;
  int Seed;

  std::string SimulName;

  SGFBase Container;
  Measurable MeasuredOperators;

  std::string outputConfiguration;


  Simulation(const SGF::Parameters& p) {


    int Seed = p.Seed;

    RNG::Initialize(Seed);

    Warm.Time = p.WarmTime;
    Warm.Iterations = p.WarmIterations;

    Meas.Time = p.MeasTime;
    Meas.Iterations = p.MeasIterations;

    NumBins = p.NBins;

    SimulName = p.modelID;

    BrokenHistorgram.resize(MAXNUMBROKENLINES);

    p.MakeContainer(Container);
    p.MakeMeasurables(Container, MeasuredOperators);

    outputConfiguration = p.outputConfiguration;

  }

  template<class OperatorStringType>
  void run() {

    OperatorStringType OpString(Container);

    ProgressBar warmProgress(SimulName + std::string(": Thermalizing "), &Warm.Updates, Warm.Iterations, Warm.Time);
    Timer time(&Warm.Updates, Warm.Iterations, Warm.Time);

    do {

      do {
        Warm.increment(OpString.directed_update());
      } while (OpString.NBrokenLines() != 0);

      warmProgress.Update();
    } while ( time.Continue() );

    Warm.ActualTime = time.ElapsedSeconds();

    warmProgress.reset_cout();
    ProgressBar measProgress(SimulName + std::string(": Measuring    "), &Meas.Updates, Meas.Iterations, Meas.Time);

    OperatorStringType OpString2(Container);

    SGF::BrokenLines BrokenLineTracer(OpString2);  // It traces the list of broken lines



    for (unsigned int i = 0; i < NumBins; ++i) {

      Timer time(&Meas.Updates, Meas.Iterations / NumBins, Meas.Time / NumBins);

      do {

        do {
          Meas.increment(OpString2.directed_update());
          MeasuredOperators.measure(OpString2, BrokenLineTracer());         // Perform measurements.
          BrokenHistorgram[OpString2.NBrokenLines()] += 1;

        } while (OpString2.NBrokenLines() != 0);

        measProgress.Update();

      } while ( time.Continue() );

      MeasuredOperators.flush();                               // Bin the data.

    }

    Meas.ActualTime = measProgress.ElapsedSeconds();

    measProgress.reset_cout();

  }


  void run() {

    if (Container.Ensemble == Canonical)
      run<CanonicalOperatorString>();
    else
      run<GrandOperatorString>();


  }

  void Results(std::ostream& o) {



    o << "\nResults:\n";

    unsigned int w = 30;
    unsigned int w2 = 13;
    o << "\n";
    o << "  Thermalization:\n";
    Warm.print(o, 2);
    o << "\n";
    o << "  Measurements:\n";
    Meas.print(o, 2);
    o << "    Measurement Count:              " << BrokenHistorgram[0] << "\n";
    o << "    Iterations per measurement:     " << std::fixed << double(Meas.Updates) / BrokenHistorgram[0] << "\n";

#ifdef USEMPI
    o << "    Number of Processors:           " << NumProcessors << "\n";
#endif

    o << "\n";
    o << "  Broken worldlines:\n\n";

    double Normalization = 0;
    for (int i = 0; i < BrokenHistorgram.size(); ++i)
      Normalization += BrokenHistorgram[i];

    for (int i = 0; i < BrokenHistorgram.size(); ++i)
      if (BrokenHistorgram[i] != 0)
        o << "    - [ " << std::right << std::fixed << std::setw(3) << i << "," << std::setw(12) << BrokenHistorgram[i] << std::setprecision(9) << std::left << "," << std::setw(13) << std::right << BrokenHistorgram[i] / Normalization << " ]\n";

    o << std::endl;

    MeasuredOperators.print(o);


    if (outputConfiguration != "")
      Container.write(outputConfiguration);


  }

};



}

int main(int argc, char* argv[]) {

  SGF::InitializeEnvironment(argc, argv);

  if (argc != 2) {
    std::cerr << "No input file" << std::endl;
    exit(1);
  }

  std::string json = readInputFile(argv[1]);

  SGF::Parameters p(json);

  SGF::Simulation S(p);
  S.run();
  p.print(std::cout);
  S.Results(std::cout);
  SGF::FinalizeEnvironment();

}
