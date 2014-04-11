#include "Parameters.hh"
#include <iostream>

#include <fstream>
#include <ctime>

#include "OperatorString.hh"
#include "Measurable.hh"
#include "Utilities.hh"

namespace SGF {


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



class Simulation {
public:

  SimCounters Warm, Meas;

  unsigned long NumBins;
  int Seed;

  std::string SimulName;
  std::string inputfname;

  SGFBase Container;
  Measurable MeasuredOperators;

  std::string outputConfiguration;

  unsigned long Measurements;

  Simulation(const SGF::Parameters& p) {

    Measurements = 0;

    int Seed = p.Seed;

    RNG::Initialize(Seed);

    inputfname=p.inputfname;

    Warm.Time = p.WarmTime;
    Warm.Iterations = p.WarmIterations;

    Meas.Time = p.MeasTime;
    Meas.Iterations = p.MeasIterations;

    NumBins = p.NBins;

    SimulName = p.inputfname;

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

    OpString2.insert_update(&BrokenLineTracer);

    for (unsigned int i = 0; i < NumBins; ++i) {

      Timer time(&Meas.Updates, Meas.Iterations / NumBins, Meas.Time / NumBins);

      do {

        do {
          Meas.increment(OpString2.directed_update());
          MeasuredOperators.measure(OpString2, BrokenLineTracer());         // Perform measurements.

        } while (OpString2.NBrokenLines() != 0);

        ++Measurements;

        measProgress.Update();

      } while ( time.Continue() );

      MeasuredOperators.flush();                               // Bin the data.
      Meas.ActualTime += time.ElapsedSeconds();
    }


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
    o << "\n";
    o << "  Thermalization:\n";
    Warm.print(o, 2);
    o << "\n";
    o << "  Measurements:\n";
    Meas.print(o, 2);
    o << "    Measurement Count:              " << Measurements << "\n";
    o << "    Iterations per measurement:     " << std::fixed << double(Meas.Updates) / Measurements << "\n";

#ifdef USEMPI
    o << "    Number of Processors:           " << NumProcessors << "\n";
#endif


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

  SGF::Parameters p(argv[1]);

  SGF::Simulation S(p);
  S.run();
  p.print(std::cout);
  S.Results(std::cout);
  SGF::FinalizeEnvironment();

}
