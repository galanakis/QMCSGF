#include <fstream>
#include <ctime>

#include "OperatorString.hh"
#include "Measurable.hh"
#include "Utilities.hh"

// Maximum number of broken lines is only used for the BrokenLineHistogram.
#define MAXNUMBROKENLINES 100

class Simulation {
public:

  typedef  std::vector<unsigned long> BrokenHistogramType;
  BrokenHistogramType BrokenHistorgram;

  unsigned long NumWarmUpdates,NumMeasUpdates;
  unsigned long NumDirectedWarmUpdates,NumDirectedMeasUpdates;
  double ActualWarmTime,ActualMeasTime;
  std::string SimulName;


  Simulation(const std::string &Name) : SimulName(Name) {

    NumWarmUpdates=0;
    NumDirectedWarmUpdates=0;
    NumMeasUpdates=0;
    NumDirectedMeasUpdates=0;
    ActualWarmTime=0;
    ActualMeasTime=0;

    BrokenHistorgram.resize(MAXNUMBROKENLINES);

  }


  void Thermalize(SGF::OperatorStringType &OpString,unsigned long WarmIterations,unsigned long WarmTime) {

    NumWarmUpdates=0;
    NumDirectedWarmUpdates=0;

    ProgressBar pbar(SimulName+std::string(": Thermalizing "),&NumWarmUpdates,WarmIterations,WarmTime);
    Timer time(&NumWarmUpdates,WarmIterations,WarmTime);

    do {

      do {
        NumWarmUpdates+=OpString.directed_update();
        ++NumDirectedWarmUpdates;
      } while(OpString.NBrokenLines()!=0);

      pbar.Update();

    } while ( time.Continue() );

    ActualWarmTime=time.ElapsedSeconds();

#ifdef USEMPI

    unsigned long send_NumWarmUpdates(NumWarmUpdates);
    unsigned long send_NumDirectedWarmUpdates(NumDirectedWarmUpdates);
    double send_ActualWarmTime(ActualWarmTime);

    MPI_Reduce(&send_NumWarmUpdates,        &NumWarmUpdates,        1,MPI_UNSIGNED_LONG,MPI_SUM,Master,MPI_COMM_WORLD);
    MPI_Reduce(&send_NumDirectedWarmUpdates,&NumDirectedWarmUpdates,1,MPI_UNSIGNED_LONG,MPI_SUM,Master,MPI_COMM_WORLD);
    MPI_Reduce(&send_ActualWarmTime,        &ActualWarmTime,        1,MPI_DOUBLE,       MPI_SUM,Master,MPI_COMM_WORLD);

#endif

  }

  void Measure(SGF::OperatorStringType &OpString,SGF::Measurable &MeasuredOp,unsigned long NumBins,unsigned long MeasIterations,unsigned long MeasTime) {

    ProgressBar pbar(SimulName+std::string(": Measuring    "),&NumMeasUpdates,MeasIterations,MeasTime);

    NumMeasUpdates=0;
    NumDirectedMeasUpdates=0;


    for (unsigned int i=0; i<NumBins; ++i) {

      Timer time(&NumMeasUpdates,MeasIterations/NumBins,MeasTime/NumBins);

      do {

        do {
          NumMeasUpdates+=OpString.directed_update();   // Perform an update.
          ++NumDirectedMeasUpdates;
          MeasuredOp.measure(OpString);          // Perform measurements.
          BrokenHistorgram[OpString.NBrokenLines()]+=1;

        } while(OpString.NBrokenLines()!=0);

        pbar.Update();

      } while ( time.Continue() );

      MeasuredOp.flush();                               // Bin the data.

    }

    ActualMeasTime=pbar.ElapsedSeconds();

#ifdef USEMPI

    unsigned long send_NumMeasUpdates(NumMeasUpdates);
    unsigned long send_NumDirectedMeasUpdates(NumDirectedMeasUpdates);
    double send_ActualMeasTime(ActualMeasTime);

    MPI_Reduce(&send_NumMeasUpdates,        &NumMeasUpdates,        1,MPI_UNSIGNED_LONG,MPI_SUM,Master,MPI_COMM_WORLD);
    MPI_Reduce(&send_NumDirectedMeasUpdates,&NumDirectedMeasUpdates,1,MPI_UNSIGNED_LONG,MPI_SUM,Master,MPI_COMM_WORLD);
    MPI_Reduce(&send_ActualMeasTime,        &ActualMeasTime,        1,MPI_DOUBLE,       MPI_SUM,Master,MPI_COMM_WORLD);

#endif


  }

  void Results(SGF::Measurable &MeasuredOp) {

    std::cout<< "*************************\n";
    std::cout<< "* Results of simulation *\n";
    std::cout<< "*************************\n\n";
    std::cout<< "  ******************************\n";
    std::cout<< "  * Operator string statistics *\n";
    std::cout<< "  ******************************\n\n";


    std::cout<< "    == Thermalization ==\n\n";
    std::cout<<std::setw(60)<< "    Number of creations/annihilations:" << NumWarmUpdates << std::endl;
    std::cout<<std::setw(60)<< "    Number of creations/annihilations per second per node:" << NumWarmUpdates/ActualWarmTime << "\n";
    std::cout<<std::setw(60)<< "    Number of directed updates:" << NumDirectedWarmUpdates << "\n";
    std::cout<<std::setw(60)<< "    Directed update length:"<< double(NumWarmUpdates)/NumDirectedWarmUpdates<<std::endl;
    std::cout<< "\n    == Measurements   ==\n\n";
    std::cout<<std::setw(60)<< "    Number of creations/annihilations:" << NumMeasUpdates << "\n";
    std::cout<<std::setw(60)<< "    Number of creations/annihilations per second per node:" << NumMeasUpdates/ActualMeasTime << "\n";
    std::cout<<std::setw(60)<< "    Number of directed updates:" << NumDirectedMeasUpdates << "\n";
    std::cout<<std::setw(60)<< "    Directed update length:"<< double(NumMeasUpdates)/NumDirectedMeasUpdates<<std::endl;

#ifdef USEMPI
    std::cout<<std::setw(60)<< "    Number of Processors used "<<NumProcessors<<"\n\n";
#endif

    std::cout<<std::setw(60)<< "    Number of measurements: " << BrokenHistorgram[0] <<std::endl;
    std::cout<<std::setw(60)<<std::fixed<< "    Number of updates per measurements: " << double(NumMeasUpdates)/BrokenHistorgram[0] << "\n\n";
    std::cout<<::std::endl;
    std::cout<<"  *******************************\n";
    std::cout<<"  * Broken worldlines histogram *\n";
    std::cout<<"  *******************************\n\n";
    std::cout<<"    "<<std::setw(10)<<std::left<<"N lines"<<std::setw(10)<<"Count"<<std::setw(10)<<"Probability\n\n";


    double Normalization=0;
    for(int i=0; i<BrokenHistorgram.size(); ++i)
      Normalization+=BrokenHistorgram[i];

    for(int i=0; i<BrokenHistorgram.size(); ++i)
      if(BrokenHistorgram[i]!=0)
        std::cout<<"    "<<std::fixed<<std::setprecision(9)<<std::setw(10)<<std::left<<i<<std::setw(10)<<BrokenHistorgram[i]<<std::setw(10)<<std::right<<BrokenHistorgram[i]/Normalization<<std::endl;

    std::cout<< std::endl;

    std::cout<< "  ***********************************************************************************\n";
    std::cout<< "  * Energies (obtained from operator string length and Green operator state energy) *\n";
    std::cout<< "  ***********************************************************************************\n\n";
    std::cout<< "    " << MeasuredOp.TotalEnergy() << "\n";
    std::cout<< "    " << MeasuredOp.Potential() << "\n";
    std::cout<< "    " << MeasuredOp.Kinetic() << "\n\n";
    std::cout<< "  ******************************\n";
    std::cout<< "  * User's defined measurables *\n";
    std::cout<< "  ******************************\n\n";

    for(int i=0; i<MeasuredOp.size(); ++i)
      std::cout<<"    "<<MeasuredOp.Quantity(i)<<std::endl;

  }
};

