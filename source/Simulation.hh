/******************************************************************\

Copyright 2015 Dimitrios Galanakis

This file is part of QMCSGF

QMCSGF is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 3.

QMCSGF is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with QMCSGF.  If not, see <http://www.gnu.org/licenses/>.

\*******************************************************************/

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

    SGF::BrokenLines BrokenLineTracer(OpString);  // It traces the list of broken lines



    for (unsigned int i=0; i<NumBins; ++i) {

      Timer time(&NumMeasUpdates,MeasIterations/NumBins,MeasTime/NumBins);

      do {

        do {
          NumMeasUpdates+=OpString.directed_update();   // Perform an update.
          ++NumDirectedMeasUpdates;
          MeasuredOp.measure(OpString,BrokenLineTracer());          // Perform measurements.
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

    std::cout<<"\nResults:\n";

    unsigned int w=30;
    unsigned int w2=13;
    std::cout<<"\n";
    std::cout<<"  Thermalization:\n";
    std::cout<<"    Iterations:                     "<<NumWarmUpdates<<"\n";
    std::cout<<"    Time:                           "<<ActualWarmTime<<std::endl;
    std::cout<<"    Iterations per sec:             "<<std::setprecision(2)<<std::fixed<<NumWarmUpdates/ActualWarmTime<<"\n";
    std::cout<<"    Directed Updates:               "<<NumDirectedWarmUpdates<<"\n";
    std::cout<<"    Update Length:                  "<<double(NumWarmUpdates)/NumDirectedWarmUpdates<<"\n";
    std::cout<<"\n";
    std::cout<<"  Measurements:\n";
    std::cout<<"    Iterations:                     "<<NumMeasUpdates<<"\n";
    std::cout<<"    Time:                           "<<ActualMeasTime<<"\n";
    std::cout<<"    Iterations per sec:             "<<std::setprecision(2)<<std::fixed<<NumMeasUpdates/ActualMeasTime<<"\n";
    std::cout<<"    Directed Updates:               "<<NumDirectedMeasUpdates<<"\n";
    std::cout<<"    Update Length:                  "<<double(NumMeasUpdates)/NumDirectedMeasUpdates<<"\n";
    std::cout<<"    Measurement Count:              "<< BrokenHistorgram[0]<<"\n";
    std::cout<<"    Iterations per measurement:     "<<std::fixed<< double(NumMeasUpdates)/BrokenHistorgram[0]<<"\n";

#ifdef USEMPI
    std::cout<<"    Number of Processors:           "<<NumProcessors<<"\n";
#endif

    std::cout<<"\n";
    std::cout<<"  Broken worldlines:\n\n";

    double Normalization=0;
    for(int i=0; i<BrokenHistorgram.size(); ++i)
      Normalization+=BrokenHistorgram[i];

    for(int i=0; i<BrokenHistorgram.size(); ++i)
      if(BrokenHistorgram[i]!=0)
        std::cout<<"    - [ "<<std::right<<std::fixed<<std::setw(3)<<i<<","<<std::setw(12)<<BrokenHistorgram[i]<<std::setprecision(9)<<std::left<<","<<std::setw(13)<<std::right<<BrokenHistorgram[i]/Normalization<<" ]\n";

    std::cout<< std::endl;

    std::cout<< MeasuredOp.TotalEnergy() << "\n";
    std::cout<< MeasuredOp.Potential() << "\n";
    std::cout<< MeasuredOp.Kinetic() << "\n\n";

    for(int i=0; i<MeasuredOp.size(); ++i)
      std::cout<<MeasuredOp.Quantity(i)<<std::endl;

  }
};

