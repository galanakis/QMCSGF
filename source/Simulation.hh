#include <fstream>
#include <ctime>

#include "OperatorString.hh"
#include "Measurable.hh"
#include "Utilities.hh"


class Simulation {
public:


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
               MeasuredOp.measure();          // Perform measurements.
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




      std::cout<< "    == Thermalization ==\n";
      std::cout<< "    Number of creations/annihilations: \t" << NumWarmUpdates << "\t("<<ActualWarmTime*1000000000/NumWarmUpdates << " seconds per billion updates per node)\n";
      std::cout<< "    Number of directed updates:       \t" << NumDirectedWarmUpdates << "\t("<<ActualWarmTime*1000000000/NumDirectedWarmUpdates << " seconds per billion updates per node)\n";
      std::cout<< "    Directed update length:           \t"<< double(NumWarmUpdates)/NumDirectedWarmUpdates<<std::endl;
      std::cout<< "    == Measurements   ==\n";
      std::cout<< "    Number of creations/annihilations: \t" << NumMeasUpdates << "\t("<<ActualMeasTime*1000000000/NumMeasUpdates << " seconds per billion updates per node)\n";
      std::cout<< "    Number of directed updates:       \t" << NumDirectedMeasUpdates << "\t("<<ActualMeasTime*1000000000/NumDirectedMeasUpdates << " seconds per billion updates per node)\n";
      std::cout<< "    Directed update length:           \t"<< double(NumMeasUpdates)/NumDirectedMeasUpdates<<std::endl;

#ifdef USEMPI
      std::cout<< "    Number of Processors used "<<NumProcessors<<"\n\n";
#endif

      const std::vector<unsigned long> &BrokenHistorgram=MeasuredOp.Histogram();

      std::cout<< "    Number of measurements: " << BrokenHistorgram[0] << "\n\n";
      std::cout<<::std::endl;
      std::cout<<"  *******************************\n";
      std::cout<<"  * Broken worldlines histogram *\n";
      std::cout<<"  *******************************\n\n";
      std::cout<<"    N lines\tCount\tProbability\n\n";


      double Normalization=0;
      for(int i=0; i<BrokenHistorgram.size(); ++i)
         Normalization+=BrokenHistorgram[i];

      for(int i=0; i<BrokenHistorgram.size(); ++i)
         if(BrokenHistorgram[i]!=0)
            std::cout<<"    "<<i<<"\t\t"<<BrokenHistorgram[i]<<"\t"<<BrokenHistorgram[i]/Normalization<<std::endl;

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

