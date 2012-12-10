#include <fstream>
#include <ctime>

#include "OperatorString.hh"
#include "Measurable.hh"

template<class T> T Min(const T &a,const T &b) {
   return (a<b)?a:b;
}
template<class T> T Max(const T &a,const T &b) {
   return (a>b)?a:b;
}


/*

   class Timer

   Timer(&counter, MaxIter, Time)

   The class provides the method Continue()
   which checks to see if the counter is still smaller than MaxIter
   and also the less time than "Time" has elapsed. In this case
   it returns true.

   Also it provides Progress which returns a number between 0 and 1.

*/

class Timer {
   clock_t StartTime;
   clock_t EndTime;
   unsigned long StartIter;
   unsigned long MaxIter;
   bool ignore_time;
protected:
   unsigned long *counter;
public:
   Timer(unsigned long *c,unsigned long Iter,unsigned int Time) : counter(c) {

      StartTime=clock();
      EndTime=StartTime+Time*CLOCKS_PER_SEC;
      StartIter=*counter;
      MaxIter=StartIter+Iter;

      // If requested time is more than a year ignore time.
      ignore_time=false;
      if((EndTime-StartTime)/CLOCKS_PER_SEC>31556940)
         ignore_time=true;

   }
   inline double Progress() const {
      return Max(static_cast<double>(*counter-StartIter)/MaxIter,static_cast<double>(clock()-StartTime)/(EndTime-StartTime));
   }
   inline bool Continue() {
      return *counter<MaxIter && (ignore_time || clock()<EndTime);
   }
   double ElapsedSeconds() {
      return double(clock()-StartTime)/CLOCKS_PER_SEC;
   }
   double Iterations() {
      return *counter-StartIter;
   }
   virtual ~Timer() {}
};


class ProgressBar : public Timer {
   static const int StringLength=300;
   char Status[StringLength],*Ptr;
   double Progress;
   unsigned long NumUpdates;
   clock_t Time;
public:
   ProgressBar(const std::string &Name,unsigned long *c,unsigned long Iter,unsigned int Time) : Timer(c,Iter,Time), Progress(0), NumUpdates(0) {
      Time=clock();
      strcpy(Status,Name.c_str());
      Ptr=Status+strlen(Status);
      sprintf(Ptr," - 000 %% - 000h 00m 00s -       0 updates per second");
#ifdef CMDLINEPROGRESS
      std::cerr<<Status<<std::flush;
#else
      std::fstream File;
      File.open(Status,std::ios::out);
      File.close();
#endif

   }
   ~ProgressBar() {
#ifdef CMDLINEPROGRESS
      reset_cout();
#else
      remove(Status);
#endif
   }

   inline void reset_cout() {
      for(unsigned int i=0; i<strlen(Status); ++i)
         std::cerr<<'\b';
   }
   inline void Update() {
      unsigned long NewNumUpdates=*counter;
      double prog=Timer::Progress();
      if(prog<0) prog=0.0;
      if(prog>1) prog=1.0;
      if(static_cast<int>(100*(prog-Progress))!=0) {
         clock_t Now=clock();
         double DeltaSeconds=double(Now-Time)/CLOCKS_PER_SEC;
         int Speed=static_cast<unsigned int>((NewNumUpdates-NumUpdates)/DeltaSeconds);
         int SecondsLeft=(1-prog)*DeltaSeconds/(prog-Progress);
         Time=Now;
         Progress=prog;
         NumUpdates=NewNumUpdates;

         int HoursLeft=SecondsLeft/3600;
         SecondsLeft%=3600;
         int MinutesLeft=SecondsLeft/60;
         SecondsLeft%=60;

         char OldStatus[StringLength];
         strcpy(OldStatus,Status);
         int percent=static_cast<int>(100*Progress);


         sprintf(Ptr," - %.3d %% - %.3dh %.2dm %.2ds - %6d updates per second",percent,HoursLeft,MinutesLeft,SecondsLeft,Speed);

#ifdef CMDLINEPROGRESS
         reset_cout();
         std::cerr<<Status<<std::flush;
#else
         rename(OldStatus,Status);
#endif


      }
   }

};





class Simulation {
   unsigned long NumWarmUpdates,NumMeasUpdates;
   unsigned long NumDirectedWarmUpdates,NumDirectedMeasUpdates;

   double ActualWarmTime,ActualMeasTime;

   std::string SimulName;

public:


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

      std::cout<< "    Number of Processors used "<<NumProcessors<<"\n\n";

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

      std::cout<< "    Number of Processors used "<<NumProcessors<<"\n\n";

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

