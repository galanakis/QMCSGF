#ifndef __UTILITIES__
#define __UTILITIES__

#include <ctime>
#include <string>
#include <iostream>
#include <fstream>

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

   inline double ProgressIter() const {
      return static_cast<double>(*counter-StartIter)/MaxIter;
   }
   inline double ProgressTime() const {
      return (ignore_time) ? 0 : static_cast<double>(clock()-StartTime)/(EndTime-StartTime);
   }
   inline double Progress() const {
      return Max(ProgressIter(),ProgressTime());
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
   static const int StringLength=500;
   char Status[StringLength],*Ptr;
   double Progress;
   unsigned long NumUpdates;
   clock_t Time;
public:
   ProgressBar(const std::string &Name,unsigned long *c,unsigned long Iter,unsigned int TotalTime) : Timer(c,Iter,TotalTime), Progress(0), NumUpdates(0) {
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

         int percent=static_cast<int>(100*Progress);



#ifdef CMDLINEPROGRESS
         reset_cout();
         sprintf(Ptr," - %.3d %% - %.3dh %.2dm %.2ds - %6d updates per second",percent,HoursLeft,MinutesLeft,SecondsLeft,Speed);
         std::cerr<<Status<<std::flush;
#else
         char OldStatus[StringLength];
         strcpy(OldStatus,Status);
         sprintf(Ptr," - %.3d %% - %.3dh %.2dm %.2ds - %6d updates per second",percent,HoursLeft,MinutesLeft,SecondsLeft,Speed);
         rename(OldStatus,Status);
#endif


      }
   }

};

#endif
