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

#ifndef __UTILITIES__
#define __UTILITIES__

#include <ctime>
#include <string>
#include <iostream>
#include <algorithm>


template<typename T>
struct Counter {
  const T start;
  const T iter;
  virtual T count() const = 0;
  Counter(T _start,T _iter) : start(_start), iter(_iter) {}
  inline T elapsed() const { 
    return count()-start;
  }
  inline bool running() const { 
    return elapsed() < iter;
  }
  inline double progress() const {
    return double(elapsed())/iter;
  }
};

struct IterCounter : Counter<unsigned long> {
  unsigned long * const counter;
  unsigned long count() const { 
    return *counter;
  }
  IterCounter(unsigned long *c,unsigned long Iter) : Counter<unsigned long>(*c,Iter), counter(c) {}
};

struct TimeCounter : public Counter<clock_t> {
  clock_t count() const { return clock(); }
  TimeCounter(unsigned int Seconds) : Counter<clock_t>(clock(),Seconds*CLOCKS_PER_SEC) {}
  double ElapsedSeconds() const { 
    return double(elapsed())/CLOCKS_PER_SEC;
  }
};

/*

   class Timer

   Timer(&counter, MaxIter, Time)

   The class provides the method Continue()
   which checks to see if the counter is still smaller than MaxIter
   and also the less time than "Time" has elapsed. In this case
   it returns true.

*/

struct Timer : public IterCounter, public TimeCounter {
  Timer(unsigned long *c,unsigned long Iter,unsigned long Time) : IterCounter(c,Iter), TimeCounter(Time) {}
  inline bool running() const {
    return IterCounter::running() && TimeCounter::running();
  }
};

class ProgressBar : public Timer {
  char Status[500],*Ptr;
  unsigned long Iter;
  clock_t Time;
public:
  ProgressBar(unsigned long *c,unsigned long _iter,unsigned long _time,const std::string &Name="") : Timer(c,_iter,_time), Ptr(Status) {
    Iter=IterCounter::count();
    Time=TimeCounter::count();
    prefix(Name);
  }

  ~ProgressBar() {
    clear();
  }

  inline void clear() {
    std::clog<<'\r'<<std::flush;
  }

  void prefix(const std::string &Name) {
    clear();
    strcpy(Status,Name.c_str());
    Ptr=Status+strlen(Status);
  }

  inline bool running() {
    unsigned long NewIter=IterCounter::count();
    clock_t NewTime=TimeCounter::count();

    // Update the progress bar every second
    if(NewTime>Time+CLOCKS_PER_SEC) {
  
      double progress=std::max(IterCounter::progress(),TimeCounter::progress());
      int sec=ElapsedSeconds()*(1.0-progress)/progress;
      int speed=CLOCKS_PER_SEC*(NewIter-Iter)/(NewTime-Time);

      sprintf(Ptr,"\x1b[31;1m%3.0f%%\x1b[0m  \x1b[33;1m%.3dh %.2dm %.2ds\x1b[36;1m  %7d updates/second\x1b[0m",100*progress,sec/3600,sec%3600/60,sec%60,speed);
      clear();
      std::clog<<Status<<std::flush;

      Time=NewTime;
      Iter=NewIter;

    }

    return Timer::running();
  }

};

#endif
