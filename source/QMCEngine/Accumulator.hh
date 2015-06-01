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

#ifndef __ACCUMULATOR__
#define __ACCUMULATOR__

namespace SGF {


/*

  class Accumulator

  It takes in a stream of numbers and calculates their mean and standard deviation.

  ======== Usage ==========
  Accumulator<double> acc; // a binned accumulator of doubles
  acc.push(data);                // push a value (is meaned to be called frequently).
  double average=acc.average();  // get the average of the bins plus the constant value.
  double sigma=acc.sigma();      // get the standard deviation.

*/


template<class T>
class Accumulator {
  T n;      // Number of data
  T mean;   // mean
  T Moment; // second moment
public:
  Accumulator() : n(0), mean(0), Moment(0) {}
  inline void push(const T& x) {
    ++n;
    T delta = x - mean;
    T R = delta / n;
    mean = mean + R;
    Moment = Moment + delta * (x - mean);
  }
  inline T sigma() const {
    return sqrt(Moment / (n - 1));
  }
  inline T average() const {
    return mean;
  }

};

}

#endif
