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

#ifndef __RANDOMNUMBERGENERATOR__
#define __RANDOMNUMBERGENERATOR__


#include <iostream>
#include <cmath>
#include <cstdlib>


/*

* Note this documentation is not updated
* The ExponentialHelper returns -log(1-r*(1-exp(-L)))/L
* but the exponential returns log(1-r*(1-exp(L)))/L (so it differs in the sign).

  It returns a number x in the interval [0,1) which includes zero and excludes 1.
  The probability distribution on this interval is exp(-Lx).
  Formally it returns
  -log(1-r*(1-exp(-L)))/L
  where is a random number chosen uniformly from the interval [0,1).
  The algorithm is based on bisecting the interval [0,1)
  repeatedly until it's width is smaller than a cutoff
  and the use a fast polynomial approximation.


  At first the interval is bisected repeatedly. In the first iteration
  the relative weigth of the leftmost interval (the one that contains zero)
  is 1 and the weight of the rightmost is exp(-L/2). If (1+exp(-L/2))*r<1
  then the left interval is selected. The index i is the index of the
  interval and takes values between 0 and 2^k-1 where k is the number
  of bisections. The index n=2^k. The bisections continue until the
  width L<T.

  At that point the interval which starts at i/n has been selected
  and we need to proceed by selecting one number in this interval.
  Next we need to pick a number in the interval [0,L]
  We can use the formula X=-log(1-r(1-exp(-L)))/L, but we can Taylor expand it in a smart way.
  Let u=r(1-exp(-L)) and s=u/(2-u) such that 1-u=(1-s)/(1+s)
  Then
  -log(1-u)=log(1+s)-log(1-s)=2*s(1+s^2/3+s^4/5...)=2*s+s*R
  where R has only even powers. Then I follow http://www.netlib.org/fdlibm/e_log.c
  where R is approximated by a 14 order polynomial such that in the s interval
  [0,0.1716] the accuracy is better than 2**-58.45. This corresponds to
  L<0.346629.

  Finally if we let a=(1-exp(-L))/L;
  X=r*a/L(2+R)/(2-u)
  where all the quantities are finite for small values of L.

  To calculate "a" without floating point errors we use the identity
  a=2*exp(-L/2)sinh(L/2)/L
  We then only need to check if L==0.


  Note. In the above derivation we assume that L>0.
  In the case of L<0 we have the distribution
  -log(1-r*(1-exp(-L)))/L = 1+log(1-(1-r)(1-exp(L)))/L
  = 1-log(1-(1-r)(1-exp(-|L|)))/|L|
  So all we need to do is flip the result and the
  uniform random number generator.


*/

template<typename RNGBase>
class _RNG : public RNGBase {

  // This variables are used by Exponential(double)
  static const double Lg1 = 6.666666666666735130e-01;
  static const double Lg2 = 3.999999999940941908e-01;
  static const double Lg3 = 2.857142874366239149e-01;
  static const double Lg4 = 2.222219843214978396e-01;
  static const double Lg5 = 1.818357216161805012e-01;
  static const double Lg6 = 1.531383769920937332e-01;
  static const double Lg7 = 1.479819860511658591e-01;

  /*
    In this function r is the uniform random generator
    and 0<L<0.346629.
    It calculates accurately the expression
    -log(1-r*(1-exp(-L)))/L
  */
  static inline double ExponentialHelper(double r,double L) {
    double u=r*(1-exp(L));
    double s=u/(2-u);
    double ss=s*s;
    double R=((((((ss*Lg7+Lg6)*ss+Lg5)*ss+Lg4)*ss+Lg3)*ss+Lg2)*ss+Lg1)*ss;
    double a= (L==0) ? exp(L/2) : 2*exp(L/2)*sinh(L/2)/L;
    return r*a*(2+R)/(2-u);
  }

public:
// Initialize the random number generator with a non zero seed.
  static void Initialize(int S) {RNGBase::_Initialize(S);};
// Returns a uniformly distributed random number in the interval [0,1[
  inline static double Uniform() {return RNGBase::Uniform();};
// Returns a uniformly distributed random number in the interval ]0,1[
  inline static double NZUniform() {return RNGBase::NZUniform();};
// Returns one integer between 0 and N-1
  static inline int UniformInteger(unsigned int N) {return Uniform()*N;}
// Returns a number in the interval [0,1[ with distribution exp(A x)
  static inline double Exponential(double L) {
  // Do not dare to touch this before you read the documentation above.
    const double T=0.346629;
    unsigned long n=1;
    unsigned long i=0;

    while(fabs(L)>=T) {
      L/=2;
      i<<=1;
      double r=Uniform();
      if(exp(-L)<=r/(1-r))
        i+=1;
      n<<=1;
    }
    double r=Uniform();
    if(L>0)
      return (i+ExponentialHelper(r,L))/n;
    else
      return (i+1-ExponentialHelper(1-r,-L))/n;
  }

// Returns a number in the interval ]0,1[ with distribution exp(A x)
  static inline double NZExponential(double A) {
    double result=0;
    do { result=Exponential(A); }
    while(result==0);
    return result;
  }

};


#if defined(RNG_STLMT)

// STL Mersenne Twister
#include <random>
class STLRNGMT {
  static std::mt19937 rand;
public:
  inline static void _Initialize(int s) {
    rand.seed(s);
  }
  inline static double Uniform() {return double(rand())/(rand.max()+1.0);}
  inline static double NZUniform() {return (double(rand()) + 0.5 )/(rand.max()+1.0);}
};

std::mt19937 STLRNGMT::rand;
typedef _RNG<STLRNGMT> RNG;

#elif defined(RNG_DSFMT)

#include "dSFMT.c"

class RNGDFSMT {
  static dsfmt_t dsfmt;
public:
  inline static void _Initialize(int s) {
    dsfmt_init_gen_rand(&dsfmt, s);
  }
  inline static double Uniform() {return dsfmt_genrand_close_open(&dsfmt);}
  inline static double NZUniform() {return dsfmt_genrand_open_open(&dsfmt);}
};

dsfmt_t RNGDFSMT::dsfmt;
typedef _RNG<RNGDFSMT> RNG;

#elif defined(RNG_WELL)

#include "WELL44497.h"

class RNGWELL {
  static well44497 rand;
public:
  inline static void _Initialize(int s) {
    unsigned int seed[well44497::size];
    for(int i=0;i<well44497::size;++i)
      seed[i]=(s*=16807)%2147483647;
    rand.init(seed);
  }
  inline static double Uniform() {return rand()/4294967296.0;}
  inline static double NZUniform() {return (double(rand())+0.5)/4294967296.0;}
};

well44497 RNGWELL::rand;
typedef _RNG<RNGWELL> RNG;


#elif defined (RNG_TINYMT32)

#include "tinymt32.c"
class RNGTINYMT32 {
  static tinymt32_t tinymt;
public:
  inline static void _Initialize(int s) { tinymt32_init(&tinymt, s); }
  inline static double Uniform() {return tinymt32_generate_uint32(&tinymt)/4294967296.0;}
  inline static double NZUniform() {return double(tinymt32_generate_uint32(&tinymt)+0.5)/4294967296.0;}

};

tinymt32_t RNGTINYMT32::tinymt;
typedef _RNG<RNGTINYMT32> RNG;

#elif defined (RNG_TINYMT64)

#include "tinymt64.c"
class RNGTINYMT64 {
  static tinymt64_t tinymt;
public:
  inline static void _Initialize(int s) { tinymt64_init(&tinymt, s); }
  inline static double Uniform() {return tinymt64_generate_uint64(&tinymt)/18446744073709551616.0;}
  inline static double NZUniform() {return double(tinymt64_generate_uint64(&tinymt)+0.5)/18446744073709551616.0;}

};

tinymt64_t RNGTINYMT64::tinymt;
typedef _RNG<RNGTINYMT64> RNG;

#elif defined (RNG_LC)

// Linear Congruence method
class RNGLC {
  static std::minstd_rand0 rand;
public:
  inline static void _Initialize(int s) {
    rand.seed(s);
  }
// Returns a uniformly distributed random number in the interval [0,1[
  inline static double Uniform() { return double(rand()-1)/rand.max(); }
// Returns a uniformly distributed random number in the interval ]0,1[
  inline static double NZUniform() {return (double(rand()) - 0.5 )/(rand.max());}
};

std::minstd_rand0 RNGLC::rand;

typedef _RNG<RNGLC> RNG;

#else

#include "MersenneTwister.h"
class RNGMT {
  static MTRand rand;
public:
  inline static void _Initialize(int s) { rand.seed(s); }
  inline static double Uniform() {return rand.randExc();}
  inline static double NZUniform() {return rand.randDblExc();}
};

MTRand RNGMT::rand;
typedef _RNG<RNGMT> RNG;


#endif







#endif
