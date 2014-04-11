#ifndef __RANDOMNUMBERGENERATOR__
#define __RANDOMNUMBERGENERATOR__


#include <iostream>
#include <cmath>
#include <cstdlib>

#ifdef USEMPI
extern int Rank;
#endif


#if defined(CPPRNG_MT)

// C++ Mersenne Twister
#include <random>
class CPPRNGMT {
  static std::mt19937 rand;
public:
  inline static void Initialize(int s) {
//    cout<<"Using the Mersenne Twister random number generator."<<std::endl;
#ifdef USEMPI
    s=s*(Rank+1);
//    std::cout<<"Initialize processor "<<Rank<<" with seed "<<s<<std::endl;
#endif
    rand.seed(s);

//    std::cout<<"Seed:          "<<s<<std::endl;
//    std::cout<<"Minimum value: "<<rand.min()<<std::endl;
//    std::cout<<"Maximum value: "<<rand.max()<<std::endl;

  }
  inline static double Uniform() {return double(rand())/rand.max();}
  inline static double NZUniform() {return (double(rand()) + 0.5 ) * (1.0/rand.max());}
};

std::mt19937 CPPRNGMT::rand;
typedef CPPRNGMT RNGBase;

#elif defined(RNG_DSFMT)

#include "dSFMT.c"

class RNGDFSMT {
  static dsfmt_t dsfmt;
public:
  inline static void Initialize(int s) {
//    std::cout<<"Using the DSFMT random number generator."<<std::endl;
#ifdef USEMPI
    s=s*(Rank+1);
//    std::cout<<"Initialize processor "<<Rank<<" with seed "<<s<<std::endl;
#endif
    dsfmt_init_gen_rand(&dsfmt, s);
  }
  inline static double Uniform() {return dsfmt_genrand_close_open(&dsfmt);}
  inline static double NZUniform() {return dsfmt_genrand_open_open(&dsfmt);}
};

dsfmt_t RNGDFSMT::dsfmt;
typedef RNGDFSMT RNGBase;

#elif defined(RNG_WELL)

#include "WELL44497a_new.c"

class RNGWELL {
  static unsigned int _seed;
public:
  inline static void Initialize(int s) {
//    std::cout<<"Using the DSFMT random number generator."<<std::endl;
#ifdef USEMPI
    s=s*(Rank+1);
//    std::cout<<"Initialize processor "<<Rank<<" with seed "<<s<<std::endl;
#endif
    _seed=s;
    InitWELLRNG44497(&_seed);
  }
  inline static double Uniform() {return WELLRNG44497()/4294967296.0;}
  inline static double NZUniform() {return (double(WELLRNG44497())+0.5)/4294967296.0;}
};

unsigned int RNGWELL::_seed;
typedef RNGWELL RNGBase;


#elif defined(RNG_MT)

#include "MersenneTwister.h"
class RNGMT {
  static MTRand rand;
public:
  inline static void Initialize(int s) {
//    std::cout<<"Using the Mersenne Twister random number generator."<<std::endl;
#ifdef USEMPI
    s=s*(Rank+1);
//    std::cout<<"Initialize processor "<<Rank<<" with seed "<<s<<std::endl;
#endif
    rand.seed(s);
  }
  inline static double Uniform() {return rand.randExc();}
  inline static double NZUniform() {return rand.randDblExc();}
};

MTRand RNGMT::rand;
typedef RNGMT RNGBase;

#elif defined (RNG_TINYMT32)

#include "tinymt32.c"
class RNGTINYMT32 {
  static tinymt32_t tinymt;
public:
  inline static void Initialize(int s) {
//    std::cout<<"Using the tinymt32 Mersenne Twister random number generator."<<std::endl;
#ifdef USEMPI
    s=s*(Rank+1);
//    std::cout<<"Initialize processor "<<Rank<<" with seed "<<s<<std::endl;
#endif
    tinymt32_init(&tinymt, s);
  }
  inline static double Uniform() {return tinymt32_generate_uint32(&tinymt)/4294967296.0;}
  inline static double NZUniform() {return double(tinymt32_generate_uint32(&tinymt)+0.5)/4294967296.0;}

};

tinymt32_t RNGTINYMT32::tinymt;
typedef RNGTINYMT32 RNGBase;

#elif defined (RNG_TINYMT64)

#include "tinymt64.c"
class RNGTINYMT64 {
  static tinymt64_t tinymt;
public:
  inline static void Initialize(int s) {
//    std::cout<<"Using the tinymt64 Mersenne Twister random number generator."<<std::endl;
#ifdef USEMPI
    s=s*(Rank+1);
//    std::cout<<"Initialize processor "<<Rank<<" with seed "<<s<<std::endl;
#endif
    tinymt64_init(&tinymt, s);
  }
  inline static double Uniform() {return tinymt64_generate_uint64(&tinymt)/18446744073709551616.0;}
  inline static double NZUniform() {return double(tinymt64_generate_uint64(&tinymt)+0.5)/18446744073709551616.0;}

};

tinymt64_t RNGTINYMT64::tinymt;
typedef RNGTINYMT64 RNGBase;

#elif defined(RNG_SIMPLE_SPRNG)

#include <sprng.h>

// Sprng: parallel random number generator
class RNGSPRNG {
  inline static void Initialize(int Seed) {
//    cout<<"Using the SPRNG random number generator."<<std::endl;
    init_sprng(SPRNG_CMRG,Seed,SPRNG_DEFAULT);
  }
  inline static double Uniform() {return sprng();}
  // Returns a uniformly distributed random number in the interval ]0,1[
  inline static double NZUniform() {
    double result=0;
    do { result=Uniform(); }
    while(result==0);
    return result;
  }
};

typedef RNGSPRNG RNGBase;

#else

// Linear Congruence method (Val's method)
class RNGLC {
  static int Seed;
public:
// Initialize the random number generator with a non zero seed.
  inline static void Initialize(int S) {
//    cout<<"Using the Linear Congruence random number generator."<<std::endl;
#ifdef USEMPI
    S=S*(Rank+1);
#endif
    Seed=S;
    if(S==0) {
      std::cerr<<"Error in RNG::Initialize: The initialization seed cannot be zero."<<std::endl;
      exit(1);
    }
    for (int i=0; i<1000; ++i)
      Seed*=16807;
  }
// Returns a uniformly distributed random number in the interval [0,1[
  inline static double Uniform() { return (Seed*=16807)/4294967296.0+0.5; }
// Returns a uniformly distributed random number in the interval ]0,1[
  inline static double NZUniform() {
    double result=0;
    do { result=Uniform(); }
    while(result==0);
    return result;
  }
};

int RNGLC::Seed;

typedef RNGLC RNGBase;

#endif


class RNG : public RNGBase {

// This variables are used by Exponential(double)
  static const double Lg1,Lg2,Lg3,Lg4,Lg5,Lg6,Lg7;
  inline static double ExponentialHelper(double r,double L);
public:
// Initialize the random number generator with a non zero seed.
  static void Initialize(int S) {RNGBase::Initialize(S);};
// Returns a uniformly distributed random number in the interval [0,1[
  inline static double Uniform() {return RNGBase::Uniform();};
// Returns a uniformly distributed random number in the interval ]0,1[
  inline static double NZUniform() {return RNGBase::NZUniform();};
// Returns 0 or 1 with relative probabilities a and b respectively
  static inline int CoinFlip(double a,double b) { return (a+b)*Uniform()<b;}
// Returns one integer between 0 and N-1
  static inline int UniformInteger(unsigned int N) {return Uniform()*N;}
// Returns a number with relative distribution exp(-x)
  static inline double Exponential() { return -log(1.0-Uniform()); }
// Returns a number in the interval [0,1[ with distribution exp(A x)
  static inline double Exponential(double A);
// Returns a number in the interval ]0,1[ with distribution exp(A x)
  static inline double NZExponential(double A) {
    double result=0;
    do { result=Exponential(A); }
    while(result==0);
    return result;
  }

};


/*

* Note this documentation is not updated
*	The ExponentialHelper returns -log(1-r*(1-exp(-L)))/L
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


const double RNG::Lg1 = 6.666666666666735130e-01;
const double RNG::Lg2 = 3.999999999940941908e-01;
const double RNG::Lg3 = 2.857142874366239149e-01;
const double RNG::Lg4 = 2.222219843214978396e-01;
const double RNG::Lg5 = 1.818357216161805012e-01;
const double RNG::Lg6 = 1.531383769920937332e-01;
const double RNG::Lg7 = 1.479819860511658591e-01;


/*
	In this function r is the uniform random generator
	and 0<L<0.346629.
	It calculates accurately the expression
	-log(1-r*(1-exp(-L)))/L
*/
inline double RNG::ExponentialHelper(double r,double L) {
  double u=r*(1-exp(L));
  double s=u/(2-u);
  double ss=s*s;
  double R=((((((ss*Lg7+Lg6)*ss+Lg5)*ss+Lg4)*ss+Lg3)*ss+Lg2)*ss+Lg1)*ss;
  double a= (L==0) ? exp(L/2) : 2*exp(L/2)*sinh(L/2)/L;
  return r*a*(2+R)/(2-u);
}

inline double RNG::Exponential(double L) {

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




#endif
