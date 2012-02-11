#ifndef __RANDOMNUMBERGENERATOR__
#define __RANDOMNUMBERGENERATOR__

// ********************************************************************************
// * This class is a high speed random number generator.                          *
// *                                                                              *
// * Initialization: Before first use, the function RNG::Initialize(Seed) must    *
// * be called, where Seed is any non zero integer value.                         *
// *                                                                              *
// * RNG::Uniform() returns a random number in [0;1[ with a uniform distribution. *
// *                                                                              *
// * RNG::Exponential(A) returns a random number in [0;+inf[ with an exponential  *
// * distribution of the form A*Exp(-A*x). A must be positive.                    *
// *                                                                              *
// * RNG::Exponential(A,T) returns a random number in [0;T[ with an exponential   *
// * distribution of the form A*Exp(-A*x)/(1-Exp(-A*T)). A can be any real value. *
// *                                                                              *
// * Valy Rousseau - October 2007                                                 *
// *                                                                              *
// * Modified by Dimitris Galanakis to include sprng and MersenneTwister          *
// ********************************************************************************

#include <iostream>
#include <cmath>
#include <cstdlib>
 
#if defined(RNG_MT)

// Mersenne Twister
#include "MersenneTwister.h" 
class RNGMT {
  static MTRand rand; 
public:
  inline static void Initialize(int s) {
//    cout<<"Using the Mersenne Twister random number generator."<<std::endl;
#ifdef USEMPI
    s=s*(Rank+1);
//		std::cout<<"Initialize processor "<<Rank<<" with seed "<<s<<std::endl;
#endif
    rand.seed(s);
  }
  inline static double Uniform() {return rand.randExc();}
  inline static double NZUniform() {return rand.randDblExc();}
};

MTRand RNGMT::rand;
typedef RNGMT RNGBase;

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
    do { result=Uniform(); } while(result==0);
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
      cout<<"Error in RNG::Initialize: The initialization seed cannot be zero."<<std::endl;
      exit(1);
    }
    for (int i=0;i<1000;++i)
      Seed*=16807;
  }
// Returns a uniformly distributed random number in the interval [0,1[
  inline static double Uniform() { return (Seed*=16807)/4294967296.0+0.5; }
// Returns a uniformly distributed random number in the interval ]0,1[
  inline static double NZUniform() {
    double result=0;
    do { result=Uniform(); } while(result==0);
    return result;
  }      
};

int RNGLC::Seed;

typedef RNGLC RNGBase;

#endif


class RNG : public RNGBase {
public:
// Initialize the random number generator with a non zero seed. 
  static void Initialize(int S) {RNGBase::Initialize(S);};
// Returns a uniformly distributed random number in the interval [0,1[
  inline static double Uniform() {return RNGBase::Uniform();};
// Returns a uniformly distributed random number in the interval ]0,1[
  inline static double NZUniform() {return RNGBase::NZUniform();};
// Returns 0 or 1 with relative probabilities a and b respectively 
  static inline int CoinFlip(double a,double b) { return (a+b)*Uniform()<b;}
// Returns a number with relative distribution exp(-A x)
  static inline double Exponential(double A) { return -log(1.0-Uniform())/A; }
  static inline double Exponential(double A,double T) {

    if (fabs(A*T)<0.000001)
      return T*Uniform();

    if (A*T<-20.0)
      return T*Uniform();

    double Temp=exp(-A*T);
    return -log((1.0-Uniform())*(1.0-Temp)+Temp)/A;
  }
};




#endif
