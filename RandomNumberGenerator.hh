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
// * Modified by Dimitris Galanakis to include sprng.                             *
// ********************************************************************************
#ifdef SIMPLE_SPRNG
#include <sprng.h>
#endif

#include <iostream>
#include <cmath>
#include <cstdlib>


class RNG {
public:

  static int Seed;
// Initialize the random number generator with a non zero seed. 
  static void Initialize(int);
// Returns a uniformly distributed random number in the interval [0,1[
  inline static double Uniform();
// Returns a uniformly distributed random number in the interval ]0,1[
  inline static double NZUniform();
// Returns 0 or 1 with relative probabilities a and b respectively 
  inline static int CoinFlip(double,double);
// Returns a number with relative distribution exp(-A x)
  inline double Exponential(double A);
  inline double Exponential(double A,double T);

};

int RNG::Seed;

void RNG::Initialize(int S) {

  Seed=S;

  if(S==0) {
    std::cout<<"Error in RNG::Initialize: The initialization seed cannot be zero."<<std::endl;
    exit(1);
  }

#ifdef SIMPLE_SPRNG
  init_sprng(SPRNG_CMRG,Seed,SPRNG_DEFAULT);
#else
  for (int i=0;i<1000;++i)
    Seed*=16807;
#endif
}

inline double RNG::Uniform() {
#ifdef SIMPLE_SPRNG
  return sprng();
#else
  return (Seed*=16807)/4294967294.0+0.5;
#endif
}

inline int RNG::CoinFlip(double a,double b) { return (a+b)*Uniform()<b;}

inline double RNG::NZUniform() {
  double result=0;
  do { result=Uniform(); } while(result==0);
  return result;
}

inline double RNG::Exponential(double A) {
  return -log(1.0-Uniform())/A;
}

inline double RNG::Exponential(double A,double T) {
  static double Temp;

  if (fabs(A*T)<0.000001)
    return T*Uniform();

  if (A*T<-20.0)
    return T*Uniform();

  Temp=exp(-A*T);
  return -log((1.0-Uniform())*(1.0-Temp)+Temp)/A;
}




#endif
