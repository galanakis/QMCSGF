#ifndef __OPERATORSTRING__
#define __OPERATORSTRING__

#include <map>
#include <cmath>

#include "Probabilities.hh"
#include "RandomNumberGenerator.hh"
#include "SGFBase.hh"

#include <iomanip>

namespace SGF {





/*
  class OperatorStringBase
  defines a list appropriate to use as an kinetic term string.
  It is a CircDList (circular double linked list) for kinetic
  operators. Each operator is a HamiltonianTerm (product of
  creation/annihilation operators) and a time in the imaginary
  axis (CircularTime). The time of the Green operator (Root of
  the list) is an additional member of the OperatorStringBase.

  Notable methods:

  push(direction,term): pushes a term in the string in the right
                        or left of the Green operator. The time
                        of this term will be the same as the
                        Green operator's
  pop(direction,term):  remove a term from the right or left of
                        the Green operator. The time of this operator
                        will be come the new Green time.
*/





enum {LATER,EARLIER};

template<typename T>
class OperatorStringBase : public T {



  CircDList<Operator> &CDList;

  double _Beta;             // Inverse temperature.
  double _Alpha;            // directionality parameter

  // Some running parameters
  double _Renormalization;  // We need the renormalization for the measurements
  double _DeltaRenormalization; // Holds the difference between left minus right normalition.
  double WCLEFT,WCRIGHT,SumWD,DeltaWD,WSum,WDiff;

  _float_accumulator _diagonal_energy; // Keeps track of the diagonal energy

  // This is crazy
  template<int direction>
  inline _float_accumulator delta_diagonal() const {
    return (CDList.length()>=2) ? CDList.top<direction>(direction).Energy*(CDList.top<direction>(!direction).Time-CDList.top<direction>(direction).Time).time() : 0;
  }
  inline double ER_Dtau() const {
    return CDList.empty() ? Energy<RIGHT>() : Energy<RIGHT>()*(CDList.top<LEFT>().Time-CDList.top<RIGHT>().Time).time();
  }

  template<int rl>
  inline const _float_accumulator &Energy() const {
    return T::template Energy<rl>();
  }

  template<int rl>
  inline double weight() const {
    return T::template weight<rl>();
  }

  inline MatrixElement G(int offset=0) const {
    return T::G();
  }


  inline _float_accumulator DiagonalEnergyIntegral() const {
    _float_accumulator result=0;
    for(unsigned int i=0; i+1<CDList.length(); ++i) {
      result+=CDList[i+1].Energy*(CDList[i].Time-CDList[i+1].Time).time();
    }
    return result;
  }

  // Calculates the diagonal energy by doing the time integral explicitly.
  inline double GetDiagonalEnergy() const {
    return ER_Dtau()+DiagonalEnergyIntegral();
  }

  inline double DeltaV() const {
    return Energy<LEFT>()-Energy<RIGHT>();  // The difference between the energies
  }
  inline double DeltaTau() const {
    return Beta()*(CDList.top<LEFT>().Time-CDList.top<RIGHT>().Time).time();
  }

  inline CircularTime newtime() const {
    if(CDList.empty())
      return 1;
    else {
      CircularTime TL=CDList.top<LEFT>().Time;
      CircularTime TR=CDList.top<RIGHT>().Time;

      double dt=(TL-TR).time();
      double L=dt*Beta()*DeltaV();

      return TR+dt*RNG::NZExponential(L);


    }

  }

  /*
  	Executes a single directed update and returns its length.
  	In this method direction is the direction of motion of the Green operator.

  	Lets define

  	dV = VL-VR
  	dt = tL-tR

  	Creation Weight LEFT			WCLEFT  = <GT>/G()
  	Creation Weight RIGHT		WCRIGHT = <TG>/G()
  	Destruction Weight LEFT		WD[LEFT]  = dV/(1-exp(-dV dt))
  	Destruction Weight RIGHT 	WD[RIGHT] = -dV/(1-exp(dV dt))

  	And from this
  	SumWD 	= WD[LEFT] + WD[RIGHT]	= dV Coth[dt dV/2]
  	DeltaWD	= WD[LEFT] - WD[RIGHT]	=	dV
  	Note here that the destruction weights will be zero for an empty string.

  	We can parametrize
  	WD[LEFT]	= (SumWD + DeltaWD)/2;
  	WD[RIGHT]	= (SumWD - DeltaWD)/2;

  	If we define Sign[LEFT] = -1 and Sign[RIGHT]=1 we can write
  	WD[direction] = (SumWD - Sign[direction]*DeltaWD)/2;

  	The total creation destruction weight per direction is
     WT[LEFT]	   =	WCLEFT  +WD[LEFT]
     WT[RIGHT]	=	WCRIGHT +WD[RIGHT]

  	and corresponds to Eq. 35 and 36.
  	Similarly I can define
  	WSum	=WT[LEFT]+WT[RIGHT]		=	WCLEFT+WCRIGHT+SumWD;
  	WDiff	=WT[RIGHT]+WT[RIGHT]	=	WCLEFT-WCRIGHT+DeltaWD;

  	Again we can parametrize WD[LEFT or RIGHT] by these sum and difference
  	WT[direction] = (WSum - Sign[direction]*WDiff)/2;

  	The loop probability from Eq. 42 is
  	PLoop[direction] = alpha min(1,WT[direction]/WT[opposite direction])
  	                 = alpha min(WT[LEFT],WT[RIGHT])/WT[opposite direction])
  	                 = alpha min(WSum-WDiff,WSum+WDiff)/(WSum-Sign[opposite direction]*WDiff)
  		              = alpha (WSum-fabs(WDiff))/(WSum+Sign[direction]*WDiff)

   	where in the last line I used Sign[opposite direction]=-Sign[direction] and min(a+b,a-b)=a+|b|

  	The normalization per direction is (see Eq. 44)
  	Q[direction]= WT[direction] (1-alpha min(1,WR[opposite direction]/WT[direction]))
  	 				= WT[direction]-alpha min(WT[LEFT],WR[RIGHT])
  					= (WSum-Sign[direction]*WDiff)/2 - alpha min(WSum-WDiff,WSum+WDiff)/2
  					= (  WSum-Sign[direction]*WDiff-alpha (WSum-abs(WDiff)) )/2

  	From this we can calculate
  	Renormalization 			= Q[LEFT]+Q[RIGHT] = (1-alpha) WSum + alpha abs(WDiff)
  	DeltaRenormalization 	= Q[LEFT]-Q[RIGHT] = WDiff

  	and again we can write
  	Q[LEFT] 	= (Renormalization+DeltaRenormalization)/2
  	Q[RIGHT]  = (Renormalization-DeltaRenormalization)/2

  	Let r be a random number uniformly distributed in [0,1[ to decide if the direction is RIGHT
  	(Q[LEFT]+Q[RIGHT])*r<Q[RIGHT]

   	This can be written as
  	2*Renormalization*r<Renormalization-DeltaRenormalization
  	or
  	DeltaRenormalization<Renormalization*(1-2*r)

  	The Boltzmann weight is 1/Renormalization and because we need to call
  	it often for measurements we store Renormalization it in a variable
  	and update it after each operator insertion removal.

  	When the string is initialized it is empty and by symmetry <GT>=<TG>
  	Therefore the initialization is
  	_Renormalization=(1-_Alpha)*(WCLEFT+WCRIGHT)
  	_DeltaRenormalization=0;


  */


  /*
    It calculates <GT>/<G> is direction=LEFT and <TG>/<G> is direction is right
    Here the direction is the direction of motion of the Green operator.
    If the operator moves to the left operators are added to their right
    and vice versa.

  template<int rl>
  inline double WeightCreation() const {
    return weight<!rl>()/G();
  }
  */


  /*
    This returns the sum WeightDestruction(LEFT)+WeightDestruction(RIGHT)
    which is equal to dV/(1-exp(-dV dt))-dV/(1-exp(dV dt)) = dV Coth[dt dV/2]

  inline double SumWeightDestuction() const {

    if(CDList.empty())
      return 0;
    else {
      double dV=DeltaV();
      double dt=DeltaTau();
      return 2*xcoth(dV*dt/2)/dt;

    }
  }

  */
  /*
    This returns the difference WeightDestruction(LEFT)+WeightDestruction(RIGHT
    which is equal to dV/(1-exp(-dV dt))+dV/(1-exp(dV dt)) = dV

  inline double DeltaWeightDestuction() const {
    return CDList.empty() ? 0 : DeltaV();
  }

  */


  inline void EvalueRenormalizations() {
    _Renormalization=(1-_Alpha)*WSum+_Alpha*fabs(WDiff);
    _DeltaRenormalization=WCLEFT-WCRIGHT+DeltaWD;

  }

  inline void EvaluateRunningPars() {
    double g=G();
    WCLEFT=weight<RIGHT>()/g;
    WCRIGHT=weight<LEFT>()/g;

    if(CDList.empty()) {
      SumWD= 0;
      DeltaWD=0;
    } else {
      double a=DeltaTau()/2;
      DeltaWD = DeltaV();
      SumWD= (DeltaWD*a==0) ? 1/a : DeltaWD/tanh(DeltaWD*a);
    }

    WSum=WCLEFT+WCRIGHT+SumWD;
    WDiff=WCLEFT-WCRIGHT+DeltaWD;
  }

// Destroy the operator in the front of the Green operator
  template<int rl>
  void destroy() {
    T::template update<rl,REMOVE>(CDList.top<rl>().Term);
    _diagonal_energy -= delta_diagonal<rl>();
    CDList.pop<rl>();
  }

// Create an operator to the right of the Green operator
// Note that the order of "update", "push" and _diagonal_energy increment matters
  void create_RIGHT() {
    const HamiltonianTerm* term=T::template choose<RIGHT>();
    CircularTime tau=newtime();
    T::template update<RIGHT,ADD>(term);
    CDList.push<RIGHT>(Operator(tau,term,Energy<RIGHT>()));
    _diagonal_energy += delta_diagonal<RIGHT>();
  }

// Create an operator to the left of the Green operator
  void create_LEFT() {
    const HamiltonianTerm* term=T::template choose<LEFT>();
    CircularTime tau=newtime();
    CDList.push<LEFT>(Operator(tau,term,Energy<LEFT>()));
    _diagonal_energy += delta_diagonal<LEFT>();
    T::template update<LEFT,ADD>(term);
  }

public:


  OperatorStringBase(SGFBase &base) :  T(base), CDList(base.OperatorCDL), _Beta(base.Beta), _Alpha(base.Alpha) {

    _diagonal_energy=DiagonalEnergyIntegral();

    EvaluateRunningPars();
    EvalueRenormalizations();

    //print_internal();

  }


  void print_internal() {

    T::print_probabilities();
    std::cout<<"Ostring Data"<<std::endl;
    std::cout<<_Beta<<std::endl;
    std::cout<<_Alpha<<std::endl;
    // Some running parameters
    std::cout<< _Renormalization<<std::endl;  // We need the renormalization for the measurements
    std::cout<< _DeltaRenormalization<<std::endl; // Holds the difference between left minus right normalition.
    std::cout<< WCLEFT <<std::endl;
    std::cout<<WCRIGHT<<std::endl;
    std::cout<<SumWD<<std::endl;
    std::cout<<DeltaWD<<std::endl;
    std::cout<<WSum<<std::endl;
    std::cout<<WDiff<<std::endl;

    std::cout<<_diagonal_energy<<std::endl; // Keeps track of the diagonal energy

    std::cout<<"Operator String (size= "<<length()<<")"<<std::endl;

    Operator *o;

    unsigned int precision=17;

    for(unsigned long i=0; i<length(); ++i) {

      o=&CDList[i];

      std::cout<<std::right<<"    ["<<std::setw(6)<<T::KinHash(o->Term)<<", "<<std::fixed<<std::setprecision(precision)<<std::setw(precision+5)<<std::left<<o->Time.time()<<", "<<std::setw(21)<<o->Energy<<"],"<<std::endl;

    }



  }

  ~OperatorStringBase() {

    //print_internal();

  }

  unsigned long length() const {
    return CDList.length();
  }

// Calculates the diagonal energy fast using the updated variable.
  double DiagonalEnergy() const {
    return ER_Dtau()+_diagonal_energy;
  }

// Returns the Kinetic energy from the length of the string
  double KineticEnergy() const {
    return -static_cast<_float_accumulator>(length())/Beta();
  }

  double alpha() const {
    return _Alpha;
  }

  double Beta() const {
    return _Beta;
  }

// The weight for the measurements as shown in Eq. 62 of the paper
  double BoltzmannWeight() const {
    return 1.0/_Renormalization/G();
  }

  bool ChooseDirectionRight() const {

#ifdef DEBUG
    if(_Renormalization==0) {
      std::cerr<<std::endl<<"Error in OperatorString::ChooseDirectionRight()"<<std::endl<<"The renormalization is zero, no meaningful choice can be made."<<std::endl;
      exit(30);
    }
#endif

    return _DeltaRenormalization < _Renormalization*(1-2*RNG::Uniform());
  }

  bool ChooseCreateRight() const {

#ifdef DEBUG
    if(WCRIGHT==0) {
      std::cerr<<std::endl<<"Error in OperatorString::ChooseCreateRight()"<<std::endl<<"The total probabilities are both zero, no meaningful choice can be made."<<std::endl;
      exit(31);
    }
#endif

    return (WCRIGHT+(SumWD-DeltaWD)/2)*RNG::Uniform()<WCRIGHT;
  }

  bool ChooseCreateLeft() const {

#ifdef DEBUG
    if(WCLEFT==0) {
      std::cerr<<std::endl<<"Error in OperatorString::ChooseCreateLeft()"<<std::endl<<"The total probabilities are both zero, no meaningful choice can be made."<<std::endl;
      exit(32);
    }
#endif

    return (WCLEFT +(SumWD+DeltaWD)/2)*RNG::Uniform()<WCLEFT;
  }

  bool KeepGoingRight() const {
    return RNG::Uniform()*(WSum+WDiff) < _Alpha*(WSum-fabs(WDiff));
  }

  bool KeepGoingLeft() const {
    return RNG::Uniform()*(WSum-WDiff) < _Alpha*(WSum-fabs(WDiff));
  }

  unsigned int directed_update() {

    unsigned int UpdateLength=0;

    if(ChooseDirectionRight())
      do {
        if(ChooseCreateRight())
          create_LEFT();
        else
          destroy<RIGHT>();

        EvaluateRunningPars();

        ++UpdateLength;

      } while(KeepGoingRight());
    else
      do {
        if(ChooseCreateLeft())
          create_RIGHT();
        else
          destroy<LEFT>();

        EvaluateRunningPars();

        ++UpdateLength;

      } while(KeepGoingLeft());


    EvalueRenormalizations();

    return UpdateLength;

  }


};

typedef OperatorStringBase<CanonicalProbabilities> CanonicalOperatorString;
typedef OperatorStringBase<GrandProbabilities> GrandOperatorString;
//typedef CanonicalOperatorString OperatorStringType;

}

#endif
