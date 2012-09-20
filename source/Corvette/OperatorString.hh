#ifndef __OPERATORSTRING__
#define __OPERATORSTRING__

#include <map>
#include <deque>
#include <cmath>

#include "CircularTime.hh"
#include "HamiltonianTerm.hh"
#include "Probabilities.hh"
#include "RandomNumberGenerator.hh"
#include "Accumulator.hh"

namespace SGF {

/*
  class CircDList
  This is a wrapper to deque, which replaces push_back, push_front, etc
  with push(direction,data). 
*/

template<class T>
class CircDList {
protected:
  std::deque<T> que;
public:
	typedef typename std::deque<T>::size_type string_size_type;
  CircDList() : que() {};
  inline void pop(int direction) { direction==RIGHT ? que.pop_front() : que.pop_back(); }
  inline void push(int direction,const T& data) { direction==RIGHT ?  que.push_front(data) : que.push_back(data); }
  inline const T& top(int direction) const { return direction==RIGHT ? que.front() : que.back(); }
	inline const T& top(int direction,int depth) const {return direction==RIGHT ? que[depth] : que[que.size()-depth-1];}
  inline string_size_type length() const {return que.size();}
  inline bool empty() const {return que.empty();}
	inline const T &operator[](string_size_type i) const {return que[i];}
}; 






/*
  class OperatorStringType
  defines a list appropriate to use as an kinetic term string.
  It is a CircDList (circular double linked list) for kinetic
  operators. Each operator is a HamiltonianTerm (product of
  creation/annihilation operators) and a time in the imaginary
  axis (CircularTime). The time of the Green operator (Root of 
  the list) is an additional member of the OperatorStringType.
  
  Notable methods:
  
  push(direction,term): pushes a term in the string in the right
                        or left of the Green operator. The time
                        of this term will be the same as the
                        Green operator's
  pop(direction,term):  remove a term from the right or left of
                        the Green operator. The time of this operator
                        will be come the new Green time.
*/
  


  
  
struct Operator {
  CircularTime Time;
  const HamiltonianTerm* Term;
	_float_accumulator Energy;
  
  Operator(const Operator &o) : Time(o.Time), Term(o.Term), Energy(o.Energy) {}
	Operator &operator=(const Operator &o) {Time=o.Time; Term=o.Term; Energy=o.Energy; return *this;}
  Operator(const CircularTime &_time,const HamiltonianTerm *_term,const _float_accumulator &_energy) : Time(_time), Term(_term), Energy(_energy) {}
};                       

enum {LATER,EARLIER};

inline double xcoth(double x) { return (x==0) ? 1 : x/tanh(x); }


class OperatorStringType : public CircDList<Operator>, public Probabilities {
  
  double _Beta;             // Inverse temperature.
	double _Alpha;            // directionality parameter
	double _Renormalization;  // We need the renormalization for the measurements
	double _DeltaRenormalization; // Holds the difference between left minus right normalition.
  
	_float_accumulator _diagonal_energy; // Keeps track of the diagonal energy

	inline _float_accumulator delta_diagonal(const int &direction) const { return (length()>=2) ? top(direction,direction).Energy*(top(direction,!direction).Time-top(direction,direction).Time).time() : 0;  }
	inline double ER_Dtau() const {return empty() ? 0 : Energy(RIGHT)*(top(LEFT).Time-top(RIGHT).Time).time();}
	

	// Calculates the diagonal energy by doing the time integral explicitly.
	inline double GetDiagonalEnergy() const {
		if(empty()) 
			return 0;
		
		_float_accumulator result=0;
		for(unsigned int i=0;i<length()-1;++i) 
			result+=que[i+1].Energy*(que[i].Time-que[i+1].Time).time();
			
		return ER_Dtau()+result;
	}  
	
  inline double DeltaV() const {return Energy(LEFT)-Energy(RIGHT);} // The difference between the energies
  inline double DeltaTau() const { return Beta()*(top(LEFT).Time-top(RIGHT).Time).time(); }

	// This returns log(1-x*(1-exp(A)))/A which is a monotonic function between 0 and 1 if x is in the same interval
	inline static long double logexponential(long double A,long double x) { return fabs(A)<0.000000001?x:log(1.0L-x*(1.0L-exp(A)))/A; }
  
 	inline CircularTime newtime() const {
		if(empty())
			return 1;
		else {
			CircularTime TL=top(LEFT).Time;
			CircularTime TR=top(RIGHT).Time; 

			double dt=(TL-TR).time();
			double L=dt*Beta()*DeltaV();
			double r=RNG::NZUniform();
     
			return TR+dt*RNG::NZExponential(L);
							
			 
    } 
			
	}
  
	/*
		It calculates <GT>/<G> is direction=LEFT and <TG>/<G> is direction is right
		Here the direction is the direction of motion of the Green operator.
		If the operator moves to the left operators are added to their right
		and vice versa.
	*/
  inline double WeightCreation(int direction) const { return weight(!direction)/G(); }
  


  /* 
		This returns the sum WeightDestruction(LEFT)+WeightDestruction(RIGHT)
		which is equal to dV/(1-exp(-dV dt))-dV/(1-exp(dV dt)) = dV Coth[dt dV/2]
	*/
	inline double SumWeightDestuction() const {

		if(empty()) 
			return 0;
		else {
			double dV=DeltaV();
			double dt=DeltaTau();
			return 2*xcoth(dV*dt/2)/dt;
			
		}  		
	}
 
 	/*
		This returns the difference WeightDestruction(LEFT)+WeightDestruction(RIGHT
		which is equal to dV/(1-exp(-dV dt))+dV/(1-exp(dV dt)) = dV
	*/ 
	inline double DeltaWeightDestuction() const {return empty() ? 0 : DeltaV(); }

 	inline double Renormalization() const { return _Renormalization; }


  /* 
		Executes a single directed update and returns its length. 
		In this method direction is the direction of motion of the Green operator.
		
		Lets define 
		
		dV = VL-VR
		dt = tL-tR
		
		Creation Weight LEFT			WC[LEFT]  = <GT>/G()  
		Creation Weight RIGHT			WC[RIGHT] = <TG>/G()
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
	  WT[LEFT]	=	WC[LEFT]	+WD[LEFT]  
	  WT[RIGHT]	=	WC[RIGHT]	+WD[RIGHT]  
	  
		and corresponds to Eq. 35 and 36.
		Similarly I can define 
		WSum	=WT[LEFT]+WT[RIGHT]		=	WC[LEFT]+WC[RIGHT]+SumWD;
		WDiff	=WT[RIGHT]+WT[RIGHT]	=	WC[LEFT]-WC[RIGHT]+DeltaWD;

		Again we can parametrize WD[LEFT or RIGHT] by these sum and difference
		WT[direction] = (WSum - Sign[direction]*WDiff)/2;
		
		The loop probability from Eq. 42 is 
		PLoop[direction] = alpha min(1,WT[direction]/WT[opposite direction])
		                 = alpha min(WT[LEFT],WT[RIGHT])/WT[opposite direction])
		                 = alpha min(WSum-WDiff,WSum+WDiff)/(WSum-Sign[opposite direction]*WDiff)
										 = alpha (WSum-fabs(WDiff))/(WSum+Sign[direction]*WDiff)
										
	 	where in the last line I used Sign[opposite direction]=-Sign[direction] and min(a+b,a-b)=a+|b|
	
		The normalization per direction is (see Eq. 44)
		Q[direction] = WT[direction] (1-alpha min(1,WR[opposite direction]/WT[direction]))
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
		_Renormalization=(1-_Alpha)*(WC[LEFT]+WC[RIGHT])
		_DeltaRenormalization=0; 
		

	*/ 
	
public:


	void print_times() const {

		//for(unsigned int i=0;i<length();++i) 
		//	std::cout<<que[i].Time<<std::endl;
		
		
	}

  OperatorStringType(const Hamiltonian &T,const Hamiltonian &V,double beta,GreenOperator<long double> &g,double alpha) : CircDList<Operator>(), Probabilities(T,V,g), _Beta(beta), _Alpha(alpha) {
		_diagonal_energy=0;

		double WL=WeightCreation(LEFT);
		double WR=WeightCreation(RIGHT);

		_Renormalization=(1-_Alpha)*(WL+WR);
		_DeltaRenormalization=0; 

  }

  // Calculates the diagonal energy fast using the updated variable.
	inline double DiagonalEnergy() const { return ER_Dtau()+_diagonal_energy; }

	inline double alpha() const {return _Alpha;}
  inline double Beta() const {return _Beta;} 

	// The weight for the measurements as shown in Eq. 62 of the paper
	inline double BoltzmannWeight() const { return 1.0/Renormalization()/G(); }
  
	   
  inline uint directed_update() {
   
		double WC[2],WSum,WDiff,SumWD,DeltaWD;

		unsigned int direction,opposite_direction;
		int sign=1; 
		
		if( _DeltaRenormalization < _Renormalization*(1-2*RNG::Uniform()) ) {
			direction=RIGHT;
			opposite_direction=LEFT;
			sign=1;
		}
		else {
			direction=LEFT;
			opposite_direction=RIGHT;
			sign=-1;
		}
		
    uint UpdateLength=0;

		WC[direction]=WeightCreation(direction);
		SumWD=SumWeightDestuction();
		DeltaWD=DeltaWeightDestuction();

		do {
       
			double TotalWeight=WC[direction]+(SumWD-sign*DeltaWD)/2;
  		if( TotalWeight*RNG::Uniform()<WC[direction]) {
				// Create an operator behind the Green operator
		    const HamiltonianTerm *term=choose(opposite_direction);
		    push(opposite_direction,Operator(newtime(),term,Energy(opposite_direction)));
				_diagonal_energy += delta_diagonal(opposite_direction);
		    update(term,opposite_direction,ADD);		
			}
			else {
        // Destroy the operator in the front of the Green operator
		    update(top(direction).Term,direction,REMOVE);
				_diagonal_energy -= delta_diagonal(direction);
		    pop(direction);		
			}
       
			++UpdateLength;

			WC[LEFT]=WeightCreation(LEFT);
			WC[RIGHT]=WeightCreation(RIGHT);
			SumWD=SumWeightDestuction();
			DeltaWD=DeltaWeightDestuction();
			WSum=WC[LEFT]+WC[RIGHT]+SumWD;
			WDiff=WC[LEFT]-WC[RIGHT]+DeltaWD;     
			
      
		} while(RNG::Uniform()*(WSum+sign*WDiff) < _Alpha*(WSum-fabs(WDiff)) );

		_Renormalization=(1-_Alpha)*WSum+_Alpha*fabs(WDiff);
		_DeltaRenormalization=WDiff;


		return UpdateLength;

  }

  
};

}

#endif
