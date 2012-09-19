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
  Time(direction):      Return the time of the operator in the Right or
                        left of the Green operator.
  GreenTime(direction): Green operator time.
  Term(direction):      The kinetic term to the right or left of
                        The Green operator.
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



class OperatorStringType : public CircDList<Operator>, public Probabilities {
  CircularTime _GreenTime;  // The time of the Green operator
  
	double _Alpha;            // directionality parameter
  double _Beta;             // Inverse temperature.
	double _Renormalization;  // We need the renormalization for the measurements
	double _DeltaRenormalization; // Holds the difference between left minus right normalition.
  
	_float_accumulator _diagonal_energy; // Keeps track of the diagonal energy

	inline _float_accumulator delta_diagonal(const int &direction) const { return (length()>=2) ? top(direction,direction).Energy*(top(direction,!direction).Time-top(direction,direction).Time).time() : 0;  }
	inline double ER_Dtau() const {return Energy(RIGHT)*(top(LEFT).Time-top(RIGHT).Time).time();}
	
public:
  OperatorStringType(const Hamiltonian &T,const Hamiltonian &V,double beta,GreenOperator<long double> &g,double alpha) : CircDList<Operator>(), Probabilities(T,V,g), _Beta(beta), _Alpha(alpha) {
		_diagonal_energy=0;

		double WL=WeightCreation(LEFT);
		double WR=WeightCreation(RIGHT);

		_Renormalization=(1-_Alpha)*(WL+WR)+_Alpha*fabs(WL-WR);
		_DeltaRenormalization=WL-WR; 


  }

  // Calculates the diagonal energy fast using the updated variable.
	inline double DiagonalEnergy() const { return ER_Dtau()+_diagonal_energy; }
   
	// Calculates the diagonal energy by doing the time integral explicitly.
	inline double GetDiagonalEnergy() const {
		_float_accumulator result=0;
		for(unsigned int i=0;i<length()-1;++i) 
			result+=que[i+1].Energy*(que[i].Time-que[i+1].Time).time();
			
		return ER_Dtau()+result;
	}


  inline double DeltaV() const {return Energy(LEFT)-Energy(RIGHT);} // The difference between the energies
  inline double DeltaTau() const { return Beta()*(top(LEFT).Time-top(RIGHT).Time).time(); }
  inline const double &Beta() const {return _Beta;}
  inline double& Beta() {return _Beta;} 
  
	inline double &alpha() {return _Alpha;}

  
 	inline CircularTime newtime() const {
		if(empty())
			return RNG::NZUniform();
		CircularTime TL=top(LEFT).Time;
		CircularTime TR=top(RIGHT).Time;
	  double dV = DeltaV();
		double dt = (TL-TR).time(); 
		double r=RNG::NZExponential(fabs(dV)*dt);
		return (dV<0) ? TR+dt*r : TL-dt*r ;
	}
 
  inline double WeightCreation(int direction) const { return weight(!direction)/G(); }
  

	inline static double xcoth(double x) { return (x==0) ? 1 : x/tanh(x); }


	inline double SumWeightDestuction() const {
		if(empty()) 
			return 0;
		else {
			double dV=DeltaV();
			double dt=DeltaTau();
			double x=DeltaTau()*DeltaV()/2; 
			return 2*xcoth(x)/dt;
			
		}  		
	}
  
	inline double DeltaWeightDestuction() const {return empty() ? 0 : DeltaV(); }

	inline double WeightDestruction(int direction) const { return 0.5*(SumWeightDestuction()-Sign[direction]*DeltaWeightDestuction()); } 
	
  
 	inline double Renormalization() const { return _Renormalization; }
  inline double BoltzmannWeight() const { return 1.0/Renormalization(); }


  // Executes a single directed update and returns its length. In this method direction is the direction of motion of the Green operator.  
  inline uint directed_update() {
   
	  
		double WC[2],WSum,WDiff,SumWD,DeltaWD;

		unsigned int direction,opposite_direction; 
		
		if(_DeltaRenormalization/_Renormalization < (1-2*RNG::Uniform())) {
			direction=RIGHT;
			opposite_direction=LEFT;
		}
		else {
			direction=LEFT;
			opposite_direction=RIGHT;
		}
		
    uint UpdateLength=0;

		WC[direction]=WeightCreation(direction);
		SumWD=SumWeightDestuction();
		DeltaWD=DeltaWeightDestuction();

		do {
      
  		if((0.5*(SumWD-Sign[direction]*DeltaWD)+WC[direction])*RNG::Uniform()<WC[direction]) {
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
      
      
		} while(RNG::Uniform()*(WSum+Sign[direction]*WDiff) < _Alpha*(WSum-fabs(WDiff)) );

		//_GreenTime=newtime(); 

		_Renormalization=(1-_Alpha)*WSum+_Alpha*fabs(WDiff);
		_DeltaRenormalization=WDiff;


		return UpdateLength;

  }

  
};

}

#endif
