#ifndef __OPERATORSTRING__
#define __OPERATORSTRING__

#include <map>
#include <deque>

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




// This returns x/(1-exp(-x))
inline long double expweight(long double x) { long double d=1.0L-exp(-x); return (d!=0.0L)?x/d:1.0L; }

// This returns log(1-x*(1-exp(A)))/A which is a monotonic function between 0 and 1 if x is in the same interval
inline long double logexponential(long double A,long double x) { return fabs(A)<0.000000001?x:log(1.0L-x*(1.0L-exp(A)))/A; }


/*
	class AlphaData
	
	This class contains and keeps track of the alpha parameters and it
	knows how to update them.
	
	There is no need for high accuracy since the signs that appear
	in the averaging are all the same and the accumulators are
	reset relatively often.
	
*/	
	
class AlphaData {

	double Alpha[2],AlphaParameter;
	double _accumulator[2];
	unsigned int _count[2];

public:
	AlphaData() {
		Alpha[0]=Alpha[1]=AlphaParameter=0;
		reset();
	}

	inline const double &operator()(int ar) const {return Alpha[ar];}
	inline double &operator()(int ar) {return Alpha[ar];}
	inline const double &operator()() const {return AlphaParameter;}
	inline double &operator()() {return AlphaParameter;}
 
	inline void push(int ar,const double &w) {
		_accumulator[ar]+=w;
		++_count[ar];
	}

	inline void reset() {
		_accumulator[0]=_accumulator[1]=0;
		_count[0]=_count[1]=0;		
	}

	inline void flush() {
		double Acc[2];

		Acc[0]=_accumulator[0]/_count[0];
		Acc[1]=_accumulator[1]/_count[1];
		
 		Alpha[0]=AlphaParameter*Acc[1]/Max(Acc[1],Acc[0]);
		Alpha[1]=AlphaParameter*Acc[0]/Max(Acc[1],Acc[0]);

		reset();
		
	}
	
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
  
	AlphaData Alpha;          // Container for the alpha parameters.
  double _Beta;             // Inverse temperature.

	_float_accumulator _diagonal_energy; // Keeps track of the diagonal energy

	inline _float_accumulator delta_diagonal(const int &direction) const { return (length()>=2) ? top(direction,direction).Energy*(top(direction,!direction).Time-top(direction,direction).Time).time() : 0;  }
	inline double ER_Dtau() const {return Energy(RIGHT)*(top(LEFT).Time-top(RIGHT).Time).time();}
	
public:
  OperatorStringType(const Hamiltonian &T,const Hamiltonian &V,double _beta) : CircDList<Operator>(), Probabilities(T,V), _Beta(_beta) {
		_diagonal_energy=0;
  }

	inline void AlphaUpdate() { Alpha.flush(); }

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
  
	inline double &alpha() {return Alpha();}
  inline double &alpha(int action) {return Alpha(action);}
  inline const CircularTime& GreenTime() const {return _GreenTime;}

  /* Chose an operator randomly and push it in the string. direction is the direction opposite to the direction of motion of the green operator */
  inline bool create(int direction) {
    const HamiltonianTerm *term=choose(direction); 
    push(direction,Operator(_GreenTime,term,Energy(direction)));
		_diagonal_energy += delta_diagonal(direction);
    update(term,direction,ADD); 
    double KeepDestroy=KeepDestroying(direction);
    Alpha.push(REMOVE,KeepDestroy);
    return RNG::Uniform()<Alpha(REMOVE)*KeepDestroy;
  }
  
  // Destroy an operator along the direction of motion of the green operator
  inline bool destroy(int direction) {
    _GreenTime=top(direction).Time;
		const HamiltonianTerm *term=top(direction).Term;
    update(term,direction,REMOVE);
		_diagonal_energy -= delta_diagonal(direction);
    pop(direction);
    double KeepCreate=KeepCreating(!direction);
    Alpha.push(ADD,KeepCreate);
    return RNG::Uniform()<Alpha(ADD)*KeepCreate;
  }

/* 
  The following lines contain all the complexity of the SGF algorithm.
  KeepCreating and KeepDestroying are just definitions.
  In CreationWeight and DestructionWeight the direction is the direction
  of motion of the Green operator.

  KeepCreating(LEFT) =weight(LEFT)/Max(weight(LEFT),weight(RIGHT))
  KeepCreating(RIGHT)=weight(RIGHT)/Max(weight(LEFT),weight(RIGHT))
  KeepDestroying(LEFT) =Min(1.0,exp(-DV DT))
  KeepDestroying(RIGHT)=Min(1.0,exp(+DV DT))
  
  CreationWeight(LEFT)  =(1-Alpha[ADD]*weight(LEFT)/Max(weight(LEFT),weight(RIGHT)))*weight(RIGHT)/G();
  CreationWeight(RIGHT) =(1-Alpha[ADD]*weight(RIGHT)/Max(weight(LEFT),weight(RIGHT)))*weight(LEFT )/G();
  DestructionWeight(LEFT) =empty() ? 0 : (1-Alpha[REMOVE]*Min(1.0,exp(-DeltaV()*DeltaTau())))*expweight( DeltaV()*DeltaTau())/DeltaTau();
  DestructionWeight(RIGHT)=empty() ? 0 : (1-Alpha[REMOVE]*Min(1.0,exp( DeltaV()*DeltaTau())))*expweight(-DeltaV()*DeltaTau())/DeltaTau();
*/

  inline double KeepCreating(int direction) const { return Min(1.0,weight(direction)/weight(!direction)); }
  inline double KeepDestroying(int direction) const { return Min(1.0,exp(Sign[direction]*DeltaV()*DeltaTau())); }
  inline double CreationWeight(int direction) const { return ((1-Alpha(ADD)*KeepCreating(direction))*weight(!direction))/G(); }
  inline double DestructionWeight(int direction) const { return empty() ? 0 : (1-Alpha(REMOVE)*KeepDestroying(direction))*expweight(-Sign[direction]*DeltaV()*DeltaTau())/DeltaTau(); }
  inline double BoltzmannWeight() const {return 1.0/G()/(CreationWeight(LEFT)+CreationWeight(RIGHT)+DestructionWeight(LEFT)+DestructionWeight(RIGHT));}

/* 
   Returns the timeshift of the Green operator towards the direction of motion.
   Left motion:  _GreenTime+=log(1-R*(1-exp( DeltaTau*DeltaV)))/DeltaV
   Right motion: _GreenTime+=log(1-R*(1-exp(-DeltaTau*DeltaV)))/DeltaV
   where R=RNG::NZUniform() is a random number in the interval ]0,1[ (never gets zero).
   If we define le(A,R)=log(1-R*(1-exp(A)))/A then both cases can be parametrized by
   _GreenTime+=-Sign[direction]*DeltaTau*le(-Sign[direction]*DeltaTau*DeltaV)  
*/   
    
  inline CircularTime TimeShift(int direction) const { return TimeShift(direction,RNG::NZUniform()); }
  inline CircularTime TimeShift(int direction,double random) const {return -Sign[direction]*(DeltaTau()/Beta())*logexponential(-DeltaTau()*DeltaV()*Sign[direction],random);}  

  

/*

This is the main algorithm of the directed updates. It starts by choosing a direction
and an initial action. Then it alternates creating and destroying operators and evaluating
the probability for the next step. More explicitly after each creation it calculates the 
probability for destroying (rather than stopping) and vice versa.

The two actions are included in the routines create_operator and destroy_operator which
apart from updating the operator string accordingly, they also calculate the probability
for taking the next step and they decide randomly if the next step should be taken. 
There are two possible sequenses of creation (C) and destruction (D)

1) initial action = C: CDCDCDCDCDC ...
2) initial action = D: DCDCDCDCDCD ...

The sequence stops when one of the probability to continue evaluates to false.

This looks like the expression
C && D && C && D && ...
where the logical operators are evaluated from left to right and execuation stops 
when a false is encountered.

We can exploit this by grouping the sequence as (CD)(CD)(CD)(CD)..., which is easier 
to loop, because we can simply loop over the statement
 
bool flag=true;
while(flag=flag && C && D);

A small complication arises because we want to know what is the final action that was 
executed. We can do this by using a counting variable as follows

bool flag=true;
uint count=0;
while(flag=flag && ++count && C && ++count && D);

Because "count" is incremented right before every action it will be equal to the total
number of actions executed. If count is odd (count&1==1) the last action was the same 
as the original one and if it is even the other way around. We can tabulate the cases:

action      count&1==0  count&1==1
Destroy=0        1          0
Create =1        0          1

It is easy to see that the final action was a creation if (action==(count&1)). Beware 
the operator precedence: (action==count&1) is not the same.

A second complication is because there are two kinds of sequenses that differ in the 
inital action:

1) CDCDCDCD -> C(DC)(DC)(DC) ...
2) DCDCDCDC -> (DC)(DC)(DC)(DC) ...

So if the initial action happens to be a creation, we can execute it separately before 
we enter the main loop.

*/


  /* Executes a single directed update and returns its length. In this method direction 
     is the direction of motion of the Green operator.  */
  inline uint directed_update() {

    double Create[2] ={   CreationWeight(LEFT),    CreationWeight(RIGHT)};
    double Destroy[2]={DestructionWeight(LEFT), DestructionWeight(RIGHT)};
    
    // In those two lines the order of the arguments is dictated by the convensions
    unsigned int direction=RNG::CoinFlip(Create[LEFT] + Destroy[LEFT], Create[RIGHT] + Destroy[RIGHT]);
    unsigned int action=RNG::CoinFlip(Destroy[direction],Create[direction]);
    
    return directed_update(direction,action);
  }

  
   inline uint directed_update(unsigned int direction,unsigned int action) {
     uint UpdateLength=0;

     bool keepgoing=true;    
     if(action==ADD) keepgoing=keepgoing && ++UpdateLength && create(!direction);
     while((keepgoing=keepgoing && ++UpdateLength && destroy(direction) && ++UpdateLength && create(!direction)))
       ; // Nothing to do.
     if(action==(UpdateLength&1))
			_GreenTime+=TimeShift(direction);
     return UpdateLength;     
   }
  
};

}

#endif
