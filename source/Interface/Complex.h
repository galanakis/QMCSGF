#ifndef ComplexDef
#define ComplexDef

#include <cmath>
#include <iostream>
using std::ostream;

class Complex
  {
    double R;
    double I;

    public:
      Complex(void) {R=0.0; I=0.0;}
      Complex(double r) {R=r; I=0.0;}
      Complex(double r,double i) {R=r; I=i;}
      inline double Re(void) {return R;}
      inline double Im(void) {return I;}
      inline double SquareModulous(void) {return R*R+I*I;}
      inline double Modulous(void) {return sqrt(SquareModulous());}
      inline bool IsReal(void) {return I==0.0;}
      inline bool IsImaginary(void) {return R==0.0;}
      friend inline Complex operator+(const Complex &,const Complex &);
      friend inline Complex operator-(const Complex &,const Complex &);
      friend inline Complex operator+(const Complex &);
      friend inline Complex operator-(const Complex &);
      friend inline Complex operator*(const Complex &,const Complex &);
      friend inline Complex operator/(const Complex &,const Complex &);
      friend inline Complex operator%(const Complex &,const Complex &);
      friend inline Complex operator^(const Complex &,const Complex &);
      friend inline Complex sqrt(const Complex &);
      friend inline Complex exp(const Complex &);
      friend inline Complex log(const Complex &);
      friend inline Complex cos(const Complex &);
      friend inline Complex sin(const Complex &);
      friend inline Complex tan(const Complex &);
      friend inline bool operator==(const Complex &,const Complex &);
      friend inline bool operator!=(const Complex &,const Complex &);
      friend ostream &operator<<(ostream &,const Complex &);
  };
  
inline Complex operator+(const Complex &X,const Complex &Y)
  {
    return Complex(X.R+Y.R,X.I+Y.I);
  }

inline Complex operator-(const Complex &X,const Complex &Y)
  {
    return Complex(X.R-Y.R,X.I-Y.I);
  }

inline Complex operator+(const Complex &Z)
  {
    return Z;
  }

inline Complex operator-(const Complex &Z)
  {
    return -Z;
  }

inline Complex operator*(const Complex &X,const Complex &Y)
  {
    return Complex(X.R*Y.R-X.I*Y.I,X.R*Y.I+X.I*Y.R);
  }

inline Complex operator/(const Complex &X,const Complex &Y)
  {
    static double Mod;
    Mod=Y.R*Y.R+Y.I*Y.I;
    return Complex((X.R*Y.R+X.I*Y.I)/Mod,(X.I*Y.R-X.R*Y.I)/Mod);
  }

inline Complex operator%(const Complex &X,const Complex &Y)
  {
    return ((int) X.R)%((int) Y.R);
  }

inline Complex operator^(const Complex &X,const Complex &Y)
  {
    if (X==0)
      return 1.0;
    
    return exp(Y*log(X));
  }

inline Complex sqrt(const Complex &Z)
  {
    static double Mod,Theta;
    Mod=sqrt(sqrt(Z.R*Z.R+Z.I*Z.I));
    Theta=0.5*atan(Z.I/Z.R);
    return Complex(Mod*cos(Theta),Mod*sin(Theta));
  }

inline Complex exp(const Complex &Z)
  {
    static Complex i(0.0,1.0);
    return exp(Z.R)*(cos(Z.I)+i*sin(Z.I));
  }

inline Complex log(const Complex &Z)
  {
    static Complex i(0.0,1.0);
    static double Theta;
    
    if (!Z.R)
      if (Z.I>0)
	Theta=1.57079632679489661923132169163975144;
      else
	Theta=-1.57079632679489661923132169163975144;
    else
      if (Z.R>0)
	Theta=atan(Z.I/Z.R);
      else
	Theta=atan(Z.I/Z.R)+3.14159265358979323846264338327950288;
      
    return 0.5*log(Z.R*Z.R+Z.I*Z.I)+i*Theta;
  }

inline Complex cos(const Complex &Z)
  {
    static Complex i(0.0,1.0);
    return cos(Z.R)*cosh(Z.I)-i*sin(Z.R)*sinh(Z.I);
  }

inline Complex sin(const Complex &Z)
  {
    static Complex i(0.0,1.0);
    return sin(Z.R)*cosh(Z.I)+i*cos(Z.R)*sinh(Z.I);
  }

inline Complex tan(const Complex &Z)
  {
    return sin(Z)/cos(Z);
  }
  
inline bool operator==(const Complex &A,const Complex &B)
  {
    return (A.R==B.R && A.I==B.I);
  }

inline bool operator!=(const Complex &A,const Complex &B)
  {
    return (A.R!=B.R || A.I!=B.I);
  }

ostream &operator<<(ostream &Out,const Complex &C)
      {
	if (!C.R && !C.I) Out << "0";
	if (C.R) Out << C.R;
	if (!C.R && C.I>0.0) Out << "i*" << C.I;
	if (C.R && C.I>0.0) Out << "+i*" << C.I;
	if (C.I<0.0) Out << "-i*" << -C.I;
      	return Out;
      }

#endif
