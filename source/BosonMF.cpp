#include <iostream>
#include <cmath>

#include "Models.hh"


#include "mkl_lapack.h"
typedef double _doublereal;
typedef int _integer;

/*

This is a wrapper for a diagonalization
routine dsyev from lapack.

SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

Input:
SIZE: the dimension of the matrix
A: the symmetric matrix it self in row major format

OUTPUT:
val: eigen values

Note: the eigen vectors will be stored
in the input matrix

*/

//
// N: the linear size of the matrix
// jobz: 'N': compute only eigenvalues, 'V': also compute eigenvectors
// sym_a: the symmetric matrix. If jobz='V' it is overwritten by the eigenvectors (in row major format, like fortran)
// W: it stores the eigenvalues. It should have a size of N.
//

inline void eig(_integer N,char jobz,_doublereal sym_a[],_doublereal W[]) {

  char uplo='L';                 // 'U': upper triangle of A is stored, 'L': lower triangle.

  _integer lwork=-1;
  _integer liwork=-1;            // Note: for jobz='N', it has to be at least 3+5*N
  _integer info;                 // Exit status: should be 0 if no errors occur


  // The first run is only to determine the optimal value for lwork and liwork
  _doublereal lwork_query;
  _integer liwork_query;
  dsyevd_(&jobz,&uplo,&N,sym_a,&N,W,&lwork_query,&lwork,&liwork_query,&liwork,&info);

  lwork=lwork_query;
  liwork=liwork_query;

  _doublereal *work=new _doublereal[lwork];
  _integer *iwork=new _integer[liwork];
  dsyevd_(&jobz,&uplo,&N,sym_a,&N,W,work,&lwork,iwork,&liwork,&info);

  delete[] iwork;
  delete[] work;

}


//
// Solves a system of linear equation a*x=b for a square matrix a
// The solution is returned in b.
//

inline void linsolve(_integer N,_doublereal a[],_doublereal b[]) {

  char trans='N';
  _integer *ipiv=new _integer[N];
  _integer info;
  _integer NRHS=1;

  dgetrf_(&N,&N,a,&N,ipiv,&info);

  dgetrs_(&trans,&N,&NRHS,a,&N,ipiv,b,&N,&info);

  delete[] ipiv;

}

namespace SGF {

// The Boltzman-Einstein distribution
inline long double B(long double x) {
  return 1/(exp(x)-1);
}

// The derivative of Boltzman-Einstein distribution
inline long double dB(long double x) {
  return 0.5/(1.0-cosh(x));
}


const double Tolerance=1e-14;

template<class T> T Min(const T &a,const T &b) {return (a<b)?a:b;}

//
// It finds the minumum value of a v
//
inline double Min(const std::vector<double> &energies) {
  double result=energies[0];
  for(unsigned long i=0; i<energies.size(); ++i)
    result=Min(energies[i],result);
  return result;
}

//
// Uses Newton-Raphson method to find the chemical potential
// N is the number of particles
// beta is the inverse temperature
// E is a vector with all the energies
//
double findmu_(unsigned int N,double beta,const std::vector<double> &E) {

  if(N>E.size()) {
    std::cerr<< "Error:findmu: the population cannot be larger the the number of available states" <<std::endl;
    exit(33);
  }

  // The ground state energy
  long double E0=Min(E);
  unsigned long M=E.size();

  // Contains exp(-beta*(E[i]-E0))
  std::vector<long double> u(M);
  for(unsigned long i=0; i<M; ++i)
    u[i]=exp(-beta*(E[i]-E0));

  long double z=1.0+1.0/N;

  long double f=0;
  long double df=0;

  do {

    f=0;
    df=0;
    for(unsigned long a=0; a<M; ++a) {
      unsigned long i=M-1-a;
      f+=u[i]/(z-u[i]);
      df+=u[i]/(z-u[i])/(z-u[i]);
    }

    z+=(f-N)/df;

  } while(f-N>Tolerance);

  return E0-log(z)/beta;


}


double findmu(unsigned int N,double beta,const std::vector<double> &E) {

  if(N>E.size()) {
    std::cerr<< "Error:findmu: the population cannot be larger the the number of available states" <<std::endl;
    exit(33);
  }

  // The ground state energy
  long double E0=Min(E);
  unsigned long M=E.size();

  // Contains exp(-beta*(E[i]-E0))
  std::vector<long double> e(M);
  for(unsigned long i=0; i<M; ++i)
    e[i]=beta*(E[i]-E0);

  long double bm=-log(1.0+1.0/N);

  long double f=0;
  long double df=0;

  do {

    f=0;
    df=0;
    for(unsigned long a=0; a<M; ++a) {
      unsigned long i=M-1-a;
      long double u=0.5*(e[i]-bm);
      long double su=2.0*sinh(u);
      f+=exp(-u)/su;
      df+=1.0/su/su;
    }

    bm-=(f-N)/df;

  } while(fabs((f-N)/df)>Tolerance);

  return E0+bm/beta;


}


inline double sumabs(const std::vector<double> &v) {
  double result=0;
  for(std::vector<double>::const_iterator it=v.begin(); it!=v.end(); ++it)
    result+=fabs(*it);
  return result;
}

// returns (b(x)-b(y))/(x-y)
inline long double specialf(long double x,long double y) {

  long double u=0.5*(x-y);
  long double a= u!=0 ? sinh(u)/u : 1.0;
  return -a/sinh(0.5*x)/sinh(0.5*y)/4.0;

}

struct LatticeData {
  // The total number of sites
  unsigned int NSites;
  // Contains pairs of nearest neighbors
  list_links_t links;
  // It contains the external trap potantial
  std::vector<double> V0;
  // The non interacting Hamiltonian;
  std::vector<double> H0;
  // A container for the occupancy;
  std::vector<double> occupancy;


public:
  LatticeData(Dimensionality dim,unsigned int L,double t,double W) {
    switch (dim) {
    case D1:
      NSites=L;
      links=links_square(L,open);
      V0=trap(L,W);
      break;
    case D2:
      NSites=L*L;
      links=links_square(L,open,L,open);
      V0=trap(L,W,L,W);
      break;
    case D3:
      NSites=L*L*L;
      links=links_square(L,open,L,open,L,open);
      V0=trap(L,W,L,W,L,W);
      break;
    }

    // It contains the non interacting hamiltonian
    H0.resize(NSites*NSites,0);

    for(list_links_t::const_iterator it=links.begin(); it!=links.end(); ++it) {

      unsigned long i=it->first;
      unsigned long j=it->second;
      H0[i+j*NSites]=H0[j+i*NSites]=t;

    }

    for(unsigned long i=0; i<NSites; ++i)
      H0[i+i*NSites]=V0[i];

    occupancy.resize(NSites,0);

  }

};


std::vector<double> generate_occupancy(double mu,const std::vector<double> &vals,const std::vector<double> &vec) {

  unsigned int NSites=vals.size();

  std::vector<double> result(NSites,0);

  for(unsigned long j=0; j<NSites; ++j) {
    double bj=B(vals[j]-mu);
    for(unsigned long i=0; i<NSites; ++i) {
      const double &uij=vec[i+j*NSites];
      result[i]+=uij*uij*bj;
    }
  }

  return result;

}

std::vector<double> generate_jacobi(double mu,const std::vector<double> &vals,const std::vector<double> &vec) {

  unsigned int NSites=vals.size();

  std::vector<double> Jacobi(NSites*NSites,0);

  for(unsigned long n=0; n<NSites; ++n) {
    for(unsigned long m=0; m<NSites; ++m) {

      double bnm=specialf(vals[n]-mu,vals[m]-mu);

      for(unsigned long i=0; i<NSites; ++i) {
        for(unsigned long j=0; j<NSites; ++j) {

          const double ui_nm=vec[i+n*NSites]*vec[i+m*NSites];
          const double uj_mn=vec[j+m*NSites]*vec[j+n*NSites];
          Jacobi[i+j*NSites]-=bnm*ui_nm*uj_mn;

        }
      }
    }
  }

  return Jacobi;

}

/*
std::vector<double> generate_occupancy(double beta,double U,unsigned int NSites,unsigned int Population,const std::vector<double> &H0,const std::vector<double> &occupancy) {
  std::vector<double> vals(NSites,0);
  std::vector<double> vec(NSites*NSites,0);

  vec=H0;
  for(unsigned long i=0; i<NSites; ++i)
    vec[i+i*NSites]+=U*occupancy[i];
  eig(NSites,'V',&vec[0],&vals[0]);
  // Calculating the appropriate chemical potential for the population
  double mu=findmu(Population,beta,vals);
  // vals now contains beta*(E[i]-mu)>0
  for(unsigned long i=0; i<NSites; ++i)
    vals[i]=beta*(vals[i]-mu);

  //return generate_occupancy(vals,vec);

}

std::vector<double> generate_jacobi_derivative(double beta,double U,unsigned int NSites,unsigned int Population,const std::vector<double> &H0,const std::vector<double> &occupancy) {

  std::vector<double> Jacobi(NSites*NSites,0);

  double epsilon=0.001;

  std::vector<double> n=occupancy;

  std::vector<double> n0=generate_occupancy(beta,U,NSites,Population,H0,n);


  for(unsigned int j=0; j<NSites; ++j)
    std::cout<< n0[j] <<std::endl;

  for(unsigned int j=0; j<NSites; ++j) {


    n[j]+=epsilon;

    std::vector<double> temp=generate_occupancy(beta,U,NSites,Population,H0,n);

    for(unsigned int i=0; i<NSites; ++i)
      Jacobi[i+j*NSites]=-(temp[i]-n0[i])/epsilon;

    n[j]=occupancy[j];

  }

  return Jacobi;

}

*/

inline void eigensystem(double beta,const std::vector<double> &a,double u,const std::vector<double> &diag,std::vector<double> &vec,std::vector<double> &vals) {

  _integer N=diag.size();
  vec=a;
  for(unsigned long i=0; i<N; ++i)
    vec[i+i*N]+=u*diag[i];
  eig(N,'V',&vec[0],&vals[0]);

  // vals now contains beta*E[i]
  for(unsigned long i=0; i<N; ++i)
    vals[i]=beta*vals[i];

}

double findmu(unsigned int N,const std::vector<double> &E) {

  if(N>E.size()) {
    std::cerr<< "Error:findmu: the population cannot be larger the the number of available states" <<std::endl;
    exit(33);
  }

  // The ground state energy
  long double E0=Min(E);
  unsigned long M=E.size();

  // Contains exp(-beta*(E[i]-E0))
  std::vector<long double> e(M);
  for(unsigned long i=0; i<M; ++i)
    e[i]=E[i]-E0;

  long double bm=-log(1.0+1.0/N);

  long double f=0;
  long double df=0;

  do {

    f=0;
    df=0;
    for(unsigned long a=0; a<M; ++a) {
      unsigned long i=M-1-a;
      long double u=0.5*(e[i]-bm);
      long double su=2.0*sinh(u);
      f+=exp(-u)/su;
      df+=1.0/su/su;
    }

    bm-=(f-N)/df;

  } while(fabs((f-N)/df)>Tolerance);

  return E0+bm;


}


void SolveMF(double beta,unsigned int NSites,const std::vector<double> &H0,std::vector<double> &occupancy,double filling,double U) {

  unsigned int Population=NSites*filling;

  std::vector<double> vec(NSites*NSites,0);
  std::vector<double> vals(NSites,0);
  std::vector<double> docc(NSites,0);

  unsigned long iter=0;


  do {

    eigensystem(beta,H0,U,occupancy,vec,vals);
    double mu=findmu(Population,vals);

    docc=generate_occupancy(mu,vals,vec);
    for(unsigned long i=0; i<NSites; ++i)
      docc[i]=occupancy[i]-docc[i];


    std::vector<double> Jacobi=generate_jacobi(mu,vals,vec);
    for(unsigned long i=0; i<NSites; ++i) {
      for(unsigned long j=0; j<NSites; ++j)
        Jacobi[i+j*NSites]*=beta*U;
      Jacobi[i+i*NSites]+=1.0;
    }


    if(iter==10) {




    }

    linsolve(NSites,&Jacobi[0],&docc[0]);

    for(unsigned long i=0; i<NSites; ++i)
      occupancy[i]-=docc[i];

    std::cout<< "Convergence= " << sumabs(docc) <<" center occ= "<<occupancy[NSites/2]<<std::endl;

    ++iter;

  } while(sumabs(docc)>1e-5 );



  std::cout<< "occupancy" <<std::endl;
  for(unsigned long i=0; i<NSites; ++i)
    std::cout<< occupancy[i] <<std::endl;

  std::cout<< "energies" <<std::endl;
  for(unsigned long i=0; i<NSites; ++i)
    std::cout<< vals[i] <<std::endl;


  std::cout<<"beta= "<<beta<<"\tmu= "<<-vals[0]/beta<<"\tU="<<U<<"\tpopulation= "<<Population<<"\tNSites= "<<NSites<<std::endl;

  std::cout<<"Done after "<<iter<<" iterations."<<std::endl;


}


void SimpleSolveMF(double beta,unsigned int NSites,const std::vector<double> &H0,std::vector<double> &occupancy,double filling,double U) {

  unsigned int Population=NSites*filling;

  std::vector<double> vec(NSites*NSites,0);
  std::vector<double> vals(NSites,0);
  std::vector<double> occsave(NSites,0);

  unsigned long iter=0;
  double convergence=0;

  do {

    occsave=occupancy;

    eigensystem(beta,H0,U,occupancy,vec,vals);
    double mu=findmu(Population,vals);

    occupancy=generate_occupancy(mu,vals,vec);

    convergence=0;
    for(unsigned long i=0; i<NSites; ++i)
      convergence+=fabs(occupancy[i]-occsave[i]);

    std::cout<< "Convergence= " << convergence <<" center occ= "<<occupancy[NSites/2]<<std::endl;

    ++iter;

  } while(convergence>1e-5 );



  std::cout<< "occupancy" <<std::endl;
  for(unsigned long i=0; i<NSites; ++i)
    std::cout<< occupancy[i] <<std::endl;

  std::cout<< "energies" <<std::endl;
  for(unsigned long i=0; i<NSites; ++i)
    std::cout<< vals[i] <<std::endl;


  std::cout<<"beta= "<<beta<<"\tmu= "<<-vals[0]/beta<<"\tU="<<U<<"\tpopulation= "<<Population<<"\tNSites= "<<NSites<<std::endl;

  std::cout<<"Done after "<<iter<<" iterations."<<std::endl;


}



}




void testlinsolve() {

  double a[9]= {3,1,5,6,3,1,9,1,0};
  double b[3]= {0,-1,2};

  std::cout<<"Input Matrix"<<std::endl;
  for(unsigned int i=0; i<3; ++i) {
    for(unsigned int j=0; j<3; ++j)
      std::cout<<a[i+3*j]<<"\t";
    std::cout<<std::endl;
  }

  std::cout<<"\nInput Right Hand Side"<<std::endl;
  for(unsigned int i=0; i<3; ++i) {
    std::cout<<b[i]<<std::endl;
  }

  linsolve(3,a,b);

  std::cout<<"\nSolution"<<std::endl;
  for(unsigned int i=0; i<3; ++i) {
    std::cout<<b[i]<<std::endl;
  }




}

void testfindmu() {

  double beta=1000;

  std::vector<double> energies(100);
  for(unsigned int i=0; i<energies.size(); ++i)
    energies[i]=i;

  double mu=SGF::findmu(90,beta,energies);

  std::cout<< mu <<std::endl;

}

int main(int argc, char *argv[]) {


  SGF::Dimensionality dim=SGF::D1;
  unsigned int N=10;
  double filling=0.5;
  double U=0.5;
  double W=10;
  double t=1;
  double beta=10;

  SGF::LatticeData lattice(dim,N,t,W);

  /*

    std::cout<< "Simple" <<std::endl;
    for(double U=0; U<7; U+=0.01) {

      SGF::SimpleSolveMF(beta,lattice.NSites,lattice.H0,lattice.occupancy,filling,U);
      for(int i=0; i<lattice.occupancy.size(); ++i)
        std::cout<< lattice.occupancy[i] <<std::endl;
    }

    //std::cout<< "Jacobian" <<std::endl;
  */
  //SGF::SolveMF(beta,lattice.NSites,lattice.H0,lattice.occupancy,filling,U);


  SGF::SolveMF(beta,lattice.NSites,lattice.H0,lattice.occupancy,filling,U);

  /*
    std::vector<double> betas;
    for(double beta=0.2; beta<=8; beta+=0.2)
      betas.push_back(beta);

    #pragma omp parallel for
    for(int i=0; i<betas.size(); ++i)
      SGF::SolveMF(betas[i],dim,60,filling,8,8);
  */

}
