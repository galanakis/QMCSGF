#include <iostream>
#include <cstring>
#include <limits>
#include <set>
#include <map>
#include <vector>
#include <stack>
#include <string>

#include "HamiltonianTerm.hh"
#include <OperatorString.hh>
#include <Simulation.hh>
#include "Conventions.hh"
#include "SGFBase.hh"


#ifdef USEMPI
// *******************************************
// * Declaration of global variables for MPI *
// *******************************************
#include <mpi.h>
#define Master 0
int NumProcessors,Rank,NameLength;
char ProcessorName[MPI_MAX_PROCESSOR_NAME];

#endif



std::ostream cout(std::cout.rdbuf());

typedef std::pair<unsigned int,unsigned int> link_t;
typedef std::vector<link_t> list_links_t;

typedef enum {periodic,open} boundary_t;


inline unsigned int periodic_site_index(unsigned int x,unsigned int Lx) {
   return x%Lx;
}

inline unsigned int periodic_site_index(unsigned int x,unsigned int Lx,unsigned int y,unsigned int Ly) {
   return x%Lx+Lx*(y%Ly);
}

inline unsigned int periodic_site_index(unsigned int x,unsigned int Lx,unsigned int y,unsigned int Ly,unsigned int z,unsigned int Lz) {
   return x%Lx+Lx*((y%Ly)+Ly*(z%Lz));
}


list_links_t links_square(unsigned int Lx,boundary_t Bx) {

   list_links_t result;
   for(unsigned int x=0; x<Lx; ++x) {
      unsigned int i=periodic_site_index(x,Lx);
      unsigned int a=periodic_site_index(x+1,Lx);
      if(Bx==periodic || x+1!=Lx) result.push_back(link_t(i,a));
   }
   return result;

}


list_links_t links_square(unsigned int Lx,boundary_t Bx,unsigned int Ly,boundary_t By) {

   list_links_t result;
   for(unsigned int x=0; x<Lx; ++x) {
      for(unsigned int y=0; y<Ly; ++y) {
         unsigned int i=periodic_site_index(x,Lx,y,Ly);
         unsigned int a=periodic_site_index(x+1,Lx,y,Ly);
         unsigned int b=periodic_site_index(x,Lx,y+1,Ly);
         if(Bx==periodic || x+1!=Lx) result.push_back(link_t(i,a));
         if(By==periodic || y+1!=Ly) result.push_back(link_t(i,b));
      }
   }
   return result;

}


list_links_t links_square(unsigned int Lx,boundary_t Bx,unsigned int Ly,boundary_t By,unsigned int Lz,boundary_t Bz) {

   list_links_t result;
   for(unsigned int x=0; x<Lx; ++x) {
      for(unsigned int y=0; y<Ly; ++y) {
         for(unsigned int z=0; z<Lz; ++z) {
            unsigned int i=periodic_site_index(x,Lx,y,Ly,z,Lz);
            unsigned int a=periodic_site_index(x+1,Lx,y,Ly,z,Lz);
            unsigned int b=periodic_site_index(x,Lx,y+1,Ly,z,Lz);
            unsigned int c=periodic_site_index(x,Lx,y,Ly,z+1,Lz);
            if(Bx==periodic || x+1!=Lx) result.push_back(link_t(i,a));
            if(By==periodic || y+1!=Ly) result.push_back(link_t(i,b));
            if(Bz==periodic || z+1!=Lz) result.push_back(link_t(i,c));
         }
      }
   }
   return result;

}

std::vector<double> trap(unsigned int Lx,double Vx) {

  std::vector<double> result(Lx);
  for(unsigned int x=0; x<Lx; ++x) {
      double rx=(2.0*x-Lx)/Lx;
      unsigned int i=periodic_site_index(x,Lx);
      result[i]=rx*rx*Vx;
  }

  return result;

}

std::vector<double> trap(unsigned int Lx,double Vx,unsigned int Ly,double Vy) {

  std::vector<double> result(Lx*Ly);
  for(unsigned int x=0; x<Lx; ++x) {
      double rx=(2.0*x-Lx)/Lx;
      for(unsigned int y=0; y<Ly; ++y) {
         double ry=(2.0*y-Ly)/Ly;
         unsigned int i=periodic_site_index(x,Lx,y,Ly);
         result[i]=rx*rx*Vx+ry*ry*Vy;
      }
   }

   return result;

}

std::vector<double> trap(unsigned int Lx,double Vx,unsigned int Ly,double Vy,unsigned int Lz,double Vz) {

  std::vector<double> result(Lx*Ly*Lz);
  for(unsigned int x=0; x<Lx; ++x) {
      for(unsigned int y=0; y<Ly; ++y) {
         for(unsigned int z=0; z<Lz; ++z) {
            double rx=(2.0*x-Lx)/Lx;
            double ry=(2.0*y-Ly)/Ly;
            double rz=(2.0*z-Lz)/Lz;
            unsigned int i=periodic_site_index(x,Lx,y,Ly,z,Lz);
            result[i]=rx*rx*Vx+ry*ry*Vy+rz*rz*Vz;
         }
      }
   }

   return result;

}



void BoseHubbardPeriodic1D() {


   // Model specific parameters
   SGF::MatrixElement U=8.0;
   SGF::MatrixElement t=1.0;

   unsigned int Lx=20;

   unsigned int NSites=Lx;
   unsigned int Population=NSites;
   unsigned int Nmax=0;



   int Seed=34715;
   unsigned int GreenOperatorLines=4;
   unsigned long WarmTime=40000;
   unsigned long WarmIterations=100000000;
   unsigned long MeasTime=40000;
   unsigned long MeasIterations=100000000;
   unsigned long NBins=20;

   SGF::SGFBase Container;

   Container.Beta=1.0/1.02;
   Container.Alpha=0.95;
   Container.Ensemble=SGF::Canonical;


   Container.g.initialize(NSites,GreenOperatorLines);


   RNG::Initialize(Seed);


   Container.Psi.resize(NSites);

   for(unsigned int i=0; i<NSites; ++i)
      Container.Psi[i].nmax() = Nmax;

   for(int p=0; p<Population; ++p) {
      unsigned int i=p%NSites;
      Container.Psi[i].n(0)++;
      Container.Psi[i].n(1)++;
   }




   for(unsigned int i=0; i<NSites; ++i) {
      const SGF::IndexedProductElement atom(SGF::C*SGF::C*SGF::A*SGF::A,&Container.Psi[i]);
      Container.V.push_back(SGF::HamiltonianTerm(U/2.0,atom));
   }


   list_links_t links=links_square(Lx,periodic);

   for(list_links_t::const_iterator it=links.begin(); it!=links.end(); ++it) {
      unsigned int i=it->first;
      unsigned int j=it->second;
      const SGF::IndexedProductElement ci(SGF::C,&Container.Psi[i]);
      const SGF::IndexedProductElement cj(SGF::C,&Container.Psi[j]);
      const SGF::IndexedProductElement ai(SGF::A,&Container.Psi[i]);
      const SGF::IndexedProductElement aj(SGF::A,&Container.Psi[j]);

      Container.T.push_back(SGF::HamiltonianTerm(t,ci,aj));
      Container.T.push_back(SGF::HamiltonianTerm(t,cj,ai));


   }

   // Building the list of measurable operators
   std::vector<SGF::Hamiltonian> _MeasurableOperators;
   std::vector<std::string> _MeasurableNameList;

   _MeasurableOperators.push_back(Container.V);
   _MeasurableNameList.push_back("Potential Energy");

   _MeasurableOperators.push_back(Container.T);
   _MeasurableNameList.push_back("Atom Kinetic energy");





   

   SGF::OperatorStringType OperatorString(Container);

   /* Initializing the simulation. Thermalize, Measure and pring the results */
   Simulation simul("BoseHubbard",cout);




   // We start warm up iterations
   simul.Thermalize(OperatorString,WarmIterations,WarmTime);

   // This defines the measurable objects some of which delay updates even if not measured.
   // This is why I declare the measurable operators after the thermalization.
   SGF::Measurable MeasuredOperators(OperatorString);
   MeasuredOperators.insert(_MeasurableNameList,_MeasurableOperators);

   //We start measurement iterations
   simul.Measure(OperatorString,MeasuredOperators,NBins,MeasIterations,MeasTime);

   // We diplay the results of the simulation
   simul.Results(MeasuredOperators);

}

// ***************************
// * Here starts the program *
// ***************************

int main() {

#ifdef USEMPI
   // ***********************************
   // * Initialization of MPI functions *
   // ***********************************

   MPI_Init(&NumArg,&Arg);
   MPI_Comm_size(MPI_COMM_WORLD,&NumProcessors);
   MPI_Comm_rank(MPI_COMM_WORLD,&Rank);
   MPI_Get_processor_name(ProcessorName,&NameLength);

   // Silence the other nodes. Only the root node can print
   cout.rdbuf( Rank==Master ? std::cout.rdbuf() : 0 );

#endif

   BoseHubbardPeriodic1D();


#ifdef USEMPI
   return MPI_Finalize();
#endif
   return 0;
}
