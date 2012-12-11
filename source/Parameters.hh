#include <string>
#include <cstdio>
#include <iostream>
#include <vector>
#include <iomanip>
#include "rapidjson/document.h"

#include "SGF.hh"


static std::string
readInputFile( const char *path ) {
   FILE *file = fopen( path, "rb" );
   if ( !file )
      return std::string("");
   fseek( file, 0, SEEK_END );
   long size = ftell( file );
   fseek( file, 0, SEEK_SET );
   std::string text;
   char *buffer = new char[size+1];
   buffer[size] = 0;
   if ( fread( buffer, 1, size, file ) == (unsigned long)size )
      text = buffer;
   fclose( file );
   delete[] buffer;
   return text;
}


namespace SGF {


std::ostream &operator<<(std::ostream &o,boundary_t e) {
   return o << ((e==open) ? "open" : "periodic");
}


std::ostream &operator<<(std::ostream &o,ensemble_t e) {
   return o << ((e==Canonical) ? "Canonical" : "GrandCanonical");
}

class Lattice {
public:
   Lattice() {}
   virtual unsigned int size() const =0;
   virtual std::ostream &print(std::ostream &o) const =0;
   virtual const std::vector<double> potential() const =0;
   virtual ~Lattice() {};
};

class LatticeLinks : public Lattice {
public:
   virtual list_links_t links() const =0;
};

class Chain : public LatticeLinks {
   boundary_t Bx;
   unsigned int Lx;
   double Vx;
public:
   Chain(unsigned int lx,boundary_t bx) : Bx(bx), Lx(lx), Vx(0) {}
   unsigned int size() const {
      return Lx;
   }
   std::ostream &print(std::ostream &o) const {
      int w=25;
      return o<<"    "<<std::setw(w)<<std::left<<"Lattice"<<"Chain Lx="<<Lx<<" ("<<Bx<<")"<<std::endl;
   }
   ~Chain() {};
   list_links_t links() const {
      return links_square(Lx,Bx);
   }
   void set_trap(double a) {
      Vx=a;
   }
   const std::vector<double> potential() const {
      return trap(Lx,Vx);
   }

};

class Square : public LatticeLinks {
   boundary_t Bx,By;
   unsigned int Lx,Ly;
   double Vx,Vy;
public:
   Square(unsigned int lx,boundary_t bx,unsigned int ly,boundary_t by) : Bx(bx), Lx(lx), By(by), Ly(ly), Vx(0), Vy(0) {}
   unsigned int size() const {
      return Lx*Ly;
   }
   std::ostream &print(std::ostream &o) const {
      int w=25;
      return o<<"    "<<std::setw(w)<<std::left<<"Lattice"<<"Square Lx="<<Lx<<" ("<<Bx<<") Ly="<<Ly<<" ("<<By<<")"<<std::endl;
   }
   ~Square() {};
   list_links_t links() const {
      return links_square(Lx,Bx,Ly,By);
   }
   void set_trap(double a,double b) {
      Vx=a;
      Vy=b;
   }
   const std::vector<double> potential() const {
      return trap(Lx,Vx,Ly,Vy);
   }

};

class Cubic : public LatticeLinks {
   boundary_t Bx,By,Bz;
   unsigned int Lx,Ly,Lz;
   double Vx,Vy,Vz;
public:
   Cubic(unsigned int lx,boundary_t bx,unsigned int ly,boundary_t by,unsigned int lz,boundary_t bz) : Bx(bx), Lx(lx), By(by), Ly(ly), Bz(bz), Lz(lz), Vx(0), Vy(0), Vz(0) {}
   unsigned int size() const {
      return Lx*Ly*Lz;
   }
   std::ostream &print(std::ostream &o) const {
      int w=25;
      return o<<"    "<<std::setw(w)<<std::left<<"Lattice"<<"Cubic Lx="<<Lx<<" ("<<Bx<<") Ly="<<Ly<<" ("<<By<<") "<<Ly<<" ("<<By<<")"<<std::endl;
   }
   ~Cubic() {};
   list_links_t links() const {
      return links_square(Lx,Bx,Ly,By,Lz,Bz);
   }
   void set_trap(double a,double b,double c) {
      Vx=a;
      Vy=b;
      Vz=c;
   }
   const std::vector<double> potential() const {
      return trap(Lx,Vx,Ly,Vy,Lz,Vz);
   }

};


class Model {
public:
   Model() {}
   virtual ~Model() {}
   virtual std::ostream &print(std::ostream &o) const = 0;
   virtual void MakeHamiltonian(SGFBase &Container) const =0;
   virtual unsigned int nsites() const =0;
};

struct BoseHubbard : public Model {

   double Hoppingt;
   double CoulombU;
   unsigned int Population;
   ensemble_t ensemble;
   double mu;
   LatticeLinks *lattice;
   std::vector<double> ParabolicTrapHeights;
   int Nmax;


   unsigned int nsites() const {
      return lattice->size();
   }

   void MakeHamiltonian(SGFBase &Container) const {

      Container.Ensemble=ensemble;

      int NSites=lattice->size();

      Container.Psi.resize(NSites);

      for(unsigned int i=0; i<NSites; ++i)
         Container.Psi[i].nmax() = Nmax;

      for(int p=0; p<Population; ++p) {
         unsigned int i=p%NSites;
         Container.Psi[i].nL()++;
         Container.Psi[i].nR()++;
      }


      for(unsigned int i=0; i<NSites; ++i) {
         const IndexedProductElement atom(C*C*A*A,&Container.Psi[i]);
         Container.V.push_back(HamiltonianTerm(CoulombU/2.0,atom));
      }


      list_links_t links=lattice->links();

      for(list_links_t::const_iterator it=links.begin(); it!=links.end(); ++it) {
         unsigned int i=it->first;
         unsigned int j=it->second;
         const IndexedProductElement ci(C,&Container.Psi[i]);
         const IndexedProductElement cj(C,&Container.Psi[j]);
         const IndexedProductElement ai(A,&Container.Psi[i]);
         const IndexedProductElement aj(A,&Container.Psi[j]);

         Container.T.push_back(HamiltonianTerm(Hoppingt,ci,aj));
         Container.T.push_back(HamiltonianTerm(Hoppingt,cj,ai));


      }

      if(HasTrap() || mu!=0) {
         std::vector<double> potential=lattice->potential();

         if(Container.Psi.size()!=potential.size()) {
            std::cout<<"Potential size and psi size differ!"<<std::endl;
            exit(34);
         }

         for(unsigned int i=0; i<Container.Psi.size(); ++i) {
            double TotalV=mu+potential[i];
            if(TotalV!=0) {
               const IndexedProductElement ni(C*A,&Container.Psi[i]);
               Container.V.push_back(HamiltonianTerm(TotalV,ni));
            }

         }
      }


   }

   bool HasTrap() const {
      bool result = false;
      for(int i=0; i<ParabolicTrapHeights.size(); ++i)
         result = result || (ParabolicTrapHeights[i]!=0);
      return result;
   }

   std::ostream &print(std::ostream &o) const {
      int w=25;
      o<<"    "<<std::setw(w)<<std::left<<"Model"<<"BoseHubbard"<<std::endl;
      lattice->print(o);
      o<<"    "<<std::setw(w)<<std::left<<"Hoppingt"<<Hoppingt<<std::endl;
      o<<"    "<<std::setw(w)<<std::left<<"CoulombU"<<CoulombU<<std::endl;
      o<<"    "<<std::setw(w)<<std::left<<"Population"<<Population<<std::endl;
      o<<"    "<<std::setw(w)<<std::left<<"Ensemble"<<ensemble<<std::endl;
      o<<"    "<<std::setw(w)<<std::left<<"Mu"<<mu<<std::endl;
      o<<"    "<<std::setw(w)<<std::left<<"Population"<<Population<<std::endl;
      o<<"    "<<std::setw(w)<<std::left<<"Nmax"<<Nmax<<std::endl;

      o<<"    "<<std::setw(w)<<std::left<<"ParabolicTrapHeights";
      for(int i=0; i<ParabolicTrapHeights.size(); ++i)
         o<<ParabolicTrapHeights[i]<<" ";
      o<<"    "<<std::endl;

      return o;
   }
   ~BoseHubbard() {
      delete lattice;
   }

   BoseHubbard(const rapidjson::Value &m) {

      // The ensemble can only be Canonical or GrandCanonical;
      assert(m["Ensemble"].GetString()==std::string("Canonical") || m["Ensemble"].GetString()==std::string("GrandCanonical"));
      ensemble = m["Ensemble"].GetString()==std::string("Canonical") ? Canonical : GrandCanonical;

      assert(m["Population"].IsInt());
      Population=m["Population"].GetInt();

      // If the ensemble is Canonical there should be a chemical potential
      mu=0;
      if(m.HasMember("Mu")) {
         assert(m["Mu"].IsNumber());
         mu=m["Mu"].GetDouble();
      }

      // If Nmax is not determined assume it is zero (i.e. infinite)
      Nmax = m["Nmax"].IsInt() ? m["Nmax"].GetInt() : 0;

      // If Nmax==1 then the Coulomb interaction is meaningless and can be omitted.
      if(Nmax!=1) {
         assert(m["CoulombU"].IsNumber());
         CoulombU=m["CoulombU"].GetDouble();
      }

      // The hoping matrix elment
      assert(m["Hoppingt"].IsNumber());
      Hoppingt=m["Hoppingt"].GetDouble();


      /*
         Geometry specific assignments
      */
      assert(m["Lattice"].IsObject());
      const rapidjson::Value &l=m["Lattice"];

      std::string LatticeType=l["Type"].GetString();
      assert(LatticeType==std::string("Square"));
      assert(l["Dimensions"].IsArray() && l["Boundaries"].IsArray() && l["Dimensions"].Size()==l["Boundaries"].Size());
      unsigned int dimensions=l["Dimensions"].Size();
      assert(dimensions<=3);
      std::vector<unsigned int> LSizes;
      std::vector<boundary_t> BCond;




      ParabolicTrapHeights.resize(dimensions);
      // There may be a parabolic trap. The heights of it at the boundaries are provided as an array
      // The array must have the same number of elements as the the dimensionality of the system.
      if(m.HasMember("ParabolicTrap")) {
         assert(m["ParabolicTrap"].IsArray() && m["ParabolicTrap"].Size()==dimensions);
         for(int i=0; i<dimensions; ++i) {
            ParabolicTrapHeights[i]=m["ParabolicTrap"][i].GetDouble();
         }
      }
      for(int i=0; i<dimensions; ++i) {
         LSizes.push_back(l["Dimensions"][i].GetInt());
         std::string bcstring=l["Boundaries"][i].GetString();
         assert(bcstring==std::string("open") || bcstring==std::string("periodic"));
         boundary_t bc= (bcstring==std::string("open")) ? open : periodic;
         BCond.push_back(bc);
      }

      const std::vector<double> &PV=ParabolicTrapHeights;

      switch(dimensions) {
      case 1: {
         Chain *latt=new Chain(LSizes[0],BCond[0]);
         latt->set_trap(PV[0]);
         lattice=latt;
         break;
      }
      case 2: {
         Square *latt=new Square(LSizes[0],BCond[0],LSizes[1],BCond[1]);
         latt->set_trap(PV[0],PV[1]);
         lattice=latt;
         break;
      }
      case 3: {
         Cubic *l=new Cubic(LSizes[0],BCond[0],LSizes[1],BCond[1],LSizes[2],BCond[2]);
         l->set_trap(PV[0],PV[1],PV[2]);
         lattice=l;
         break;
      }
      }

   }


};



//
// This class parses the input string and
// extracts all the parameters.
//
struct Parameters {

   std::string json;

   double Alpha;
   double Beta;
   unsigned int GreenOperatorLines;
   int Seed;
   unsigned int NBins;
   unsigned long WarmTime;
   unsigned long WarmIterations;
   unsigned long MeasTime;
   unsigned long MeasIterations;

   std::vector<std::string> Measurables;

   Model *model;

   std::string model_id;

   std::string label() const {
      return model_id;
   }

   bool HasMeasurable(const std::string &s) const {
      bool result=false;
      for(int i=0; i<Measurables.size(); ++i)
         result = result || (Measurables[i]==s);
      return result;
   }


   void MakeContainer(SGFBase &Container) const {

      Container.Alpha=Alpha;
      Container.Beta=Beta;
      model->MakeHamiltonian(Container);
      Container.g.initialize(model->nsites(),GreenOperatorLines);

   }

   std::ostream &print(std::ostream &o) const {
      int w=25;

      model->print(o);

      o<<"    "<<std::setw(w)<<std::left<<"Beta"<<Beta<<std::endl;
      o<<"    "<<std::setw(w)<<std::left<<"Alpha"<<Alpha<<std::endl;
      o<<"    "<<std::setw(w)<<std::left<<"GreenOperatorLines"<<GreenOperatorLines<<std::endl;
      o<<"    "<<std::setw(w)<<std::left<<"Seed"<<Seed<<std::endl;
      o<<"    "<<std::setw(w)<<std::left<<"NBins"<<NBins<<std::endl;
      o<<"    "<<std::setw(w)<<std::left<<"WarmTime"<<WarmTime<<std::endl;
      o<<"    "<<std::setw(w)<<std::left<<"WarmIterations"<<WarmIterations<<std::endl;
      o<<"    "<<std::setw(w)<<std::left<<"MeasTime"<<MeasTime<<std::endl;
      o<<"    "<<std::setw(w)<<std::left<<"MeasIterations"<<MeasIterations<<std::endl;

      o<<"    "<<std::setw(w)<<std::left<<"Measurables";
      for(int i=0; i<Measurables.size(); ++i)
         o<<Measurables[i]<<" ";
      o<<std::endl;

      return o;
   }

   Parameters(std::string &json_input) : json(json_input) {

      rapidjson::Document d;
      d.Parse<0>(json.c_str());


      // Only one of "Beta" or "Temperature" should be present, but not both.
      assert(d["Beta"].IsNumber() ^ d["Temperature"].IsNumber());
      Beta= d.HasMember("Beta") ? d["Beta"].GetDouble() : 1.0/d["Temperature"].GetDouble();

      // The cutoff for the Green operator lines.
      GreenOperatorLines=d["GreenOperatorLines"].GetInt();

      // The alpha parameter for the updates
      Alpha=d["Alpha"].GetDouble();

      // Initial Seed for random number generator
      Seed=d["Seed"].GetInt();

      // Number of Binds for measurements
      NBins=d["Bins"].GetInt();

      // at least one of "WarmTime" and "WarmIterations" should exist
      assert(d["WarmTime"].IsInt() || d["WarmIterations"].IsInt());

      // at least one of "MeasTime" and "WarmIterations" should exist
      assert(d["MeasTime"].IsInt() || d["MeasIterations"].IsInt());

      WarmTime=d["WarmTime"].IsInt() ? d["WarmTime"].GetInt() : -1;
      WarmIterations=d["WarmIterations"].IsInt() ? d["WarmIterations"].GetInt() : -1 ;
      MeasTime=d["MeasTime"].IsInt() ? d["MeasTime"].GetInt() : -1;
      MeasIterations=d["MeasIterations"].IsInt() ? d["MeasIterations"].GetInt() : -1;

      // Names of the operators to measure
      const rapidjson::Value &measure=d["Measure"];
      assert(measure.IsArray());
      for(rapidjson::SizeType i = 0; i < measure.Size(); i++) {
         Measurables.push_back(measure[i].GetString());
      }


      /*
         Model specific stuff.
         Those are relevant to the Hubbard model on a chain/square/cubic lattice
      */

      assert(d["Model"].IsObject());
      const rapidjson::Value &m=d["Model"];

      // The only model we know :-)
      assert(m["Label"].GetString()==std::string("BoseHubbard"));

      model_id=m["Label"].GetString();

      BoseHubbard *model_ptr=new BoseHubbard(d["Model"]);

      model=model_ptr;

   }

   ~Parameters() {
      delete model;
   }

};

}
