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


std::ostream &operator<<(std::ostream &o,ensemble_t e) {
  return o << ((e==Canonical) ? "Canonical" : "GrandCanonical");
}


inline void AppendHamiltonianTerm(Hamiltonian &H,const HamiltonianTerm &term) {
  if(term.coefficient()!=0)
    H.push_back(term);
}

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
  Lattice *lattice;
  int Nmax;


  unsigned int nsites() const {
    return lattice->size();
  }

  void MakeHamiltonian(SGFBase &Container) const {

    Container.Ensemble=ensemble;

    int NSites=lattice->size();

    if(Population>NSites*Nmax) {
      std::cerr<<"The population is larger than the capacity of the system!"<<std::endl;
      exit(33);
    }

    Container.Psi.resize(NSites);

    for(unsigned int i=0; i<NSites; ++i)
      Container.Psi[i].nmax() = Nmax;

    for(int p=0; p<Population; ++p) {
      unsigned int i=p%NSites;
      Container.Psi[i].nL()++;
      Container.Psi[i].nR()++;
    }


    for(unsigned int i=0; i<NSites; ++i) {
      const IndexedProductElement atom(CCAA,&Container.Psi[i]);
      AppendHamiltonianTerm(Container.V,HamiltonianTerm(CoulombU/2.0,atom));
    }


    LinkDataList links=lattice->links();

    for(LinkDataList::const_iterator it=links.begin(); it!=links.end(); ++it) {
      unsigned int i=it->first;
      unsigned int j=it->second;
      const IndexedProductElement ci(C,&Container.Psi[i]);
      const IndexedProductElement cj(C,&Container.Psi[j]);
      const IndexedProductElement ai(A,&Container.Psi[i]);
      const IndexedProductElement aj(A,&Container.Psi[j]);

      AppendHamiltonianTerm(Container.T,HamiltonianTerm(Hoppingt,ci,aj));
      AppendHamiltonianTerm(Container.T,HamiltonianTerm(Hoppingt,cj,ai));

    }

    std::vector<double> potential=lattice->sitedata();

    for(unsigned int i=0; i<Container.Psi.size(); ++i) {
      double TotalV=mu+potential[i];
      if(TotalV!=0) {
        const IndexedProductElement ni(CA,&Container.Psi[i]);
        AppendHamiltonianTerm(Container.V,HamiltonianTerm(TotalV,ni));
      }

    }

  }

  std::ostream &print(std::ostream &o) const {
    int w=25;
    o<<"    "<<std::setw(w)<<std::left<<"Model"<<"BoseHubbard"<<std::endl;
    o<<"    "<<std::setw(w)<<std::left<<"Hoppingt"<<Hoppingt<<std::endl;
    o<<"    "<<std::setw(w)<<std::left<<"CoulombU"<<CoulombU<<std::endl;
    o<<"    "<<std::setw(w)<<std::left<<"Population"<<Population<<std::endl;
    o<<"    "<<std::setw(w)<<std::left<<"Ensemble"<<ensemble<<std::endl;
    o<<"    "<<std::setw(w)<<std::left<<"Mu"<<mu<<std::endl;
    o<<"    "<<std::setw(w)<<std::left<<"Population"<<Population<<std::endl;
    o<<"    "<<std::setw(w)<<std::left<<"Nmax"<<Nmax<<std::endl;
    lattice->print(o);
    return o;
  }
  ~BoseHubbard() {
    delete lattice;
  }

  BoseHubbard(const rapidjson::Value &m) {

    // The ensemble can only be Canonical or GrandCanonical;
    assert(m["Ensemble"].GetString()==std::string("Canonical") || m["Ensemble"].GetString()==std::string("GrandCanonical"));
    ensemble = m["Ensemble"].GetString()==std::string("Canonical") ? Canonical : GrandCanonical;

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
    } else {
      CoulombU=0;
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
    std::vector<unsigned int> LSizes(dimensions);
    for(int i=0; i<dimensions; ++i)
      LSizes[i]=l["Dimensions"][i].GetInt();

    assert(!m.HasMember("ParabolicTrap") || (m["ParabolicTrap"].IsNumber() || m["ParabolicTrap"].IsArray() && m["ParabolicTrap"].Size()==dimensions));

    if(m.HasMember("ParabolicTrap") && m["ParabolicTrap"].IsNumber() ) {
      double ParabolicHeight=m["ParabolicTrap"].GetDouble();
      switch(dimensions) {
      case 1: {
        lattice=new Chain(LSizes[0],open,ParabolicHeight);
        break;
      }
      case 2: {
        lattice=new Ellipsis(LSizes[0],LSizes[1],ParabolicHeight);
        break;
      }
      case 3: {
        lattice=new Ellipsoid(LSizes[0],LSizes[1],LSizes[2],ParabolicHeight);
        break;
      }
      }

    } else {
      std::vector<boundary_t> BCond(dimensions);
      std::vector<double> PV(dimensions,0);
      for(int i=0; i<dimensions; ++i) {
        std::string bcstring=l["Boundaries"][i].GetString();
        assert(bcstring==std::string("open") || bcstring==std::string("periodic"));
        BCond[i]= (bcstring==std::string("open")) ? open : periodic;
      }
      if(m.HasMember("ParabolicTrap")) {
        for(int i=0; i<dimensions; ++i) {
          PV[i]=m["ParabolicTrap"][i].GetDouble();
        }
      }
      switch(dimensions) {
      case 1: {
        lattice=new Chain(LSizes[0],BCond[0],PV[0]);
        break;
      }
      case 2: {
        lattice=new Square(LSizes[0],BCond[0],PV[0],LSizes[1],BCond[1],PV[1]);
        break;
      }
      case 3: {
        lattice=new Cubic(LSizes[0],BCond[0],PV[0],LSizes[1],BCond[1],PV[1],LSizes[2],BCond[2],PV[2]);
        break;
      }
      }
    }




    unsigned int NSites=lattice->size();

    assert(m["Population"].IsInt() || m["Filling"].IsNumber());

    Population=m["Population"].IsInt() ? m["Population"].GetInt() : m["Filling"].GetDouble()*NSites;

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

    WarmTime=d["WarmTime"].IsInt64() ? d["WarmTime"].GetInt64() : -1;
    WarmIterations=d["WarmIterations"].IsInt64() ? d["WarmIterations"].GetInt64() : -1 ;
    MeasTime=d["MeasTime"].IsInt64() ? d["MeasTime"].GetInt64() : -1;
    MeasIterations=d["MeasIterations"].IsInt64() ? d["MeasIterations"].GetInt64() : -1;

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
