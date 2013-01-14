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


class Model {
public:
  Model() {}
  virtual ~Model() {}
  virtual std::ostream &print(std::ostream &o,const std::string &) const = 0;
  virtual void MakeContainer(SGFBase &Container) const =0;
  virtual unsigned int nsites() const =0;
};

struct BoseHubbard : public Model {

  double Beta;
  double Hoppingt;
  double CoulombU;
  unsigned int Population;
  ensemble_t ensemble;
  double mu;
  Lattice *lattice;
  int Nmax;
  boundary_t boundary;

  std::vector<unsigned int> InitDistribution;


  unsigned int nsites() const {
    return lattice->size();
  }

  void MakeContainer(SGFBase &Container) const {

    Container.Beta=Beta;
    Container.Ensemble=ensemble;

    int NSites=lattice->size();

    if(Nmax!=0 && Population>NSites*Nmax) {
      std::cerr<<"The population is larger than the capacity of the system!"<<std::endl;
      exit(33);
    }

    Container.Psi.resize(NSites);

    for(unsigned int i=0; i<NSites; ++i) {
      Container.Psi[i].nmax() = Nmax;
      Container.Psi[i].nL() = InitDistribution[i];
      Container.Psi[i].nR() = InitDistribution[i];
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

  std::ostream &print(std::ostream &o,const std::string &indent) const {
    int w=25;
    o<<indent<<std::setw(w)<<std::left<<"Model:"<<"BoseHubbard\n";
    o<<indent<<std::setw(w)<<std::left<<"Beta:"<<std::setprecision(10)<<Beta<<"\n";
    o<<indent<<std::setw(w)<<std::left<<"Temperature:"<<std::setprecision(10)<<1.0/Beta<<"\n";
    o<<indent<<std::setw(w)<<std::left<<"Hoppingt:"<<Hoppingt<<"\n";
    o<<indent<<std::setw(w)<<std::left<<"CoulombU:"<<CoulombU<<"\n";
    o<<indent<<std::setw(w)<<std::left<<"Population:"<<Population<<"\n";
    o<<indent<<std::setw(w)<<std::left<<"Ensemble:"<<ensemble<<"\n";
    o<<indent<<std::setw(w)<<std::left<<"Mu:"<<mu<<"\n";
    o<<indent<<std::setw(w)<<std::left<<"Nmax:"<<Nmax<<"\n";
    o<<indent<<std::setw(w)<<std::left<<"Boundary:"<<boundary<<"\n";
    lattice->yaml_print(o,"  ");
    o<<std::endl;
    return o;
  }
  ~BoseHubbard() {
    delete lattice;
  }

  BoseHubbard(const rapidjson::Value &m) {

    // The hoping matrix element
    assert(m["Hoppingt"].IsNumber());
    Hoppingt=m["Hoppingt"].GetDouble();

    // Only one of "Beta" or "Temperature" should be present, but not both.
    assert(m["Beta"].IsNumber() ^ m["Temperature"].IsNumber());
    Beta= m.HasMember("Beta") ? m["Beta"].GetDouble() : 1.0/m["Temperature"].GetDouble();

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
      CoulombU=m["CoulombU"].GetDouble();
    } else {
      CoulombU=0;
    }

    /*
       Geometry specific assignments
    */
    assert(m["Lattice"].GetString()==std::string("Square"));

    std::map<std::string,boundary_t> bcmap;
    bcmap["open"]=open;
    bcmap["periodic"]=periodic;
    bcmap["trap"]=trap;



    std::string bcstring=m["Boundary"].GetString();
    assert(bcstring==std::string("open") || bcstring==std::string("periodic") || bcstring==std::string("trap"));
    boundary= bcmap[bcstring];

    assert(m["Dimensions"].IsArray());
    unsigned int dimensions=m["Dimensions"].Size();
    assert(dimensions<=3 && dimensions>0);
    std::vector<unsigned int> LSizes(dimensions);
    for(int i=0; i<dimensions; ++i)
      LSizes[i]=m["Dimensions"][i].GetInt();

    std::vector<double> W(dimensions,0);
    if(m.HasMember("ParabolicTrapHeight")) {
      assert(m["ParabolicTrapHeight"].IsArray() && m["ParabolicTrapHeight"].Size()==dimensions);
      for(int i=0; i<dimensions; ++i)
        W[i]=m["ParabolicTrapHeight"][i].GetDouble();
    }

    if(boundary==trap) {
      switch(dimensions) {
      case 1: {
        lattice=new Chain(LSizes[0],open,W[0]);
        break;
      }
      case 2: {
        lattice=new Ellipsis(LSizes[0],LSizes[1],W[0],W[1]);
        break;
      }
      case 3: {
        lattice=new Ellipsoid(LSizes[0],LSizes[1],LSizes[2],W[0],W[1],W[3]);
        break;
      }
      }
    } else {
      switch(dimensions) {
      case 1: {
        lattice=new Chain(LSizes[0],boundary,W[0]);
        break;
      }
      case 2: {
        lattice=new Square(LSizes[0],boundary,LSizes[1],boundary,W[0],W[1]);
        break;
      }
      case 3: {
        lattice=new Cubic(LSizes[0],boundary,LSizes[1],boundary,LSizes[2],boundary,W[0],W[1],W[2]);
        break;
      }
      }
    }



    unsigned int NSites=lattice->size();

    assert(m["Population"].IsInt() || m["Filling"].IsNumber());

    Population=m["Population"].IsInt() ? m["Population"].GetInt() : m["Filling"].GetDouble()*NSites;


    InitDistribution.resize(NSites);
    for(int p=0; p<Population; ++p) {
      unsigned int i=p%NSites;
      InitDistribution[i]++;
    }

    if(Nmax==0 && m.HasMember("InitDistribution") && m["InitDistribution"].GetString()==std::string("Gaussian")) {
      std::vector<double> g=lattice->NonInteractingGS();
      double norm=0;
      unsigned int index0=0;
      double max=0;
      for(unsigned int i=0; i<g.size(); ++i) {
        norm+=g[i];
        if(g[i]>max) {
          max=g[i];
          index0=i;
        }
      }
      unsigned int count=0;
      for(unsigned int i=0; i<g.size(); ++i) {
        unsigned int LocalPopulation=floor(Population*g[i]/norm);
        count+=LocalPopulation;
        InitDistribution[i]+=LocalPopulation;
      }
      if(count<Population) {
        InitDistribution[index0]+=Population-count;
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
    model->MakeContainer(Container);
    Container.g.initialize(model->nsites(),GreenOperatorLines);

  }

  std::ostream &print(std::ostream &o) const {
    int w=25;
    std::string indent="  ";
    o<<"Input:\n\n";

    model->print(o,indent);

    o<<indent<<std::setw(w)<<std::left<<"Alpha:"<<Alpha<<std::endl;
    o<<indent<<std::setw(w)<<std::left<<"GreenOperatorLines:"<<GreenOperatorLines<<std::endl;
    o<<indent<<std::setw(w)<<std::left<<"Seed:"<<Seed<<std::endl;
    o<<indent<<std::setw(w)<<std::left<<"NBins:"<<NBins<<std::endl;
    //o<<indent<<std::setw(w)<<std::left<<"WarmTime:"<<WarmTime<<std::endl;
    //o<<indent<<std::setw(w)<<std::left<<"WarmIterations:"<<WarmIterations<<std::endl;
    //o<<indent<<std::setw(w)<<std::left<<"MeasTime:"<<MeasTime<<std::endl;
    //o<<indent<<std::setw(w)<<std::left<<"MeasIterations:"<<MeasIterations<<std::endl;

    return o;
  }

  Parameters(std::string &json_input) : json(json_input) {

    rapidjson::Document d;
    d.Parse<0>(json.c_str());

    assert(d["SGF"].IsObject());
    const rapidjson::Value &sgf=d["SGF"];


    // The cutoff for the Green operator lines.
    GreenOperatorLines=sgf["GreenOperatorLines"].GetInt();

    // The alpha parameter for the updates
    Alpha=sgf["Alpha"].GetDouble();

    // Initial Seed for random number generator
    Seed=sgf["Seed"].GetInt();

    // Number of Binds for measurements
    NBins=sgf["Bins"].GetInt();

    // at least one of "WarmTime" and "WarmIterations" should exist
    assert(sgf["WarmTime"].IsInt() || sgf["WarmIterations"].IsInt());

    // at least one of "MeasTime" and "WarmIterations" should exist
    assert(sgf["MeasTime"].IsInt() || sgf["MeasIterations"].IsInt());

    WarmTime=sgf["WarmTime"].IsInt64() ? sgf["WarmTime"].GetInt64() : -1;
    WarmIterations=sgf["WarmIterations"].IsInt64() ? sgf["WarmIterations"].GetInt64() : -1 ;
    MeasTime=sgf["MeasTime"].IsInt64() ? sgf["MeasTime"].GetInt64() : -1;
    MeasIterations=sgf["MeasIterations"].IsInt64() ? sgf["MeasIterations"].GetInt64() : -1;

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
