/******************************************************************\

Copyright 2015 Dimitrios Galanakis

This file is part of QMCSGF

QMCSGF is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 3.

QMCSGF is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with QMCSGF.  If not, see <http://www.gnu.org/licenses/>.

\*******************************************************************/

#include <string>
#include <cstdio>
#include <iostream>
#include <vector>
#include <iomanip>
#include "rapidjson/document.h"

#include "SGF.hh"


static std::string
readInputFile( const char* path ) {
  FILE* file = fopen( path, "rb" );
  if ( !file )
    return std::string("");
  fseek( file, 0, SEEK_END );
  long size = ftell( file );
  fseek( file, 0, SEEK_SET );
  std::string text;
  char* buffer = new char[size + 1];
  buffer[size] = 0;
  if ( fread( buffer, 1, size, file ) == (unsigned long)size )
    text = buffer;
  fclose( file );
  delete[] buffer;
  return text;
}


namespace SGF {


std::ostream& operator<<(std::ostream& o, ensemble_t e) {
  return o << ((e == Canonical) ? "Canonical" : "GrandCanonical");
}


class Model {
public:
  Model() {}
  virtual ~Model() {}
  virtual std::ostream& print(std::ostream& o, const std::string&) const = 0;
  virtual void MakeContainer(SGFBase& Container) const = 0;
  virtual unsigned int nsites() const = 0;
};

struct BoseHubbard : public Model {

  double Beta;
  double Hoppingt;
  double CoulombU;
  unsigned int Population;
  ensemble_t ensemble;
  double mu;
  Lattice* lattice;
  int Nmax;

  std::vector<unsigned int> InitDistribution;

  unsigned int nsites() const {
    return lattice->size();
  }


  void MakeContainer(SGFBase& Container) const {

    Container.Beta = Beta;
    Container.Ensemble = ensemble;

    unsigned int NSites = lattice->size();

    if (Nmax != 0 && Population > NSites * Nmax) {
      std::cerr << "The population is larger than the capacity of the system!" << std::endl;
      exit(33);
    }

    Container.Psi.resize(NSites);

    for (unsigned int i = 0; i < NSites; ++i) {
      Container.Psi[i].nmax() = Nmax;
      Container.Psi[i].nL() = InitDistribution[i];
      Container.Psi[i].nR() = InitDistribution[i];
    }

    for (unsigned int i = 0; i < NSites; ++i) {
      SGFBase::AppendHamiltonianTerm(Container.V, SGFBase::CreateHamiltonianTerm(CoulombU / 2.0, CCAA, &Container.Psi[i]));
    }


    LinkDataList links = lattice->links();

    for (LinkDataList::const_iterator it = links.begin(); it != links.end(); ++it) {
      unsigned int i = it->first;
      unsigned int j = it->second;
      SGFBase::AppendHamiltonianTerm(Container.T, SGFBase::CreateHamiltonianTerm(Hoppingt, C, &Container.Psi[i], A, &Container.Psi[j]));
      SGFBase::AppendHamiltonianTerm(Container.T, SGFBase::CreateHamiltonianTerm(Hoppingt, C, &Container.Psi[j], A, &Container.Psi[i]));
    }

    std::vector<double> potential = lattice->sitedata();

    for (unsigned int i = 0; i < Container.Psi.size(); ++i) {
      SGFBase::AppendHamiltonianTerm(Container.V, SGFBase::CreateHamiltonianTerm(mu + potential[i], CA, &Container.Psi[i]));
    }

  }

  std::ostream& print(std::ostream& o, const std::string& indent) const {
    int w = 25;
    o << indent << std::setw(w) << std::left << "Model:" << "BoseHubbard\n";
    o << indent << std::setw(w) << std::left << "Beta:" << std::setprecision(10) << Beta << "\n";
    o << indent << std::setw(w) << std::left << "Temperature:" << std::setprecision(10) << 1.0 / Beta << "\n";
    o << indent << std::setw(w) << std::left << "Hoppingt:" << Hoppingt << "\n";
    o << indent << std::setw(w) << std::left << "CoulombU:" << CoulombU << "\n";
    o << indent << std::setw(w) << std::left << "Population:" << Population << "\n";
    o << indent << std::setw(w) << std::left << "Ensemble:" << ensemble << "\n";
    o << indent << std::setw(w) << std::left << "Mu:" << mu << "\n";
    o << indent << std::setw(w) << std::left << "Nmax:" << Nmax << "\n";
    //o<<indent<<std::setw(w)<<std::left<<"Boundary:"<<boundary<<"\n";
    lattice->yaml_print(o, "  ");
    o << std::endl;
    return o;
  }
  ~BoseHubbard() {
    delete lattice;
  }

  BoseHubbard(const rapidjson::Value& m) {

    // The hoping matrix element
    assert(m["Hoppingt"].IsNumber());
    Hoppingt = m["Hoppingt"].GetDouble();

    // Only one of "Beta" or "Temperature" should be present, but not both.
    assert(m["Beta"].IsNumber() ^ m["Temperature"].IsNumber());
    Beta = m.HasMember("Beta") ? m["Beta"].GetDouble() : 1.0 / m["Temperature"].GetDouble();

    // The ensemble can only be Canonical or GrandCanonical;
    assert(m["Ensemble"].GetString() == std::string("Canonical") || m["Ensemble"].GetString() == std::string("GrandCanonical"));
    ensemble = m["Ensemble"].GetString() == std::string("Canonical") ? Canonical : GrandCanonical;

    // If the ensemble is Canonical there should be a chemical potential
    mu = 0;
    if (m.HasMember("Mu")) {
      assert(m["Mu"].IsNumber());
      mu = m["Mu"].GetDouble();
    }

    // If Nmax is not determined assume it is zero (i.e. infinite)
    Nmax = m["Nmax"].IsInt() ? m["Nmax"].GetInt() : 0;

    // If Nmax==1 then the Coulomb interaction is meaningless and can be omitted.
    if (Nmax != 1) {
      CoulombU = m["CoulombU"].GetDouble();
    } else {
      CoulombU = 0;
    }


    /*
       Geometry specific assignments
    */
    assert(m["Lattice"].GetString() == std::string("Square") || m["Lattice"].GetString() == std::string("Penrose"));

    if (m["Lattice"].GetString() == std::string("Square")) {

      std::map<std::string, boundary_t> bcmap;
      bcmap["open"] = open;
      bcmap["periodic"] = periodic;
      bcmap["trap"] = trap;


      boundary_t boundary;
      std::string bcstring = m["Boundary"].GetString();
      assert(bcstring == std::string("open") || bcstring == std::string("periodic") || bcstring == std::string("trap"));
      boundary = bcmap[bcstring];

      assert(m["Dimensions"].IsArray());
      unsigned int dimensions = m["Dimensions"].Size();
      assert(dimensions <= 3 && dimensions > 0);
      std::vector<unsigned int> LSizes(dimensions);
      for (unsigned int i = 0; i < dimensions; ++i)
        LSizes[i] = m["Dimensions"][i].GetInt();

      std::vector<double> W(dimensions, 0);
      if (m.HasMember("ParabolicTrapHeight")) {
        assert(m["ParabolicTrapHeight"].IsArray() && m["ParabolicTrapHeight"].Size() == dimensions);
        for (unsigned int i = 0; i < dimensions; ++i)
          W[i] = m["ParabolicTrapHeight"][i].GetDouble();
      }

      if (boundary == trap) {
        switch (dimensions) {
        case 1: {
          lattice = new Chain(LSizes[0], open, W[0]);
          break;
        }
        case 2: {
          lattice = new Ellipsis(LSizes[0], LSizes[1], W[0], W[1]);
          break;
        }
        case 3: {
          lattice = new Ellipsoid(LSizes[0], LSizes[1], LSizes[2], W[0], W[1], W[2]);
          break;
        }
        }
      } else {
        switch (dimensions) {
        case 1: {
          lattice = new Chain(LSizes[0], boundary, W[0]);
          break;
        }
        case 2: {
          lattice = new Square(LSizes[0], boundary, LSizes[1], boundary, W[0], W[1]);
          break;
        }
        case 3: {
          lattice = new Cubic(LSizes[0], boundary, LSizes[1], boundary, LSizes[2], boundary, W[0], W[1], W[2]);
          break;
        }
        }
      }
    }
    else if (m["Lattice"].GetString() == std::string("Penrose")) {
      assert(m["Divisions"].IsInt());
      int divisions = m["Divisions"].GetInt();

      lattice = new Penrose(divisions);

    }

    unsigned int NSites = lattice->size();

    assert(m["Population"].IsInt() || m["Filling"].IsNumber());

    Population = m["Population"].IsInt() ? m["Population"].GetInt() : m["Filling"].GetDouble() * NSites;

    InitDistribution.resize(NSites, 0);

    if (Nmax == 0 && m.HasMember("InitDistribution") && m["InitDistribution"].GetString() == std::string("Gaussian")) {
      std::vector<double> g = lattice->NonInteractingGS();
      std::vector<double> gint(NSites + 1, 0);
      for (unsigned int i = 0; i < g.size(); ++i) {
        gint[i + 1] = gint[i] + g[i];
      }
      unsigned int P = Population;
      while (P > 0) {
        double rnd;
        do {
          rnd = double(rand()) / (RAND_MAX);
        } while (rnd == 1.0);

        double r = rnd * gint[NSites];

        //std::cout << rnd <<"\t"<<gint[NSites]<<"\t"<<NSites<<std::endl;

        unsigned int i = 0;
        while (r >= gint[i + 1]) {
          ++i;
        }

        ++InitDistribution[i];
        --P;

      }

    } else {
      for (unsigned int p = 0; p < Population; ++p) {
        unsigned int i = p % NSites;
        InitDistribution[i]++;
      }
    }

  }

};


struct SSL2D : public Model {

  double J1, J2, J3, J4;      // Hoping parameters
  double Delta;               // Anisotropy
  double Beta;                // Temperature
  ensemble_t ensemble;        // Allow the magnetization to fluctuate or not
  double Hz;                  // magnetic field.
  boundary_t boundary;        // periodic or open boundary conditions;
  unsigned int Lx, Ly;
  unsigned int NSites;
  unsigned int Population;

  std::vector< std::vector< unsigned int > > map;

  unsigned int nsites() const {
    return NSites;
  }


  SSL2D(const rapidjson::Value& m) {

    J1 = m["J1"].GetDouble();
    J2 = m["J2"].GetDouble();
    J3 = m["J3"].IsNumber() ? m["J3"].GetDouble() : 0;
    J4 = m["J4"].IsNumber() ? m["J4"].GetDouble() : 0;


    Delta = m["Delta"].GetDouble();

    Lx = m["Lx"].GetInt();
    Ly = m["Ly"].GetInt();

    // Only one of "Beta" or "Temperature" should be present, but not both.
    assert(m["Beta"].IsNumber() ^ m["Temperature"].IsNumber());
    Beta = m.HasMember("Beta") ? m["Beta"].GetDouble() : 1.0 / m["Temperature"].GetDouble();

    // The ensemble can only be Canonical or GrandCanonical;
    assert(m["Ensemble"].GetString() == std::string("Canonical") || m["Ensemble"].GetString() == std::string("GrandCanonical"));
    ensemble = m["Ensemble"].GetString() == std::string("Canonical") ? Canonical : GrandCanonical;

    // If the ensemble is Canonical there should be a chemical potential
    Hz = m["Hz"].IsNumber() ? m["Hz"].GetDouble() : 0;

    std::map<std::string, boundary_t> bcmap;
    bcmap["open"] = open;
    bcmap["periodic"] = periodic;

    std::string bcstring = m["Boundary"].GetString();
    assert(bcstring == std::string("open") || bcstring == std::string("periodic"));
    boundary = bcmap[bcstring];

    Population = m["Population"].GetDouble();

    NSites = Lx * Ly;

    map.resize(Lx);

    for (unsigned int x = 0; x < Lx; ++x)
      map[x].resize(Ly, nsites());

    unsigned int count = 0;

    for (unsigned int y = 0; y < Ly; ++y)
      for (unsigned int x = 0; x < Lx; ++x)
        map[x][y] = count++;


    if (Population > NSites) {
      std::cerr << "The magnetization is larger than the capacity of the system!" << std::endl;
      exit(33);
    }

    if ( (Population == 0 || Population == NSites) && ensemble == Canonical) {
      std::cerr << "Maximum/Minimum magnetization in the canonical ensemble. Unique configuration and meaningless to proceed." << std::endl;
      exit(33);
    }



  }

  std::ostream& print(std::ostream& o, const std::string& indent) const {
    int w = 25;
    o << indent << std::setw(w) << std::left << "Model:" << "HeisenbergSSL\n";
    o << indent << std::setw(w) << std::left << "Lx:" << Lx << "\n";
    o << indent << std::setw(w) << std::left << "Ly:" << Ly << "\n";
    o << indent << std::setw(w) << std::left << "NSites:" << NSites << "\n";
    o << indent << std::setw(w) << std::left << "Population:" << Population << "\n";
    o << indent << std::setw(w) << std::left << "Magnetization:" << double(Population) - 0.5 * NSites << "\n";
    o << indent << std::setw(w) << std::left << "Beta:" << std::setprecision(10) << Beta << "\n";
    o << indent << std::setw(w) << std::left << "Temperature:" << std::setprecision(10) << 1.0 / Beta << "\n";
    o << indent << std::setw(w) << std::left << "J1:" << J1 << "\n";
    o << indent << std::setw(w) << std::left << "J2:" << J2 << "\n";
    o << indent << std::setw(w) << std::left << "J3:" << J3 << "\n";
    o << indent << std::setw(w) << std::left << "J4:" << J4 << "\n";
    o << indent << std::setw(w) << std::left << "Delta:" << Delta << "\n";
    o << indent << std::setw(w) << std::left << "Ensemble:" << ensemble << "\n";
    if (ensemble == GrandCanonical)
      o << indent << std::setw(w) << std::left << "Hz:" << Hz << "\n";

    o << indent << std::setw(w) << std::left << "Boundary:" << boundary << "\n";
    o << std::endl;
    return o;
  }

  void insert_hoping(double t, unsigned int i, unsigned int j, std::vector<Boson>& psi, Hamiltonian& T, Hamiltonian& V) const {

    if (i != nsites() && j != nsites()) {

      double Hoppingt = 0.5 * Delta * t;

      SGFBase::AppendHamiltonianTerm(T, SGFBase::CreateHamiltonianTerm(Hoppingt, C, &psi[i], A, &psi[j]));
      SGFBase::AppendHamiltonianTerm(T, SGFBase::CreateHamiltonianTerm(Hoppingt, C, &psi[j], A, &psi[i]));

      SGFBase::AppendHamiltonianTerm(V, SGFBase::CreateHamiltonianTerm(t, CA, &psi[i], CA, &psi[j]));
    }
  }

  unsigned int getSite(unsigned int x, unsigned int y) const {

    if (boundary == periodic) {
      return map[x % Lx][y % Ly];
    } else {
      if (x < Lx && y < Ly) {
        return map[x][y];
      } else {
        return nsites();
      }
    }

  }

  void MakeContainer(SGFBase& Container) const {

    Container.Beta = Beta;
    Container.Ensemble = ensemble;

    Hamiltonian& T = Container.T;
    Hamiltonian& V = Container.V;
    std::vector<Boson>& psi = Container.Psi;


    psi.resize(NSites);

    std::vector<unsigned int> InitDistribution(NSites, 0);

    for (unsigned int p = 0; p < Population; ++p) {
      InitDistribution[p]++;
    }

    for (unsigned int i = 0; i < NSites; ++i) {
      psi[i].nmax() = 1;
      psi[i].nL() = InitDistribution[i];
      psi[i].nR() = InitDistribution[i];
    }

    // The nearest neighbor links
    if (J2 != 0) {
      for (unsigned int x = 0; x < Lx; ++x)
        for (unsigned int y = 0; y < Ly; ++y) {
          insert_hoping(J2, getSite(x, y), getSite(x + 1, y), psi, T, V);
          insert_hoping(J2, getSite(x, y), getSite(x, y + 1), psi, T, V);
        }
    }

    // The diagonal links
    if (J1 != 0) {
      for (unsigned int x = 0; x < Lx; x += 2)
        for (unsigned int y = 0; y < Ly; y += 2) {
          insert_hoping(J1, getSite(x + 1, y), getSite(x, y + 1), psi, T, V);
          insert_hoping(J1, getSite(x + 1, y + 1), getSite(x + 2, y + 2), psi, T, V);
        }
    }

    if (J3 != 0) {
      for (unsigned int x = 0; x < Lx; x += 2)
        for (unsigned int y = 0; y < Ly; y += 2) {
          insert_hoping(J3, getSite(x + 1, y), getSite(x + 2, y + 1), psi, T, V);
          insert_hoping(J3, getSite(x + 1, y + 1), getSite(x + 2, y), psi, T, V);
        }
    }

    if (J4 != 0) {
      for (unsigned int x = 0; x < Lx; x += 2)
        for (unsigned int y = 0; y < Ly; y += 2) {
          insert_hoping(J4, getSite(x, y), getSite(x + 2, y), psi, T, V);
          insert_hoping(J4, getSite(x + 2, y), getSite(x + 2, y + 2), psi, T, V);
        }
    }

    if (ensemble == GrandCanonical && Hz != 0) {
      for (unsigned long i = 0; i < nsites(); ++i) {
        SGFBase::AppendHamiltonianTerm(V, SGFBase::CreateHamiltonianTerm(-Hz, CA, &psi[i]));
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
  std::string inputfname;

  double Alpha;
  unsigned int GreenOperatorLines;
  int Seed;
  unsigned int NBins;
  unsigned long WarmTime;
  unsigned long WarmIterations;
  unsigned long MeasTime;
  unsigned long MeasIterations;

  std::string inputConfiguration;
  std::string outputConfiguration;

  std::vector<std::string> Measurables;

  double ExtraTermProbability;

  Model* model;

  std::string modelID;

  void read_configuration(SGFBase& Container) const {

    rapidjson::Document d;
    std::string json = readInputFile(inputConfiguration.c_str());

    if (json == "")
      return;

    d.Parse<0>(json.c_str());



    const rapidjson::Value& psi = d["configuration"];

    const rapidjson::Value& terms = d["operators"];

    unsigned long nsites = psi.Size();
    unsigned long nterms = terms.Size();

    for (unsigned long i = 0; i < nsites; ++i) {
      Container.Psi[i].nR() = Container.Psi[i].nL() = psi[i].GetInt();
    }

    enum {_Term, _Time, _Energy};

    for (unsigned long i = 0; i < nterms; ++i) {

      const rapidjson::Value& oo = terms[i];

      unsigned long iterm = oo[_Term].GetInt();
      CircularTime Time(oo[_Time].GetDouble());

      Container.OperatorCDL.push<LEFT>(Time, iterm);


    }

  }



  bool HasMeasurable(const std::string& s) const {
    bool result = false;
    for (unsigned int i = 0; i < Measurables.size(); ++i)
      result = result || (Measurables[i] == s);
    return result;
  }


  void MakeContainer(SGFBase& Container) const {

    Container.ExtraTermProbability=ExtraTermProbability;
    Container.Alpha = Alpha;
    model->MakeContainer(Container);
    Container.g.initialize(model->nsites(), GreenOperatorLines);

    read_configuration(Container);

  }

  void MakeMeasurables(SGFBase& Container, Measurable& MeasuredOperators) const {
    MeasuredOperators.InsertOperator("Potential Energy", Container.V);
    MeasuredOperators.InsertOperator("Kinetic Energy", Container.T);

    if (HasMeasurable("Number")) {
      MeasuredOperators.InsertOperator("Particle Number", SGFBase::GenerateNumberOperator(Container.Psi));
    }

    if (HasMeasurable("LocalDensity")) {
      MeasuredOperators.InsertLocalDensity("Local Density", Container.Psi);
    }

#ifdef WITHMKL
    if (HasMeasurable("DensityMatrixEigenSystem")) {
      MeasuredOperators.InsertDensityMatrixEigenSystem("Density Matrix EigenSystem", Container.Psi);
    }
#endif
    if (HasMeasurable("DensityMatrix"))
      MeasuredOperators.InsertFunnyDensityMatrix("Density Matrix", Container.Psi);

    if (HasMeasurable("DensityDensityMatrix"))
      MeasuredOperators.InsertDensityDensityMatrix("Density Density Matrix", Container.Psi);

    if (HasMeasurable("DensityMatrixSlow")) {
      MeasuredOperators.InsertDensityMatrix("Density Matrix", Container.Psi);
    }


  }

  std::ostream& print(std::ostream& o) const {
    int w = 25;
    std::string indent = "  ";
    o << "Input:\n\n";

    model->print(o, indent);

    o << indent << std::setw(w) << std::left << "Alpha:" << Alpha << std::endl;
    o << indent << std::setw(w) << std::left << "GreenOperatorLines:" << GreenOperatorLines << std::endl;
    o << indent << std::setw(w) << std::left << "Seed:" << Seed << std::endl;
    o << indent << std::setw(w) << std::left << "NBins:" << NBins << std::endl;

    return o;
  }

  Parameters(const std::string& fname) : inputfname(fname) {

    json=readInputFile(inputfname.c_str());

    rapidjson::Document d;
    d.Parse<0>(json.c_str());

    assert(d["SGF"].IsObject());
    const rapidjson::Value& sgf = d["SGF"];


    // The cutoff for the Green operator lines.
    GreenOperatorLines = sgf["GreenOperatorLines"].GetInt();

    // The alpha parameter for the updates
    Alpha = sgf["Alpha"].GetDouble();

    // Initial Seed for random number generator
    Seed = sgf["Seed"].GetInt();

    // Number of Binds for measurements
    NBins = sgf["Bins"].GetInt();

    // at least one of "WarmTime" and "WarmIterations" should exist
    assert(sgf["WarmTime"].IsInt() || sgf["WarmIterations"].IsInt());

    // at least one of "MeasTime" and "WarmIterations" should exist
    assert(sgf["MeasTime"].IsInt() || sgf["MeasIterations"].IsInt());

    WarmTime = sgf["WarmTime"].IsInt64() ? sgf["WarmTime"].GetInt64() : -1;
    WarmIterations = sgf["WarmIterations"].IsInt64() ? sgf["WarmIterations"].GetInt64() : -1 ;
    MeasTime = sgf["MeasTime"].IsInt64() ? sgf["MeasTime"].GetInt64() : -1;
    MeasIterations = sgf["MeasIterations"].IsInt64() ? sgf["MeasIterations"].GetInt64() : -1;

    // Names of the operators to measure
    const rapidjson::Value& measure = d["Measure"];
    assert(measure.IsArray());
    for (rapidjson::SizeType i = 0; i < measure.Size(); i++) {
      Measurables.push_back(measure[i].GetString());
    }

    ExtraTermProbability = sgf["ExtraTermProbability"].IsNumber() ? sgf["ExtraTermProbability"].GetDouble() : 0;
    assert(ExtraTermProbability>=0 && ExtraTermProbability<=1.0);

    inputConfiguration = "";
    if (sgf.HasMember("InitOperatorString") && sgf["InitOperatorString"].IsString())
      inputConfiguration = sgf["InitOperatorString"].GetString();

    outputConfiguration = "";
    if (sgf.HasMember("OutputOperatorString") && sgf["OutputOperatorString"].IsString())
      outputConfiguration = sgf["OutputOperatorString"].GetString();



    /*
       Model specific stuff.
       Those are relevant to the Hubbard model on a chain/square/cubic lattice
    */

    assert(d["Model"].IsObject());
    const rapidjson::Value& m = d["Model"];

    modelID = m["Label"].GetString();

    if (modelID == std::string("BoseHubbard")) {
      model = new BoseHubbard(d["Model"]);
    } else if (modelID == std::string("HeisenbergSSL")) {
      model = new SSL2D(d["Model"]);
    } else {
      std::cout << "Unknown Model" << std::endl;
      exit(33);
    }


  }


  ~Parameters() {
    delete model;
  }

};

}
