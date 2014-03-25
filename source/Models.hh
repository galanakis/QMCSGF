#ifndef __MODELS__
#define __MODELS__


#include <vector>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <cmath>

namespace SGF {

typedef enum {periodic,open,trap} boundary_t;

std::string bouncary_t_names[3]= {"periodic","open","trap"};

std::ostream &operator<<(std::ostream &o,const boundary_t b) {
  return o<< bouncary_t_names[b];
}

typedef std::pair<unsigned int,unsigned int> LinkData;
typedef std::vector<LinkData> LinkDataList;
typedef std::vector<double> SiteData;


class Lattice {
protected:
  LinkDataList _links;
  SiteData _sitedata;
  SiteData _gaussian;
  std::string label;
public:
  unsigned int size() const {return _sitedata.size();};
  virtual std::ostream &yaml_print(std::ostream &,const std::string &) const =0;
  const SiteData &sitedata() const {return _sitedata;}
  const SiteData &NonInteractingGS() const {return _gaussian;}
  LinkDataList links() const {return _links;}
  virtual ~Lattice() {};

};

template<typename T,unsigned int N>
struct Vector {
  T v[N];
  Vector(const Vector &o) {
    for(unsigned int i=0;i<N;++i)
      v[i]=o.v[i];
  }
  T operator[](unsigned int i) {return v[i];}
};

template<typename T,unsigned int N>
std::ostream & yaml_print_site(std::ostream &o,const Vector<T,N> &s) {
  o<<"[ ";
  for(int i=0;i<N-1;++i)
    o<<std::right<<std::setw(3)<<s[i]<<",";
  o<<std::right<<std::setw(3)<<s[N-1]<<" ]";
  return o;
}


struct Site1D {
  unsigned int x;
  Site1D() : x(0) {}
  Site1D(const Site1D &o) : x(o.x) {}
  Site1D(unsigned int a) : x(a) {}
};

struct Site2D {
  unsigned int x,y;
  Site2D() : x(0), y(0) {}
  Site2D(const Site2D &o) : x(o.x), y(o.y) {}
  Site2D(unsigned int a,unsigned int b) : x(a), y(b) {}
};

struct Site3D {
  unsigned x,y,z;
  Site3D() : x(0), y(0), z(0) {}
  Site3D(const Site3D &o) : x(o.x), y(o.y), z(o.z) {}
  Site3D(unsigned int a,unsigned int b,unsigned int c) : x(a), y(b), z(c) {}
};


struct DimensionData {
  unsigned int L;
  boundary_t B;
  double V;
  DimensionData(unsigned int _L,boundary_t _B,double _V) : L(_L), B(_B), V(_V) {}
  DimensionData(const DimensionData &o) : L(o.L), B(o.B), V(o.V) {}

  double cratio(unsigned int x) {
    return (2.0*x-L+1)/(L-1);
  }

  double parabola(unsigned int x) {
    double r=(2.0*x-L+1)/(L-1);
    return V*r*r;
  }

  double gaussian(unsigned int x) {
    double k=2.0*sqrt(V)/L;
    double x0=0.5*(L-1);
    return exp(-k*(x-x0)*(x-x0));
  }

};

std::ostream & yaml_print_site(std::ostream &o,const Site1D &s) {
  o<<"[ "<<std::right<<std::setw(3)<<s.x<<" ]";
  return o;
}

std::ostream & yaml_print_site(std::ostream &o,const Site2D &s) {
  o<<"[ "<<std::right<<std::setw(3)<<s.x<<","<<std::setw(3)<<s.y<<" ]";
  return o;
}

std::ostream & yaml_print_site(std::ostream &o,const Site3D &s) {
  o<<"[ "<<std::right<<std::setw(3)<<s.x<<","<<std::setw(3)<<s.y<<","<<std::setw(3)<<s.z<<" ]";
  return o;
}

template<class Site>
std::ostream & yaml_print_data(std::ostream &o,const std::string &indent,const std::string &label,const std::vector<DimensionData> &d,const std::vector<Site> &_site_coordinates,const SiteData &_sitedata) {
  o<<"\n";
  o<<indent<<"Lattice:\n\n";
  o<<indent<<"  Type: "<<label<<"\n";
  o<<indent<<"  Size: "<<_site_coordinates.size()<<"\n";

  o<<"\n";
  o<<indent<<"  Dims:\n";
  for(std::vector<DimensionData>::const_iterator it=d.begin(); it!=d.end(); ++it) {
    o<<indent<<"    - Length:      "<<it->L<<"\n";
    o<<indent<<"      Boundary:    "<<it->B<<"\n";
    o<<indent<<"      Trap Height: "<<it->V<<"\n";

  }

  o<<"\n";
  o<<indent<<"  Sites:\n";
  for(typename std::vector<Site>::const_iterator it=_site_coordinates.begin(); it!=_site_coordinates.end(); ++it) {
    o<<indent<<"    - ";
    yaml_print_site(o,*it);
    o<<"\n";
  }

  o<<"\n";
  o<<indent<<"  Potential:\n";
  for(unsigned int i=0; i<_sitedata.size(); ++i)
    o<<indent<<"    - "<<_sitedata[i]<<"\n";

  return o;

}


class Lattice1D : public Lattice {
protected:
  std::vector<DimensionData> dim;
  std::vector<Site1D> _site_coordinates;
  double parabola(unsigned int x) {
    return dim[0].parabola(x);
  }
  double gaussian(unsigned int x) {
    return dim[0].gaussian(x);
  }
  double relativeradious(unsigned x) {
    double X=dim[0].cratio(x);
    return X*X;
  }

};

class Lattice2D : public Lattice {
protected:
  std::vector<DimensionData> dim;
  std::vector<Site2D> _site_coordinates;
  double parabola(unsigned int x,unsigned int y) {
    return dim[0].parabola(x)+dim[1].parabola(y);
  }
  double gaussian(unsigned int x,unsigned int y) {
    return dim[0].gaussian(x)*dim[1].gaussian(y);
  }
  double relativeradious(unsigned int x,unsigned int y) {
    double X=dim[0].cratio(x);
    double Y=dim[1].cratio(y);
    return X*X+Y*Y;
  }

};

class Lattice3D : public Lattice {
protected:
  std::vector<DimensionData> dim;
  std::vector<Site3D> _site_coordinates;
  double parabola(unsigned int x,unsigned int y,unsigned int z) {
    return dim[0].parabola(x)+dim[1].parabola(y)+dim[2].parabola(z);
  }
  double gaussian(unsigned int x,unsigned int y,unsigned int z) {
    return dim[0].gaussian(x)*dim[1].gaussian(y)*dim[2].gaussian(z);
  }

  double relativeradious(unsigned int x,unsigned int y,unsigned int z) {
    double X=dim[0].cratio(x);
    double Y=dim[1].cratio(y);
    double Z=dim[2].cratio(z);
    return X*X+Y*Y+Z*Z;
  }

};

class Chain : public Lattice1D {

  static unsigned int periodic_site_index(unsigned int x,unsigned int Lx) {
    return x%Lx;
  }

public:

  Chain(unsigned int Lx,boundary_t Bx,double Vx) {

    label="Chain";
    dim.push_back(DimensionData(Lx,Bx,Vx));

    _sitedata.resize(Lx);
    _gaussian.resize(Lx);

    for(unsigned int x=0; x<Lx; ++x) {
      unsigned int i=periodic_site_index(x,Lx);
      unsigned int a=periodic_site_index(x+1,Lx);
      if(Bx==periodic || x+1!=Lx) _links.push_back(LinkData(i,a));
      _sitedata[i]=parabola(x);
      _gaussian[i]=gaussian(x);
      _site_coordinates.push_back(Site1D(x));
    }
  }

  std::ostream &yaml_print(std::ostream &o,const std::string &indent) const {
    return yaml_print_data(o,indent,label,dim,_site_coordinates,_sitedata);
  }

};


class Square : public Lattice2D {

  inline unsigned int periodic_site_index(unsigned int x,unsigned int Lx,unsigned int y,unsigned int Ly) {
    return x%Lx+Lx*(y%Ly);
  }

public:
  Square(unsigned int Lx,boundary_t Bx,unsigned int Ly,boundary_t By,double Vx,double Vy) {

    label="Square Lattice";
    dim.push_back(DimensionData(Lx,Bx,Vx));
    dim.push_back(DimensionData(Ly,By,Vy));

    _sitedata.resize(Lx*Ly);
    _gaussian.resize(Lx*Ly);

    for(unsigned int x=0; x<Lx; ++x) {
      for(unsigned int y=0; y<Ly; ++y) {
        unsigned int i=periodic_site_index(x,Lx,y,Ly);
        unsigned int a=periodic_site_index(x+1,Lx,y,Ly);
        unsigned int b=periodic_site_index(x,Lx,y+1,Ly);
        if(Bx==periodic || x+1!=Lx) _links.push_back(LinkData(i,a));
        if(By==periodic || y+1!=Ly) _links.push_back(LinkData(i,b));
        _sitedata[i]=parabola(x,y);
        _gaussian[i]=gaussian(x,y);
        _site_coordinates.push_back(Site2D(x,y));
      }
    }
  }

  std::ostream &yaml_print(std::ostream &o,const std::string &indent) const {
    return yaml_print_data(o,indent,label,dim,_site_coordinates,_sitedata);
  }

};


class Cubic : public Lattice3D {

  inline unsigned int periodic_site_index(unsigned int x,unsigned int Lx,unsigned int y,unsigned int Ly,unsigned int z,unsigned int Lz) {
    return x%Lx+Lx*((y%Ly)+Ly*(z%Lz));
  }

public:
  Cubic(unsigned int Lx,boundary_t Bx,unsigned int Ly,boundary_t By,unsigned int Lz,boundary_t Bz,double Vx,double Vy,double Vz) {

    label="Cubic Lattice";
    dim.push_back(DimensionData(Lx,Bx,Vx));
    dim.push_back(DimensionData(Ly,By,Vy));
    dim.push_back(DimensionData(Lz,Bz,Vz));

    _sitedata.resize(Lx*Ly*Lz);
    _gaussian.resize(Lx*Ly*Lz);

    for(unsigned int x=0; x<Lx; ++x) {
      for(unsigned int y=0; y<Ly; ++y) {
        for(unsigned int z=0; z<Lz; ++z) {
          unsigned int i=periodic_site_index(x,Lx,y,Ly,z,Lz);
          unsigned int a=periodic_site_index(x+1,Lx,y,Ly,z,Lz);
          unsigned int b=periodic_site_index(x,Lx,y+1,Ly,z,Lz);
          unsigned int c=periodic_site_index(x,Lx,y,Ly,z+1,Lz);
          if(Bx==periodic || x+1!=Lx) _links.push_back(LinkData(i,a));
          if(By==periodic || y+1!=Ly) _links.push_back(LinkData(i,b));
          if(Bz==periodic || z+1!=Lz) _links.push_back(LinkData(i,c));
          _sitedata[i]=parabola(x,y,z);
          _gaussian[i]=gaussian(x,y,z);
          _site_coordinates.push_back(Site3D(x,y,z));
        }
      }
    }
  }

  std::ostream &yaml_print(std::ostream &o,const std::string &indent) const {
    return yaml_print_data(o,indent,label,dim,_site_coordinates,_sitedata);
  }

};


class Ellipsis : public Lattice2D {

  void push_back(unsigned int i,unsigned int j) {
    _links.push_back(LinkData(i,j));
  }

public:
  Ellipsis(unsigned int Lx,unsigned int Ly,double Vx,double Vy) {

    label="Square Lattice, Elliptical Domain";
    dim.push_back(DimensionData(Lx,open,Vx));
    dim.push_back(DimensionData(Ly,open,Vy));

    unsigned int Invalid=Lx*Ly;
    std::vector<std::vector<unsigned int> > index_map(Lx,std::vector<unsigned int>(Ly,Invalid));

    unsigned int count=0;

    for(unsigned int x=0; x<Lx; ++x) {
      for(unsigned int y=0; y<Ly; ++y) {
        if(relativeradious(x,y)<=1) {
          _site_coordinates.push_back(Site2D(x,y));
          _sitedata.push_back(parabola(x,y));
          _gaussian.push_back(gaussian(x,y));
          index_map[x][y]=count;
          ++count;
        }
      }
    }

    for(std::vector<Site2D>::const_iterator it=_site_coordinates.begin(); it!=_site_coordinates.end(); ++it) {
      unsigned int x=it->x;
      unsigned int y=it->y;
      unsigned int i=index_map[x][y];
      if(x+1!=Lx && index_map[x+1][y]!=Invalid)
        push_back(i,index_map[x+1][y]);

      if(y+1!=Ly && index_map[x][y+1]!=Invalid)
        push_back(i,index_map[x][y+1]);

    }

  }

  std::ostream &yaml_print(std::ostream &o,const std::string &indent) const {
    return yaml_print_data(o,indent,label,dim,_site_coordinates,_sitedata);
  }


};


class Ellipsoid : public Lattice3D {

  void push_back(unsigned int i,unsigned int j) {
    _links.push_back(LinkData(i,j));
  }

public:
  Ellipsoid(unsigned int Lx,unsigned int Ly,unsigned int Lz,double Vx,double Vy,double Vz) {

    label="Cubic Lattice, Ellipsoid Domain";
    dim.push_back(DimensionData(Lx,open,Vx));
    dim.push_back(DimensionData(Ly,open,Vy));
    dim.push_back(DimensionData(Lz,open,Vz));

    unsigned int Invalid=Lx*Ly*Lz;
    std::vector<std::vector<std::vector<unsigned int> > > index_map(Lx,std::vector<std::vector<unsigned int> >(Ly,std::vector<unsigned int>(Lz,Invalid)));
    unsigned int count=0;


    for(unsigned int x=0; x<Lx; ++x) {
      for(unsigned int y=0; y<Ly; ++y) {
        for(unsigned int z=0; z<Lz; ++z) {
          if(relativeradious(x,y,z)<=1) {
            _site_coordinates.push_back(Site3D(x,y,z));
            _sitedata.push_back(parabola(x,y,z));
            _gaussian.push_back(gaussian(x,y,z));
            index_map[x][y][z]=count;
            ++count;
          }
        }
      }
    }

    for(std::vector<Site3D>::const_iterator it=_site_coordinates.begin(); it!=_site_coordinates.end(); ++it) {
      unsigned int x=it->x;
      unsigned int y=it->y;
      unsigned int z=it->z;

      unsigned int i=index_map[x][y][z];

      if(x+1!=Lx && index_map[x+1][y][z]!=Invalid)
        push_back(i,index_map[x+1][y][z]);

      if(y+1!=Ly && index_map[x][y+1][z]!=Invalid)
        push_back(i,index_map[x][y+1][z]);

      if(z+1!=Lz && index_map[x][y][z+1]!=Invalid)
        push_back(i,index_map[x][y][z+1]);

    }

  }

  std::ostream &yaml_print(std::ostream &o,const std::string &indent) const {
    return yaml_print_data(o,indent,label,dim,_site_coordinates,_sitedata);
  }


};


const long double GoldenRatio=1.61803398874989484820458683437L;


class Penrose : public Lattice {

  struct Polar {
    long double r;
    long double phi;
    Polar(long double _r,long double _phi) : r(_r), phi(_phi) {}
  };

  struct Site {
    long double x;
    long double y;
  public:
    Site(const Site &o) : x(o.x), y(o.y) {}
    Site(const Polar &p) : x(p.r*cos(p.phi)), y(p.r*sin(p.phi)) {}
    Site(long double _x,long double _y) : x(_x), y(_y) {}

    Site operator+(const Site &o) const {return Site(x+o.x,y+o.y);}
    Site operator-(const Site &o) const {return Site(x-o.x,y-o.y);}
    Site operator/(long double r) const {return Site(x/r,y/r);}

  };

  
  struct triangle {
    bool color;
    Site A,B,C;
    triangle(const triangle &o) : color(o.color), A(o.A), B(o.B), C(o.C) {}
    triangle(bool _color,Site _A,Site _B,Site _C) : color(_color), A(_A), B(_B), C(_C) {}

  };


  long double scale() const {
    return pow(1/GoldenRatio,numdivisions);
  }

  Site chop(const Site &s) const {
    Site result(s);
    result.x = fabs(s.x)<0.1*scale() ? 0 : s.x;
    result.y = fabs(s.y)<0.1*scale() ? 0 : s.y;
    return result;  
  }

  static std::list<triangle> subdivide(const std::list<triangle> &triangles) {
    std::list<triangle> result;
    for(std::list<triangle>::const_iterator it=triangles.begin(); it!=triangles.end(); ++it) {
      Site A=it->A;
      Site B=it->B;
      Site C=it->C;
      // Subdivide the red triangle
      if(it->color==0) {
        Site P=A+(B-A)/GoldenRatio;
        result.push_back(triangle(0,C,P,B));
        result.push_back(triangle(1,P,C,A));
      }
      // subdivide the blue triangle
      else {
        Site Q=B+(A-B)/GoldenRatio;
        Site R=B+(C-B)/GoldenRatio;
        result.push_back(triangle(1,R,C,A));
        result.push_back(triangle(1,Q,R,B));
        result.push_back(triangle(0,R,Q,A));
      }
    }

    return result;

  }

  unsigned int numdivisions;
  std::vector<Site> Coordinates;
  friend std::ostream& operator<< (std::ostream& stream, const Site &s);
  friend std::ostream& operator<< (std::ostream& stream, const triangle &s);
  friend bool operator<(const Site &a,const Site &b);

  unsigned long SiteHash(const Site &o) const {
    unsigned long X=100000*o.x;
    unsigned long Y=100000*o.y;
    return X*100000+Y;
  }

public:
  Penrose(unsigned int subdivisions) : numdivisions(subdivisions) {

    label=std::string("Penrose");

    std::list<triangle> triangles;
    // Create wheel of red triangles around the origin
    for(int i=0;i<10;++i) {
      int s=(i%2==0)? 1 : -1;
      Site B(Polar(1,(2*i+s)*M_PI/10));
      Site C(Polar(1,(2*i-s)*M_PI/10));
      triangles.push_back(triangle(0,Site(0,0),B,C));
    }

    for(int i=0;i<subdivisions;++i)
      triangles=subdivide(triangles);


    std::map<unsigned long,unsigned int> CoordinateMap;

    unsigned int index=0;
    for(std::list<triangle>::const_iterator it=triangles.begin(); it!=triangles.end(); ++it) {
      Site A=it->A;
      Site B=it->B;
      Site C=it->C;
      if(CoordinateMap.count(SiteHash(A))==0) {
        CoordinateMap[SiteHash(A)]=index;
        Coordinates.push_back(A);
        ++index;
      }     
      if(CoordinateMap.count(SiteHash(B))==0) {
        CoordinateMap[SiteHash(B)]=index;
        Coordinates.push_back(B);
        ++index;
      }     
      if(CoordinateMap.count(SiteHash(C))==0) {
        CoordinateMap[SiteHash(C)]=index;
        Coordinates.push_back(C);
        ++index;
      }     

    }

    std::set<LinkData> edges;

    for(std::list<triangle>::const_iterator it=triangles.begin(); it!=triangles.end(); ++it) {
      unsigned int iA=CoordinateMap[SiteHash(it->A)];
      unsigned int iB=CoordinateMap[SiteHash(it->B)];
      unsigned int iC=CoordinateMap[SiteHash(it->C)];

      edges.insert(LinkData(iA,iB));
      edges.insert(LinkData(iA,iC));

    }

    for(std::set<LinkData>::const_iterator it=edges.begin(); it!=edges.end(); ++it) {
      _links.push_back(*it);
    }


    unsigned int nsites=Coordinates.size();

    _sitedata.resize(nsites);
    _gaussian.resize(nsites);

  }


  // Tremmendous effort to simply align the decimal points.
  std::ostream & yaml_print_site(std::ostream &o,const Site &s) const {
    if(s.x>=0) {
      o<<" "<<std::setprecision(12)<<std::left<<std::setw(16)<<s.x;
    }
    else {
      o<<std::setprecision(12)<<std::left<<std::setw(17)<<s.x;
    }

    if(s.y>=0) {
      o<<" "<<std::setprecision(12)<<std::left<<std::setw(16)<<s.y<<" ]";
    }
    else {
      o<<std::setprecision(12)<<std::left<<std::setw(17)<<s.y<<" ]";
    }
    return o;
  }


  std::ostream & yaml_print(std::ostream &o,const std::string &indent) const {
      o<<"\n";
      o<<indent<<"Lattice:\n\n";
      o<<indent<<"  Type: "<<label<<"\n";
      o<<indent<<"  Size: "<<size()<<"\n";
      o<<indent<<"  Divisions: "<<numdivisions<<"\n";

      o<<"\n";
      o<<indent<<"  Sites:\n";
      for(typename std::vector<Site>::const_iterator it=Coordinates.begin(); it!=Coordinates.end(); ++it) {
        o<<indent<<"    - ";
        yaml_print_site(o,chop(*it));
        o<<"\n";
      }

      return o;
  }

};

std::ostream& operator<< (std::ostream& stream, const Penrose::Site &s) {
  return stream<<"("<<s.x<<","<<s.y<<")";
}

std::ostream& operator<< (std::ostream& stream, const Penrose::triangle &s) {
  return stream<<"["<<s.color<<","<<s.A<<","<<s.B<<","<<s.C<<"]";
}

bool operator<(const Penrose::Site &a,const Penrose::Site &b) {
  return (a.x==b.x)? (a.x<b.x) : (a.y<b.y);

}



}


#endif