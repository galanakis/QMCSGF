#ifndef __MODELS__
#define __MODELS__


#include <vector>
#include <sstream>
#include <iomanip>


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
public:
  unsigned int size() const {return _sitedata.size();};
  virtual std::ostream &yaml_print(std::ostream &,const std::string &) const =0;
  const SiteData &sitedata() const {return _sitedata;}
  LinkDataList links() const {return _links;}
  virtual ~Lattice() {};

};

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
    o<<indent<<"    - Length: "<<it->L<<"\n";
    o<<indent<<"      Boundary: "<<it->B<<"\n";
    o<<indent<<"      Trap: "<<it->V<<"\n";

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

class Chain : public Lattice {

  static unsigned int periodic_site_index(unsigned int x,unsigned int Lx) {
    return x%Lx;
  }

  std::string label;
  std::vector<DimensionData> dim;
  std::vector<Site1D> _site_coordinates;

public:

  Chain(unsigned int Lx,boundary_t Bx,double Vx) {

    label="Chain";
    dim.push_back(DimensionData(Lx,Bx,Vx));

    _sitedata.resize(Lx);

    for(unsigned int x=0; x<Lx; ++x) {
      unsigned int i=periodic_site_index(x,Lx);
      unsigned int a=periodic_site_index(x+1,Lx);
      if(Bx==periodic || x+1!=Lx) _links.push_back(LinkData(i,a));
      double rx=(2.0*x-Lx+1)/(Lx-1);
      _sitedata[i]=rx*rx*Vx;
      _site_coordinates.push_back(Site1D(x));
    }
  }

  std::ostream &yaml_print(std::ostream &o,const std::string &indent) const {
    return yaml_print_data(o,indent,label,dim,_site_coordinates,_sitedata);
  }

};


class Square : public Lattice {

  inline unsigned int periodic_site_index(unsigned int x,unsigned int Lx,unsigned int y,unsigned int Ly) {
    return x%Lx+Lx*(y%Ly);
  }

  std::string label;
  std::vector<DimensionData> dim;
  std::vector<Site2D> _site_coordinates;

public:
  Square(unsigned int Lx,boundary_t Bx,double Vx,unsigned int Ly,boundary_t By,double Vy) {

    label="Square Lattice";
    dim.push_back(DimensionData(Lx,Bx,Vx));
    dim.push_back(DimensionData(Ly,By,Vy));

    _sitedata.resize(Lx*Ly);

    for(unsigned int x=0; x<Lx; ++x) {
      for(unsigned int y=0; y<Ly; ++y) {
        unsigned int i=periodic_site_index(x,Lx,y,Ly);
        unsigned int a=periodic_site_index(x+1,Lx,y,Ly);
        unsigned int b=periodic_site_index(x,Lx,y+1,Ly);
        if(Bx==periodic || x+1!=Lx) _links.push_back(LinkData(i,a));
        if(By==periodic || y+1!=Ly) _links.push_back(LinkData(i,b));
        double rx=(2.0*x-Lx+1)/(Lx-1);
        double ry=(2.0*y-Ly+1)/(Ly-1);
        _sitedata[i]=rx*rx*Vx+ry*ry*Vy;
        _site_coordinates.push_back(Site2D(x,y));
      }
    }
  }

  std::ostream &yaml_print(std::ostream &o,const std::string &indent) const {
    return yaml_print_data(o,indent,label,dim,_site_coordinates,_sitedata);
  }

};


class Cubic : public Lattice {

  inline unsigned int periodic_site_index(unsigned int x,unsigned int Lx,unsigned int y,unsigned int Ly,unsigned int z,unsigned int Lz) {
    return x%Lx+Lx*((y%Ly)+Ly*(z%Lz));
  }

  std::string label;
  std::vector<DimensionData> dim;
  std::vector<Site3D> _site_coordinates;


public:
  Cubic(unsigned int Lx,boundary_t Bx,double Vx,unsigned int Ly,boundary_t By,double Vy,unsigned int Lz,boundary_t Bz,double Vz) {

    label="Cubic Lattice";
    dim.push_back(DimensionData(Lx,Bx,Vx));
    dim.push_back(DimensionData(Ly,By,Vy));
    dim.push_back(DimensionData(Lz,Bz,Vz));

    _sitedata.resize(Lx*Ly*Lz);

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
          double rx=(2.0*x-Lx+1)/(Lx-1);
          double ry=(2.0*y-Ly+1)/(Ly-1);
          double rz=(2.0*z-Lz+1)/(Lz-1);
          _sitedata[i]=rx*rx*Vx+ry*ry*Vy+rz*rz*Vz;
          _site_coordinates.push_back(Site3D(x,y,z));
        }
      }
    }
  }

  std::ostream &yaml_print(std::ostream &o,const std::string &indent) const {
    return yaml_print_data(o,indent,label,dim,_site_coordinates,_sitedata);
  }

};


class Ellipsis : public Lattice {

  std::string label;
  std::vector<DimensionData> dim;
  std::vector<Site2D> _site_coordinates;

  void push_back(unsigned int i,unsigned int j) {
    _links.push_back(LinkData(i,j));
  }

public:
  Ellipsis(unsigned int Lx,unsigned int Ly,double V) {

    label="Square Lattice, Elliptical Domain";
    dim.push_back(DimensionData(Lx,open,V));
    dim.push_back(DimensionData(Ly,open,V));

    unsigned int Invalid=Lx*Ly;
    std::vector<std::vector<unsigned int> > index_map(Lx,std::vector<unsigned int>(Ly,Invalid));

    unsigned int count=0;

    for(unsigned int x=0; x<Lx; ++x) {
      for(unsigned int y=0; y<Ly; ++y) {
        double rx=(2.0*x-Lx+1)/(Lx-1);
        double ry=(2.0*y-Ly+1)/(Ly-1);
        double rr=rx*rx+ry*ry;
        if(rr<=1) {
          _site_coordinates.push_back(Site2D(x,y));
          _sitedata.push_back(rr*V);
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


class Ellipsoid : public Lattice {


  std::string label;
  std::vector<DimensionData> dim;
  std::vector<Site3D> _site_coordinates;

  void push_back(unsigned int i,unsigned int j) {
    _links.push_back(LinkData(i,j));
  }

public:
  Ellipsoid(unsigned int Lx,unsigned int Ly,unsigned int Lz,double V) {

    label="Cubic Lattice, Ellipsoid Domain";
    dim.push_back(DimensionData(Lx,open,V));
    dim.push_back(DimensionData(Ly,open,V));
    dim.push_back(DimensionData(Lz,open,V));

    unsigned int Invalid=Lx*Ly*Lz;
    std::vector<std::vector<std::vector<unsigned int> > > index_map(Lx,std::vector<std::vector<unsigned int> >(Ly,std::vector<unsigned int>(Lz,Invalid)));
    unsigned int count=0;


    for(unsigned int x=0; x<Lx; ++x) {
      for(unsigned int y=0; y<Ly; ++y) {
        for(unsigned int z=0; z<Lz; ++z) {
          double rx=(2.0*x-Lx+1)/(Lx-1);
          double ry=(2.0*y-Ly+1)/(Ly-1);
          double rz=(2.0*z-Lz+1)/(Lz-1);
          double rr=rx*rx+ry*ry+rz*rz;
          if(rr<=1) {
            _site_coordinates.push_back(Site3D(x,y,z));
            _sitedata.push_back(rr*V);
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


}


#endif