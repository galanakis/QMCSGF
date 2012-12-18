#ifndef __MODELS__
#define __MODELS__


#include <vector>


namespace SGF {


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
    double rx=(2.0*x-Lx+1)/(Lx-1);
    unsigned int i=periodic_site_index(x,Lx);
    result[i]=rx*rx*Vx;
  }

  return result;

}

std::vector<double> trap(unsigned int Lx,double Vx,unsigned int Ly,double Vy) {

  std::vector<double> result(Lx*Ly);
  for(unsigned int x=0; x<Lx; ++x) {
    double rx=(2.0*x-Lx+1)/(Lx-1);
    for(unsigned int y=0; y<Ly; ++y) {
      double ry=(2.0*y-Ly+1)/(Ly-1);
      unsigned int i=periodic_site_index(x,Lx,y,Ly);
      result[i]=rx*rx*Vx+ry*ry*Vy;
    }
  }

  return result;

}

std::vector<double> trap(unsigned int Lx,double Vx,unsigned int Ly,double Vy,unsigned int Lz,double Vz) {

  std::vector<double> result(Lx*Ly*Lz);
  for(unsigned int x=0; x<Lx; ++x) {
    double rx=(2.0*x-Lx+1)/(Lx-1);
    for(unsigned int y=0; y<Ly; ++y) {
      double ry=(2.0*y-Ly+1)/(Ly-1);
      for(unsigned int z=0; z<Lz; ++z) {
        double rz=(2.0*z-Lz+1)/(Lz-1);
        unsigned int i=periodic_site_index(x,Lx,y,Ly,z,Lz);
        result[i]=rx*rx*Vx+ry*ry*Vy+rz*rz*Vz;
      }
    }
  }

  return result;

}



}


#endif