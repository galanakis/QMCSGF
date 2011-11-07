#ifndef __ADJACENCYLIST__
#define __ADJACENCYLIST__
 
#include <vector>
#include <map>
#include <set>

#include "HamiltonianTerm.hh"
#include "Boson.hh"


/***********************************************************************\

In this header we define an adjacency list for Hamilonian terms. The class
AdjacencyList is the core structure, whereas AdjListDataElement and 
AdjListData are helper classes. The AdjacencyList is initialized from
a single or a pair of Hamiltonians:
AdjacencyList a(T,V);
AdjacencyLust a(T);
where H and V are of the type Hamiltonian (vectors of HamiltonianTerm).
The second syntax is the first one for V=T. The first syntact, builds 
an array of lists. The i^th list contains all the V terms that have one
index (site, particle flavor, band) in common with the i^th T term.
Each list member is a structure of type AdjListData, which contains
information necessary to evaluate the offset and the matrix element
of the V term if a T term is added or removed. For instance

a[i][j].offset(ADD);
a[i][j].me(REMOVE);

The first will calculate the offset of the j^th "neighbor" of the i^th
T term, if the latter is added. Similarly for the second line.

An estimate for the memory consumption is 32 bytes (256 bits) per entry.

\***********************************************************************/


/*
  class AdjListDataElement
  
  The purpose of this structure is to allow for an efficient way of 
  calculating the offset and matrix element of a product of creation
  annihilation operators A, after another such product B has been added 
  or removed to the right or left. To do this we record A together with 
  the change of occupancy, dn, caused by B. Whenever we want the amplitude
  of A we just add to the occupancy of the state the occupancy of B.
  
  Note that we keep a pointer to the IndexedProductElement rather than
  a copy to save memory (this structure will be used a lot!).
*/

namespace SGF {

class AdjListDataElement {
  short dn;
  const IndexedProductElement * ipeptr;
  inline int n(int direction,int action) const {return ipeptr->n(direction)+dn*Sign[direction==action];}
  inline Boson* particle_id() const {return ipeptr->particle_id();}
public:
  AdjListDataElement() {}
  AdjListDataElement(const IndexedProductElement *h,const IndexedProductElement *o) : dn(h->delta()), ipeptr(o) {}
  inline int offset(int action) const {return ipeptr->ProductElement::offset(particle_id()->delta()+dn*Sign[action]);}
  inline uint amplitude(int direction,int action) const {return ipeptr->ProductElement::amplitude(n(direction,action),particle_id()->nmax());}
};

/*
  class AdjListData
   
  is initialized with two HamiltonianTerm pointers h and o. It separates
  their particle indices to common and different ones. The common ones
  are stored in an array of AdjListDataElement (able to calculate offsets
  and matrix elements) and the different ones to an array of IndexedProductElement.
  I used a simple array rather than a vector, because the memory needs to remina
  as little as possible. I store a pointer of the "neighbor" o. From this I can
  lookup the coeffient along with the the original offsets and matrix elements.
  
  Note that I don't have to keep _ndifferent=neighbor->size()-_common.
  However due to structure alignment, it does not cost any extra memory.
  Do not change the order of the declarations!
*/

class AdjListData {

  typedef AdjListDataElement common_type;
  typedef const IndexedProductElement* diff_type;
  
  ushort _ndifferent;                  // The number of different indices
  ushort _ncommon;                     // the number of common indices
  common_type  *_common;               // Holds the list of common indices
  diff_type *_different;               // Holds the different indices
  const HamiltonianTerm *neighbour;    // Pointer to the term that is affected

  inline int ncommon() const {return _ncommon;}
  inline int ndifferent() const {return _ndifferent;}

public:

  AdjListData(const HamiltonianTerm *h,const HamiltonianTerm *o) : neighbour(o) {
    
      
    std::map<Boson*,const IndexedProductElement*> hmap,omap;
    std::map<Boson*,const IndexedProductElement*>::const_iterator hmit,omit;
    
    std::vector<IndexedProductElement>::const_iterator hit, oit;
    for(oit=o->product().begin();oit!=o->product().end();++oit)
      omap[oit->particle_id()]= &(*oit);
    for(hit=h->product().begin();hit!=h->product().end();++hit)
      hmap[hit->particle_id()]= &(*hit);
    
    std::vector<common_type> v_common;
    std::vector<diff_type> v_different;

    for(omit=omap.begin();omit!=omap.end();++omit) {
      hmit=hmap.find(omit->first);
      if(hmit!=hmap.end())
        v_common.push_back(AdjListDataElement(hmit->second,omit->second));
      else
        v_different.push_back(omit->second);
    }
    
    _ncommon=v_common.size();
    _ndifferent=v_different.size();
    
    _common=new common_type[ncommon()];
    for(int i=0;i<ncommon();++i)
      _common[i]=v_common[i];
      
    _different=new diff_type[ndifferent()];
    for(int i=0;i<ndifferent();++i)
      _different[i]=v_different[i];
    
  }
  
  AdjListData(const AdjListData &o) :_ndifferent(o._ndifferent),  _ncommon(o._ncommon), neighbour(o.neighbour) {
    
    _common=new common_type[ncommon()];
    _different=new diff_type[ndifferent()];

    for(int i=0;i<ncommon();++i)
      _common[i]=o._common[i];
    for(int i=0;i<ndifferent();++i)
      _different[i]=o._different[i];
    
  }
 
  ~AdjListData() {
    delete [] _common;
    delete [] _different;
  }
      
  inline const HamiltonianTerm* term() const {return neighbour;}
  
  // This is the offset of the "neighbor" term after the addition/removal of "this" term.
  inline int offset(int action) const {
    int result=0;
    for(int i=0;i<ncommon();++i)
      result+=_common[i].offset(action);
    for(int i=0;i<ndifferent();++i)
      result+=_different[i]->offset();
    return result;
  }

  // This is the original offset of the "neighbour".
  inline int offset() const {return neighbour->offset();}

  // This is the original matrix element of the "neighbor".
  inline MatrixElement me(int direction) const {return neighbour->me(direction);}

  // This is the matrix element of the "neighbour" after the addition/removal of "this" term.
  inline MatrixElement me(int direction,int action) const {
    int result=1;
    for(int i=0;i<ncommon();++i) 
      result*=_common[i].amplitude(direction,action); 
    for(int i=0;i<ndifferent();++i) 
      result*=_different[i]->amplitude(direction);
    return sqrt(result)*neighbour->coefficient();
  }  
    
};

typedef std::vector<AdjListData> adjacency_list_t;  // This is the type of the adjacency list

/*
  class AdjacencyList
  
  This is the desired class. It holds the data needed to calculate
  offsets and matrix elements.
*/

class AdjacencyList {
public:
  void initialize(const Hamiltonian &,const Hamiltonian &);
protected:
  std::vector<adjacency_list_t> _adjacency; // The adjacency list is stored here
public:  
  /* 
    Trow is the row terms and Tcol the column terms with which Trow have common indices.
  */  
  AdjacencyList(const Hamiltonian &Trow,const Hamiltonian &Tcol) { initialize(Trow,Tcol); }
  AdjacencyList(const Hamiltonian &T) { initialize(T,T); }
  
  // Total number of row terms
  std::vector<adjacency_list_t>::size_type size() const {return _adjacency.size();}
  const adjacency_list_t &operator[](Hamiltonian::size_type i) const {return _adjacency[i];}
   
};

/* 
   AdjacencyList initialization

   At first we categorize the Tcol terms by index and store them 
   in the map map_to_set. Then for each term in Trow, we merge 
   the sets of Tcol that correspond to its indices and store the 
   resulting set to "merged". Then for each element in merged we 
   construct the corresponding AdjListData and push them to the 
   corresponding vector in the _adjacency. 
*/

void AdjacencyList::initialize(const Hamiltonian &Trow,const Hamiltonian &Tcol)  {

  // Categorize the Tcol terms by index.
  std::map<Boson*,std::set<Hamiltonian::size_type> > map_to_set;
  for(Hamiltonian::size_type i=0;i<Tcol.size();++i) 
    for(Hamiltonian::size_type j=0;j<Tcol[i].product().size();++j)
      map_to_set[Tcol[i].product()[j].particle_id()].insert(i);
  
  _adjacency.clear();
  _adjacency.resize(Trow.size());
 
  /* For each term in Trow, merge the sets corresponding to it's indices
     Then, copy the set elements to a vector. */
  for(Hamiltonian::size_type i=0;i<Trow.size();++i) {
    std::set<Hamiltonian::size_type> merged;
    for(Hamiltonian::size_type j=0;j<Trow[i].product().size();++j) {
      Boson* pid=Trow[i].product()[j].particle_id();
      std::set<Hamiltonian::size_type> &s=map_to_set[pid];
      merged.insert(s.begin(),s.end());
    }
    
    _adjacency[i].reserve(merged.size());
    for(std::set<Hamiltonian::size_type>::const_iterator sit=merged.begin();sit!=merged.end();++sit)
      _adjacency[i].push_back(AdjListData(&Trow[i],&Tcol[*sit]));
    
  }
  
}   

}

#endif
