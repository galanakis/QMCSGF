#ifndef __ADJACENCYLIST__
#define __ADJACENCYLIST__
 
#include <vector>
#include <map>
#include <set>

#include "HamiltonianTerm.hh"


namespace SGF {


/*
  class AdjacencyList
  
  This is the desired class. It holds the data needed to calculate
  offsets and matrix elements.
*/

class AdjacencyList {
  void initialize(const Hamiltonian &,const Hamiltonian &);
	typedef std::vector<const HamiltonianTerm *> adjacency_list_t;  // This is the type of the adjacency list
  std::vector<adjacency_list_t> _adjacency; // The adjacency list is stored here
public:  


  /* 
    Trow is the row terms and Tcol the column terms with which Trow have common indices.
  */  
  AdjacencyList(const Hamiltonian &Trow,const Hamiltonian &Tcol) { initialize(Trow,Tcol); }
  AdjacencyList(const Hamiltonian &T) { initialize(T,T); }
	AdjacencyList(const AdjacencyList &o) : _adjacency(o._adjacency), _ranges(o._ranges) {}
  
	typedef adjacency_list_t::const_iterator const_iterator;
	typedef std::pair<const_iterator,const_iterator> range_type;

	std::vector<range_type> _ranges;
  
	inline const range_type &range(Hamiltonian::size_type i) const {return _ranges[i];}

   
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
			if(Trow[i].product()[j].delta()!=0) {
      	Boson* pid=Trow[i].product()[j].particle_id();
      	std::set<Hamiltonian::size_type> &s=map_to_set[pid];
      	merged.insert(s.begin(),s.end());
			}
    }
    
    _adjacency[i].reserve(merged.size());
    for(std::set<Hamiltonian::size_type>::const_iterator sit=merged.begin();sit!=merged.end();++sit)
      _adjacency[i].push_back(&Tcol[*sit]);
    
  }

	for(std::vector<adjacency_list_t>::const_iterator it=_adjacency.begin();it!=_adjacency.end();++it) {
		_ranges.push_back(range_type(it->begin(),it->end()));
	}

  
}   

}

#endif
