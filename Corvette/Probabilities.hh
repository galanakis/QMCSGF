#ifndef __PROBABILITIES__
#define __PROBABILITIES__

#include "RandomNumberGenerator.hh"
#include "HamiltonianTerm.hh"
#include "GreenOperator.hh"
#include "AdjacencyList.hh"
#include <vector>
#include <list>
#include <set>

namespace SGF {

uint RebuildFrequency=100;

/* Crop small doubles. Use them when you expect a finite double
   in which case you can disregard the small ones as numerical
   errors */
const double Tolerance=1e-6;  // This is the square root of the smallet double.
long double Crop(long double x) { return (fabs(x)>Tolerance)?x:0; }


/* 
  class OffsetMap
  Helper class. It defines a map that can convert offsets into consecutive
  integers and vice versa.
  
  For an operator product at the same index, the offset is defined as
  Abs(dn+NR-NL)-Abs(NR-NL) and for a given dn it takes values -dn,-dn+2,...,dn-2,dn.
  If we have products at different indices, the offset is additive. This reminds a 
  bit addition of spins. The maximum offset is the Abs(delta()) 
  (# of creation -# of annihilation operators) of the product.
  The step of the offset is always 2. Now if we have a lot of terms, some of them
  will have an even and some of them will have an odd and some others and even
  offset. If all of them are even then we can set the offset step=2, otherwise
  it will be 1.
  
  All this does not assume the case of hard core bosons, in which the number
  of offsets is limitted.
  In this case we need to consider Abs(dn+DN)-Abs(DN), -Nmax<=(DN=NR-NL)<=Nmax.
  It is easy tos show that given the constraint:
  |dn|-2*Min(|dn|,Nmax) <= Abs(|dn|+DN)-Abs(DN) <= |dn|
  
  This is handled by IndexedProductElement, which can return the minimum and
  maximum offset.
  Given that we can scan through all kinetic terms and generate a list of 
  offsets (between min and max with step 2). Then we can map them to
  sequential integers.

*/

class OffsetMap {
  std::vector<uint> _index;    // Map from an offset to consecutive integers
  std::vector<int> _offsets;   // Keeps a sorted list of the offsets

  inline int MinOffset() const {return _offsets[0];}    // Minimum possible offset
  inline int MaxOffset() const {return _offsets[_offsets.size()-1];}    // Maximum possible offset

  inline uint nindex() const {return MaxOffset()-MinOffset()+1;}

    
  void initialize(const Hamiltonian &T) {
    // Scan all kinetic terms to find all the lengths
    std::set<int> set_offsets;
    for(int i=0;i<T.size();++i)
      for(int offset=T[i].minoffset();offset<=T[i].maxoffset();offset+=2)
        set_offsets.insert(offset);

    _index.clear();
    _offsets.clear();
    _offsets.insert(_offsets.begin(),set_offsets.begin(),set_offsets.end());
     
    _index.reserve(nindex());
    for(int i=0;i<nindex();++i)
      _index.push_back(Infinity); // The default value is set to infinity

    std::vector<int>::const_iterator sit;
    int count=0;
    for(sit=_offsets.begin();sit!=_offsets.end();++sit)
      _index[(*sit)-MinOffset()]=count++;    

  }
  
public:

  OffsetMap(const Hamiltonian &T) {initialize(T);}
    
  // Number of different offsets
  inline uint size() const {return _offsets.size();}
  // OffsetMap(int offset) will return a consecutive integer
  inline uint operator()(int offset) const {return _index[offset-MinOffset()];}
  // OffsetMap[int i] will return the i^th smallest offset
  inline int operator[](int i) const {return _offsets[i];}
  
};


/*
  class Probabilities
  This is the core class of the algorithm. It is able to evalute
  the probabilities that are necessary for the SGF algorithm
  with directed updates.
  
  To use the class one needs to declare the variable Kinetic
  which is of type Hamiltonian and stores all Kinetic energy terms.
  Then it can be called to provide all sorts of probabilities.
  
  The class provices
  int choose(): randomly chooses a term
  update(int,int): which adds or removes a term
  weight(): returns the total probability of the tree (the weight of the root)

  The class is itinialized only by an array of type Hamiltonian.

  Typical Usage of the class:
  
  static Probabilites::Kinetic=Kinetic_Terms;   // Gives the path to the Kinetic energy Terms;
  Probabilities p();                            // Initializes the internal structures
  unsigned long p.choose(RIGHT);                // pick a term for RIGHT addition. Similarly for LEFT addition.
  p.update(term,RIGHT or LEFT,ADD or REMOVE);   // update the probabilities, BEFORE updating the Psi's!
  p.weight(RIGHT or LEFT);                      // gives the probability normalization constant
  p.rebuild();                                  // Fixes accumulated errors, by starting from scratch. This is slow
  
*/

class Probabilities {
protected:
  const Hamiltonian &Kinetic;      // Local copy of the kinetic operators
  const Hamiltonian &Potential;    // Local copy of the potential operators    
private:
  const AdjacencyList kin_adjacency;// For each term it holds a list of other kinetic terms with one common site
  const AdjacencyList pot_adjacency;// For each term it holds a list of other kinetic terms with one common site

  GreenOperator GF;                // Defines the Green operator function
  int _NBWL;                       // The number of broken world lines

  double Energies[2];              // energy of right and left state
  std::vector<Boson*> _indices;    // Vector containing all indices appearing in the row terms  
  std::set<Boson*> _broken_lines;  // A set of the boson indices that are broken.
  
  // The fast index should be the term
  MatrixElement *TSum;             // We will record the partial sums here.
  const OffsetMap offsets;               // This will map the offsets to consecutive integers

  long long NUpdates;              // Number of updates since last rebuild


  inline uint tsum_index(uint direction,int offset,uint iterm) const {return (direction*noffsets()+offsets(offset))*nterms()+iterm;}
  inline uint tsum_rawindex(uint direction,int indoffset,uint iterm) const {return (direction*noffsets()+indoffset)*nterms()+iterm;}
  
  inline uint TSumLength() const {return 2*nterms()*noffsets();}
  inline uint noffsets() const {return offsets.size();}
  inline uint nterms() const {return Kinetic.size();}
  inline double tsum(int direction,int offset,int iterm) const {return (iterm>=nterms())?0:Crop(TSum[tsum_index(direction,offset,iterm)]);}  
  inline double tsum(int direction,int offset) const { return Crop(TSum[tsum_index(direction,offset,0)]); }

/* Optimization hack: It assumes that the term index is the fastest index 
  This way you don't have to evaluate an expression to get the index. This
  hack alone increases speed by 25%.*/
  inline void tsum_update(int rl,int index,int offset,MatrixElement me) {
    if(me!=0) {
      index++;
      MatrixElement *ptr=&TSum[tsum_index(rl,offset,0)];
      while(index>0) {
        /*
          After setting the value of the sum, check if the sum is close to zero.
          If it is within the Tolerance then set if exactly to zero. This costs 5%
          but it may result in avoiding numerical crashes.
          //ptr[index-1]+=me;
        */
        ptr[index-1]=Crop(ptr[index-1]+=me);
        index>>=1;
      }
    }
  }


public:
  Probabilities(const Hamiltonian &T,const Hamiltonian &V) : Kinetic(T), Potential(V), kin_adjacency(Kinetic), pot_adjacency(Kinetic,Potential), offsets(Kinetic) {
     
    // Scan all kinetic terms to find all the indices
    std::set<Boson*> indexset;
    for(int i=0;i<Kinetic.size();++i)
      for(int j=0;j<Kinetic[i].product().size();++j)
        indexset.insert(Kinetic[i].product()[j].particle_id());
    _indices.insert(_indices.begin(),indexset.begin(),indexset.end());

    /* Initialize the Green operator function.
       The number of sites is just the number of different
       indices appearing in the Kinetic operators */
    GF.initialize(_indices.size());


    NUpdates=0;
    /* Finding the minimum of all Kinetic operator coefficients */
    MatrixElement MinCoefficient=T[0].coefficient();
    for(int term=0;term<T.size();++term)
      MinCoefficient=Min(MinCoefficient,T[term].coefficient());
    
    TSum=new MatrixElement[TSumLength()];
     
    rebuild();

  }
  
  ~Probabilities() { delete [] TSum; }

  inline void GreenInit(int nsites) {GF.initialize(nsites);}

  inline double Energy(int direction) const {return Energies[direction];}
  
  inline double G() const {return GF(_NBWL);};                  // The value of the Green Operator for the given broken lines
  inline double G(int offset) const {return GF(_NBWL+offset);}  // The value of the Green Operator given the total broken lines and the offset.
  inline int NBrokenLines() const {return _NBWL;}
  const std::set<Boson*> &ListBrokenLines() const {return _broken_lines;};  // A set of the boson indices that are broken.
  
  
  inline double weight(int rl) const {
    double s=0.0;
    for(int i=0;i<noffsets();++i) 
      s+=G(offsets[i])*Crop(TSum[tsum_rawindex(rl,i,0)]);
    return s;
  }

  /* Choose the offset first. This version
     calls TSum rather than weight and so it's faster */
  const HamiltonianTerm* choose(int rl) const {

    double R=RNG::Uniform()*weight(rl);
    int i=0;
    while((R-=G(offsets[i])*Crop(TSum[tsum_rawindex(rl,i,0)]))>=0)
      ++i;

    MatrixElement *ptr=&TSum[tsum_index(rl,offsets[i],0)]; 
    uint _nterms=nterms();

    int index=0;
    while(index<_nterms) {
      int indr=(index+1)<<1;
      int indl=indr-1;
      
      double w =Crop(ptr[index]);
      double wr=(indr<_nterms) ? Crop(ptr[indr]) : 0;
      double wl=(indl<_nterms) ? Crop(ptr[indl]) : 0;
      
      if(w*RNG::Uniform()>=wl+wr)
        return &Kinetic[index]; 
      else
        index=((wl+wr)*RNG::Uniform()>=wr)?indl:indr;
    }
    std::cout<<"Probabilities::choose has reached a dead end. Cannot choose term"<<std::endl;
    exit(13);
    return NULL;
  }
    
  // Evaluates the matrix elements and populates the trees
  void rebuild() {
     
     _NBWL=0;
     _broken_lines.clear();
     for(int i=0;i<_indices.size();++i) {
      int delta=Abs(_indices[i]->delta());
      if(delta!=0) {
        _NBWL+=delta; 
        _broken_lines.insert(_indices[i]);
       }
     }

    for(int i=0;i<TSumLength();++i) TSum[i]=0;

    for(int rl=0;rl<2;++rl) 
      for(int i=0;i<nterms();++i) 
        tsum_update(rl,i,Kinetic[i].offset(),Kinetic[i].me(rl));
        
    Energies[LEFT]=Energies[RIGHT]=0;  
    for(int direction=0;direction<2;++direction)
      for(uint i=0;i<Potential.size();++i)
        Energies[direction]+=Potential[i].me(direction);
 
    NUpdates=0;

  }
  
  // Copy the TSum, rebuild and compare
  bool verify() {
    
    double Tolerance=1e-10;
    
    MatrixElement *TSumCopy=new MatrixElement[TSumLength()];
    for(int i=0;i<TSumLength();++i)
      TSumCopy[i]=TSum[i];
      
    double EnergiesCopy[2];
    EnergiesCopy[LEFT]=Energies[LEFT];
    EnergiesCopy[RIGHT]=Energies[RIGHT];
    
    int _NBWLCopy=_NBWL;
    std::set<Boson*> _broken_lines_copy=_broken_lines;
    
    rebuild();
    MatrixElement result=0;
    for(int i=0;i<TSumLength();++i)
      result+=Abs(TSumCopy[i]-TSum[i]);

    double ediff=Abs(EnergiesCopy[LEFT]-Energies[LEFT])+Abs(EnergiesCopy[RIGHT]-Energies[RIGHT]);
    
    bool tree_success=(result>Tolerance);
    bool energy_success=(ediff>Tolerance);
    bool worldline_success=(_NBWL!=_NBWLCopy); 
    bool broken_line_success=(_broken_lines==_broken_lines_copy);
     
    if(tree_success) std::cout<<"Probabilities: Tree verification failed"<<std::endl;
    else std::cout<<"Probabilities: Tree verification succeeded"<<std::endl;
    
    if(energy_success) std::cout<<"Probabilities: Energy verification failed"<<std::endl;
    else std::cout<<"Probabilities: Energy verification succeeded"<<std::endl;

    if(worldline_success) std::cout<<"Probabilities: Worldline verification failed "<<_NBWL<<", "<<_NBWLCopy<<std::endl;
    else std::cout<<"Probabilities: Worldline verification succeeded"<<std::endl;

    if(broken_line_success) std::cout<<"Broken line verification succeeded"<<std::endl;
    else std::cout<<"Broken line verification failed"<<std::endl;

    
    for(int i=0;i<offsets.size();++i)
      std::cout<<tsum(LEFT,offsets[i])<<", "<<tsum(RIGHT,offsets[i]);
    std::cout<<std::endl;

    delete [] TSumCopy;
    
    return tree_success && energy_success && worldline_success;
  }
 
  inline int term_index(const HamiltonianTerm* term) const {return term-&Kinetic[0];}

  inline void update(const HamiltonianTerm* term,int rl,int arflag) {update(term_index(term),rl,arflag);}
  
  inline void update(int index,int rl,int arflag) { 
    const HamiltonianTerm* term=&Kinetic[index];
    _NBWL+=term->offset(arflag);

    adjacency_list_t::const_iterator nbr;
    for(nbr=kin_adjacency[index].begin();nbr!=kin_adjacency[index].end();++nbr) {
      int fndex= term_index(nbr->term());
      int ioffset = nbr->offset();
      int foffset = nbr->offset(arflag);
      MatrixElement ime = nbr->me(rl);
      MatrixElement fme = nbr->me(rl,arflag);
      // The if statements decrease the time by 30%
      if(ioffset!=foffset) {
        MatrixElement jme =nbr->me(!rl);
        // Note that those statements are parallelizable
        tsum_update( rl,fndex,ioffset,-ime);
        tsum_update( rl,fndex,foffset,+fme);
        tsum_update(!rl,fndex,ioffset,-jme);
        tsum_update(!rl,fndex,foffset,+jme);
      }
      else 
        tsum_update( rl,fndex,ioffset,fme-ime);

    }
    
    for(nbr=pot_adjacency[index].begin();nbr!=pot_adjacency[index].end();++nbr) 
      Energies[rl]+=nbr->me(rl,arflag)-nbr->me(rl);
      
    term->update_psi(rl,arflag);   
    
    for(int i=0;i<term->product().size();++i) {
        Boson *pind=term->product()[i].particle_id();
        if(pind->delta()!=0)
           _broken_lines.insert(pind);
        else
           _broken_lines.erase(pind);
    }
    

    /* This needs to be confirmed: the tree needs to be rebuilt from time
      to time to fix accumulated floating points errors */
    if((++NUpdates % RebuildFrequency)==0)
      rebuild();

  }
 

};        

}

#endif
