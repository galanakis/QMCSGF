#ifndef __MEASURABLE__
#define __MEASURABLE__

#include <map>
#include <vector>
#include <set>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>


#include "HamiltonianTerm.hh"
#include "Accumulator.hh"
#include "OperatorString.hh"


#ifdef USEMPI
// *******************************************
// * Declaration of global variables for MPI *
// *******************************************
#include <mpi.h>
#define Master 0
int NumProcessors,Rank,NameLength;
char ProcessorName[MPI_MAX_PROCESSOR_NAME];

#endif


namespace SGF {
void InitializeEnvironment(int NumArg,char *Arg[]) {
#ifdef USEMPI
  // ***********************************
  // * Initialization of MPI functions *
  // ***********************************

  MPI_Init(&NumArg,&Arg);
  MPI_Comm_size(MPI_COMM_WORLD,&NumProcessors);
  MPI_Comm_rank(MPI_COMM_WORLD,&Rank);
  MPI_Get_processor_name(ProcessorName,&NameLength);

  // Silence the other nodes. Only the root node can print
  std::cout.rdbuf( Rank==Master ? std::cout.rdbuf() : 0 );
  std::cerr.rdbuf( Rank==Master ? std::cout.rdbuf() : 0 );

#endif
}

void FinalizeEnvironment() {
#ifdef USEMPI
  MPI_Finalize();
#endif

}

}


namespace SGF {

// How to print the Bins
template<class T> inline std::ostream& operator<<(std::ostream& o, const BinnedAccumulator<T> &b) {
  return o<<"[ "<<std::fixed<<std::setprecision(12)<<std::setw(15)<<std::right<<b.average()<<", "<<std::setw(20)<<b.sigma()<<" ]";
}

typedef BinnedAccumulator<_float_accumulator> BinnedAccumulatorME;

class MeasurableObject {
protected:
  std::string _tag;
public:
  MeasurableObject(const std::string &tag) :_tag(tag) {}
  virtual ~MeasurableObject() {}
  virtual void push(const _float_accumulator &Weight) =0;
  virtual std::ostream& print(std::ostream& o) const = 0;
};

inline std::ostream& operator<<(std::ostream& o,const MeasurableObject &m) {
  return m.print(o);
}

class MeasurableFunction : public MeasurableObject {
protected:
  BinnedAccumulatorME _bin;
  virtual _float_accumulator evaluate() const = 0;
public:
  MeasurableFunction(const std::string &tag) : MeasurableObject(tag) {}
  void push(const _float_accumulator &Weight) {
#ifdef USEMPI
    if(Rank==Master)
#endif
      _bin.push(evaluate()/Weight);
  }

  std::ostream &print(std::ostream &o) const {
    return o<<"  "<<std::setw(25)<<std::left<<_tag+": "<<_bin;
  }


};

class MeasurableSum : public MeasurableFunction {
  typedef std::vector<std::pair<_float_accumulator*,MatrixElement> > vector_pair_t;
  vector_pair_t data;

public:
  MeasurableSum(const std::string &tag) : MeasurableFunction(tag) {}
  inline _float_accumulator evaluate() const {
    _float_accumulator result=0;
    for(vector_pair_t::const_iterator it=data.begin(); it!=data.end(); ++it)
      result += *it->first * it->second;
    return result;
  }
  void push_back(_float_accumulator *a,MatrixElement m) {
    data.push_back(std::pair<_float_accumulator*,MatrixElement>(a,m));
  }

};

class MeasurableNumber : public MeasurableFunction {
  _float_accumulator* data;
  MatrixElement coefficient;
public:
  MeasurableNumber(const std::string &tag,_float_accumulator* a,MatrixElement c=1.0) : MeasurableFunction(tag), data(a), coefficient(c) {}
  inline _float_accumulator evaluate() const {
    return *data * coefficient;
  }
};


class ExtraMeasurables : public MeasurableObject {
public:
  ExtraMeasurables(std::string tag) : MeasurableObject(tag) {}
  typedef std::vector<std::pair<Boson*,int> > KeyType;
  virtual void reset() =0;
  virtual void measure(_float_accumulator Weight,const KeyType &key)=0;
};

class MeasurableDensityMatrix : public ExtraMeasurables {

  std::vector<BinnedAccumulator<double> > BinsElements;

  const Boson *Psi0;
  typedef std::vector<Boson>::size_type size_type;
  size_type N;
  std::string tag;
  std::vector<double> data;
  double _number;
  BinnedAccumulator<long double> _numberbin;

  _float_accumulator me_ca(size_type i,size_type j) const {
    return sqrt( (Psi0[i].nR()+1)*Psi0[j].nR() );
  }

  size_type index(size_type i,size_type j) const {
    return (i<j) ? (2*N-i+1)*i/2+j-i : (2*N-j+1)*j/2+i-j;
  }

public:
  MeasurableDensityMatrix(const std::string &tag,std::vector<Boson> &Psi) : ExtraMeasurables(tag), _numberbin() {
    Psi0=&Psi[0];
    N=Psi.size();
    size_type ndata=N*(N+1)/2;
    data.resize(ndata,0);
    BinsElements.resize(ndata);
    _number=0;
  }

  void measure(_float_accumulator Weight,const KeyType &key) {

    if(key.size()==0) {
      double count=0;
      for(unsigned long i=0; i<N; ++i) {
        unsigned long ind = index(i,i);
        data[ind]+=Psi0[i].nR()*Weight;
        count+=Psi0[i].nR();
      }
      _number+=count*Weight;
    } else if(key.size()==2) {

      if(key[0].second==1 && key[1].second==-1) {
        size_type j=key[0].first-Psi0;
        size_type i=key[1].first-Psi0;
        data[index(i,j)]+=me_ca(i,j)*Weight/2.0;
      }

      if(key[1].second==1 && key[0].second==-1) {
        size_type j=key[1].first-Psi0;
        size_type i=key[0].first-Psi0;
        data[index(i,j)]+=me_ca(i,j)*Weight/2.0;
      }

    }

  }

  void reset() {
    std::fill(data.begin(),data.end(),0);
    _number=0;
  }


  void reduce() {

#ifdef USEMPI

    std::vector<double> send_data(data);
    MPI_Reduce(&send_data[0],&data[0],data.size(),MPI_LONG_DOUBLE,MPI_SUM,Master,MPI_COMM_WORLD);

#endif

  }

  void push(const _float_accumulator &Weight) {

    reduce();

    for(unsigned long i=0; i<data.size(); ++i)
      BinsElements[i].push(data[i]/Weight);

    _numberbin.push(_number/Weight);

  }



  std::ostream& print(std::ostream& o) const  {

    o<<"\n";
    o<<"  Number: "<<_numberbin<<std::endl;
    o<<"\n";
    o<<"  "<<_tag<<":";
    o<<"\n";

    for(unsigned long i=0; i<N; ++i)
      for(unsigned long j=i; j<N; ++j) {
        unsigned long ind=index(i,j);
        o<<"    - "<<BinsElements[ind]<<std::endl;
      }

    return o;
  }



};


/*

	class Measurable

	It deals with the measurements. It stores the buffers of all the quantities
	that need to be tracked. It also holds the bins of the quantities that are
	being measured. Note that the quantities being tracked are not the same
	as the quantities being measured. For example to measure the total
	potential energy we need to track all the potential energy terms, but
	we only need one bin.

	Notes for MPI parallelization.

	The BinnedAccumulatorME variables should only be defined in the Root node.
 	Also the print and flush methods should only be called by the root.
	During the flush operation the mpi reduce should be called for every bin.

*/


class Measurable  {

  struct TermBuffer {
    _float_accumulator buffer;
    std::vector<IndexedProductElement> product;
    TermBuffer(const std::vector<IndexedProductElement> &_p) : buffer(0), product(_p) {}
    TermBuffer(const TermBuffer &o) : buffer(o.buffer), product(o.product) {}
    MatrixElement me() const {
      return MultiplyMe<IndexedProductElement,RIGHT>(product);
    }
  };

  typedef BrokenLines::BosonDeltaMapType KeyType;
  typedef std::multimap<KeyType,TermBuffer> multimap_type;
  typedef std::pair<multimap_type::iterator,multimap_type::iterator> equal_range_type;

//  const OperatorStringType &OperatorString; // Holds a reference for the operator string being measured
  BrokenLines BrokenLineTracer;  // It traces the list of broken lines


  _float_accumulator buffer_BoltzmannWeight;
  _float_accumulator buffer_kinetic;
  _float_accumulator buffer_potential;

  multimap_type _multimap;       // It stores the buffers of all measurable terms

  // It stores the addresses of all _multimap's buffers together with the kinetic and potential buffers
  std::vector<_float_accumulator*> _buffers;

  MeasurableNumber _Kinetic;
  MeasurableNumber _Potential;
  MeasurableSum _TotalEnergy;

  std::vector<MeasurableObject*> _Meas_Ptr;


  std::vector<ExtraMeasurables*> _ExtraMeas_Ptr;


  //
  // For a given HamiltonianTerm, it will search to find it in the map
  // If it does not find it, it will try to create it and insert it in the map.
  // In either case it will return an iterator to the structure holding
  // the operator together with its buffer. Note that the coefficient
  // of the term is ignored.
  //


  multimap_type::iterator insert_term(const KeyType &key,const std::vector<IndexedProductElement> &term_product) {

    equal_range_type equal_range=_multimap.equal_range(key);

    // Try to look up the item
    multimap_type::iterator it=equal_range.first;
    while(it!=equal_range.second && it->second.product != term_product )
      ++it;

    // If it is not found create it.
    if(it==equal_range.second) {
      it=_multimap.insert(std::pair<KeyType,TermBuffer>( key, TermBuffer(term_product) ) );
    }

    return it;

  }

public:
  Measurable(OperatorStringType &OS) : BrokenLineTracer(OS.GetListBrokenLines(),OS), _Kinetic("Non-diagonal energy",&buffer_kinetic), _Potential("Diagonal energy",&buffer_potential), _TotalEnergy("Total energy") {
    reset();
    _buffers.push_back(&buffer_BoltzmannWeight);
    _buffers.push_back(&buffer_kinetic);
    _buffers.push_back(&buffer_potential);
    _TotalEnergy.push_back(&buffer_kinetic,1.0);
    _TotalEnergy.push_back(&buffer_potential,1.0);
  }

  ~Measurable() {
    for(std::vector<MeasurableObject*>::iterator it=_Meas_Ptr.begin(); it!=_Meas_Ptr.end(); ++it)
      delete *it;
  }

  // Looks up a term and returns an address to its buffer
  _float_accumulator *new_buffer_pointer(const HamiltonianTerm &term) {
    _float_accumulator *ptr=&insert_term(BrokenLines::map(term.product()),term.product())->second.buffer;
    _buffers.push_back(ptr);
    return ptr;
  }

  // Inserts an arbitrary measurable function to the list of measurables
  void insert(MeasurableObject *_meas_ptr) {
    _Meas_Ptr.push_back(_meas_ptr);
  }

  void insert_extra(ExtraMeasurables *ptr) {
    _ExtraMeas_Ptr.push_back(ptr);
    _Meas_Ptr.push_back(ptr);
  }


  inline void measure(const OperatorStringType &OperatorString) {

    const _float_accumulator Weight=OperatorString.BoltzmannWeight();

    const KeyType &key=BrokenLineTracer();

    equal_range_type equal_range=_multimap.equal_range(key);
    for(multimap_type::iterator it=equal_range.first; it!=equal_range.second; ++it) {
      it->second.buffer += it->second.me()*Weight;
    }

    for(std::vector<ExtraMeasurables*>::iterator it=_ExtraMeas_Ptr.begin(); it!=_ExtraMeas_Ptr.end(); ++it)
      (*it)->measure(Weight,key);

    if(OperatorString.NBrokenLines()==0) {
      buffer_BoltzmannWeight+=Weight;
      _float_accumulator Kinetic=OperatorString.KineticEnergy();
      _float_accumulator Potential=OperatorString.DiagonalEnergy();
      buffer_kinetic+=Kinetic*Weight;
      buffer_potential+=Potential*Weight;
    }

  }


  inline void flush() {

#ifdef USEMPI

    //
    // When using MPI we need to sum all the buffers of all the nodes
    // *before* we run the evaluate functions.
    //
    std::vector<_float_accumulator*>::size_type buffer_size=_buffers.size();
    std::vector<_float_accumulator> send_buffers(buffer_size),receive_buffers(buffer_size);
    for(std::vector<_float_accumulator*>::size_type i=0; i<buffer_size; ++i)
      send_buffers[i]=*_buffers[i];

    MPI_Reduce(&send_buffers[0],&receive_buffers[0],buffer_size,MPI_LONG_DOUBLE,MPI_SUM,Master,MPI_COMM_WORLD);

    for(std::vector<_float_accumulator*>::size_type i=0; i<buffer_size; ++i)
      *_buffers[i]=receive_buffers[i];

#endif

    for(std::vector<MeasurableObject*>::size_type i=0; i<_Meas_Ptr.size(); ++i)
      _Meas_Ptr[i]->push(buffer_BoltzmannWeight);

    _Kinetic.push(buffer_BoltzmannWeight);
    _Potential.push(buffer_BoltzmannWeight);
    _TotalEnergy.push(buffer_BoltzmannWeight);


    reset();

  }

  // resets all the buffers
  inline void reset() {
    for(std::vector<_float_accumulator*>::size_type i=0; i<_buffers.size(); ++i)
      *_buffers[i]=0;

    for(std::vector<ExtraMeasurables*>::iterator it=_ExtraMeas_Ptr.begin(); it!=_ExtraMeas_Ptr.end(); ++it)
      (*it)->reset();


    buffer_BoltzmannWeight=0;
    buffer_kinetic=0;
    buffer_potential=0;
  }

  std::vector<MeasurableObject*>::size_type size() const {
    return _Meas_Ptr.size();
  }
  const MeasurableObject &TotalEnergy() const {
    return _TotalEnergy;
  }
  const MeasurableObject &Potential() const {
    return _Potential;
  }
  const MeasurableObject &Kinetic() const {
    return _Kinetic;
  }
  const MeasurableObject &Quantity(std::vector<MeasurableObject*>::size_type i) const {
    return *_Meas_Ptr[i];
  }

};

void InsertOperator(const std::string &tag,const Hamiltonian &Operator,Measurable &MeasuredOperators) {
  MeasurableSum *_meas_ptr=new MeasurableSum(tag);
  for(Hamiltonian::const_iterator term_ptr=Operator.begin(); term_ptr!=Operator.end(); ++term_ptr) {
    _meas_ptr->push_back( MeasuredOperators.new_buffer_pointer(*term_ptr), term_ptr->coefficient() );
  }
  MeasuredOperators.insert(_meas_ptr);
}

// Measure an individual term
void InsertTerm(const std::string &tag,const HamiltonianTerm &term,Measurable &MeasuredOperators) {
  MeasurableNumber *_meas_ptr=new MeasurableNumber(tag, MeasuredOperators.new_buffer_pointer(term), term.coefficient());
  MeasuredOperators.insert(_meas_ptr);
}



}


#endif
