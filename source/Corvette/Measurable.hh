#ifndef __MEASURABLE__
#define __MEASURABLE__

#include <map>
#include <vector>
#include <set>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>

#include "HamiltonianTerm.hh"
#include "Accumulator.hh"
#include "OperatorString.hh"


#include "mkl_lapack.h"
typedef double _doublereal;
typedef int _integer;

/*

This is a wrapper for a diagonalization
routine dsyev from lapack.

SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

Input:
SIZE: the dimension of the matrix
A: the symmetric matrix it self in row major format

OUTPUT:
val: eigen values

Note: the eigen vectors will be stored
in the input matrix

*/

//
// N: the linear size of the matrix
// jobz: 'N': compute only eigenvalues, 'V': also compute eigenvectors
// sym_a: the symmetric matrix. If jobz='V' it is overwritten by the eigenvectors (in row major format, like fortran)
// W: it stores the eigenvalues. It should have a size of N.
//

void eig(_integer N, char jobz, _doublereal sym_a[], _doublereal W[]) {

  char uplo = 'L';               // 'U': upper triangle of A is stored, 'L': lower triangle.
  _integer order = N;
  _integer lda = N;              // Number of rows

  _integer lwork = -1;
  _integer liwork = -1;          // Note: for jobz='N', it has to be at least 3+5*N
  _integer info;                 // Exit status: should be 0 if no errors occur


  // The first run is only to determine the optimal value for lwork
  _doublereal lwork_query;
  _integer liwork_query;
  dsyevd_(&jobz, &uplo, &order, sym_a, &lda, W, &lwork_query, &lwork, &liwork_query, &liwork, &info);

  lwork = lwork_query;
  liwork = liwork_query;

  _doublereal* work = new _doublereal[lwork];
  _integer* iwork = new _integer[liwork];
  dsyevd_(&jobz, &uplo, &order, sym_a, &lda, W, work, &lwork, iwork, &liwork, &info);

  delete [] iwork;
  delete [] work;

};


#ifdef USEMPI
// *******************************************
// * Declaration of global variables for MPI *
// *******************************************
#include <mpi.h>
#define Master 0
int NumProcessors, Rank, NameLength;
char ProcessorName[MPI_MAX_PROCESSOR_NAME];

#endif


namespace SGF {
void InitializeEnvironment(int NumArg, char* Arg[]) {
#ifdef USEMPI
  // ***********************************
  // * Initialization of MPI functions *
  // ***********************************

  MPI_Init(&NumArg, &Arg);
  MPI_Comm_size(MPI_COMM_WORLD, &NumProcessors);
  MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
  MPI_Get_processor_name(ProcessorName, &NameLength);

  // Silence the other nodes. Only the root node can print
  std::cout.rdbuf( Rank == Master ? std::cout.rdbuf() : 0 );
  std::cerr.rdbuf( Rank == Master ? std::cout.rdbuf() : 0 );

#endif
}

void FinalizeEnvironment() {
#ifdef USEMPI
  MPI_Finalize();
#endif

}

}


namespace SGF {


typedef BinnedAccumulator<_float_accumulator> BinnedAccumulatorME;

template<typename T>
std::ostream& yaml_print_inline(std::ostream& o, const BinnedAccumulator<T>& b) {
  return o << "[ " << std::fixed << std::setprecision(12) << std::setw(15) << std::right << b.average() << ", " << std::setw(20) << b.sigma() << " ]";
}

template<typename T>
std::ostream& yaml_print_inline(std::ostream& o, const std::string& tag, const BinnedAccumulator<T>& b) {
  o << std::setw(25) << std::left << tag + ": ";
  yaml_print_inline(o, b);
  return o;
}


template<typename T>
std::ostream& yaml_print_inline(std::ostream& o, unsigned int depth, const std::string& tag, const BinnedAccumulator<T>& b) {
  std::string indent(2 * depth, ' ');
  o << indent;
  yaml_print_inline(o, tag, b);
  return o;
}

// How to print a sequence of bins
template<typename T>
std::ostream& yaml_print(std::ostream& o, unsigned int depth, const std::vector<T>& v) {

  std::string indent(2 * depth, ' ');

  for (typename std::vector<T>::const_iterator it = v.begin(); it != v.end(); ++it) {
    o << indent << "- ";
    yaml_print_inline(o, *it);
    o << std::endl;
  }

  return o;

}


template<typename T>
std::ostream& yaml_print(std::ostream& o, unsigned int depth, const std::string& tag, const std::vector<T>& v) {

  o << std::endl;

  std::string indent(2 * depth, ' ');
  o << indent << tag << ":";
  o << std::endl;
  o << std::endl;
  yaml_print(o, depth + 1, v);
  o << std::endl;

  return o;

}



class MeasurableObject {
protected:
  std::string _tag;
public:
  MeasurableObject(const std::string& tag) : _tag(tag) {}
  virtual ~MeasurableObject() {}
  virtual void push(const _float_accumulator& Weight) = 0 ;
  virtual std::ostream& print(std::ostream& o, unsigned int depth) const = 0;
};


class MeasurableScalar : public MeasurableObject {
protected:
  BinnedAccumulatorME _bin;
  virtual _float_accumulator value() const = 0;
public:
  MeasurableScalar(const std::string& tag) : MeasurableObject(tag) {}
  void push(const _float_accumulator& Weight) {
#ifdef USEMPI
    if (Rank == Master)
#endif
      _bin.push(value() / Weight);
  }

  std::ostream& print(std::ostream& o, unsigned int depth = 1) const {
    yaml_print_inline(o, depth, _tag, _bin);
    return o;
  }

};

class MeasurableNumber : public MeasurableScalar {
  _float_accumulator* data;
  MatrixElement coefficient;
public:
  MeasurableNumber(const std::string& tag, _float_accumulator* a, MatrixElement c = 1.0) : MeasurableScalar(tag), data(a), coefficient(c) {}
  inline _float_accumulator value() const {
    return *data * coefficient;
  }
};


class MeasurableSum : public MeasurableScalar {
  typedef std::pair<_float_accumulator*, MatrixElement> Pair;
  typedef std::vector<Pair> vector_pair_t;
  vector_pair_t data;

public:
  MeasurableSum(const std::string& tag) : MeasurableScalar(tag) {}
  inline _float_accumulator value() const {
    _float_accumulator result = 0;
    for (vector_pair_t::const_iterator it = data.begin(); it != data.end(); ++it)
      result += *it->first * it->second;
    return result;
  }
  void push_back(_float_accumulator* a, MatrixElement m) {
    data.push_back(Pair(a, m));
  }

};

class MeasurableSequence : public MeasurableObject {
protected:
  std::vector<BinnedAccumulatorME> Bins;
public:
  MeasurableSequence(const std::string& s, unsigned long nn) : MeasurableObject(s), Bins(nn) {}

  std::ostream& print(std::ostream& o, unsigned int depth) const  {
    yaml_print(o, depth, _tag, Bins);
    return o;
  }

  virtual void push(const _float_accumulator& Weight) {};

  void push(const std::vector<MatrixElement>& data) {
    assert(data.size() == Bins.size());
    for (unsigned long i = 0; i < data.size(); ++i)
      Bins[i].push(data[i]);
  }

};


class MeasurableVector : public MeasurableSequence {
protected:
  std::vector<_float_accumulator*> data_ptr;
public:
  MeasurableVector(const std::string& tag, unsigned long nn) : MeasurableSequence(tag, nn), data_ptr(nn) {}
  virtual void push(const _float_accumulator& Weight) {
    for (unsigned long i = 0; i < data_ptr.size(); ++i)
      Bins[i].push(*data_ptr[i] / Weight);
  }
  _float_accumulator*& set(unsigned long i) {
    return data_ptr[i];
  }


};


class MeasurableEigenSystem : public MeasurableVector {
protected:

  MeasurableSequence Elements;
  MeasurableSequence EigenVectors;
  MeasurableSequence EigenValues;

  unsigned long N;
public:
  MeasurableEigenSystem(const std::string& s, unsigned long nn) :
    MeasurableVector(s, nn* nn),
    Elements("Elements", nn* nn),
    EigenVectors("Eigenvectors", nn* nn),
    EigenValues("Eigenvalues", nn),
    N(nn)  {}

  unsigned long index(unsigned long i, unsigned long j) const {
    return i + N * j;
  }

  void push(const _float_accumulator& Weight) {

    std::vector<_doublereal> data_copy(N * N);
    std::vector<_doublereal> vals(N);

    for (unsigned long i = 0; i < data_copy.size(); ++i) {
      data_copy[i] = *data_ptr[i] / Weight;
    }

    // Symmetrize the matrix
    for (unsigned long i = 0; i < N; ++i)
      for (unsigned long j = 0; j < N; ++j) {
        _doublereal temp = data_copy[index(i, j)] + data_copy[index(j, i)];
        data_copy[index(i, j)] = data_copy[index(j, i)] = 0.5 * temp;
      }

    Elements.push(data_copy);
    eig(N, 'V', &data_copy[0], &vals[0]);
    EigenValues.push(vals);
    EigenVectors.push(data_copy);
  }



  std::ostream& print(std::ostream& o,unsigned int depth) const  {

    std::string indent(2*depth,' ');
    o << std::endl;
    o << indent << _tag << ":" << std::endl;
    o << std::endl;
    o << std::endl;
    Elements.print(o, 1+depth);
    EigenValues.print(o, 1+depth);
    EigenVectors.print(o, 1+depth);
    return o;
  }


};



class MeasurableEigenvalues : public MeasurableVector {
  std::vector<_float_accumulator*> data_ptr;
  unsigned long N;
public:
  MeasurableEigenvalues(const std::string& tag, unsigned long nn) : MeasurableVector(tag, nn), data_ptr(nn* nn), N(nn) {}

  void push(const _float_accumulator& Weight) {

    std::vector<_doublereal> data_copy(N * N);
    std::vector<_doublereal> vals(N);

    for (unsigned long i = 0; i < data_copy.size(); ++i) {
      data_copy[i] = *data_ptr[i] / Weight;
    }

    for (unsigned long i = 0; i < N; ++i)
      for (unsigned long j = 0; j < N; ++j) {
        _doublereal temp = data_copy[index(i, j)] + data_copy[index(j, i)];
        data_copy[index(i, j)] = data_copy[index(j, i)] = 0.5 * temp;
      }


    eig(N, 'N', &data_copy[0], &vals[0]);

    for (unsigned long i = 0; i < N; ++i)
      Bins[i].push(vals[i]);
  }

  unsigned long index(unsigned long i, unsigned long j) const {
    return i + N * j;
  }

  _float_accumulator*& set(unsigned long i) {
    return data_ptr[i];
  }

};


class ExtraMeasurables : public MeasurableObject {
public:
  ExtraMeasurables(std::string tag) : MeasurableObject(tag) {}
  typedef std::vector<std::pair<Boson*, int> > KeyType;
  virtual void reset() = 0;
  virtual void measure(_float_accumulator Weight, const KeyType& key) = 0;
};

class MeasurableDensityMatrix : public ExtraMeasurables {

  std::vector<BinnedAccumulator<double> > BinsElements;

  const Boson* Psi0;
  typedef std::vector<Boson>::size_type size_type;
  size_type N;
  std::string tag;
  std::vector<double> data;
  double _number;
  BinnedAccumulator<long double> _numberbin;

  _float_accumulator me_ca(size_type i, size_type j) const {
    return sqrt( (Psi0[i].nR() + 1) * Psi0[j].nR() );
  }

  size_type index(size_type i, size_type j) const {
    return (i < j) ? (2 * N - i + 1) * i / 2 + j - i : (2 * N - j + 1) * j / 2 + i - j;
  }

public:
  MeasurableDensityMatrix(const std::string& tag, std::vector<Boson>& Psi) : ExtraMeasurables(tag), _numberbin() {
    Psi0 = &Psi[0];
    N = Psi.size();
    size_type ndata = N * (N + 1) / 2;
    data.resize(ndata, 0);
    BinsElements.resize(ndata);
    _number = 0;
  }

  void measure(_float_accumulator Weight, const KeyType& key) {

    if (key.size() == 0) {
      double count = 0;
      for (unsigned long i = 0; i < N; ++i) {
        unsigned long ind = index(i, i);
        data[ind] += Psi0[i].nR() * Weight;
        count += Psi0[i].nR();
      }
      _number += count * Weight;
    } else if (key.size() == 2) {

      if (key[0].second == 1 && key[1].second == -1) {
        size_type j = key[0].first - Psi0;
        size_type i = key[1].first - Psi0;
        data[index(i, j)] += me_ca(i, j) * Weight / 2.0;
      }

      if (key[1].second == 1 && key[0].second == -1) {
        size_type j = key[1].first - Psi0;
        size_type i = key[0].first - Psi0;
        data[index(i, j)] += me_ca(i, j) * Weight / 2.0;
      }

    }

  }

  void reset() {
    std::fill(data.begin(), data.end(), 0);
    _number = 0;
  }


  void reduce() {

#ifdef USEMPI

    std::vector<double> send_data(data);
    MPI_Reduce(&send_data[0], &data[0], data.size(), MPI_LONG_DOUBLE, MPI_SUM, Master, MPI_COMM_WORLD);

#endif

  }

  void push(const _float_accumulator& Weight) {

    reduce();

    for (unsigned long i = 0; i < data.size(); ++i)
      BinsElements[i].push(data[i] / Weight);

    _numberbin.push(_number / Weight);

  }



  std::ostream& print(std::ostream& o, unsigned int depth = 1) const  {

    o << std::endl;

    yaml_print_inline(o, depth, "Number", _numberbin);
    yaml_print(o, depth, _tag, BinsElements);

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

// Maximum number of broken lines is only used for the BrokenLineHistogram.
#define MAXNUMBROKENLINES 100


class Measurable {

  typedef  std::vector<unsigned long> BrokenHistogramType;
  BrokenHistogramType BrokenHistorgram;

  struct TermBuffer {
    _float_accumulator buffer;
    std::vector<IndexedProductElement> product;
    TermBuffer(const std::vector<IndexedProductElement>& _p) : buffer(0), product(_p) {}
    TermBuffer(const TermBuffer& o) : buffer(o.buffer), product(o.product) {}
    MatrixElement me() const {
      return MultiplyMe<MatrixElement, IndexedProductElement, RIGHT>(product);
    }
  };

  typedef BosonDeltaMapType KeyType;
  typedef std::multimap<KeyType, TermBuffer> multimap_type;
  typedef std::pair<multimap_type::iterator, multimap_type::iterator> equal_range_type;

  //  const OperatorStringType &OperatorString; // Holds a reference for the operator string being measured

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


  multimap_type::iterator insert_term(const KeyType& key, const std::vector<IndexedProductElement>& term_product) {

    equal_range_type equal_range = _multimap.equal_range(key);

    // Try to look up the item
    multimap_type::iterator it = equal_range.first;
    while (it != equal_range.second && it->second.product != term_product )
      ++it;

    // If it is not found create it.
    if (it == equal_range.second) {
      it = _multimap.insert(std::pair<KeyType, TermBuffer>( key, TermBuffer(term_product) ) );
    }

    return it;

  }

public:
  Measurable() : _Kinetic("Non-diagonal energy", &buffer_kinetic), _Potential("Diagonal energy", &buffer_potential), _TotalEnergy("Total energy") {
    reset();
    _buffers.push_back(&buffer_BoltzmannWeight);
    _buffers.push_back(&buffer_kinetic);
    _buffers.push_back(&buffer_potential);
    _TotalEnergy.push_back(&buffer_kinetic, 1.0);
    _TotalEnergy.push_back(&buffer_potential, 1.0);
    BrokenHistorgram.resize(MAXNUMBROKENLINES);
  }

  ~Measurable() {
    for (std::vector<MeasurableObject*>::iterator it = _Meas_Ptr.begin(); it != _Meas_Ptr.end(); ++it)
      delete *it;
  }

  // Looks up a term and returns an address to its buffer
  _float_accumulator* new_buffer_pointer(const HamiltonianTerm& term) {
    _float_accumulator* ptr = &insert_term(GetTermHash(term.product()), term.product())->second.buffer;
    _buffers.push_back(ptr);
    return ptr;
  }

  // Inserts an arbitrary measurable function to the list of measurables
  void insert(MeasurableObject* _meas_ptr) {
    _Meas_Ptr.push_back(_meas_ptr);
  }

  void insert(ExtraMeasurables* ptr) {
    _ExtraMeas_Ptr.push_back(ptr);
    _Meas_Ptr.push_back(ptr);
  }

  inline void measure() {}

  template<class OperatorStringType>
  inline void measure(const OperatorStringType& OperatorString, const KeyType& key) {

    const _float_accumulator Weight = OperatorString.BoltzmannWeight();

    equal_range_type equal_range = _multimap.equal_range(key);
    for (multimap_type::iterator it = equal_range.first; it != equal_range.second; ++it) {
      it->second.buffer += it->second.me() * Weight;
    }

    for (std::vector<ExtraMeasurables*>::iterator it = _ExtraMeas_Ptr.begin(); it != _ExtraMeas_Ptr.end(); ++it)
      (*it)->measure(Weight, key);

    if (OperatorString.NBrokenLines() == 0) {
      buffer_BoltzmannWeight += Weight;
      _float_accumulator Kinetic = OperatorString.KineticEnergy();
      _float_accumulator Potential = OperatorString.DiagonalEnergy();
      buffer_kinetic += Kinetic * Weight;
      buffer_potential += Potential * Weight;
    }

    BrokenHistorgram[OperatorString.NBrokenLines()] += 1;


  }


  inline void flush() {

#ifdef USEMPI

    //
    // When using MPI we need to sum all the buffers of all the nodes
    // *before* we run the evaluate functions.
    //
    std::vector<_float_accumulator*>::size_type buffer_size = _buffers.size();
    std::vector<_float_accumulator> send_buffers(buffer_size), receive_buffers(buffer_size);
    for (std::vector<_float_accumulator*>::size_type i = 0; i < buffer_size; ++i)
      send_buffers[i] = *_buffers[i];

    MPI_Reduce(&send_buffers[0], &receive_buffers[0], buffer_size, MPI_LONG_DOUBLE, MPI_SUM, Master, MPI_COMM_WORLD);

    for (std::vector<_float_accumulator*>::size_type i = 0; i < buffer_size; ++i)
      *_buffers[i] = receive_buffers[i];

#endif

    for (std::vector<MeasurableObject*>::size_type i = 0; i < _Meas_Ptr.size(); ++i)
      _Meas_Ptr[i]->push(buffer_BoltzmannWeight);

    _Kinetic.push(buffer_BoltzmannWeight);
    _Potential.push(buffer_BoltzmannWeight);
    _TotalEnergy.push(buffer_BoltzmannWeight);


    reset();

  }

  // resets all the buffers
  inline void reset() {
    for (std::vector<_float_accumulator*>::size_type i = 0; i < _buffers.size(); ++i)
      *_buffers[i] = 0;

    for (std::vector<ExtraMeasurables*>::iterator it = _ExtraMeas_Ptr.begin(); it != _ExtraMeas_Ptr.end(); ++it)
      (*it)->reset();


    buffer_BoltzmannWeight = 0;
    buffer_kinetic = 0;
    buffer_potential = 0;
  }

  std::ostream&  print(std::ostream& o, unsigned int depth = 1) {

    std::string indent(2 * depth, ' ');
    o << std::endl;
    o << indent << "Broken worldlines";
    o << std::endl;
    o << std::endl;

    double Normalization = 0;
    for (unsigned int i = 0; i < BrokenHistorgram.size(); ++i)
      Normalization += BrokenHistorgram[i];

    for (unsigned int i = 0; i < BrokenHistorgram.size(); ++i)
      if (BrokenHistorgram[i] != 0)
        o << indent << "  " << "- [ " << std::right << std::fixed << std::setw(3) << i << "," << std::setw(12) << BrokenHistorgram[i] << std::setprecision(9) << std::left << "," << std::setw(13) << std::right << BrokenHistorgram[i] / Normalization << " ]\n";

    o << std::endl;


    _TotalEnergy.print(o, depth);
    o << std::endl;
    _Potential.print(o, depth);
    o << std::endl;
    _Kinetic.print(o, depth);
    o << std::endl;
    o << std::endl;


    for (unsigned int i = 0; i < _Meas_Ptr.size(); ++i) {
      _Meas_Ptr[i]->print(o, depth);
      o << std::endl;
    }
    return o;

  }

  void InsertOperator(const std::string& tag, const Hamiltonian& Operator) {
    MeasurableSum* _meas_ptr = new MeasurableSum(tag);
    for (Hamiltonian::const_iterator term_ptr = Operator.begin(); term_ptr != Operator.end(); ++term_ptr) {
      _meas_ptr->push_back( new_buffer_pointer(*term_ptr), term_ptr->coefficient() );
    }
    insert(_meas_ptr);
  }

  void InsertTerm(const std::string& tag, const HamiltonianTerm& term) {
    MeasurableNumber* _meas_ptr = new MeasurableNumber(tag, new_buffer_pointer(term), term.coefficient());
    insert(_meas_ptr);
  }

  void InsertFunnyDensityMatrix(const std::string& tag, std::vector<Boson>& Psi) {
    insert(new MeasurableDensityMatrix(tag, Psi) );
  }


  void InsertVector(const std::string& tag, const SGF::Hamiltonian& T, MeasurableVector* MVec) {
    for (SGF::Hamiltonian::size_type i = 0; i != T.size(); ++i) {
      _float_accumulator* term_ptr = new_buffer_pointer(T[i]);
      MVec->set(i) = term_ptr;
    }
    insert(MVec);
  }

  void InsertOperatorTerms(const std::string& tag, const SGF::Hamiltonian& T) {
    InsertVector(tag, T, new MeasurableVector(tag, T.size()));
  }

  void InsertDensityMatrixEigenvalues(const std::string& tag, std::vector<Boson>& Psi) {
    InsertVector(tag, SGFBase::GenerateDensityMatrix(Psi), new MeasurableEigenvalues(tag, Psi.size()));
  }


  void InsertDensityMatrixEigenSystem(const std::string& tag, std::vector<SGF::Boson>& Psi) {
    InsertVector(tag, SGFBase::GenerateDensityMatrix(Psi), new MeasurableEigenSystem(tag, Psi.size()));
  }

  void InsertDensityMatrix(const std::string& tag, std::vector<SGF::Boson>& Psi) {
    InsertOperatorTerms(tag, SGFBase::GenerateDensityMatrix(Psi));
  }

  void InsertDensityDensityMatrix(const std::string& tag, std::vector<SGF::Boson>& Psi) {
    InsertOperatorTerms(tag, SGFBase::GenerateDensityDensityMatrix(Psi));
  }

  void InsertLocalDensity(const std::string& tag, std::vector<SGF::Boson>& Psi) {
    InsertOperatorTerms(tag, SGFBase::GenerateNumberOperator(Psi));
  }

};




}


#endif
