#ifndef __MEASURABLE__
#define __MEASURABLE__

#include <map>
#include <vector>
#include <set>
#include <string>
#include <cmath>

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

// How to print the Bins
template<class T> inline std::ostream& operator<<(std::ostream& output, const BinnedAccumulator<T> &o) { 
   return output<<o.average()<<" +/- "<<o.sigma(); 
} 


typedef BinnedAccumulator<_float_accumulator> BinnedAccumulatorME;

class MeasurableFunction {
   std::string _tag;
   BinnedAccumulatorME _bin;
public:
   MeasurableFunction(const std::string &tag) :_tag(tag) {}
   virtual inline _float_accumulator evaluate() const = 0;
   virtual ~MeasurableFunction() {}
   void push(const _float_accumulator &Weight) {
#ifdef USEMPI
      if(Rank==Master)
#endif
         _bin.push(evaluate()/Weight);
   }
   std::ostream &print(std::ostream &o) const {
      return o<<_tag<<": "<<_bin;
   }
};

inline std::ostream& operator<<(std::ostream& o,const MeasurableFunction &m) {
   return m.print(o);
}

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

   const OperatorStringType &OperatorString; // Holds a reference for the operator string being measured


   _float_accumulator buffer_BoltzmannWeight;
   _float_accumulator buffer_kinetic;
   _float_accumulator buffer_potential;

   BrokenLines BrokenLineTracer;  // It traces the list of broken lines
   multimap_type _multimap;       // It stores the buffers of all measurable terms

   // It stores the addresses of all _multimap's buffers together with the kinetic and potential buffers
   std::vector<_float_accumulator*> _buffers;

   MeasurableNumber _Kinetic;
   MeasurableNumber _Potential;
   MeasurableSum _TotalEnergy;

   std::vector<MeasurableFunction*> _Meas_Ptr;

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
   Measurable(OperatorStringType &OS) : OperatorString(OS), BrokenLineTracer(OS.GetListBrokenLines(),OS), _Kinetic("Non-diagonal energy",&buffer_kinetic), _Potential("Diagonal energy",&buffer_potential), _TotalEnergy("Total energy") {
      reset();
      _buffers.push_back(&buffer_BoltzmannWeight);
      _buffers.push_back(&buffer_kinetic);
      _buffers.push_back(&buffer_potential);
      _TotalEnergy.push_back(&buffer_kinetic,1.0);
      _TotalEnergy.push_back(&buffer_potential,1.0);
   }

   ~Measurable() {
      for(std::vector<MeasurableFunction*>::iterator it=_Meas_Ptr.begin(); it!=_Meas_Ptr.end(); ++it)
         delete *it;
   }

   // Looks up a term and returns an address to its buffer
   _float_accumulator *new_buffer_pointer(const HamiltonianTerm &term) {
      _float_accumulator *ptr=&insert_term(BrokenLines::map(term.product()),term.product())->second.buffer;
      _buffers.push_back(ptr);
      return ptr;
   }

   // Inserts an arbitrary measurable function to the list of measurables
   void insert(MeasurableFunction *_meas_ptr) {
      _Meas_Ptr.push_back(_meas_ptr);
   }

   // Measure an operator (sum of terms)
   void insert(const std::string &tag,const Hamiltonian &Operator) {
      MeasurableSum *_meas_ptr=new MeasurableSum(tag);
      for(Hamiltonian::const_iterator term_ptr=Operator.begin(); term_ptr!=Operator.end(); ++term_ptr) {
         _meas_ptr->push_back( new_buffer_pointer(*term_ptr), term_ptr->coefficient() );
      }
      _Meas_Ptr.push_back(_meas_ptr);
   }

   // Measure an individual term
   void insert(const std::string &tag,const HamiltonianTerm &term) {
      MeasurableNumber *_meas_ptr=new MeasurableNumber(tag, new_buffer_pointer(term), term.coefficient());
      _Meas_Ptr.push_back(_meas_ptr);
   }

   inline void measure() {

      const _float_accumulator Weight=OperatorString.BoltzmannWeight();

      equal_range_type equal_range=_multimap.equal_range(BrokenLineTracer());
      for(multimap_type::iterator it=equal_range.first; it!=equal_range.second; ++it) {
         it->second.buffer += it->second.me()*Weight;
      }

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

      for(std::vector<MeasurableFunction*>::size_type i=0; i<_Meas_Ptr.size(); ++i)
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
      buffer_BoltzmannWeight=0;
      buffer_kinetic=0;
      buffer_potential=0;
   }

   std::vector<MeasurableFunction*>::size_type size() const {
      return _Meas_Ptr.size();
   }
   const MeasurableFunction &TotalEnergy() const {
      return _TotalEnergy;
   }
   const MeasurableFunction &Potential() const {
      return _Potential;
   }
   const MeasurableFunction &Kinetic() const {
      return _Kinetic;
   }
   const MeasurableFunction &Quantity(std::vector<MeasurableFunction*>::size_type i) const {
      return *_Meas_Ptr[i];
   }

};

}


#endif
