#include "Measurable.hh"
#include <sstream>

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

template<class T>
inline void transpose(_integer N,T *sym_a) {

   T temp;

   for(_integer i=0; i<N; i++)
      for(_integer j=i+1; j<N; j++) {
         temp=sym_a[i*N+j];
         sym_a[i*N+j]=sym_a[j*N+i];
         sym_a[j*N+i]=temp;
      }

}

template<class A,class B>
inline void copy(_integer N,A *from,B *to) {
   for(_integer i=0; i<N; i++)
      to[i]=static_cast<B>(from[i]);

}


//
// N: the linear size of the matrix
// jobz: 'N': compute only eigenvalues, 'V': also compute eigenvectors
// sym_a: the symmetric matrix. If jobz='V' it is overwritten by the eigenvectors (in row major format, like fortran)
// W: it stores the eigenvalues. It should have a size of N.
//

void eig(_integer N,char jobz,_doublereal sym_a[],_doublereal W[]) {

   char uplo='L';						 // 'U': upper triangle of A is stored, 'L': lower triangle.
   _integer order=N;
   _integer lda=N;			       // Number of rows

   _integer lwork=-1;
   _integer liwork=-1; 			 // Note: for jobz='N', it has to be at least 3+5*N
   _integer info;			 			 // Exit status: should be 0 if no errors occur


   // The first run is only to determine the optimal value for lwork
   _doublereal lwork_query;
   _integer liwork_query;
   dsyevd_(&jobz,&uplo,&order,sym_a,&lda,W,&lwork_query,&lwork,&liwork_query,&liwork,&info);

   lwork=lwork_query;
   liwork=liwork_query;

   _doublereal *work=new _doublereal[lwork];
   _integer *iwork=new _integer[liwork];
   dsyevd_(&jobz,&uplo,&order,sym_a,&lda,W,work,&lwork,iwork,&liwork,&info);

   delete [] iwork;
   delete [] work;

};

namespace SGF {

class MeasurableEigenvalues : public MeasurableFunction {

   std::vector<_float_accumulator*> data_ptr;
   unsigned int N;

public:
   MeasurableEigenvalues(const std::string &tag,unsigned int nn) : MeasurableFunction(tag), N(nn), data_ptr(nn*nn) {}

   inline _float_accumulator evaluate() const {

      std::vector<_doublereal> data_copy(N*N);
      std::vector<_doublereal> vals(N);

      for(unsigned int i=0; i<data_copy.size(); ++i) {
         data_copy[i]=*data_ptr[i];
      }
      eig(N,'N',&data_copy[0],&vals[0]);

      return vals[N-1]-vals[N-2];
   }

   unsigned int index(unsigned int i,unsigned int j) const {
      return i+N*j;
   }

   _float_accumulator* &set(unsigned int i,unsigned int j) {
      return data_ptr[index(i,j)];
   }

};

void InsertCondensateFraction(std::vector<SGF::Boson> &Psi,Measurable &MeasuredOperators) {
   
   // Insert the calculation of the condensate fraction
   
	MeasurableEigenvalues *MEig=new MeasurableEigenvalues("Condensate Fraction",Psi.size());

   // Inserts the density matrix
   for(unsigned int i=0;i<Psi.size();++i) {
      for(unsigned int j=0;j<Psi.size();++j) {
         const IndexedProductElement ci(C,&Psi[i]);
         const IndexedProductElement aj(A,&Psi[j]);
         std::stringstream ss;
         ss<<"C["<<i<<"]A["<<j<<"]";
         HamiltonianTerm term(1.0,ci,aj);
         _float_accumulator *term_ptr=MeasuredOperators.new_buffer_pointer(term);
         MEig->set(i,j)=term_ptr;

      }
   }

   MeasuredOperators.insert(MEig);


}


}




