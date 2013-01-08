#include "Measurable.hh"
#include <sstream>
#include <iomanip>

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

void eig(_integer N,char jobz,_doublereal sym_a[],_doublereal W[]) {

  char uplo='L';                 // 'U': upper triangle of A is stored, 'L': lower triangle.
  _integer order=N;
  _integer lda=N;                // Number of rows

  _integer lwork=-1;
  _integer liwork=-1;            // Note: for jobz='N', it has to be at least 3+5*N
  _integer info;                 // Exit status: should be 0 if no errors occur


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



class MeasurableMultiFunction : public MeasurableObject {
protected:
  std::vector<BinnedAccumulatorME> Bins;
public:
  MeasurableMultiFunction(const std::string &s,unsigned int nn) : Bins(nn), MeasurableObject(s) {}

  std::ostream& print(std::ostream& o) const  {
    int N=Bins.size();
    o<<std::endl;
    o<<"    "<<_tag<<" size= "<<N<<std::endl;
    o<<std::endl;

    for(unsigned int i=0; i<N; ++i) {
      o<<"     "<<std::setw(6)<<std::left<<i<<std::setprecision(9)<<Bins[i]<<std::endl;
    }

    o<<std::endl;
    return o;
  }

};

class MeasurableEigenSystem : public MeasurableObject {
protected:
  std::vector<BinnedAccumulatorME> BinsElements;
  std::vector<BinnedAccumulatorME> BinsEigenVectors;
  std::vector<BinnedAccumulatorME> BinsEigenValues;
  std::vector<_float_accumulator*> data_ptr;
  unsigned long N;
public:
  MeasurableEigenSystem(const std::string &s,unsigned int nn) : MeasurableObject(s), N(nn), BinsElements(nn*nn), BinsEigenVectors(nn*nn), BinsEigenValues(nn), data_ptr(nn*nn) {}

  unsigned int index(unsigned int i,unsigned int j) const {
    return i+N*j;
  }

  _float_accumulator* &set(unsigned int i,unsigned int j) {
    return data_ptr[index(i,j)];
  }

  void push(const _float_accumulator &Weight) {

    std::vector<_doublereal> data_copy(N*N);
    std::vector<_doublereal> vals(N);

    for(unsigned int i=0; i<data_copy.size(); ++i) {
      data_copy[i]=*data_ptr[i]/Weight;
    }

    // Symmetrize the matrix
    for(unsigned int i=0; i<N; ++i)
      for(unsigned int j=0; j<N; ++j) {
        _doublereal temp=data_copy[index(i,j)]+data_copy[index(j,i)];
        data_copy[index(i,j)]=data_copy[index(j,i)]=0.5*temp;
      }

    for(unsigned int i=0; i<N*N; ++i)
      BinsElements[i].push(data_copy[i]);

    eig(N,'V',&data_copy[0],&vals[0]);

    for(unsigned int i=0; i<N; ++i)
      BinsEigenValues[i].push(vals[i]);

    for(unsigned int i=0; i<N*N; ++i)
      BinsEigenVectors[i].push(data_copy[i]);

  }



  std::ostream& print(std::ostream& o) const  {

    o<<std::endl;
    o<<"    "<<_tag<<" (LDA= "<<N<<" )"<<std::endl;
    o<<std::endl<<std::endl;

    o<<"    Elements\n\n";
    for(unsigned int i=0; i<N*N; ++i) {
      o<<"     "<<std::setw(6)<<std::left<<i<<std::setprecision(9)<<BinsElements[i]<<std::endl;
    }
    o<<std::endl;

    o<<"    Eigenvalues\n\n";
    for(unsigned int i=0; i<N; ++i) {
      o<<"     "<<std::setw(6)<<std::left<<i<<std::setprecision(9)<<BinsEigenValues[i]<<std::endl;
    }
    o<<std::endl;

    o<<"    Eigenvectors\n\n";
    for(unsigned int i=0; i<N*N; ++i) {
      o<<"     "<<std::setw(6)<<std::left<<i<<std::setprecision(9)<<BinsEigenVectors[i]<<std::endl;
    }
    o<<std::endl;

    return o;
  }


};


class MeasurableVector : public MeasurableMultiFunction {
  std::vector<_float_accumulator*> data_ptr;
  unsigned int N;
public:
  MeasurableVector(const std::string &tag,unsigned int nn) : MeasurableMultiFunction(tag,nn), N(nn), data_ptr() {}
  void push(const _float_accumulator &Weight) {

    for(unsigned int i=0; i<N; ++i)
      Bins[i].push(*data_ptr[i]/Weight);
  }
  void push_back(_float_accumulator* o) {
    data_ptr.push_back(o);
  }

};

class MeasurableEigenvalues : public MeasurableMultiFunction {
  std::vector<_float_accumulator*> data_ptr;
  unsigned int N;

public:
  MeasurableEigenvalues(const std::string &tag,unsigned int nn) : MeasurableMultiFunction(tag,nn), N(nn), data_ptr(nn*nn) {}

  void push(const _float_accumulator &Weight) {

    std::vector<_doublereal> data_copy(N*N);
    std::vector<_doublereal> vals(N);

    for(unsigned int i=0; i<data_copy.size(); ++i) {
      data_copy[i]=*data_ptr[i]/Weight;
    }

    for(unsigned int i=0; i<N; ++i)
      for(unsigned int j=0; j<N; ++j) {
        _doublereal temp=data_copy[index(i,j)]+data_copy[index(j,i)];
        data_copy[index(i,j)]=data_copy[index(j,i)]=0.5*temp;
      }


    eig(N,'N',&data_copy[0],&vals[0]);

    for(unsigned int i=0; i<N; ++i)
      Bins[i].push(vals[i]);
  }

  unsigned int index(unsigned int i,unsigned int j) const {
    return i+N*j;
  }

  _float_accumulator* &set(unsigned int i,unsigned int j) {
    return data_ptr[index(i,j)];
  }


};




void InsertFunnyDensityMatrix(const std::string &tag,const std::vector<Boson> &Psi,Measurable &MeasuredOperators) {
  MeasuredOperators.insert_extra(new MeasurableDensityMatrix(tag,Psi) );
}

void InsertDensityMatrixEigenvalues(const std::string &tag,std::vector<Boson> &Psi,Measurable &MeasuredOperators) {

  // Insert the calculation of the condensate fraction

  MeasurableEigenvalues *MEig=new MeasurableEigenvalues(tag,Psi.size());

  // Inserts the density matrix
  for(unsigned int i=0; i<Psi.size(); ++i) {
    for(unsigned int j=0; j<Psi.size(); ++j) {
      if(i!=j) {
        const IndexedProductElement ci(C,&Psi[i]);
        const IndexedProductElement aj(A,&Psi[j]);
        HamiltonianTerm term(1.0,ci,aj);
        _float_accumulator *term_ptr=MeasuredOperators.new_buffer_pointer(term);
        MEig->set(i,j)=term_ptr;
      } else {
        const IndexedProductElement ni(C*A,&Psi[i]);
        HamiltonianTerm term(1.0,ni);
        _float_accumulator *term_ptr=MeasuredOperators.new_buffer_pointer(term);
        MEig->set(i,j)=term_ptr;
      }
    }
  }

  MeasuredOperators.insert(MEig);


}


void InsertDensityMatrixEigenSystem(const std::string &tag,std::vector<SGF::Boson> &Psi,Measurable &MeasuredOperators) {

  // Insert the calculation of the condensate fraction

  MeasurableEigenSystem *MEig=new MeasurableEigenSystem(tag,Psi.size());

  // Inserts the density matrix
  for(unsigned int i=0; i<Psi.size(); ++i) {
    for(unsigned int j=0; j<Psi.size(); ++j) {
      if(i!=j) {
        const IndexedProductElement ci(C,&Psi[i]);
        const IndexedProductElement aj(A,&Psi[j]);
        HamiltonianTerm term(1.0,ci,aj);
        _float_accumulator *term_ptr=MeasuredOperators.new_buffer_pointer(term);
        MEig->set(i,j)=term_ptr;
      } else {
        const IndexedProductElement ni(C*A,&Psi[i]);
        HamiltonianTerm term(1.0,ni);
        _float_accumulator *term_ptr=MeasuredOperators.new_buffer_pointer(term);
        MEig->set(i,j)=term_ptr;
      }
    }
  }

  MeasuredOperators.insert(MEig);


}


// We have to be carefull here
// We need to separately build the diagonal and the off diagonal terms
// Because the code assumes that all elements of IndexedProductElement
// correspond to different bosons.
void InsertDensityMatrix(const std::string &tag,std::vector<SGF::Boson> &Psi,Measurable &MeasuredOperators) {

  MeasurableVector *MVec=new MeasurableVector("Density Matrix",Psi.size()*Psi.size());
  // Inserts the density matrix
  for(unsigned int i=0; i<Psi.size(); ++i) {
    for(unsigned int j=0; j<Psi.size(); ++j) {
      if(i!=j) {
        const IndexedProductElement ci(C,&Psi[i]);
        const IndexedProductElement aj(A,&Psi[j]);
        HamiltonianTerm term(1.0,ci,aj);
        _float_accumulator *term_ptr=MeasuredOperators.new_buffer_pointer(term);
        MVec->push_back(term_ptr);
      } else {
        const IndexedProductElement ni(C*A,&Psi[i]);
        HamiltonianTerm term(1.0,ni);
        _float_accumulator *term_ptr=MeasuredOperators.new_buffer_pointer(term);
        MVec->push_back(term_ptr);
      }

    }
  }

  MeasuredOperators.insert(MVec);

}




void InsertLocalDensity(const std::string &tag,std::vector<SGF::Boson> &Psi,Measurable &MeasuredOperators) {
  MeasurableVector *MVec=new MeasurableVector(tag,Psi.size());
  // Inserts the density matrix
  for(unsigned int i=0; i<Psi.size(); ++i) {
    const IndexedProductElement ni(C*A,&Psi[i]);
    std::stringstream ss;
    ss<<"n"<<std::setw(6)<<i<<" ";
    HamiltonianTerm term(1.0,ni);
    _float_accumulator *term_ptr=MeasuredOperators.new_buffer_pointer(term);
    MVec->push_back(term_ptr);
  }
  MeasuredOperators.insert(MVec);

}



}




