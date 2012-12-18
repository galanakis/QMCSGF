#ifndef __UNORDEREDSET__
#define __UNORDEREDSET__


namespace SGF {

class UnorderedSet {
  typedef unsigned int size_type;

  size_type _capacity;
  size_type *data_begin, *data_end;
  size_type **map;

public:

  typedef size_type* iterator;

  inline iterator begin() const {return data_begin;}
  inline iterator end() const {return data_end;}

  inline size_type capacity() const {return _capacity;}
  inline size_type size() const {return data_end-data_begin;}

  UnorderedSet() : _capacity(0), data_begin(0), data_end(0), map(0) {}

  UnorderedSet(size_type size) {initialize(size);}

  void initialize(size_type size) {

    _capacity=size;
    map=new size_type*[_capacity];
    data_end=data_begin=new size_type[_capacity];

    for(size_type i=0; i<_capacity; ++i) {
      data_begin[i]=0;
      map[i]=0;
    }

  }

  ~UnorderedSet() {
    delete [] data_begin;
    delete [] map;
  }


  inline size_type element(size_type i) const {
    return *(data_begin+i);
  }

  inline void insert(size_type i) {

    if(!map[i]) {

      *data_end=i;
      map[i]=data_end;
      ++data_end;

    }
  }

  inline void erase(size_type i) {

    if(map[i]) {

      --data_end;
      *map[i]=*data_end;
      map[*data_end]=map[i];
      map[i]=0;

    }

  }


};

}

#endif