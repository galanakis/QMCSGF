/* ***************************************************************************** */
/* Copyright:      Francois Panneton and Pierre L'Ecuyer, University of Montreal */
/*                 Makoto Matsumoto, Hiroshima University                        */
/* Notice:         This code can be used freely for personal, academic,          */
/*                 or non-commercial purposes. For commercial purposes,          */
/*                 please contact P. L'Ecuyer at: lecuyer@iro.UMontreal.ca       */
/*                 A modified "maximally equidistributed" implementations        */
/*                 by Shin Harase, Hiroshima University.                         */
/* ***************************************************************************** */



class well44497 {

  static const unsigned int R = 1391;
  static const unsigned int W = 32;
  static const unsigned int DISCARD = 15;
  static const unsigned int MASKU = 0xffffffffU>>(W-DISCARD);
  static const unsigned int MASKL = ~MASKU;
  static const unsigned int M1 = 23;
  static const unsigned int M3 = 229;
  static const unsigned int M2 = 481;
  static const unsigned int BITMASK = 0x48000000;

  static unsigned int  MAT0POS(unsigned int t,unsigned int v) { return v^(v>>t); }
  static unsigned int  MAT0NEG(unsigned int t,unsigned int v) { return v^(v<<t); }

  static unsigned int MAT3(unsigned int t,unsigned int v) { return v<<t; }

  static unsigned int MAT5(unsigned int r,unsigned int a,unsigned int ds,unsigned int dt,unsigned int v) {
    return ((v & dt)?((((v<<r)^(v>>(W-r)))&ds)^a):(((v<<r)^(v>>(W-r)))&ds));
  }

  static unsigned int mod(int i) {
    if (i<0)
      return i+R;
    else if (i<R)
      return i;
    else 
      return i-R;
  }

  unsigned int STATE[R];
  int state_i;

  unsigned int &V(int i) {
    return STATE[mod(i)];
  }


public: 
  well44497() {}

  static const unsigned int size = R;

  void init(unsigned int *init ){
    state_i=0;
    for(int j=0;j<R;j++)
      STATE[j]=init[j];
  }


  unsigned int operator()() {

    unsigned int z0 = (V(state_i-1) & MASKL) | (V(state_i-2) & MASKU);
    unsigned int z1 = MAT0NEG(24,V(state_i)) ^ MAT0POS(30,V(state_i+M1));
    unsigned int z2 = MAT0NEG(10,V(state_i+M2)) ^ MAT3(26,V(state_i+M3));
    V(state_i) = z1 ^ z2;
    V(state_i-1) = z0 ^ MAT0POS(20,z1) ^ MAT5(9,0xb729fcecU,0xfbffffffU,0x00020000U,z2) ^ z1 ^ z2;
    state_i=mod(state_i-1);
    return V(state_i) ^ (V(state_i+M2+1) & BITMASK);

  }

};

