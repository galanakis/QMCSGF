/* ***************************************************************************** */
/* Copyright:      Francois Panneton and Pierre L'Ecuyer, University of Montreal */
/*                 Makoto Matsumoto, Hiroshima University                        */
/* Notice:         This code can be used freely for personal, academic,          */
/*                 or non-commercial purposes. For commercial purposes,          */
/*                 please contact P. L'Ecuyer at: lecuyer@iro.UMontreal.ca       */
/*                 A modified "maximally equidistributed" implementations        */
/*                 by Shin Harase, Hiroshima University.                         */
/* ***************************************************************************** */

#define _W_ 32
#define _R_ 1391
#define DISCARD 15
#define MASKU (0xffffffffU>>(_W_-DISCARD))
#define MASKL (~MASKU)

#define _M1_ 23
#define _M2_ 481
#define _M3_ 229

#define MAT0POS(t,v) (v^(v>>t))
#define MAT0NEG(t,v) (v^(v<<(-(t))))
#define MAT1(v) v
#define MAT2(a,v) ((v & 1U)?((v>>1)^a):(v>>1))
#define MAT3POS(t,v) (v>>t)
#define MAT3NEG(t,v) (v<<(-(t)))
#define MAT4POS(t,b,v) (v ^ ((v>>  t ) & b))
#define MAT4NEG(t,b,v) (v ^ ((v<<(-(t))) & b))
#define MAT5(r,a,ds,dt,v) ((v & dt)?((((v<<r)^(v>>(_W_-r)))&ds)^a):(((v<<r)^(v>>(_W_-r)))&ds))
#define MAT7(v) 0

#define V0            STATE[state_i]
#define VM1Over       STATE[state_i+_M1_-_R_]
#define VM1           STATE[state_i+_M1_]
#define VM2Over       STATE[state_i+_M2_-_R_]
#define VM2           STATE[state_i+_M2_]
#define VM3Over       STATE[state_i+_M3_-_R_]
#define VM3           STATE[state_i+_M3_]
#define Vrm1          STATE[state_i-1]
#define Vrm1Under     STATE[state_i+_R_-1]
#define Vrm2          STATE[state_i-2]
#define Vrm2Under     STATE[state_i+_R_-2]

#define newV0         STATE[state_i-1]
#define newV0Under    STATE[state_i-1+_R_]
#define newV1         STATE[state_i]
#define newVRm1       STATE[state_i-2]
#define newVRm1Under  STATE[state_i-2+R]

/*output transformation parameter*/
#define newVM2Over    STATE[state_i+_M2_-_R_+1]
#define newVM2        STATE[state_i+_M2_+1]
#define BITMASK 0x48000000

static unsigned int STATE[_R_];
static unsigned int z0,z1,z2;
static int state_i=0;

static unsigned int case_1(void);
static unsigned int case_2(void);
static unsigned int case_3(void);
static unsigned int case_4(void);
static unsigned int case_5(void);
static unsigned int case_6(void);

unsigned int (*WELLRNG44497)(void);

void InitWELLRNG44497(unsigned int *init ){
   int j;
   state_i=0;
   WELLRNG44497 = case_1;
   for(j=0;j<_R_;j++)
      STATE[j]=init[j];
}

unsigned int case_1(void){
  // state_i == 0
  z0 = (Vrm1Under & MASKL) | (Vrm2Under & MASKU);
  z1 = MAT0NEG(-24,V0) ^ MAT0POS(30,VM1);
  z2 = MAT0NEG(-10,VM2) ^ MAT3NEG(-26,VM3);
  newV1  = z1 ^ z2;
  newV0Under = MAT1(z0) ^ MAT0POS(20,z1) ^ MAT5(9,0xb729fcecU,0xfbffffffU,0x00020000U,z2) ^ MAT1(newV1);
  state_i = _R_-1;
  WELLRNG44497 = case_3;
  
  return (STATE[state_i] ^ (newVM2Over & BITMASK));
}

static unsigned int case_2(void){
  // state_i == 1
  z0 = (Vrm1 & MASKL) | (Vrm2Under & MASKU);
  z1 = MAT0NEG(-24,V0) ^ MAT0POS(30,VM1);
  z2 = MAT0NEG(-10,VM2) ^ MAT3NEG(-26,VM3);
  newV1 = z1 ^ z2;
  newV0 = MAT1(z0) ^ MAT0POS(20,z1) ^ MAT5(9,0xb729fcecU,0xfbffffffU,0x00020000U,z2) ^ MAT1(newV1);
  state_i=0;
  WELLRNG44497 = case_1;
  return (STATE[state_i] ^ (newVM2 & BITMASK)); 
}
static unsigned int case_3(void){
  // state_i+M1 >= R
  z0 = (Vrm1 & MASKL) | (Vrm2 & MASKU);
  z1 = MAT0NEG(-24,V0) ^ MAT0POS(30,VM1Over);
  z2 = MAT0NEG(-10,VM2Over) ^ MAT3NEG(-26,VM3Over);
  newV1 = z1 ^ z2;
  newV0 = MAT1(z0) ^ MAT0POS(20,z1) ^ MAT5(9,0xb729fcecU,0xfbffffffU,0x00020000U,z2) ^ MAT1(newV1);
  state_i--;
  if(state_i+_M1_<_R_)
    WELLRNG44497 = case_4;
    return (STATE[state_i] ^ (newVM2Over & BITMASK));
}

static unsigned int case_4(void){
  // state_i+M3 >= R
  z0 = (Vrm1 & MASKL) | (Vrm2 & MASKU);
  z1 = MAT0NEG(-24,V0) ^ MAT0POS(30,VM1);
  z2 = MAT0NEG(-10,VM2Over) ^ MAT3NEG(-26,VM3Over);
  newV1 = z1 ^ z2;
  newV0 = MAT1(z0) ^ MAT0POS(20,z1) ^ MAT5(9,0xb729fcecU,0xfbffffffU,0x00020000U,z2) ^ MAT1(newV1);
  state_i--;
  if (state_i+_M3_ < _R_)
    WELLRNG44497 = case_5;
	return (STATE[state_i] ^ (newVM2Over & BITMASK));
}

static unsigned int case_5(void){
  //state_i+M2 >= R
  z0 = (Vrm1 & MASKL) | (Vrm2 & MASKU);
  z1 = MAT0NEG(-24,V0) ^ MAT0POS(30,VM1);
  z2 = MAT0NEG(-10,VM2Over) ^ MAT3NEG(-26,VM3);
  newV1 = z1 ^ z2;
  newV0 = MAT1(z0) ^ MAT0POS(20,z1) ^ MAT5(9,0xb729fcecU,0xfbffffffU,0x00020000U,z2) ^ MAT1(newV1);
  state_i--;
  if(state_i+_M2_ < _R_)
    WELLRNG44497 = case_6;
	return (STATE[state_i] ^ (newVM2Over & BITMASK));
}

static unsigned int case_6(void){
  // 2 <= state_i <= R-M2-1
  z0 = (Vrm1 & MASKL) | (Vrm2 & MASKU);
  z1 = MAT0NEG(-24,V0) ^ MAT0POS(30,VM1);
  z2 = MAT0NEG(-10,VM2) ^ MAT3NEG(-26,VM3);
  newV1 = z1 ^ z2;
  newV0 = MAT1(z0) ^ MAT0POS(20,z1) ^ MAT5(9,0xb729fcecU,0xfbffffffU,0x00020000U,z2) ^ MAT1(newV1);
  state_i--;
  if(state_i == 1 )
    WELLRNG44497 = case_2;
    return (STATE[state_i] ^ (newVM2 & BITMASK));
}
