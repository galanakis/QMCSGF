#ifndef __CONVENTIONS__
#define __CONVENTIONS__

#include <climits>
#include <limits>

#include <complex>

namespace SGF {

typedef unsigned int uint;
typedef unsigned long ulong;
typedef double MatrixElement;
typedef long double _float_accumulator;
typedef unsigned long long _integer_counter;

typedef uint occupancy_t;
const uint Infinity=UINT_MAX; // The maximum possible integer is infinity


const int Sign[]={-1,1};
enum {LEFT,RIGHT};            // Directions
enum {ANNIHILATION,CREATION}; // Operator type
enum {REMOVE,ADD};            // Action
enum {DESTROY,CREATE};        // another name for action
typedef enum {Canonical,GrandCanonical} ensemble_t;

// Absolute value
template<class T> inline T Abs(T x) {return x<0?-x:x;}
template<class T> inline T Max(T a,T b) {return (a>b)?a:b;}
template<class T> inline T Min(T a,T b) {return (a>b)?b:a;}

}

#endif
