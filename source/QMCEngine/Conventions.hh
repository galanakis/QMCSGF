/******************************************************************\

Copyright 2015 Dimitrios Galanakis

This file is part of QMCSGF

QMCSGF is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 3.

QMCSGF is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with QMCSGF.  If not, see <http://www.gnu.org/licenses/>.

\*******************************************************************/

#ifndef __CONVENTIONS__
#define __CONVENTIONS__

#include <climits>
#include <limits>

namespace SGF {

typedef unsigned int uint;
typedef double MatrixElement;
typedef long double _float_accumulator;
typedef unsigned long long _integer_counter;

typedef uint occupancy_t;
const uint Infinity=UINT_MAX; // The maximum possible integer is infinity


const int Sign[]= {-1,1};
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
