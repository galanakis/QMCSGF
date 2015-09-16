#!/bin/sh

########################################################################
# Copyright 2015 Dimitrios Galanakis
#
# This file is part of QMCSGF
#
# QMCSGF is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation version 3.
# 
# QMCSGF is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with QMCSGF.  If not, see <http://www.gnu.org/licenses/>.
# 
#######################################################################


EXECUTABLE="qmcsgf"
SOURCE="source/qmcsgf.cpp"
INCLUDES="-Isource -Isource/QMCEngine -Isource/Library"
FLAGS="-std=c++11 -Wall"
LIBS=""

RNGFLAG="-DRNG_MT"


# checks if the argument is an executable in the path
function isexec() {
	command -v $1 >/dev/null 2>&1
}

function seticc() {
	COMPILER="icc -m64 -fast -fp-model precise -xhost -wd11021"	
}

function setclang() {
	COMPILER="clang++ -Ofast"	
}

function setgpp() {
	COMPILER="g++ -O3"
}

#It selects a compiler by searching in the path for icc,clang and g++ in that order.
function autocpp() {
	if isexec icc
	then
		seticc
	elif isexec clang
	then
		setclang
	else isexec g++
		setgpp
	fi
}

#Use the MKL library for lapack
function setmkl() {

	if [ $MKLROOT ]
	then	
		if [ `uname` == "Darwin" ]
		then
			LIBS="$LIBS $MKLROOT/lib/libmkl_intel_lp64.a $MKLROOT/lib/libmkl_sequential.a $MKLROOT/lib/libmkl_core.a -lpthread -lm"
		elif [ `uname` == "Linux" ]
		then
			LIBS="$LIBS -openmp -Wl,--start-group  $MKLROOT/lib/intel64/libmkl_intel_lp64.a $MKLROOT/lib/intel64/libmkl_intel_thread.a $MKLROOT/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm"
		fi
		FLAGS="$FLAGS -DWITHMKL"
		INCLUDES="$INCLUDES -I$MKLROOT/include"
	fi

}

#use the open mpi library
function setmpi() {

	MPIROOT="/opt/local"
	LIBS="$LIBS -L#{mpiroot}/lib -lmpi_cxx -lmpi -lm"
	INCLUDES="$INCLUDES -I#{mpiroot}/include/openmpi"
	FLAGS="$FLAGS -DUSEMPI"

}

function setdebug() {
	FLAGS="$FLAGS -g -DDEBUG"
}


WITHMKL=false
SHOWCMD=false

# Parse the command line keys
while [[ $# > 0 ]]
do
key="$1"

case $key in
    icc)
	seticc
    ;;
    gcc)
    setgpp
    ;;
    clang)
	setclang
    ;;
    mkl)        # Enable MKL
    setmkl
    ;;
    mpi)        # Enable MPI
	setmpi
	;;
	debug)      # Enable debug options
	setdebug
	;;
	show)       # set a flag that determines if the compilation command is shown
	SHOWCMD=true
	;;
	rng_mt)     # Linear congruences random number generator
	RNGFLAG="-DRNG_LC"
	;;
	rng_well)   # WELL44497 random number generator
	RNGFLAG="-DRNG_WELL"
	;;
	rng_dsfmt)  # DSMT random number generator
	RNGFLAG="-DRNG_DSFMT -DHAVE_SSE2 -DDSFMT_MEXP=216091"
	;;
	:rng_stlmt)  # STL Mersenne Twister random number generator
	RNGFLAG="-DRNG_STLMT"
	;;
	rng_tinymt32) # TINYMT32 Mersenne Twister random number generator
	RNGFLAG="-DRNG_TINYMT32"
	;;
	rng_tinymt64) # TINYMT64 Mersenne Twister random number generator
	RNGFLAG="-DRNG_TINYMT64"
	;;
    *)
	echo "Invalid key $key"
	exit
    ;;
esac
shift
done


# If the compiler was not set, use the default
if [ -z "$COMPILER" ]
then
	autocpp
fi

# This is the compilation command
CMD="$COMPILER $FLAGS $RNGFLAG $INCLUDES $LIBS $SOURCE -o $EXECUTABLE"

# Show the command if requested
if $SHOWCMD
then
	echo $CMD
fi

# cd into the directory of the script and run the compilation
cd "$( dirname "${BASH_SOURCE[0]}" )" && $CMD

