# QMCSGF #

## Quantum Monte Carlo for bosons based on the Stockastic Green Function Algorithm.

### &copy; 2015 Dimitrios Galanakis 

This program is an implementation of the Stochastic Green Function algorithm for the simulation of bosonic Hamiltonians. 
The Stochastic Green Function (SGF) algorithm is able to simulate any Hamiltonian that does not suffer from the so-called “sign problem”.

For more information of the specific implementation you can read 
http://lanl.arxiv.org/pdf/1209.0946.pdf
which is also found under doc/1209.0946.pdf.

### Compilation
To compile the program you need intel's Math Kernel library and ruby's rake. 
If the MKLROOT variable is defined in the shell and the icc compiler is present,
simply typing rake should be able to compile the program. If there are errors try to 
start by tuning the following line (works on OS X):

    /opt/intel/composerxe/bin/icc -m64 -fast -fp-model precise -std=c++11 -xhost -wd11021 source/SGFBoson.cpp -DRNG_MT -Wall -DCMDLINEPROGRESS -Isource -Isource/Interface -Isource/Corvette -Isource/Library -I/opt/intel/composer_xe_2015.3.187/mkl/include /opt/intel/composer_xe_2015.3.187/mkl/lib/libmkl_intel_lp64.a /opt/intel/composer_xe_2015.3.187/mkl/lib/libmkl_sequential.a /opt/intel/composer_xe_2015.3.187/mkl/lib/libmkl_core.a -lpthread -lm -o qmcsgf

### Compilation options

The rake file contains a list of tasks which implement configuration options. It is possible to compile in a different random number generator and also include MPI support and enable debugging options and compiler warnings. 

    rake -T
shows all possible tasks. Consult the rake file more more options. A few examples about how to invoke rake

    rake icc release

will compile using the Intel Compiler for a release. This is also the default task.

    rake gcc test

will compile using g++ for test (it will include debugging options).

    rake clang rng_well debug standard

will compile a minimal executable without extra options, but will use the WELL random number generator and include debug options.

### Execution
Once the program is compiled it can be run as

    ./qmcsgf Input.json > Output.yaml

Input.json is a JSON formated file which contains running parameters which define the model, the running parameters the the desired types of measurements.
The program returns a YAML formated output which includes the input parameters and the results of the measurements.

A few examples for the input files are included under the Test directory.

### Copying
&copy; 2015 Dimitrios Galanakis

QMCSGF is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 3.

QMCSGF is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with QMCSGF.  If not, see <http://www.gnu.org/licenses/>.
