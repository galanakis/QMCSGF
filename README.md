# QMCSGF #

## Quantum Monte Carlo for bosons based on the Stochastic Green Function Algorithm.

### &copy; 2015 Dimitrios Galanakis 

This program is an implementation of the Stochastic Green Function algorithm for the simulation of Bosonic Hamiltonians. 
The Stochastic Green Function (SGF) algorithm is able to simulate any Hamiltonian that does not suffer from the so-called “sign problem”.

For more information of the specific implementation you can read 
http://lanl.arxiv.org/pdf/1209.0946.pdf
which is also found under doc/1209.0946.pdf.

### Compilation

The easiest way to compile the program is to use ruby's rake and the included rakefile.
If rake is not available you can try to edit the following line:

    clang++ -Ofast -std=c++11 -Wall  -Isource -Isource/QMCEngine -Isource/Library  source/qmcsgf.cpp -o qmcsgf

This will compile the program with minimal (but sufficient) options.

### Compilation options

The rake file contains a list of tasks which implement configuration options. It is possible to compile in a different random number generator and also include MPI support and enable debugging options and compiler warnings. 

    rake -T
shows all possible tasks. Consult the rake file more more options. A few examples about how to invoke rake

    rake icc

will compile using the Intel Compiler for a release. This is also the default task.

    rake gcc debug

will compile using g++ for test (it will include debugging options).

    rake clang rng_well mkl debug

will compile a minimal executable using clang, it will use the WELL random number generator and include debug options.

### Execution
Once the program is compiled it can be run as

    ./qmcsgf Input.json > Output.yaml

Input.json is a JSON formated file which contains running parameters which define the model, the running parameters the the desired types of measurements.
The program returns a YAML formated output which includes the input parameters and the results of the measurements. 

A few examples for the input files are included under the Test directory.

Note that during execution the program will print a progress bar in stderr. If you want to suppress it you can just redirect it to /dev/null: 

    ./qmcsgf Input.json 2>/dev/null > Output.yaml


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
