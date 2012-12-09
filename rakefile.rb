#
# Compilation options
# There are two main source files:
# - Corvette-SGF.cpp, uses Val's ancient interface
# - Boson-SGF.cpp, is the new version
#
# For both of them there are the following compilation options
# - Random Number Generator can be 
# 		a) MerseneTwister: -DRNG_MT
# 		b) Sprng: -DRNGSPRNG
#     c) Linear Congruence (Val's original): -DRNGLC
# - LAPACK provided by MKL can be included to provide extra functionality (diagonalizations)
# - OPENMPI can also be compiled for parallel computation
# - Progress bar on the terminal: -DCMDLINEPROGRESS
# - Debug (checks for negative occupancy): -DDEBUG
#

# MKL Library
MKLROOT = '/opt/intel/composerxe/mkl'
MKLLIB = "#{MKLROOT}/lib/libmkl_intel_lp64.a #{MKLROOT}/lib/libmkl_sequential.a #{MKLROOT}/lib/libmkl_core.a -lpthread -lm"
MKLINCLUDE = "-I#{MKLROOT}/include"

# OPENMPI library
MPIINCLUDE = '-I/opt/local/include/openmpi'
MPILIB = '-L/opt/local/lib -lmpi_cxx -lmpi -lm'

INCLUDE="-Isource/Interface -Isource/Corvette -Isource"
EXECUTABLE="corvette"+"."+%x{uname}
SOURCE="source/Corvette-SGF.cpp"
FLAGS="-DRNG_MT -DCMDLINEPROGRESS"
#FLAGS="-DRNG_MT"
LIBS=""

# Compiler selections
ICC='/opt/intel/composerxe/bin/icc -fast -fp-model precise -Wall'
GCC='g++-mp-4.7-O3 -Wall'
CLANG='clang++ -O3 -Wall'

require 'fileutils'

task :default do
 puts cmd="#{ICC} #{FLAGS} #{INCLUDE} #{LIBS} #{SOURCE} -o #{EXECUTABLE}"
 puts %x{#{cmd}}
end
task :debug do
  puts cmd="#{ICC} -DDEBUG #{FLAGS} #{INCLUDE} #{LIBS} #{SOURCE} -o #{EXECUTABLE}"
  puts %x{#{cmd}}
end 
task :mpi do
# It is a pain to change the openmpi C++ compiler. It does not accept command line arguments. You need to export OMPI_CXX=icc.
# The simplest way is to avoid using the compiler wrappers.
# puts cmd="openmpicxx -fast -Wall -DDEBUG -DUSEMPI #{FLAGS} #{INCLUDE} #{SOURCE} -o #{EXECUTABLE}"
 puts cmd="#{ICC} #{FLAGS} #{INCLUDE} #{MPIINCLUDE} #{LIBS} #{SOURCE} #{MPILIB} -o #{EXECUTABLE}"
 puts %x{#{cmd}}
end

task :clang do
 git_version=%x{git rev-parse HEAD}
 cmd="#{CLANG} -DDEBUG #{FLAGS} #{INCLUDE} #{LIBS} #{SOURCE} -o #{EXECUTABLE}"
 puts cmd
 puts %x{#{cmd}}
end
task :gcc do
 git_version=%x{git rev-parse HEAD}
 cmd="#{GCC} -DDEBUG #{FLAGS} #{INCLUDE} #{LIBS} #{SOURCE} -o #{EXECUTABLE}"
 puts cmd
 puts %x{#{cmd}}
end

task :example do
 puts cmd="#{ICC} -DDEBUG #{FLAGS} #{INCLUDE} #{MKLINCLUDE} #{LIBS} #{MKLLIB} source/Boson-SGF.cpp -o Example"
 puts %x{#{cmd}}
end

