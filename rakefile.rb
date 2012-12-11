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
#     c) Linear Congruence (Val's original): -DRNGLC (this is the default)
# - LAPACK provided by MKL can be included to provide extra functionality (diagonalizations)
# - OPENMPI can also be compiled for parallel computation
# - Progress bar on the terminal: -DCMDLINEPROGRESS
# - Debug (checks for negative occupancy): -DDEBUG
#

executable="corv"
source="source/Corvette-SGF.cpp"

# Compiler selections
icc='/opt/intel/composerxe/bin/icc -fast -fp-model precise -fvisibility=protected'
gcc='g++-mp-4.7 -O3'
clang='clang++ -O3'

flags=""
libs=""
include="-Isource -Isource/Interface -Isource/Corvette -Isource/Library"
compiler=icc

# Use the Mersenne Twister random number generator
task :rng_mt do
	flags+=" -DRNG_MT"
end

# Print a progress bar in cerr
task :cmdlineprogress do
	flags+=" -DCMDLINEPROGRESS"
end

# enable debugging options
task :debug do
	flags+=" -DDEBUG"
end

# print compiler warnings
task :wall do
	flags+=" -Wall"
end

# use the MKL library for lapack
task :mkl do
	mklroot = '/opt/intel/composerxe/mkl'
	mkllib = "#{mklroot}/lib/libmkl_intel_lp64.a #{mklroot}/lib/libmkl_sequential.a #{mklroot}/lib/libmkl_core.a -lpthread -lm"
	mklinclude = "-I#{mklroot}/include"
	libs+=" #{mkllib}"
	include+=" #{mklinclude}"
end

# use the open mpi library
task :mpi do
	mpiinclude = '-I/opt/local/include/openmpi'
	mpilib = '-L/opt/local/lib -lmpi_cxx -lmpi -lm'
	libs="#{libs} #{mpilib}"
	include+=" #{mpiinclude}"
	flags+=" -DUSEMPI"
end

# Intel compiler
task :icc do
	compiler=icc
end

# clang compiler
task :clang do
	compiler=clang
end

# gcc compiler
task :gcc do
	compiler=gcc
end

task :all => [:rng_mt,:cmdlineprogress,:wall,:debug,:mkl,:mpi] do
end

task :std => [:rng_mt,:cmdlineprogress,:debug] do
end

task :corv => [:rng_mt,:cmdlineprogress] do
	executable="corv"
	source="source/Corvette-SGF.cpp"
	puts cmd="#{compiler} #{flags} #{include} #{libs} #{source} -o #{executable}"
	puts %x{#{cmd}}	
end

task :corv_boson => [:rng_mt,:cmdlineprogress,:mkl] do
	executable="corv_boson"
	source="source/SGFBoson.cpp"
	puts cmd="#{compiler} #{flags} #{include} #{libs} #{source} -o #{executable}"
	puts %x{#{cmd}}
end

task :corv_example => [:rng_mt,:cmdlineprogress,:mkl] do
	executable="corv_example"
	source="source/SGFBosonExample.cpp"
	puts cmd="#{compiler} #{flags} #{include} #{libs} #{source} -o #{executable}"
	puts %x{#{cmd}}	
end

task :default => [:std,:corv,:icc,:corv] do
end

task :everything => [:corv,:corv_boson,:corv_example] do
end

