#
# Compilation options
# There is one source, Boson-SGF.cpp.
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

# Compiler selections
icc='/opt/intel/composerxe/bin/icc -m64 -fast -fp-model precise -std=c++11 -xhost -wd11021'
gcc='g++-mp-4.9 -O3 -std=c++11'
clang='clang++ -O3 -std=c++11'

compiler=""
flags   =[]
libs    =[]
includes=["-Isource -Isource/Interface -Isource/Corvette -Isource/Library"]

message=[]
compiler_options=[]

task :cpprng_mt do
	message << "+C++Mersenne-Twister"
	flags << "-DCPPRNG_MT"
end

# Use the Mersenne Twister random number generator
task :rng_mt do
	message << "+Mersenne-Twister"
	flags << "-DRNG_MT"
end

# Print a progress bar in cerr
task :cmdlineprogress do
	message << "+ProgressBar"
	flags << "-DCMDLINEPROGRESS"
end

# enable debugging options
task :debug do
	message << "+debug"
	flags << "-g -DDEBUG"
end

# print compiler warnings
task :wall do
	compiler_options << "+warnings"
	flags << "-Wall"
end

# use the MKL library for lapack
task :mkl do
	message << "+MKL"

	# Get the location of ICC
	mklroot=ENV['MKLROOT']
	
	if RUBY_PLATFORM.downcase.include?("darwin") # Is it a mac?
		libs     << "#{mklroot}/lib/libmkl_intel_lp64.a #{mklroot}/lib/libmkl_sequential.a #{mklroot}/lib/libmkl_core.a -lpthread -lm"
	elsif RUBY_PLATFORM.downcase.include?("linux") 	# if this a linux?
		libs     << "-openmp -Wl,--start-group  #{mklroot}/lib/intel64/libmkl_intel_lp64.a #{mklroot}/lib/intel64/libmkl_intel_thread.a #{mklroot}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm"
	end
	
	includes <<  "-I#{mklroot}/include"
end

# use the open mpi library
task :mpi do
	message << "+MPI"
	mpiroot="/opt/local"
	libs     << "-L#{mpiroot}/lib -lmpi_cxx -lmpi -lm"
	includes << "-I#{mpiroot}/include/openmpi"
	flags    << "-DUSEMPI"
end

# Intel compiler
task :icc do
	compiler_options << "intel"
	compiler=icc
end

# clang compiler
task :clang do
	compiler_options << "clang"
	compiler=clang
end

# gcc compiler
task :gcc do
	compiler_options << "gcc"
	compiler=gcc
end

task :all     => [:rng_mt,:cmdlineprogress,:mkl,:debug,:mpi]

task :release => [:rng_mt,:cmdlineprogress,:mkl]

task :debug =>   [:rng_mt,:cmdlineprogress,:mkl,:debug]

task :corv_boson => [:release] do
	puts "Compiler: "+compiler_options*" "
	puts "Options:  "+message*" "
	executable="corv_boson"
	source="source/SGFBoson.cpp"
	cmd="#{compiler} #{source} #{flags*" "} #{includes*" "} #{libs*" "} -o #{executable}"
	puts %x{#{cmd}}
end

task :default => [:icc,:wall,:corv_boson] do
end
