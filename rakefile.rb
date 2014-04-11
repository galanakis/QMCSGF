#
# Compilation options
# There is one source file, SGFBoson.cpp.
#
# For both of them there are the following compilation options
# - Random Number Generator can be 
# 		a) MerseneTwister: -DRNG_MT
# 		b) Sprng: -DRNGSPRNG
#       c) Linear Congruence (Val's original): -DRNGLC (this is the default)
# - LAPACK provided by MKL can be included to provide extra functionality (diagonalizations)
# - OPENMPI can also be compiled for parallel computation
# - Progress bar on the terminal: -DCMDLINEPROGRESS
# - Debug (checks for negative occupancy): -DDEBUG
#


compiler=""
flags   =[]
libs    =[]
includes=["-Isource -Isource/Interface -Isource/Corvette -Isource/Library"]

message=[]
compiler_options=[]
compiler_name=""

desc "Set intel compiler"
task :icc do
	compiler_name="icc"
	compiler='/opt/intel/composerxe/bin/icc -m64 -fast -fp-model precise -std=c++11 -xhost -wd11021'
end

desc "Set clang compiler"
task :clang do
	compiler_name="clang"
	compiler='clang++ -O3 -std=c++11'
end

desc "Set gcc compiler"
task :gcc do
	compiler_name="gcc"
	compiler='g++-mp-4.9 -O3 -std=c++11'
end

#desc "Use the Mersenne Twister random number generator"
task :rng_mt do
	message << "+Mersenne-Twister"
	flags << "-DRNG_MT"
end

#desc "Use the WELL44497 random number generator"
task :rng_well do
	message << "+RNG_WELL"
	flags << "-DRNG_WELL"
end


#desc "Use the DSMT random number generator"
task :rng_dsfmt do
	message << "+DSFMT"
	flags << "-DRNG_DSFMT -DHAVE_SSE2 -DDSFMT_MEXP=216091"
end


#desc "Use the STL Mersenne Twister random number generator"
task :cpprng_mt do
	message << "+C++Mersenne-Twister"
	flags << "-DCPPRNG_MT"
end

#desc "Use the TINYMT32 Mersenne Twister random number generator"
task :rng_tinymt32 do
	message << "+Tiny-Mersenne-Twister(32bit)"
	flags << "-DRNG_TINYMT32"
end

#desc "Use the TINYMT64 Mersenne Twister random number generator"
task :rng_tinymt64 do
	message << "+Tiny-Mersenne-Twister(64bit)"
	flags << "-DRNG_TINYMT64"
end

desc "Print a progress bar in cerr"
task :cmdlineprogress do
	message << "+ProgressBar"
	flags << "-DCMDLINEPROGRESS"
end

desc "Enable debugging"
task :debug do
	message << "+debug"
	flags << "-g -DDEBUG"
end

desc "Print compiler warnings"
task :wall do
	compiler_options << "+warnings"
	flags << "-Wall"
end

#desc "Use the MKL library for lapack"
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

desc "use the open mpi library"
task :mpi do
	message << "+MPI"
	mpiroot="/opt/local"
	libs     << "-L#{mpiroot}/lib -lmpi_cxx -lmpi -lm"
	includes << "-I#{mpiroot}/include/openmpi"
	flags    << "-DUSEMPI"
end

task :compile do
	executable="corv_boson"
	puts [compiler_name,compiler_options,message,executable]*" "
	source="source/SGFBoson.cpp"
	cmd="#{compiler} #{source} #{flags*" "} #{includes*" "} #{libs*" "} -o #{executable}"
	%x{#{cmd}}
end

task :random => [:rng_mt]

task :icc   => [:mkl,:random]
task :gcc   => [:mkl,:random]
task :clang => [:mkl,:random]

desc "cmdlineprogress+debug+mpi"
task :all     =>   [:icc,:cmdlineprogress,:debug,:mpi,:compile]

desc "cmdlineprogress"
task :release =>   [:icc,:cmdlineprogress,:compile]

desc "intel compiler + wall + cmdlineprogress"
task :default => [:icc,:wall,:release,:compile]


