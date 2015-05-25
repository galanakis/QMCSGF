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

compiler=""
rng=""
rng_name=""
flags   =[]
libs    =[]
includes=["-Isource -Isource/Interface -Isource/QMCEngine -Isource/Library"]
executable="qmcsgf"

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

desc "Set g++ compiler"
task :gcc do
	compiler_name="g++"
	compiler='g++ -O3 -std=c++11'
end

#desc "Use the Mersenne Twister random number generator"
task :rng_mt do
	rng_name = "+Mersenne-Twister"
	rng = "-DRNG_MT"
end

#desc "Use the WELL44497 random number generator"
task :rng_well do
	rng_name = "+RNG_WELL"
	rng = "-DRNG_WELL"
end


#desc "Use the DSMT random number generator"
task :rng_dsfmt do
	rng_name = "+DSFMT"
	rng = "-DRNG_DSFMT -DHAVE_SSE2 -DDSFMT_MEXP=216091"
end


#desc "Use the STL Mersenne Twister random number generator"
task :cpprng_mt do
	rng_name = "+C++Mersenne-Twister"
	rng = "-DCPPRNG_MT"
end

#desc "Use the TINYMT32 Mersenne Twister random number generator"
task :rng_tinymt32 do
	rng_name = "+Tiny-Mersenne-Twister(32bit)"
	rng = "-DRNG_TINYMT32"
end

#desc "Use the TINYMT64 Mersenne Twister random number generator"
task :rng_tinymt64 do
	rng_name = "+Tiny-Mersenne-Twister(64bit)"
	rng = "-DRNG_TINYMT64"
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

#desc "use the open mpi library"
task :mpi do
	message << "+MPI"
	mpiroot="/opt/local"
	libs     << "-L#{mpiroot}/lib -lmpi_cxx -lmpi -lm"
	includes << "-I#{mpiroot}/include/openmpi"
	flags    << "-DUSEMPI"
end

# This builds the code
task :compile do
	flags << rng
	message << rng_name
	puts [compiler_name,compiler_options,message,executable]*" "
	source="source/qmcsgf.cpp"
	cmd="#{compiler} #{source} #{flags*" "} #{includes*" "} #{libs*" "} -o #{executable}"
	%x{#{cmd}}
end

task :icc   => [:mkl,:rng_mt,:wall]
task :gcc   => [:mkl,:rng_mt,:wall]
task :clang => [:mkl,:rng_mt,:wall]

desc "No extra features"
task :standard =>   [:compile]

desc "Standard features (cmdlineprogress)"
task :release =>   [:cmdlineprogress,:compile]

desc "Debugging options ("
task :test =>   [:cmdlineprogress,:debug,:compile]

desc "All extras (cmdlineprogress, debug and mpi)"
task :all =>   [:cmdlineprogress,:debug,:mpi,:compile]

desc "intel compiler + wall + cmdlineprogress"
task :default => [:icc,:release]


