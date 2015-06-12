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


compiler_name=""
compiler_exec=""

rng_flag=""
rng_name=""

flags   =["-std=c++11","-Wall"]
libs    =[]
includes=["-Isource -Isource/QMCEngine -Isource/Library"]
executable="qmcsgf"
show_command=false

message=[]

desc "Set intel compiler"
task :icc do
	compiler_name="icc"
	compiler_exec='icc -m64 -fast -fp-model precise -xhost -wd11021'
end

desc "Set g++ compiler"
task :gcc do
	compiler_name="g++"
	compiler_exec='g++ -O3'
end

desc "Set clang compiler"
task :clang do
	compiler_name="clang++"
	compiler_exec='clang++ -Ofast'
end

#desc "Use the Mersenne Twister random number generator"
task :rng_mt do
	rng_name = "+Mersenne-Twister"
	rng_flag = "-DRNG_MT"
end

#desc "Use the WELL44497 random number generator"
task :rng_well do
	rng_name = "+RNG_WELL"
	rng_flag = "-DRNG_WELL"
end


#desc "Use the DSMT random number generator"
task :rng_dsfmt do
	rng_name = "+DSFMT"
	rng_flag = "-DRNG_DSFMT -DHAVE_SSE2 -DDSFMT_MEXP=216091"
end


#desc "Use the STL Mersenne Twister random number generator"
task :cpprng_mt do
	rng_name = "+C++Mersenne-Twister"
	rng_flag = "-DCPPRNG_MT"
end

#desc "Use the TINYMT32 Mersenne Twister random number generator"
task :rng_tinymt32 do
	rng_name = "+Tiny-Mersenne-Twister(32bit)"
	rng_flag = "-DRNG_TINYMT32"
end

#desc "Use the TINYMT64 Mersenne Twister random number generator"
task :rng_tinymt64 do
	rng_name = "+Tiny-Mersenne-Twister(64bit)"
	rng_flag = "-DRNG_TINYMT64"
end

desc "Enable debugging"
task :debug do
	message << "+debug"
	flags << "-g -DDEBUG"
end

#desc "Use the MKL library for lapack"
task :mkl do

	# Get the location of ICC
	mklroot=ENV['MKLROOT']
	
	if(mklroot!=nil)

		if RUBY_PLATFORM.downcase.include?("darwin") # Is it a mac?
			libs     << "#{mklroot}/lib/libmkl_intel_lp64.a #{mklroot}/lib/libmkl_sequential.a #{mklroot}/lib/libmkl_core.a -lpthread -lm"
		elsif RUBY_PLATFORM.downcase.include?("linux") 	# if this a linux?
			libs     << "-openmp -Wl,--start-group  #{mklroot}/lib/intel64/libmkl_intel_lp64.a #{mklroot}/lib/intel64/libmkl_intel_thread.a #{mklroot}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm"
		end
	
		message << "+MKL"
		flags << "-DWITHMKL"
		includes <<  "-I#{mklroot}/include"

	end
end

#desc "use the open mpi library"
task :mpi do
	message << "+MPI"
	mpiroot="/opt/local"
	libs     << "-L#{mpiroot}/lib -lmpi_cxx -lmpi -lm"
	includes << "-I#{mpiroot}/include/openmpi"
	flags    << "-DUSEMPI"
end

# This task builds the code
task :compile do
	flags << rng_flag
	message << rng_name
	puts [compiler_name,message,executable]*" "
	source="source/qmcsgf.cpp"
	cmd="#{compiler_exec} #{flags*" "} #{includes*" "} #{libs*" "} #{source} -o #{executable}"
	puts cmd if show_command
	%x{#{cmd}}
end

desc "Show the compilation line"
task :show do
	show_command = true
end

task :default

# Set the default compiler
Rake::Task["icc"].execute

# Set the default randon number generator
Rake::Task["rng_mt"].execute

# After all tasks are done do a compilation
Rake.application.top_level_tasks << :compile

