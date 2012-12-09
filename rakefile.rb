Include="-Isource/Interface -Isource/Corvette -Isource"
Executable="corvette"+"."+%x{uname}
Source="source/Corvette-SGF.cpp"
Flags="-DRNG_MT -DCMDLINEPROGRESS"
#Flags="-DRNG_MT"
BinPrefix=File.expand_path('~/')+"/.bin"

require 'fileutils'

task :default do
 puts cmd="/opt/intel/composerxe/bin/icc -fast -fp-model precise -Wall #{Flags} #{Include} #{Source} -o #{Executable}"
 puts %x{#{cmd}}
end
task :debug do
  puts cmd="/opt/intel/composerxe/bin/icc -fast -fp-model precise -Wall -DDEBUG #{Flags} #{Include} #{Source} -o #{Executable}"
  puts %x{#{cmd}}
end 
task :mpi do
# It is a pain to change the openmpi C++ compiler. It does not accept command line arguments. You need to export OMPI_CXX=icc.
# The simplest way is to avoid using the compiler wrappers.
# puts cmd="openmpicxx -fast -Wall -DDEBUG -DUSEMPI #{Flags} #{Include} #{Source} -o #{Executable}"
 puts cmd="/opt/intel/composerxe/bin/icc -fast -fp-model precise -Wall #{Flags} #{Include} -I/opt/local/include/openmpi #{Source} -o #{Executable} -L/opt/local/lib -lmpi_cxx -lmpi -lm"
 puts %x{#{cmd}}
end

task :clang do
 git_version=%x{git rev-parse HEAD}
 cmd="clang++ -O3 -Wall -DDEBUG #{Flags} #{Include} #{Source} -o #{Executable}"
 puts cmd
 puts %x{#{cmd}}
end
task :gcc do
 git_version=%x{git rev-parse HEAD}
 cmd="g++-mp-4.6 -O3 -Wall -DDEBUG #{Flags} #{Include} #{Source} -o #{Executable}"
 puts cmd
 puts %x{#{cmd}}
end

task :example do
 puts cmd="/opt/intel/composerxe/bin/icc -fast -fp-model precise -Wall -DDEBUG #{Flags} #{Include} source/Boson-SGF.cpp -o Example"
 puts %x{#{cmd}}
end

