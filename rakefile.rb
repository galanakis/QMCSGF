Include="-Isource/Interface -Isource/Corvette -Isource"
Executable="corvette"+"."+%x{uname}
Source="source/Corvette-SGF.cpp"
Flags="-DRNG_MT -DCMDLINEPROGRESS"
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
 puts cmd="mpicxx -cxx=icc -fast -Wall -DDEBUG -DUSEMPI #{Flags} #{Include} #{Source} -o #{Executable}"
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

