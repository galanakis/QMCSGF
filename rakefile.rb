Include="-I./Interface -I./Corvette -I."
Executable="corvette"+"."+%x{uname}
Source="Corvette-SGF.cpp"
Flags="-DRNG_MT"

task :default do
 puts cmd="/opt/intel/composerxe/bin/icc -fast -Wall #{Flags} #{Include} #{Source} -o #{Executable}"
 puts %x{#{cmd}}
end
task :clang do
 git_version=%x{git rev-parse HEAD}
 %x{clang++ -O3 -Wall #{Include} #{Source} -o #{Executable}}
end
task :gcc do
 git_version=%x{git rev-parse HEAD}
 %x{g++-mp-4.5 -O3 -Wall #{Include} #{Source} -o #{Executable}}
end

