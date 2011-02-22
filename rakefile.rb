Include="-I./Interface -I./Corvette -I."
Executable="Exec/corvette"+"."+%x{uname}
Source="Corvette-SGF.cpp"

task :default do
 puts %x{/opt/intel/composerxe/bin/icc -fast -Wall #{Include} #{Source} -o #{Executable}}
end
task :clang do
 git_version=%x{git rev-parse HEAD}
 %x{clang++ -O3 -Wall #{Include} #{Source} -o #{Executable}}
end
task :gcc do
 git_version=%x{git rev-parse HEAD}
 %x{g++-4.2 -O3 #{Include} #{Source} -o #{Executable}}
end

