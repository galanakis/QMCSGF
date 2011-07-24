name=Corvette_`git rev-parse --short HEAD`.`uname`
g++ -O3 source/Corvette-SGF.cpp -Isource -Isource/Corvette -Isource/Interface/ -o $name
ln -sf $name Corvette.`uname`
