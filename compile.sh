name=Corvette_`git rev-parse --short HEAD`.`uname`
icc -fast -DCMDLINEPROGRESS source/Corvette-SGF.cpp -Isource -Isource/Corvette -Isource/Interface/ -o $name
ln -sf $name Corvette.`uname`
