#!/bin/bash 

export TOP=`pwd`/build

module load boost/1.55.0/gcc/4.9.2/base
module load gcc/4.9.2/openmpi/1.6.5
module load cmake/2.8.12
if [ ! -d $TOP ]; then mkdir -p $TOP; fi

cd $TOP

if [ -d Trilinos ]; then 
    cd Trilinos; git pull; cd ..
else 
    git clone software:/git/Trilinos
fi

if [ ! -d mini-PIC ]; then 
    git clone software:/git/mini-PIC
fi
cd mini-PIC
git checkout tet-mesh
git pull
cd ..



if [ ! -d build-trilinos ]; then mkdir build-trilinos; fi

cd  build-trilinos
cp ../mini-PIC/extra/testing/face/* .
./do-configure.sh

make -j16 install

cd ..
cd mini-PIC
if [ ! -d build ]; then mkdir build; fi
cd build
cmake -DTrilinos_PREFIX=$TOP/install -DENABLE_PIC_TESTS=ON ..
make -j16
make test
exit $?


