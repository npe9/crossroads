#!/bin/bash -xe

if [ $# -gt 1 ]; then
    echo $0 " -debug | -opt )"
    exit -1
fi

DEBUG=0
SUBPATH="opt"
BUILD_FLAG="RELEASE"
if [ "$#" -gt 0 ] && [ "$1" == '-debug' ]; then
    DEBUG=1
    SUBPATH="debug"
    BUILD_FLAG="DEBUG"
fi
export BUILD_FLAG
export TOP=`pwd`

module load boost/1.55.0/gcc/4.9.2/base
module load gcc/4.9.2/openmpi/1.6.5
module load cmake/2.8.12

export BUILD_DIR=$TOP/build/build-trilinos-${SUBPATH}
export INSTALL_DIR=$TOP/build/install-$SUBPATH
if [ ! -d $BUILD_DIR ]; then mkdir -p $BUILD_DIR; fi

cd $BUILD_DIR
cp $TOP/mini-PIC/jenkins-test-scripts/* .

./do-configure.sh

make -j16 install

cd $TOP/mini-PIC
MINIPIC_BUILD_DIR=build-${SUBPATH}
if [ ! -d $MINIPIC_BUILD_DIR ]; then mkdir $MINIPIC_BUILD_DIR; fi
cd $MINIPIC_BUILD_DIR
cmake -DTrilinos_PREFIX=$INSTALL_DIR -DENABLE_PIC_TESTS=ON ..
make -j16
make test
exit $?


