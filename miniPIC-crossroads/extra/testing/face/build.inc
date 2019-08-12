MPIHOME=/projects/install/rhel6-x86_64/sems/compiler/gcc/4.9.2/openmpi/1.6.5
INSTALL_DIR=$INSTALL_DIR
BUILD_FLAG=$BUILD_FLAG
BUILD_OPENMP=ON
BUILD_PTHREADS=OFF
if [ $BUILD_OPENMP == "ON" ]; then 
  COMPILE_FLAGS="-fopenmp"
fi
