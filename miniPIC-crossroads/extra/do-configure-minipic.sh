#!/bin/bash

cmake 	\
	-DTrilinos_PREFIX=$PWD/../install \
	-DCMAKE_CXX_COMPILER=`which CC` \
	../../
