#!/bin/bash

#MSUB -S bash
#MSUB -l nodes=4:ppn=16
#MSUB -l walltime=02:00:00
#MSUB -q pbatch
#MSUB -m be
#MSUB -V
#MSUB -A coral
#MSUB    -j oe
#MSUB    -o /nfs/tmp2/ghosh4/ADC/final/UMT2_SW_Review/Teton/4nodes_16rpn_1omp.out

./run3x3x4_g200_p9a10_D1T1.bash
./run3x3x4_g200_p9a10_D1T16.bash
./run3x3x4_g200_p9a10_D1T4.bash
./run3x3x4_g200_p9a10_D64T1.bash
./run3x6x8_g200_p9a10_D16T4.bash
