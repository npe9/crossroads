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

Px=1
Py=1
Pz=1

Groups=200

quadType=2
Order=16
Polar=9
Azim=10

Xzones=3
Yzones=3
Zzones=4
Ranks=$((${Px}*${Py}*${Pz}))
#UMTHOME=${HOME}/UMT2013
#EXE=$UMTHOME/Install/bin/pyMPI
#LDPATH=$UMTHOME/Teton:$UMTHOME/cmg2Kull/sources:$UMTHOME/CMG_CLEAN/src:$UMTHOME/python
#export PYTHONPATH=$UMTHOME/python
#
# Input: SuOlsonTest $gridFileName $Groups $quadType $Order $Polar $Azim
# Allowed values: 1 <= quadType <= 2; 1 <= Polar <= 18; 1 <= Azim <= 22 
#
#=======================================================================================
#=======================================================================================

gridFileName=grid1MPI_3x3x4.cmg
for T in    16
do
export OMP_NUM_THREADS=$T     
srun   -n 1    -N 1   ../Teton/SuOlsonTest $gridFileName $Groups $quadType $Order $Polar $Azim
done

