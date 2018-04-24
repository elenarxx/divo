#!/bin/bash
#SBATCH
module add clastro
module add openmpi/gcc/1.8.4

export OMP_NUM_THREADS=$1
export OMP_STACKSIZE=32000000

#time ./voids_v2.0.x
mpirun  --mca btl ^openib ./voids_v4.1.x




