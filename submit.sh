#!/bin/bash
#SBATCH -p i64m512u
#SBATCH -J Ag2Na6Se12_AB3C6_182
#SBATCH -n 128
#SBATCH -o job.%j.out
#SBATCH -e job.%j.err

ulimit -s unlimited

module load vasp/6.4.2

MPIEXEC=`which mpirun`
VASP_STD=`which vasp_std`
###The i64M512U queue has 64 cores per computer

$MPIEXEC -np 128 $VASP_STD > vasp.out
