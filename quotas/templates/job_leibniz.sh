#!/bin/bash
#PBS -N {{name}}
#PBS -q batch
#PBS -L tasks={{nodes}}:lprocs=28
##PBS -l pmem=3GB
#PBS -l walltime=72:00:00
#PBS -j eo

export LD_BIND_NOW=1

module purge
module load VASP/5.4.4-intel-2016b

cd $PBS_O_WORKDIR

echo "Job started:" `/bin/date` >> out 2>&1
mpirun -genv LD_BIND_NOW=1 vasp_std >> out 2>&1
echo "Job finished:" `/bin/date` >> out 2>&1

grep Elapsed OUTCAR >> out