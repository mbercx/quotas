export LD_BIND_NOW=1

module load VASP/5.4.4-intel-2018a

cd $PBS_O_WORKDIR

echo "Job started:" `/bin/date` >> out 2>&1
mpirun -genv LD_BIND_NOW=1 vasp_std >> out 2>&1
echo "Job finished:" `/bin/date` >> out 2>&1

grep Elapsed OUTCAR >> out