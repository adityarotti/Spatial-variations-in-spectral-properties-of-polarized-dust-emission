#PBS -q regular
#PBS -l mppwidth=24
#PBS -l walltime=04:00:00
#PBS -j oe 

cd $PBS_O_WORKDIR
#PBS -N yr1
#PBS -e ./sysout/yr1.$PBS_JOBID.err
#PBS -o ./sysout/yr1.$PBS_JOBID.out
export OMP_NUM_THREADS=8
./run-gen-sim < submit/yr1param

