#!/bin/bash
#
#SBATCH --exclusive
#SBATCH -J himeno_gpu
##SBATCH -o %x-j%j.out
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --hint=nomultithread

echo "Starting job $SLURM_JOB_ID at `date`"

ulimit -c unlimited
ulimit -s unlimited

gpu_bind=select_gpu.sh
cpu_bind="--cpu-bind=map_cpu:48,56,16,24,0,8,32,40"

source gpu_env.sh
export MPICH_OFI_NIC_POLICY=GPU #NUMA
export OMP_NUM_THREADS=1
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export NO_STOP_MESSAGE=1

exe1=../src/himeno_omp.exe
exe2=../src/himeno_acc.exe

time srun -N ${SLURM_NNODES} -n ${SLURM_NTASKS} -c ${OMP_NUM_THREADS} ${cpu_bind} ${gpu_bind} ${exe1}
wait
time srun -N ${SLURM_NNODES} -n ${SLURM_NTASKS} -c ${OMP_NUM_THREADS} ${cpu_bind} ${gpu_bind} ${exe2}

echo "... finished job $SLURM_JOB_ID at `date`"
