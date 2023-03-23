#!/bin/bash

#SBATCH -p standard
#SBATCH -q standard
#SBATCH -J himeno
#SBATCH -t 10
#SBATCH -N 1

#NOW=$(date +"%S-%M-%H-%m-%d-%Y")

# Change the name of the perftools experiment file/directory
#export PAT_RT_EXPDIR_NAME=expfile.${CRAYPAT_LITE}.${SLURM_JOBID}
export NO_STOP_MESSAGE=1
time srun -n 8 himeno.exe #> my_output.${CRAYPAT_LITE}.${SLURM_JOBID} 2>&1


