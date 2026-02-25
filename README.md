# himeno
Fortran 90 version of himeno benchmark with ports for OpenACC and OpenMP offload to GPUs. This is a simple example and the use OpenMP target offload just focuses on the Jacobi solver, the MPI is not GPU-aware.

# Compile

Firstly one needs to load the appropriate modules for whichever architecture is being targeted. In the run directory `gpu_env.sh` gives an example of this for using the Cray compiler and targeting AMD MI205X hardware.

Then on most systems `make` will be sufficient to compile the benchmark successfully.

To setup a differenet model one can change the header file by using the provided script paramset.sh e.g XL model with a decomposition of 2x2x1 `bash paramset.sh XL 2 2 1`

# Run

Example scripts for environment, binding and job submission are given for running on AMD GPUs.

To run on Nvidia GPUs we would require the appropriate CUDA and accelerator modules e.g for Nvidia A100 GPUs

```
>$ cat gpu_env_nvidia.sh
#!/bin/bash
module load PrgEnv-cray
module load craype-x86-milan
module load cudatoolkit
module load craype-accel-nvidia80

```
and appropriate binding script e.g.
```
>$ cat select_gpus_nvidia.sh
#!/bin/bash
export CUDA_VISIBLE_DEVICES=$SLURM_LOCALID
exec $*

```
