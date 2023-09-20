# himeno
Fortran 90 version of himeno benchmark with ports for OpenACC and OpenMP offload to GPUs. The OpenMP target offload just focuses on the Jacobi solver.

# Compile

Firstly one needs to load the appropriate modules for whichever architecture is being targeted. In the run directory `gpu_env.sh` gives an example of this for using the Cray compiler and targeting AMD MI205X hardware.

Then on most systems `make` will be sufficient to compile the benchmark successfully.

# Run

Submit slurm jobscript with `sbatch job.slurm`.
