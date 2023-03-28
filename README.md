# himeno
Fortran 90 version of himeno benchmark. Target offload just focuses on the Jacobi solver.

# Compile

Firstly one needs to load the appropriate modules for whichever architecture which is being targeted.
Then on most systems `make` will be sufficient to compile the benchmark successfully.

# Run

Submit slurm jobscript with `sbatch job.slurm`.
