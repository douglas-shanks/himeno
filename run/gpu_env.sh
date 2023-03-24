#!/bin/bash

ulimit -c unlimited
ulimit -s unlimited

module load PrgEnv-cray
module load craype-x86-trento
module load craype-accel-amd-gfx90a
module load rocm
# Required Environment settings
export MPICH_GPU_SUPPORT_ENABLED=1
export LD_LIBRARY_PATH=${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH}
