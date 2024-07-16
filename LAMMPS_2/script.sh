#!/bin/bash
#SBATCH --job-name=HPCCTest
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --output=%N.%j.out
#SBATCH --error=%N.%j.err

module load gcc/10.1.0 openmpi/4.1.4 lammps/20230802-mpi-openmp

# The most efficient way of running LAMMPS on Nocona nodes with AMD architecture
# is to leverage the multi-processing and multi-threading feature in LAMPPS.
# In order to do that, we go with 32 MPI process per node (out of 128) and after
# mapping each process to L3cache, we let them to swan 4 OMP threads across four
# CPU cores which are sharing the same L3cache:

export OMP_NUM_THREADS=4

mpirun --mca btl self,vader --map-by l3cache lmp -in input.relaxation
