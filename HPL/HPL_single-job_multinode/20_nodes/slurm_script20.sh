#!/bin/bash
#SBATCH --job-name=HPL_bench_num_nodes_20
#SBATCH --nodes 20
#SBATCH --ntasks-per-node=36
#SBATCH -p quanah
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
#SBATCH --account=hpcc

# Load required modules
module purge
module load gcc/14.2.0 openmpi/4.1.6 hpl/2.3

# The xhpl executable expects the HPL.dat file to be present on the current directory
#mpirun xhpl
hostname
time -p mpirun --mca btl self,vader --map-by l3cache --report-bindings -np ${SLURM_NTASKS} xhpl
