#!/bin/bash
#SBATCH --job-name=PIPERe180
#SBATCH --output=%x.o%j
#SBATCH --error=%x.e%j
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --partition=nocona
##SBATCH --reservation=nocona_test
#SBATCH -t 10:00:00
##SBATCH --constraint="[(rack23|rack24)*16]"
##SBATCH --mail-user=sagnik.singha@ttu.edu
#SBATCH --mail-type=all   
## SBATCH --exclude=cpu-26-22
#SBATCH --exclude=cpu-24-53,cpu-23-5

#ml restore OPEN3
echo $SLURM_JOB_NODELIST
export OMP_NUM_THREADS=1

ml gcc/10.1.0 openmpi/3.1.6 openblas/0.3.10 netcdf-c/4.7.3-mpi netcdf-fortran/4.5.2-mpi fftw/3.3.8-mpi-openmp netlib-lapack/3.8.0 gnuplot/5.2.8 matlab/R2020b

#ml  gcc/10.1.0 openmpi/3.1.6  fftw/3.3.8-mpi-openmp  netlib-lapack/3.8.0
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/apps/nfs/spack/opt/spack/linux-centos8-zen2/gcc-10.1.0/netcdf-c-4.7.3-ljb3z4dr2m6fke4ck6pkitprbpfu6esf/lib/
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/apps/nfs/spack/opt/spack/linux-centos8-zen2/gcc-10.1.0/netcdf-fortran-4.5.2-y6v6qtrpb2ynw6e2bemgkvak6uflggpp/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/jieyao/CODE_TEST/CHANNEL/fftw-3.3.7-build/lib

#ml gcc/10.1.0 openmpi/3.1.6 hdf5/1.10.6-mpi fftw/3.3.8-mpi-openmp gsl/2.5 
#ml  gcc/10.1.0 openmpi/3.1.6 netcdf-fortran/4.5.2-mpi fftw/3.3.8-mpi-openmp netcdf-c/4.7.3-mpi netlib-lapack/3.8.0
module list
pwd
date
#echo $SLURM_JOB_NODELIST

mpirun  --map-by l3cache -mca pml ucx -x UCX_NET_DEVICES=mlx5_2:1 -mca coll_hcoll_enable 1 -x HCOLL_MAIN_IB=mlx5_2:1 ./main.out
#mpirun -n 4096 main.out
