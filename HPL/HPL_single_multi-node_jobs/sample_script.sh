#!/bin/bash

#if [ -z "$1" ]; then
#    echo "Usage: sbatch $0 <nodelist>"
#    exit 1
#fi

NODELIST=$1
echo $NODELIST

sbatch --nodelist=$NODELIST <<EOF
#!/bin/bash
#SBATCH -J HPL_bench_intel
#SBATCH -o ./results/%x.%j.out
#SBATCH -e ./results/%x.%j.err
#SBATCH -p quanah
#SBATCH --account=hpcc
#SBATCH -N 1
#SBATCH --ntasks-per-node=36

sleep 10
echo "This is hostnode: \$(hostname)"

module purge
module load intel-oneapi/2024.2.1 oneapi-mpi/2021.13.1 oneapi-mkl/2024.2.1-omp hpl/2.3-mkl
hostname

export MKL_VERBOSE=1
mpirun xhpl

EOF
