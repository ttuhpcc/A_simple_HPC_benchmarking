#!/bin/bash
#SBATCH --job-name=matlab-test
#SBATCH --output=%x.o%j
#SBATCH --error=%x.e%j
#SBATCH --partition nocona
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=10:00:00     ##10 hour runtime, you may change or remove this option. Default time is 48 hours
#SBATCH --reservation=nocona_test
####SBATCH --mem=512000MB  ##You can always request memory up to the max memory per node even when running serial cores, which is ~512000MB in this case. If you remove memory option, default memory per core (~3.9G/core) will be allocated. Alternatively, use --mem-per-core to request memory per core. --mem is for memory per node.

module load matlab
matlab -batch example_4
