#!/bin/bash

for g in conservative ondemand userspace powersave performance; do
	echo "Setting governor to: $g"
	cpupower frequency-set -g $g	

	echo "Running MATLAB job"
	su - ssingha /lustre/work/ssingha/A_simple_HPC_benchmarking_Quanah/matlab/Example_4/matlab_test.sh
done
