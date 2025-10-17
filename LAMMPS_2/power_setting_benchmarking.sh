#!/bin/bash

for g in conservative ondemand userspace powersave performance; do
	echo "Setting governor to: $g"
	cpupower frequency-set -g $g	

	echo "Running LAMMPS job"
	su - ssingha -c 'cd /lustre/work/ssingha/A_simple_HPC_benchmarking_Quanah/LAMMPS && ./lammps_test.sh' >> "$g-output.txt" 2>&1
done
