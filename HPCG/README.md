# HPCG Benchmarking

## Introduction

High-Performance Conjugate Gradients (HPCG) is a benchmark designed to complement HPL by focusing on computational and data access patterns more representative of real-world applications.

## Step-by-Step Guide

### 1. Install Prerequisites

**Ensure GCC, other required libraries are available on your environment:**

``` bash
Operations was conducted on the Nocona Cluster
These libraries were available to be loaded as modules, they also can be downloaded from: 
```
```OpenMPI - 4.0.4 ``` - [https://www.open-mpi.org/software/ompi/v5.0/](https://www.open-mpi.org/software/ompi/v5.0/)  
```GCC - 9.2.0 ``` - [https://gcc.gnu.org/](https://gcc.gnu.org/)  
```OpenBLAS ``` - [https://github.com/OpenMathLib/OpenBLAS)](https://github.com/OpenMathLib/OpenBLAS)

### 2. Download and Extract HPCG

**Download HPCG:**

``` bash
wget http://www.hpcg-benchmark.org/downloads/hpcg-3.1.tar.gz
tar -xvzf hpcg-3.1.tar.gz
cd hpcg-3.1
```

### 3. Configure HPCG

**Edit the Makefile to set the required libraries and compilers:**

``` bash
cd setup
cp Makefile_GCC_MPI Makefile
nano Makefile
```

Set the MPI and compiler paths.

``` bash
echo $PATH
```
To find the paths to these libraries/modules


### 4. Compile HPCG

**Build HPCG:**

``` bash
make
```

### 5. Configure the Input File

**Edit the hpcg.dat file:**

``` bash
vi bin/hpcg.dat
```

Set the problem size and runtime.  
**Example:**

    128 128 128
    30

This sets a grid size of 128x128x128 and a runtime of 30 minutes. For optimal results, recommendations suggest runtime of atleast 60

### 6. Run the Benchmark

**Execute the benchmark:**

``` bash
cd bin
mpirun -np 4 ./xhpcg
```

This runs the HPCG benchmark using 4 MPI processes.

### Results

![Sample HPCG Result](/Images/hpcg.png)

## Conclusion

After running the benchmark, the output will provide a detailed performance report, including the achieved GFLOP rate and details on computational efficiency. This helps in evaluating the system’s performance for real-world applications.

## References

- [HPCG](https://www.hpcg-benchmark.org/faq/index.html)
- [HPCG Documentation](https://www.hpcg-benchmark.org/documentation/index.html)
- [HPCG GitHub Repository](https://github.com/hpcg-benchmark/hpcg)
- [HPCG - University of Luxembourg Tutorial Documentation](https://ulhpc-tutorials.readthedocs.io/en/latest/parallel/hybrid/HPCG/)