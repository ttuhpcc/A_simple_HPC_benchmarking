
# FFTW Benchmarking

## Introduction

FFTW (Fastest Fourier Transform in the West) is a highly efficient library for computing discrete Fourier transforms (DFTs) in one or more dimensions. Benchmarking FFTW can help assess its performance and efficiency in various applications.

## Step-by-Step Guide

### 1. Install Prerequisites

**Ensure GCC, other required libraries are available on your environment:**

``` bash
Operations was conducted on the Nocona Cluster 
```
```OpenMPI - 4.0.4 ``` - [https://www.open-mpi.org/software/ompi/v5.0/](https://www.open-mpi.org/software/ompi/v5.0/)  
```GCC - 9.2.0 ``` - [https://gcc.gnu.org/](https://gcc.gnu.org/)

### 2. Download and Extract FFTW

**Download FFTW:**

``` bash
wget https://www.fftw.org/fftw-3.3.10.tar.gz
tar -xvzf fftw-3.3.10.tar.gz
cd fftw-3.3.10
```

### 3. Configure FFTW

Set the MPI and compiler paths.

``` bash
echo $PATH
```
To find the paths to these libraries/modules

#### To Build Just Serial FFTW
**Edit the Makefile to set the required libraries and compilers:**

``` bash
./configure --prefix=/your/install/path --enable-shared=yes
make
make install
```

#### To Build MPI FFTW
**Edit the Makefile to set the required libraries and compilers:**

``` bash
./configure --enable-mpi --prefix=/your/install/path --enable-shared=yes
make
make install
```


### 4. Compile FFTW

**Compile Serial FFTW Program:**

``` bash
gcc -I/your/install/pathinclude sample.c -lm -L/your/install/path/lib -lfftw3 -o outputprogram
```
**Compile MPI FFTW Program:**

``` bash
mpicc -I/your/install/mpi-path/include sample.c -lm -L/your/mpi-install/path/lib -lfftw3_mpi -lfftw3 -o outputprogram
```

> **Note:** Replace `/your/install/path` and `/your/mpi-install/path` with your specific FFTW installation directory.


### 5. Configure the Input File

**Edit the benchmark.dat file:**

``` bash
vi benchmark.dat
```

Set the problem size and runtime.  
**Example:**

    1024

This specifies a transform size

### 6. Run the Benchmark

**Execute the benchmark:**

``` bash
export LD_LIBRARY_PATH=$HOME/path/to/fftw-installation/lib:$LD_LIBRARY_PATH
mpirun -np 4 ./program inputvalue/file
```

**Sample**

``` bash
export LD_LIBRARY_PATH=$HOME/path/to/fftw-installation/lib:$LD_LIBRARY_PATH
mpirun -np 4 ./program 1024 || mpirun -np 4 ./program benchmark.dat
```
This runs the FFTW benchmark using 4 MPI processes.

### Results

| ![Sample FFTW Result](/Images/mpi-result.png) | ![Sample FFTW Result](/Images/mpi-end.png) |
|:--:|:--:|

## Conclusion

After running the benchmark, the output will provide detailed performance metrics, including the time taken for the transforms and efficiency data. This information is crucial for evaluating FFTW's performance in your specific applications.

## References

- [FFTW Documentation](https://www.fftw.org/fftw3_doc/Introduction.html)
- [FFTW Benchmark](https://www.fftw.org/fftw3_doc/Benchmarks)
- [Stanford FFTW Installation Documentation](http://micro.stanford.edu/wiki/Install_FFTW3)
- [FFTW Paper](https://fftw.org/pub/fftw/fftw-paper.pdf)
- [Sample Codes](https://github.com/pkestene/simpleFFTW)