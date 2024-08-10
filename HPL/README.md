# HPL Benchmarking

## Introduction

High-Performance Linpack (HPL) is a benchmark used to measure the
floating-point computing power of a system. It’s a standard benchmark
for evaluating the performance of supercomputers.

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

### 2. Download and Extract HPL

**Download HPL:**

``` bash
wget http://www.netlib.org/benchmark/hpl/hpl-2.3.tar.gz
tar -xvzf hpl-2.3.tar.gz
cd hpl-2.3
```

### 3. Configure HPL

**run the shell script to generate a makefile for your environment:**

``` bash
. setup/make_generic
```

**Copy the Makefile for your environment:**

``` bash
cp setup/Make.UNKNOWN Makefile
```

**Edit the Makefile to set the required libraries and compilers:**

``` bash
vi Makefile
```

Set `MPdir`, `MPinc`, and `MPlib` to the OpenMPI paths.  
Set `LAdir`, `LAinc`, and `LAlib` to the OpenBLAS paths.

``` bash
echo $PATH
```
To find the paths to these libraries/modules

**Example makefile**
``` bash

# ----------------------------------------------------------------------
# - Platform identifier ------------------------------------------------
# ----------------------------------------------------------------------
#
ARCH         = $(arch)
#
# ----------------------------------------------------------------------
# - HPL Directory Structure / HPL library ------------------------------
# ----------------------------------------------------------------------
#
TOPdir       = $(HOME)/hpl-2.3
INCdir       = $(TOPdir)/include
BINdir       = $(TOPdir)/bin/$(ARCH)
LIBdir       = $(TOPdir)/lib/$(ARCH)
#
HPLlib       = $(LIBdir)/libhpl.a 
#
# ----------------------------------------------------------------------
# - Message Passing library (MPI) --------------------------------------
# ----------------------------------------------------------------------
# MPinc tells the  C  compiler where to find the Message Passing library
# header files,  MPlib  is defined  to be the name of  the library to be 
# used. The variable MPdir is only used for defining MPinc and MPlib.
#
MPdir        = /opt/apps/nfs/spack/opt/spack/linux-centos8-zen2/gcc-9.2.0/openmpi-4.0.4-7s5s4ctaquh6we2nzffv7frfpy4qvqyw
MPinc        = -I$(MPdir)/include
MPlib        = $(MPdir)/lib/libmpi.so
#
```

### 4. Compile HPL

**Build HPL:**

``` bash
make arch=Linux
```
### 5. Configure the HPL.dat File

``` bash
cd bin/<arch>
```
**Example:**

``` bash
cd bin/Linux
```

**Edit the HPL.dat file:**

``` bash
vi HPL.dat
```

Set the problem size (`N`) and block size (`NB`).  
**Example:**

    HPLinpack benchmark input file
    Innovative Computing Laboratory, University of Tennessee
    HPL.out      output file name (if any) 
    6            device out (6=stdout,7=stderr,file)
    1            # of problems sizes (N)
    40960        Ns
    1            # of NBs
    192          NBs
    0            PMAP process mapping (0=Row-,1=Column-major)
    1            # of process grids (P x Q)
    2            Ps
    2            Qs
    16.0         threshold

**Recommended problem sizes, block Sizes and HPL.dat can be derived using either of the following websites:**
https://www.advancedclustering.com/act_kb/tune-hpl-dat-file/
https://hpl-calculator.sourceforge.net/

### 6. Run the Benchmark

**Execute the benchmark:**

``` bash
mpirun -np 4 ./xhpl
```

This runs the HPL benchmark using 4 MPI processes.

## Conclusion

After running the benchmark, you’ll see performance results in terms of
GFLOPS (Giga Floating Point Operations per Second). These results can be
used to gauge the performance of your system.
