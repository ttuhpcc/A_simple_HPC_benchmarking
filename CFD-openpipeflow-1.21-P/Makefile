INSTDIR		= ./install/
PROGDIR		= ./program/
UTILDIR		= ./utils/
UTIL		= energyspectrum
#UTIL		= statfluct
#UTIL		= prim2matlab

TRANSFORM	= fftw3
MODSOBJ		= io.o meshs.o mpi.o parameters.o \
		  timestep.o transform.o variables.o velocity.o

#COMPILER	= g95 -C  
#COMPFLAGS	= -cpp -c -O3
#COMPILER	= ifort -i-dynamic #-C #-static 
#COMPFLAGS	= -cpp -c -O3 -heap-arrays 1024 -mcmodel=medium
#COMPILER	= pgf90 #-C
#COMPFLAGS	= -Mpreprocess -c -fast #-mcmodel=medium
#COMPILER	= pathf90 #-pg -C
#COMPFLAGS	= -cpp -c -O3 -OPT:Ofast -march=opteron -fno-second-underscore


COMPILER	= mpif90
COMPFLAGS	= -cpp -c -O3 -ffree-line-length-none -mcmodel=medium\
		  -I/usr/include \
		-I/home/jieyao/CODE_TEST/CHANNEL/fftw-3.3.7-build/include\
		-I/opt/apps/nfs/spack/opt/spack/linux-centos8-zen2/gcc-10.1.0/openblas-0.3.10-x5doqjyvllsykuqxm7j3xzh6gsdex3ws/include\
		-I/opt/apps/nfs/spack/opt/spack/linux-centos8-zen2/gcc-10.1.0/netlib-lapack-3.8.0-sfyc3hooop3fxr3rjasueqqkpirjh437/include\
		-I/opt/apps/nfs/spack/opt/spack/linux-centos8-zen2/gcc-10.1.0/netcdf-c-4.7.3-ljb3z4dr2m6fke4ck6pkitprbpfu6esf//include\
		-I/opt/apps/nfs/spack/opt/spack/linux-centos8-zen2/gcc-10.1.0/netcdf-fortran-4.5.2-y6v6qtrpb2ynw6e2bemgkvak6uflggpp//include\
		#-I/opt/apps/gcc9_1/impi19_0/parallel-netcdf/4.6.2/x86_64/include/\
		  -I/opt/apps/gcc9_1/impi19_0/fftw3/3.3.8/include\
 		 #-I/home1/02890/yaojie/OPENPIPE/netcdf-fortran-4.4.4-build/include\
		  -I/opt/apps/gcc9_1/impi19_0/fftw3/3.3.8/include\
		#-I/home1/02890/yaojie/OPENPIPE/netcdf-4.6.1-build/include\
                  #-C #-pg
LIBS		= -L/home/jieyao/CODE_TEST/CHANNEL/fftw-3.3.7-build/lib\
		  -L/opt/apps/nfs/spack/opt/spack/linux-centos8-zen2/gcc-10.1.0/openblas-0.3.10-x5doqjyvllsykuqxm7j3xzh6gsdex3ws/lib\
		  -L/opt/apps/nfs/spack/opt/spack/linux-centos8-zen2/gcc-10.1.0/netcdf-c-4.7.3-ljb3z4dr2m6fke4ck6pkitprbpfu6esf/lib\
		-L/opt/apps/nfs/spack/opt/spack/linux-centos8-zen2/gcc-10.1.0/netcdf-fortran-4.5.2-y6v6qtrpb2ynw6e2bemgkvak6uflggpp//lib\
		-L/opt/apps/nfs/spack/opt/spack/linux-centos8-zen2/gcc-10.1.0/netlib-lapack-3.8.0-sfyc3hooop3fxr3rjasueqqkpirjh437/lib64\
		  cheby.o -lfftw3  -llapack -lnetcdff -lnetcdf 
		  # -lblas -lcurl
		#-L/home1/02890/yaojie/OPENPIPE/netcdf-4.6.1-build/lib\
		#-L/home1/02890/yaojie/OPENPIPE/netcdf-fortran-4.4.4-build/lib\
		-L/opt/apps/gcc9_1/impi19_0/parallel-netcdf/4.6.2/x86_64/lib\
                  -L/opt/apps/gcc9_1/impi19_0/fftw3/3.3.8/lib\

#------------------------------------------------------------------------
all : 	$(MODSOBJ) $(PROGDIR)main.f90
	$(COMPILER) $(COMPFLAGS) $(PROGDIR)main.f90
	$(COMPILER) -o ./main.out main.o $(MODSOBJ) $(LIBS)

install : main.out
	if test ! -d $(INSTDIR); then mkdir -p $(INSTDIR); fi
	mv ./main.out $(INSTDIR)
	date > $(INSTDIR)/main.info
	echo $(HOSTNAME) >> $(INSTDIR)/main.info
	pwd >> $(INSTDIR)/main.info
	echo $(COMPILER) $(COMPFLAGS) >> $(INSTDIR)/main.info
	grep "define _N" parallel.h >> $(INSTDIR)/main.info
	cut -d! -f1 $(PROGDIR)parameters.f90 | grep = | \
	   cut -d: -f3  >> $(INSTDIR)/main.info

util : 	$(MODSOBJ) $(UTILDIR)/$(UTIL).f90
	$(COMPILER) $(COMPFLAGS) $(UTILDIR)/$(UTIL).f90
	$(COMPILER) -o ./$(UTIL).out $(UTIL).o $(MODSOBJ) $(LIBS)

#------------------------------------------------------------------------
clean :
	rm -f *.o *.mod *.d *.il core *.out -r install

#------------------------------------------------------------------------
io.o : $(PROGDIR)io.f90 velocity.o 
	$(COMPILER) $(COMPFLAGS) $(PROGDIR)io.f90

meshs.o : $(PROGDIR)meshs.f90 parameters.o mpi.o
	$(COMPILER) $(COMPFLAGS) $(PROGDIR)cheby.f
	$(COMPILER) $(COMPFLAGS) $(PROGDIR)meshs.f90

mpi.o : $(PROGDIR)mpi.f90 parallel.h
	$(COMPILER) $(COMPFLAGS) $(PROGDIR)mpi.f90

parameters.o : $(PROGDIR)parameters.f90 parallel.h
	$(COMPILER) $(COMPFLAGS) $(PROGDIR)parameters.f90

timestep.o : $(PROGDIR)timestep.f90 variables.o
	$(COMPILER) $(COMPFLAGS) $(PROGDIR)timestep.f90

transform.o : $(PROGDIR)transform.$(TRANSFORM).f90 variables.o
	$(COMPILER) $(COMPFLAGS) -o transform.o \
	$(PROGDIR)transform.$(TRANSFORM).f90

variables.o : $(PROGDIR)variables.f90 meshs.o
	$(COMPILER) $(COMPFLAGS) $(PROGDIR)variables.f90

velocity.o : $(PROGDIR)velocity.f90 timestep.o transform.o
	$(COMPILER) $(COMPFLAGS) $(PROGDIR)velocity.f90

