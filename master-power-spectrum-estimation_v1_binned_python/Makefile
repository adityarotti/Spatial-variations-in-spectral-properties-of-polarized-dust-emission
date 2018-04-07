# Choose which file you want to compile.
filename=main

exexfilename = run-master-ps
# Path to the moduales
modules=./modules
subroutines=./subroutines

LDFLAGS	= -llapack -lblas -lm -lc

# PATHS ON MY LAPTOP
include=/Users/adityarotti/Documents/Work/Software/Healpix_3.20/include
healpixlib=/Users/adityarotti/Documents/Work/Software/Healpix_3.20/lib
cfitsiolib=/Users/adityarotti/Documents/Work/Software/cfitsio/lib

######################## GFORTRAN ##########################
FC=gfortran -fbounds-check -ffixed-line-length-none #-Wall -Wextra -Wconversion
F77flags=-c -fopenmp
F90flags= -DGFORTRAN -fno-second-underscore -fopenmp
############################################################

all:bips
bips:bipsobj
	#f2py -h bmaster.pyf -m master  $(modules)/master.f90
	f2py -c --fcompiler=gfortran -I$(include) -L$(healpixlib) -lhealpix  -lhpxgif -L$(cfitsiolib) -lcfitsio $(LDFLAGS) *.o -m bmaster $(modules)/master.f90
#	$(FC) -O3 *o $(filename).f90 -o nrun

bipsobj:
	$(FC) -c $(F90flags) $(subroutines)/wigner/*.f 	
	$(FC) -c -I$(include) $(F90flags) $(modules)/global.f90 
	$(FC) -c -I$(include) $(F90flags) $(modules)/data_io.f90 
	#$(FC) -O3 -I$(include) $(F90flags) $(LDFLAGS) -c $(filename).f90 -o $(filename).o

clean:
	$(RM) -r *.o *.mod *~ *dSYM

cleanall:
	$(RM) *.o *.mod *.so *.pyf *~
