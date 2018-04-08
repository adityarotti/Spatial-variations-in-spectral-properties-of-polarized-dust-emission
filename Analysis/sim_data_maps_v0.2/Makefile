# Choose which file you want to compile.
filename=main

execfilename = run-gen-sim
# Path to the moduales
modules=./modules


# PATHS ON MY LAPTOP
include=/Users/adityarotti/Documents/Work/Software/Healpix_3.20/include
healpixlib=/Users/adityarotti/Documents/Work/Software/Healpix_3.20/lib
cfitsiolib=/Users/adityarotti/Documents/Work/Software/cfitsio/lib

######################## GFORTRAN ##########################
FC=gfortran -fbounds-check -ffixed-line-length-none -Wall -Wextra -Wconversion
F77flags=-c -fopenmp
F90flags= -DGFORTRAN -fno-second-underscore -fopenmp
############################################################

all:bips
bips:bipsobj
	$(FC) -O  -I$(include) $(F90flags) -o $(execfilename) *.o -L$(healpixlib) -lhealpix  -lhpxgif -L$(cfitsiolib) -lcfitsio
#	$(FC) -O3 *o $(filename).f90 -o nrun

bipsobj:
	$(FC) -c -I$(include) $(F90flags) $(modules)/global.f90 
	$(FC) -c -I$(include) $(F90flags) $(modules)/gen_noise_maps.f90 
	$(FC) -c -I$(include) $(F90flags) $(modules)/gen_data_maps.f90 
	$(FC) -O3 -I$(include) $(F90flags) -c $(filename).f90 -o $(filename).o

clean:
	$(RM) *.o *.mod *~

# PATH ON CARVER
#include=/global/u2/t/tarun/CMBSOFT/Healpix_2.12a/include
#healpixlib=/global/u2/t/tarun/CMBSOFT/Healpix_2.12a/lib
#cfitsiolib=/global/u2/t/tarun/CMBSOFT/cfitsio/lib


# PATHS ON CRAY
#include=/home/aditya/CMBSOFT/Healpix_2.12a/include
#healpixlib=/home/aditya/CMBSOFT/Healpix_2.12a/lib
#cfitsiolib=/home/aditya/CMBSOFT/cfitsio/lib


## PATHS ON ATHENS
#include=/home/aditya/software/cmb//Healpix_3.20/include
#healpixlib=/home/aditya/software/cmb/Healpix_3.20/lib
#cfitsiolib=/home/aditya/software/cmb/cfitsio/lib



## PATHS ON ASTRO12
#include=/home/aditya/software/cmb/Healpix_3.20/include
#healpixlib=/home/aditya/software/cmb/Healpix_3.20/lib
#cfitsiolib=/home/aditya/software/cmb/cfitsio/lib64

# PATHS ON PERSEUS
#include=/data1/student/caditya/CMBSOFT/Healpix_3.11/include
#healpixlib=/data1/student/caditya/CMBSOFT/Healpix_3.11/lib
#cfitsiolib=/data1/student/caditya/CMBSOFT/cfitsio/lib

## PATH ON HOPPER
#include=/global/homes/a/adityar/Hoppersoft/Healpix_2.12a/include
#healpixlib=/global/homes/a/adityar/Hoppersoft/Healpix_2.12a/lib
#cfitsiolib=/global/homes/a/adityar/Hoppersoft/cfitsio/lib


