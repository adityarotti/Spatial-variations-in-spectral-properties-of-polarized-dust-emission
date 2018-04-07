# Choose which file you want to compile.
filename=master

exexfilename = master
subroutines=./subroutines

LDFLAGS	= -llapack -lblas -lm -lc

######################## GFORTRAN ##########################
FC=gfortran -fbounds-check -ffixed-line-length-none #-Wall -Wextra -Wconversion
F77flags=-c -fopenmp
F90flags= -DGFORTRAN #-fno-second-underscore -fopenmp
############################################################

#all:bips
bips:bipsobj
	f2py -h master.pyf -m master master.f90
	f2py -c --fcompiler=gfortran $(F90flags) $(LDFLAGS) *.o -m master master.f90
#	$(FC) -O3 *o $(filename).f90 -o nrun

bipsobj:
	$(FC) -c $(F90flags) $(subroutines)/wigner/*.f 	
	#$(FC) -O3 -I. $(F90flags) $(LDFLAGS) -c $(filename).f90 -o $(filename).o

clean:
	$(RM) *.o *.mod *~
