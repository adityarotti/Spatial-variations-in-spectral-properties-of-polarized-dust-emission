module global_param
use healpix_modules
use healpix_types

implicit none
integer(i4b) :: i, nstart, nrlz
integer(i4b) :: lmax, masklmax
integer(i4b) :: nside, npixtot
logical :: swSUFFIX=.false.
character(len=200) :: maskfile, data_prefix, data_suffix, datafile 
character(len=200) :: paramfile, fileout, pathout
character :: datatype*3,simnum*4
contains

!########################################################################
subroutine read_param(paramfile)

implicit none
character(LEN=200) :: paramfile

open(10,file=paramfile)
read(10,*) ! Reading in parameter common to simulations and observed maps.
read(10,*) datatype
read(10,*) lmax,masklmax,nside
read(10,*) maskfile
read(10,*) nstart,nrlz
read(10,*) data_prefix
read(10,*) data_suffix
read(10,*) pathout
close(10)

if (nrlz.gt.1) swSUFFIX=.true.
npixtot=nside2npix(nside)

if (lmax.gt.2*nside) then
   print*, "The resolution of the data is not compatible with the lmax&
& specified"
   stop
elseif (masklmax.gt.2*nside) then
   print*, "The resolution of the data is not compatible with the&
& masklmax specified"
   stop
endif

end subroutine read_param
!########################################################################


end module global_param

