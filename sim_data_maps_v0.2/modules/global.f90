module global_param
use healpix_modules
use healpix_types

implicit none
! Map read write variables
integer(i4b), parameter :: nlheader=80
real(dp) :: zbounds(1:2)=0.d0,nullval,fmissval=-1.6375e30
real(dp), allocatable, dimension(:,:) :: reccov,L,TQU_map,cmb,w8ring
complex(dpc), allocatable, dimension(:,:,:) :: TEB_alm,noise_alm,TEB_frg
real(dp), allocatable, dimension(:,:) :: bl,wl
character(len=80) :: mapheader(1:nlheader)
logical :: anynull
type(planck_rng) :: rng_handle_noise,rng_handle_cmb,rng_handle_frg


integer(i4b) :: nside,npixtot,nstart,nrlz,simnum
integer(i4b) :: seed1init,seed2init
real(dp) :: tnorm ! =1.d6
character(len=200) :: covfname,pathout,paramfile,filename,prefix,clfilename,beamfilename,wlfilename
character(len=200) :: clfrgfilename
character :: s1*4,s2*4,s3*4,whichmap1*1,whichmap2*1
logical :: swPWC,swBEAM,swRECCOV=.False.
logical :: swNOISE=.True.,swCMB=.False.,swFRG=.True.,swTESTNOISE=.False.
logical :: swSIMDATACALIB=.True.

contains

!########################################################################
subroutine read_param(paramfile)

implicit none
character(LEN=200) :: paramfile

open(10,file=paramfile)
! Parameters common to simulations and observed maps.
read(10,*) ! Reading in parameter common to simulations and observed maps.
read(10,*) nside
read(10,*) tnorm
read(10,*) seed1init,seed2init
read(10,*) covfname
read(10,*) clfilename
read(10,*) clfrgfilename
read(10,*) beamfilename
read(10,*) wlfilename
read(10,*) pathout
read(10,*) prefix
read(10,*) nstart,nrlz
close(10)

npixtot=nside2npix(nside) !; print*, npixtot

write(s3,"(i4.1)") nside

call write_minimal_header(mapheader,"MAP",nside=nside,ordering="RING",polar=.True.)
end subroutine read_param
!########################################################################

!########################################################################
subroutine allocate_data()

allocate(L(0:npixtot-1,1:6),TQU_map(0:npixtot-1,1:3),cmb(0:npixtot-1,1:3),w8ring(1:2*nside,1:3))
w8ring=1.d0
allocate(TEB_alm(1:3,0:2*nside,0:2*nside),TEB_frg(1:3,0:2*nside,0:2*nside),noise_alm(1:3,0:2*nside,0:2*nside))
allocate(bl(0:2*nside,3),wl(0:2*nside,3))
wl=1.d0

if (swRECCOV) then
  allocate(reccov(0:npixtot-1,1:6)) ; reccov=0.d0
endif

end subroutine allocate_data
!########################################################################

!########################################################################
subroutine deallocate_data()

deallocate(L,TQU_map,cmb,w8ring)
deallocate(TEB_alm,TEB_frg,noise_alm,bl)

if (swRECCOV) deallocate(reccov)

end subroutine deallocate_data
!########################################################################


end module global_param

