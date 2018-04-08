module shared_data
use global_param
use healpix_types
use healpix_modules
use process_mask

implicit none
real(dp) :: multipoles(0:3)
real(dp), allocatable, dimension(:,:) :: pwc,bl
real(dp), allocatable, dimension(:,:) :: mapin,mapout,residual
complex(dpc),allocatable,dimension(:,:,:) :: mapalm,mapinalm

contains

!########################################################################
subroutine allocate_data()

implicit none

allocate(bl(0:maplmax,1:3),pwc(0:maplmax,1:3))
allocate(mapin(0:npixtot-1,1:whichmap),mapout(0:npixtot-1,1:whichmap),residual(0:npixtot-1,1:2))
allocate(mapalm(1:whichmap,0:maplmax,0:maplmax),mapinalm(1:whichmap,0:maplmax,0:maplmax))

end subroutine allocate_data

!-------------------------------------------------------------------------

subroutine deallocate_data()

implicit none

deallocate(bl, pwc)
deallocate(mapin,mapout,residual)
deallocate(mapalm,mapinalm)

end subroutine deallocate_data
!########################################################################

!########################################################################
subroutine return_corrections()

implicit none
character(LEN=80)  :: blheader(1:80)

! PIXEL WINDOW CORRECTION
pwc=1.d0
if (swPWC) call pixel_window(pwc,nside)  ! Whether it returns the pwc for polarization depends on the dimensionality of the array passed to the subroutine.

! INSTRUMENT BEAM SMOOTHING
bl= 1.d0
if (swBEAM) call fits2cl(beamfile,bl,maplmax,whichcl,blheader)

end subroutine return_corrections
!########################################################################

!########################################################################
subroutine return_data(filename)

implicit none
character(LEN=200) :: filename

call read_bintab(filename,mapin,npixtot,whichmap,nullval,anynull)
!call remove_dipole(nside,mapin(:,1),1,2,multipoles,zbounds,fmissval,apomask(:,1))
!call remove_dipole(nside,mapin(:,2),1,2,multipoles,zbounds,fmissval,apomask(:,1))
!call remove_dipole(nside,mapin(:,3),1,2,multipoles,zbounds,fmissval,apomask(:,1))

! If the data has a pixel window, maybe you want to correct it here before passing it through the rest of the code.
!dw8=1.d0
!call map2alm(nside,2*nside,2*nside,mapin,mapinalm,zbounds,dw8)
!call map2alm_iterative(nside,2*nside,2*nside,3,mapin,mapinalm)

end subroutine return_data
!########################################################################

!########################################################################
subroutine write_data(filename)

implicit none
character(len=200) :: filename

call write_bintab(mapout,npixtot,3,mapheader,nlheader,"!"//filename)

end subroutine write_data
!########################################################################

end module shared_data
