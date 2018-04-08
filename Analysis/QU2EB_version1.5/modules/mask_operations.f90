module process_mask
use global_param
use healpix_types
use healpix_modules

implicit none

real(dp), allocatable, dimension(:,:) :: apomask
complex(dpc),allocatable,dimension(:,:,:) :: maskalm

contains

!########################################################################
subroutine allocate_mask_arrays()

implicit none

allocate(apomask(0:npixtot-1,1),maskalm(1,0:maplmax,0:maplmax))

apomask=0.d0
maskalm=dcmplx(0.d0,0.d0)

end subroutine allocate_mask_arrays

!-------------------------------------------------------------------------

subroutine deallocate_mask_arrays()

deallocate(apomask,maskalm)

end subroutine deallocate_mask_arrays
!########################################################################


!########################################################################
subroutine read_data_mask()
implicit none

call read_bintab(maskfile,apomask,npixtot,1,nullval,anynull)
!call convert_nest2ring(nside,mask(:,1))
call map2alm_iterative(nside,2*nside,2*nside,5,apomask,maskalm)

!dw8=1.d0
!call map2alm(nside,2*nside,2*nside,apomask(:,1),maskalm,zbounds,dw8)

end subroutine read_data_mask
!########################################################################

end module process_mask
