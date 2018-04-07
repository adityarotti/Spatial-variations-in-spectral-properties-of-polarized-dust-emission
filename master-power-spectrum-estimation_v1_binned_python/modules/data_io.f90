module data_io
use healpix_modules
use healpix_types
use global_param


implicit none
real(dp), allocatable,dimension(:,:) :: map
real(dp), allocatable, dimension(:,:) :: cl,wl,clm ! ,dw8
real(dp), allocatable, dimension(:,:) :: mllp
integer(i4b), allocatable, dimension(:) :: ipiv
integer(i4b) :: info
complex(dpc), allocatable, dimension(:,:,:) :: mapalm
real(dp) :: nullval,zbounds(1:2)=0.d0,multipoles(0:3),fmissval=-1.6375e30
logical :: anynull
character(LEN=80)  :: clheader(1:80)

contains
!########################################################################
subroutine allocate_io_data()

allocate(map(0:npixtot-1,1),mapalm(1,0:lmax,0:lmax))
allocate(cl(0:lmax,1),clm(0:lmax,1),wl(0:masklmax,1))!,dw8(1:lmax,1))
allocate(mllp(0:lmax,0:lmax),ipiv(lmax+1))
! dw8=1.d0

end subroutine allocate_io_data
!########################################################################

!########################################################################
subroutine deallocate_io_data()

deallocate(map,mapalm)
deallocate(cl,clm,wl,mllp,ipiv)!,dw8)

end subroutine deallocate_io_data
!########################################################################


!########################################################################
subroutine return_mask_cl()

implicit none
character(LEN=80)  :: clheader(1:80)
complex(dpc), dimension(1,0:masklmax,0:masklmax) :: maskalm
!real(dp),dimension(1:masklmax,1) :: dw8


if (datatype.eq."MAP") then
   call read_bintab(maskfile,map,npixtot,1,nullval,anynull)
!  call map2alm(nside,masklmax,masklmax,map(:,1),maskalm,zbounds,dw8)
   call map2alm_iterative(nside,masklmax,masklmax,3,map,maskalm)
   call alm2cl(masklmax,masklmax,maskalm,wl)
else
   call fits2cl(maskfile,wl,masklmax,1,clheader)
endif

end subroutine return_mask_cl
!########################################################################

!########################################################################
subroutine return_data_cl()

implicit none
character(LEN=80)  :: clheader(1:80)

if (datatype.eq."MAP") then
   call read_bintab(datafile,map,npixtot,1,nullval,anynull)
!  call map2alm(nside,lmax,lmax,map(:,1),mapalm,zbounds,dw8)
   call map2alm_iterative(nside,lmax,lmax,3,map,mapalm)
   call alm2cl(lmax,lmax,mapalm,clm)
else
   call fits2cl(datafile,clm,lmax,1,clheader)
endif

cl=clm

end subroutine return_data_cl
!########################################################################



end module data_io
