module map_operations
use healpix_types
use global_param
use shared_data
use process_mask

implicit none

contains

!########################################################################
subroutine convert_TQU2TEB()
! This routine uses map2alm_iterative

implicit none
real(dp) :: prefactor,ell
real(dp), dimension(0:npixtot-1,1:3) :: tempmap
complex(dpc),dimension(1,0:maplmax,0:maplmax) :: tempalm

if (.not.swMASK) then
  call map2alm_iterative(nside,2*nside,2*nside,maxiter,mapin,mapalm)
  !call map2alm(nside,2*nside,2*nside,mapin,mapalm,zbounds,dw8)
elseif(swMASK) then
  tempmap(:,1)=mapin(:,1)*apomask(:,1)
  tempmap(:,2)=mapin(:,2)*apomask(:,1)
  tempmap(:,3)=mapin(:,3)*apomask(:,1)
  call map2alm_iterative(nside,2*nside,2*nside,maxiter,tempmap,mapalm)
  !call map2alm(nside,2*nside,2*nside,tempmap,mapalm,zbounds,dw8)
  mapinalm=mapalm
endif

do i=0,maplmax
  ell=float(i)
  prefactor=0.d0
  if (i.gt.1) prefactor=sqrt((ell+2.d0)*(ell+1.d0)*ell*(ell-1.d0))
  do j=0,i
    mapalm(2,i,j)=mapalm(2,i,j)*prefactor
    mapalm(3,i,j)=mapalm(3,i,j)*prefactor 
  enddo
enddo

tempalm(1,:,:)=mapalm(1,:,:)
call alm2map(nside,2*nside,2*nside,tempalm,mapout(:,1))

tempalm(1,:,:)=mapalm(2,:,:)
call alm2map(nside,2*nside,2*nside,tempalm,mapout(:,2))

tempalm(1,:,:)=mapalm(3,:,:)
call alm2map(nside,2*nside,2*nside,tempalm,mapout(:,3))

end subroutine convert_TQU2TEB
!#######################################################################

!########################################################################
subroutine convert_TQU2TEB_tilde()
! Only to be used with full sky reconstruction

implicit none
real(dp) :: prefactor,ell
real(dp), dimension(0:npixtot-1,1:3) :: tempmap
complex(dpc),dimension(1,0:maplmax,0:maplmax) :: tempalm

if (.not.swMASK) then
  call map2alm_iterative(nside,2*nside,2*nside,maxiter,mapin,mapalm)
  !call map2alm(nside,2*nside,2*nside,mapin,mapalm,zbounds,dw8)
elseif(swMASK) then
  tempmap(:,1)=mapin(:,1)*apomask(:,1)
  tempmap(:,2)=mapin(:,2)*apomask(:,1)
  tempmap(:,3)=mapin(:,3)*apomask(:,1)
  call map2alm_iterative(nside,2*nside,2*nside,maxiter,tempmap,mapalm)
  !call map2alm(nside,2*nside,2*nside,tempmap,mapalm,zbounds,dw8)
  mapinalm=mapalm
endif

do i=0,maplmax
  ell=float(i)
  prefactor=1.d0
  do j=0,i
    mapalm(2,i,j)=mapalm(2,i,j)*prefactor
    mapalm(3,i,j)=mapalm(3,i,j)*prefactor 
  enddo
enddo

tempalm(1,:,:)=mapalm(1,:,:)
call alm2map(nside,2*nside,2*nside,tempalm,mapout(:,1))

tempalm(1,:,:)=mapalm(2,:,:)
call alm2map(nside,2*nside,2*nside,tempalm,mapout(:,2))

tempalm(1,:,:)=mapalm(3,:,:)
call alm2map(nside,2*nside,2*nside,tempalm,mapout(:,3))

end subroutine convert_TQU2TEB_tilde
!########################################################################

!########################################################################
subroutine calc_residual()

! Be extremely careful while make any changes to this section of the code.
! Its the heart of the machinery.
! This routine has been written assuming the following variables are set :
! 1. > Full sky Q/U maps. (mapin) --> Never taken harmonic transform of here.
! 2. > Apodized mask. (apomask)
! 3. > E/B alms from masked Q/U maps returned from convert_TQU2TEB. (mapinalm)

implicit none
real(dp) :: ell
real(dp), allocatable, dimension(:,:) :: dx,dw,d2w,tempmap
complex(dpc),allocatable,dimension(:,:,:) :: tempalm
!character(len=200) :: filename

allocate(tempalm(2,0:maplmax,0:maplmax))
allocate(dx(0:npixtot-1,1:2),dw(0:npixtot-1,1:2),d2w(0:npixtot-1,1:2),tempmap(0:npixtot-1,1))

tempalm=dcmplx(0.d0,0.d0)
dx=0.d0
dw=0.d0
d2w=0.d0

! Calculating the raising operator on the polarization field.
do i=0,maplmax
  ell=float(i)
  do j=0,i
    tempalm(1,i,j)= -sqrt((ell+2.d0)*(ell-1.d0))*mapinalm(2,i,j)
    tempalm(2,i,j)= -sqrt((ell+2.d0)*(ell-1.d0))*mapinalm(3,i,j)
  enddo
enddo
call alm2map_spin(nside,2*nside,2*nside,1,tempalm,dx)

!!--------------------------------------------------------------------
!! Writing out the spin raised/lowered polarization fields.
!filename=trim(adjustl(pathout))//"dxr.fits"
!tempmap(:,1)=dx(:,1)
!call write_bintab(tempmap,npixtot,1,mapheader,nlheader,"!"//filename)
!
!filename=trim(adjustl(pathout))//"dxi.fits"
!tempmap(:,1)=dx(:,2)
!call write_bintab(tempmap,npixtot,1,mapheader,nlheader,"!"//filename)
!!--------------------------------------------------------------------

! Calculating the raising/lowering operator on the mask.
do i=0,maplmax
  ell=float(i)
  do j=0,i
    tempalm(1,i,j)=-sqrt(ell*(ell+1.d0))*maskalm(1,i,j)
    tempalm(2,i,j)=dcmplx(0.d0,0.d0)
  enddo
enddo
call alm2map_spin(nside,2*nside,2*nside,1,tempalm,dw)

!!--------------------------------------------------------------------
!! Writing out the spin raised/lowered ONCE mask fields.
!filename=trim(adjustl(pathout))//"dwr.fits"
!tempmap(:,1)=dw(:,1)
!call write_bintab(tempmap,npixtot,1,mapheader,nlheader,"!"//filename)
!
!filename=trim(adjustl(pathout))//"dwi.fits"
!tempmap(:,1)=dw(:,2)
!call write_bintab(tempmap,npixtot,1,mapheader,nlheader,"!"//filename)
!!--------------------------------------------------------------------


! Calculating the (raising operator)^2 on the mask.
do i=0,maplmax
  ell=float(i)
  do j=0,i
    tempalm(1,i,j)=-sqrt((ell+2.d0)*(ell+1.d0)*ell*(ell-1.d0))*maskalm(1,i,j)
    tempalm(2,i,j)=dcmplx(0.d0,0.d0)
  enddo
enddo
call alm2map_spin(nside,2*nside,2*nside,2,tempalm,d2w)

d2w(:,1)=d2w(:,1)*apomask(:,1)
d2w(:,2)=d2w(:,2)*apomask(:,1)

!!--------------------------------------------------------------------
!! Writing out the spin raised/lowered TWICE mask fields.
!filename=trim(adjustl(pathout))//"d2wr.fits"
!tempmap(:,1)=d2w(:,1)
!call write_bintab(tempmap,npixtot,1,mapheader,nlheader,"!"//filename)
!
!filename=trim(adjustl(pathout))//"d2wi.fits"
!tempmap(:,1)=d2w(:,2)
!call write_bintab(tempmap,npixtot,1,mapheader,nlheader,"!"//filename)
!!--------------------------------------------------------------------

! E
residual(:,1) = (mapin(:,2)*d2w(:,1) + mapin(:,3)*d2w(:,2))
residual(:,1) = residual(:,1) + 2.d0*(dw(:,1)*dx(:,1) + dw(:,2)*dx(:,2))
residual(:,1) = residual(:,1) - 2.d0*(dw(:,1)*dw(:,1) - dw(:,2)*dw(:,2))*mapin(:,2) - 4.d0*dw(:,1)*dw(:,2)*mapin(:,3)
! B
residual(:,2) = (-mapin(:,2)*d2w(:,2) + mapin(:,3)*d2w(:,1))
residual(:,2) = residual(:,2) + 2.d0*(dw(:,1)*dx(:,2) - dw(:,2)*dx(:,1))
residual(:,2) = residual(:,2) - 2.d0*(dw(:,1)*dw(:,1) - dw(:,2)*dw(:,2))*mapin(:,3) + 4.d0*dw(:,1)*dw(:,2)*mapin(:,2)


!filename=trim(adjustl(pathout))//"E-residual.fits"
!tempmap(:,1)=residual(:,1)
!call write_bintab(tempmap,npixtot,1,mapheader,nlheader,"!"//filename)
!
!filename=trim(adjustl(pathout))//"B-residual.fits"
!tempmap(:,1)=residual(:,2)
!call write_bintab(tempmap,npixtot,1,mapheader,nlheader,"!"//filename)
!--------------------------------------------------------------------

mapout(:,2)=mapout(:,2)*apomask(:,1) + residual(:,1)
mapout(:,3)=mapout(:,3)*apomask(:,1) + residual(:,2)

!filename=trim(adjustl(pathout))//"WE.fits"
!tempmap(:,1)=mapout(:,2)
!call write_bintab(tempmap,npixtot,1,mapheader,nlheader,"!"//filename)
!
!filename=trim(adjustl(pathout))//"WB.fits"
!tempmap(:,1)=mapout(:,3)
!call write_bintab(tempmap,npixtot,1,mapheader,nlheader,"!"//filename)

!filename=trim(adjustl(pathout))//"q.fits"
!tempmap(:,1)=mapin(:,2)
!call write_bintab(tempmap,npixtot,1,mapheader,nlheader,"!"//filename)
!
!filename=trim(adjustl(pathout))//"u.fits"
!tempmap(:,1)=mapin(:,3)
!call write_bintab(tempmap,npixtot,1,mapheader,nlheader,"!"//filename)

deallocate(tempalm,dx,dw,d2w,tempmap)

end subroutine calc_residual
!########################################################################

end module map_operations
