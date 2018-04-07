module simulate_data_map
use healpix_modules
use healpix_types
use global_param

implicit none

contains
!########################################################################
subroutine return_eff_bl()
implicit none
character(LEN=80)  :: blheader(1:80)
real(dp) :: pwc(0:2*nside,1:3)
integer(i4b) :: i,j

call pixel_window(pwc,nside)
call fits2cl(beamfilename,bl,2*nside,3,blheader)

if (swSIMDATACALIB) call fits2cl(wlfilename,wl,2*nside,3,blheader)

do i = 0,2*nside
  do j=1,3
    bl(i,j)=bl(i,j)*pwc(i,j)
  enddo
enddo
!bl=1.d0

end subroutine return_eff_bl
!########################################################################

!########################################################################
subroutine gen_cmb_frg_add_noise_smooth_map()
implicit none
character(len=80) :: header(1:60)
logical :: filexist

   
write(s1,"(I4.1)")simnum


!filename=trim(adjustl(pathout))//trim(adjustl(prefix(:3)))//"_f"//&
!&trim(adjustl(prefix(4:7)))//"_nside"//trim(adjustl(s3))//"_nrlz"//trim(adjustl(s1))//".fits"
!inquire(file=filename,EXIST=filexist)


! Generating alms for CMB map.
TEB_alm=cmplx(0.d0,0.d0)
if (swCMB) call create_alm(nside,2*nside,2*nside,1,clfilename,rng_handle_cmb,0.d0,TEB_alm,header)


! Generating alms for FRG map.
TEB_frg=cmplx(0.d0,0.d0)
if (swFRG) call create_alm(nside,2*nside,2*nside,1,clfrgfilename,rng_handle_frg,0.d0,TEB_frg,header)
TEB_frg(1,:,:)=cmplx(0.d0,0.d0)

! Converting noise map to alms
noise_alm=cmplx(0.d0,0.d0)
if (swNOISE) call map2alm(nside,2*nside,2*nside,TQU_map,noise_alm,zbounds,w8ring)
! This rescales the harmonic coefficients of noise to match with data.
call alter_alm(nside,2*nside,2*nside,0.d0,noise_alm,"none",wl)

! --------------------------------------------------------------------
! CMB + NOISE. (C+N)
filename=trim(adjustl(pathout))//trim(adjustl(prefix(:3)))//"_cn"//&
&trim(adjustl(prefix(4:)))//"_nside"//trim(adjustl(s3))//"_nrlz"//trim(adjustl(s1))//".fits"

TEB_alm=TEB_alm + noise_alm
call alter_alm(nside,2*nside,2*nside,0.d0,TEB_alm,"none",bl)

call alm2map(nside,2*nside,2*nside,TEB_alm,TQU_map)
call write_bintab(TQU_map,npixtot,3,mapheader,nlheader,"!"//trim(adjustl(filename)))
! --------------------------------------------------------------------

! --------------------------------------------------------------------
! Foregrounds. (F)
filename=trim(adjustl(pathout))//trim(adjustl(prefix(:3)))//"_f"//&
&trim(adjustl(prefix(4:7)))//"_nside"//trim(adjustl(s3))//"_nrlz"//trim(adjustl(s1))//".fits"
call alter_alm(nside,2*nside,2*nside,0.d0,TEB_frg,"none",bl)
TQU_map=0.d0
call alm2map(nside,2*nside,2*nside,TEB_frg,TQU_map)
call write_bintab(TQU_map,npixtot,3,mapheader,nlheader,"!"//trim(adjustl(filename)))
! --------------------------------------------------------------------

! --------------------------------------------------------------------
! (CMB+NOISE) + Foregrounds. (C+N+F)
filename=trim(adjustl(pathout))//trim(adjustl(prefix(:3)))//"_cnf"//&
&trim(adjustl(prefix(4:)))//"_nside"//trim(adjustl(s3))//"_nrlz"//trim(adjustl(s1))//".fits"
TEB_alm=TEB_alm+TEB_frg
TQU_map=0.d0
call alm2map(nside,2*nside,2*nside,TEB_alm,TQU_map)
call write_bintab(TQU_map,npixtot,3,mapheader,nlheader,"!"//trim(adjustl(filename)))
! --------------------------------------------------------------------

end subroutine gen_cmb_frg_add_noise_smooth_map
!########################################################################

end module simulate_data_map
