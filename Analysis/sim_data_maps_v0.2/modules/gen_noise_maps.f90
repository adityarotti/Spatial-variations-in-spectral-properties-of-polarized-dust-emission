module generate_noise_map
use healpix_modules
use healpix_types
use global_param

implicit none
integer(i4b) :: i,j

contains

!########################################################################
subroutine read_TQU_cov_calc_LTM()

implicit none
real(dp), allocatable, dimension(:,:) :: cov

allocate(cov(0:npixtot-1,6))

call read_bintab(covfname,cov,npixtot,6,nullval,anynull)

! Calculate the lower triangular matrix (LTM). 
! We read the TQU covariance and calculate the LTM at each pixel.
! LTM is what is stored.
! L1 = L11, L2=L12, L3=L22, L4=L13, L5=L23, L6=L33
do i = 0,npixtot-1
   L(i,1) = sqrt(cov(i,1)) 
   L(i,2) = cov(i,2)/L(i,1)
   L(i,3) = sqrt(cov(i,4)-(L(i,2)**2.d0)) 
   L(i,4) = cov(i,3)/L(i,1)
   L(i,5) = (cov(i,5)  - L(i,4)*L(i,2))/L(i,3)
   L(i,6) = sqrt(cov(i,6)- (L(i,4)**2.d0) - (L(i,5)**2.d0))
   
   do j=1,6
      if (L(i,j).ne.L(i,j)) then
         print*, i,j,cov(i,j)
         L(i,j)=0.d0
      endif
   enddo
enddo

deallocate(cov)

Print*, "Populated the lower triangular matrix"

end subroutine read_TQU_cov_calc_LTM
!########################################################################

!!########################################################################
!subroutine read_TQU_cov_calc_LTM()
!
!implicit none
!real(dp), allocatable, dimension(:,:) :: cov
!real(dp), allocatable, dimension(:,:) :: TT, QQ, UU, TQ, TU, QU
!
!allocate(TT(0:npixtot-1,1),QQ(0:npixtot-1,1),UU(0:npixtot-1,1))
!allocate(TQ(0:npixtot-1,1),TU(0:npixtot-1,1),QU(0:npixtot-1,1))
!allocate(cov(0:npixtot-1,6))
!
!call read_bintab(covfname,cov,npixtot,6,nullval,anynull)
!TT(:,1)=cov(:,1) ; call write_bintab(TT,npixtot,1,mapheader,nlheader,"!"//"./tempdata/tt.fits")
!TQ(:,1)=cov(:,2) ; call write_bintab(TQ,npixtot,1,mapheader,nlheader,"!"//"./tempdata/tq.fits")
!TU(:,1)=cov(:,3) ; call write_bintab(TU,npixtot,1,mapheader,nlheader,"!"//"./tempdata/tu.fits")
!QQ(:,1)=cov(:,4) ; call write_bintab(QQ,npixtot,1,mapheader,nlheader,"!"//"./tempdata/qq.fits")
!QU(:,1)=cov(:,5) ; call write_bintab(QU,npixtot,1,mapheader,nlheader,"!"//"./tempdata/qu.fits")
!UU(:,1)=cov(:,6) ; call write_bintab(UU,npixtot,1,mapheader,nlheader,"!"//"./tempdata/uu.fits")
!
!! Calculate the lower triangular matrix (LTM). 
!! We read the TQU covariance and calculate the LTM at each pixel.
!! LTM is what is stored.
!! L1 = L11, L2=L12, L3=L22, L4=L13, L5=L23, L6=L33
!do i = 0,npixtot-1
!   L(i,1) = sqrt(TT(i,1))
!   L(i,2) = TQ(i,1)/L(i,1)
!   L(i,3) = sqrt(QQ(i,1)-(L(i,2)**2.d0)) 
!   L(i,4) = TU(i,1)/L(i,1)
!   L(i,5) = (QU(i,1)  - L(i,4)*L(i,2))/L(i,3)
!   L(i,6) = sqrt(UU(i,1)- (L(i,4)**2.d0) - (L(i,5)**2.d0))
!enddo
!
!!do i = 0,npixtot-1
!!   L(i,1) = sqrt(cov(i,1))
!!   L(i,2) = cov(i,2)/L(i,1)
!!   L(i,3) = sqrt(cov(i,4)-(L(i,2)**2.d0)) 
!!   L(i,4) = cov(i,3)/L(i,1)
!!   L(i,5) = (cov(i,5)  - L(i,4)*L(i,2))/L(i,3)
!!   L(i,6) = sqrt(cov(i,6)- (L(i,4)**2.d0) - (L(i,5)**2.d0))
!!enddo
!
!
!deallocate(TT,QQ,UU,TQ,TU,QU)
!deallocate(cov)
!
!Print*, "Populated the lower triangular matrix"
!
!end subroutine read_TQU_cov_calc_LTM
!!########################################################################

!########################################################################
subroutine return_noise_map()
use rngmod

implicit none
real(dp) :: grn3(1:3),noise(1:3)
character(len=200) :: fname

! Generating the random map
do i = 0,npixtot-1
   grn3(1) = rand_gauss(rng_handle_noise)
   grn3(2) = rand_gauss(rng_handle_noise)
   grn3(3) = rand_gauss(rng_handle_noise)
!   print*, simnum,grn3(1), grn3(2), grn3(3)

   call return_Lx(i,grn3,noise)

!  Sanity check.
!   if ((noise(1).ne.noise(1)).or.(noise(2).ne.noise(2)).or.(noise(3).ne.noise(3))) then
!      print*, i,(L(i,j),j=1,6)
!      pause
!   endif
 
   TQU_map(i,:)=noise
!   print*, simnum,noise(1), noise(2), noise(3)
!   print*, simnum,TQU_map(i,1), TQU_map(i,2), TQU_map(i,3)
   if (swRECCOV) call reconstruct_covariance(i)
enddo

if (swTESTNOISE) then
   fname=trim(adjustl(pathout))//"noise_nrlz"//trim(adjustl(s1))//".fits"
   call map2alm(nside,2*nside,2*nside,TQU_map,noise_alm,zbounds,w8ring) !; print*, noise_alm
   call alm2map(nside,2*nside,2*nside,noise_alm,TQU_map)
   call write_bintab(TQU_map,npixtot,3,mapheader,nlheader,"!"//trim(adjustl(fname)))
endif

end subroutine return_noise_map
!########################################################################

!########################################################################
subroutine return_Lx(pixnum,ranvecin,ranvecout)

implicit none
integer(i4b) :: pixnum
real(dp), dimension(1:3) :: ranvecin, ranvecout

ranvecout(1) = L(pixnum,1)*ranvecin(1)
ranvecout(2) = L(pixnum,2)*ranvecin(1) + L(pixnum,3)*ranvecin(2)
ranvecout(3) = L(pixnum,4)*ranvecin(1) + L(pixnum,5)*ranvecin(2) + L(pixnum,6)*ranvecin(3)

end subroutine return_Lx
!########################################################################


!########################################################################
subroutine reconstruct_covariance(pixnum)

implicit none
integer(i4b)::pixnum

reccov(pixnum,1)=reccov(pixnum,1)+(TQU_map(pixnum,1)*TQU_map(pixnum,1)) ! TT
reccov(pixnum,2)=reccov(pixnum,2)+(TQU_map(pixnum,1)*TQU_map(pixnum,2)) ! TQ 
reccov(pixnum,3)=reccov(pixnum,3)+(TQU_map(pixnum,1)*TQU_map(pixnum,3)) ! TU
reccov(pixnum,4)=reccov(pixnum,4)+(TQU_map(pixnum,2)*TQU_map(pixnum,2)) ! QQ
reccov(pixnum,5)=reccov(pixnum,5)+(TQU_map(pixnum,2)*TQU_map(pixnum,3)) ! QU
reccov(pixnum,6)=reccov(pixnum,6)+(TQU_map(pixnum,3)*TQU_map(pixnum,3)) ! UU

end subroutine reconstruct_covariance
!########################################################################

!########################################################################
subroutine write_reconstructed_cov()

implicit none
real(dp), allocatable, dimension(:,:) :: TT, QQ, UU, TQ, TU, QU

allocate(TT(0:npixtot-1,1),QQ(0:npixtot-1,1),UU(0:npixtot-1,1))
allocate(TQ(0:npixtot-1,1),TU(0:npixtot-1,1),QU(0:npixtot-1,1))

TT(:,1)=reccov(:,1) ; call write_bintab(TT,npixtot,1,mapheader,nlheader,"!"//"./tempdata/tt.fits")
TQ(:,1)=reccov(:,2) ; call write_bintab(TQ,npixtot,1,mapheader,nlheader,"!"//"./tempdata/tq.fits")
TU(:,1)=reccov(:,3) ; call write_bintab(TU,npixtot,1,mapheader,nlheader,"!"//"./tempdata/tu.fits")
QQ(:,1)=reccov(:,4) ; call write_bintab(QQ,npixtot,1,mapheader,nlheader,"!"//"./tempdata/qq.fits")
QU(:,1)=reccov(:,5) ; call write_bintab(QU,npixtot,1,mapheader,nlheader,"!"//"./tempdata/qu.fits")
UU(:,1)=reccov(:,6) ; call write_bintab(UU,npixtot,1,mapheader,nlheader,"!"//"./tempdata/uu.fits")

deallocate(TT,QQ,UU,TQ,TU,QU)

end subroutine write_reconstructed_cov
!########################################################################



end module generate_noise_map
