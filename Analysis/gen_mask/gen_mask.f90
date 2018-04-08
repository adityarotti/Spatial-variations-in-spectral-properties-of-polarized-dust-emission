module gen_mask
use healpix_types
use healpix_modules

implicit none

contains

!########################################################################
subroutine return_disc_mask(discsize,apow,nside,nsideout,pixel_index,npixtot,mask)!,whichapo)

implicit none
integer*4, intent(in) :: nside,nsideout,pixel_index,npixtot
integer*4 :: i !,whichapo (1 : Cosine apodization, 2 : Cosine squared apodization)
real*8 :: theta,phi,theta0,phi0,dtheta,dotprod
real*8,intent(in) :: discsize,apow ! degrees
real*8,intent(out),dimension(npixtot) :: mask

call pix2ang_ring(nsideout,pixel_index,theta0,phi0)
do i=0,npixtot-1
  call pix2ang_ring(nside,i,theta,phi)
  dotprod=sin(theta)*sin(theta0)*cos(phi-phi0)+cos(theta)*cos(theta0)
  dtheta=acos(dotprod)*180.d0/pi
  

  if (dtheta.le.discsize-apow) then
    mask(i)=1.d0
!  elseif ((dtheta.ge.cut).and.(dtheta.le.(cut+apobeamwd))) then
  elseif ((dtheta.gt.(discsize-apow)).and.(dtheta.le.discsize)) then
    dtheta=dtheta-(discsize-apow)
    !if (whichapo.eq.1) mask(i)=(1.d0+cos(dtheta*pi/apow))/2.d0 ! Cosine apodization
    !if (whichapo.eq.2) mask(i)=cos(dtheta*pi*0.5d0/apow)**2.d0 ! Cosine^2 apodization
    mask(i)=cos(dtheta*pi*0.5d0/apow)**2.d0 ! Cosine^2 apodization
  endif

!! Tests that the mask is centered around the pointed pixel.
!     if (dotprod.ge.0.999999d0) then
!        mask(i,1)=20.d0 
!        print*, theta*180/pi,theta0*180/pi,phi*180/pi,min(phi*180/pi,(360.d0-(phi*180/pi)))
!     endif
enddo

end subroutine return_disc_mask
!########################################################################

!!########################################################################
!subroutine return_apodized_disk_masks()
!
!implicit none
!integer(i4b) :: i,j,whichapo=2 ! (1 : Cosine apodization, 2 : Cosine squared apodization)
!real(dp) :: theta,phi,theta0,phi0,dtheta,dotprod
!!real(dp) :: apobeamwd=2.d0,cut=11.3d0 ! degrees
!real(dp) :: apobeamwd=2.d0,cut=11.3d0 ! degrees
!
!
!do j=0,npixtotout-1
!   mask=0.d0
!   call pix2ang_ring(nsideout, j, theta0, phi0)
!
!   if ((180*theta0/pi.lt.(90.d0-35.d0)).or.(180*theta0/pi.gt.(90.d0+35.d0))) then
!      do i=0,npixtot-1
!        call pix2ang_ring(nside,i,theta,phi)
!        dotprod=sin(theta)*sin(theta0)*cos(phi-phi0)+cos(theta)*cos(theta0)
!        dtheta=acos(dotprod)*180.d0/pi
!        
!        if (dtheta.lt.cut-apobeamwd) then
!          mask(i,1)=1.d0
!  !     elseif ((dtheta.ge.cut).and.(dtheta.le.(cut+apobeamwd))) then
!        elseif ((dtheta.ge.(cut-apobeamwd)).and.(dtheta.lt.cut)) then
!          dtheta=dtheta-(cut-apobeamwd)
!          if (whichapo.eq.1) mask(i,1)=(1.d0+cos(dtheta*pi/apobeamwd))/2.d0 ! Cosine apodization
!          if (whichapo.eq.2) mask(i,1)=cos(dtheta*pi*0.5d0/apobeamwd)**2.d0 ! Cosine^2 apodization
!        endif
!      
!  !    Multiply the disk mask with the original mask, especially important for pixels close to the mask edges.
!        mask(i,1)=mask(i,1)*origmask(i,1)
!      
!  !   ! Tests that the mask is centered around the pointed pixel.
!  !        if (dotprod.ge.0.999999d0) then
!  !           mask(i,1)=20.d0 
!  !           print*, theta*180/pi,theta0*180/pi,phi*180/pi,min(phi*180/pi,(360.d0-(phi*180/pi)))
!  !        endif
!      enddo
!
!      fskymap(j,1)=sum(mask(:,1)/float(npixtot))
!
!!     Writing out the generated disk mask.
!      write(s2,"(i4.1)")j
!      filename="./diskmasks/"//"mask"//trim(adjustl(s2))//"_nside"//trim(adjustl(strnside))//".fits" 
!      print*, "Writing the apodized mask",j ; print*, trim(adjustl(filename))
!      call write_bintab(mask,npixtot,1,maskheader,nlheader,"!"//filename)
!   endif
!
!enddo
!
!!	Writing out the fsky map.
!filename=trim(adjustl(pathout))//"fskymap.fits"
!call write_bintab(fskymap,npixtotout,1,mapheader,nlheader,"!"//filename)
!
!print*, "Returned all the disk masks"
!
!
!end subroutine return_apodized_disk_masks
!!########################################################################
!
!!########################################################################
!subroutine read_apodized_diskmasks()
!integer(i4b) :: j
!real(dp) :: theta0,phi0
!
!mask=0.d0
!do j=0,npixtotout-1
!   call pix2ang_ring(nsideout, j, theta0, phi0)
!   if ((180*theta0/pi.lt.(90.d0-35.d0)).or.(180*theta0/pi.gt.(90.d0+35.d0))) then
!      write(s2,"(i4.1)")j ; filename="./diskmasks/"//"mask"//trim(adjustl(s2))//"_nside"//trim(adjustl(strnside))//".fits"
!      call read_bintab(filename,tempmap,npixtot,1,nullval,anynull)
!      mask(:,1)=tempmap(:,1)
!      fskymap(j,1)=sum(mask(:,j)/float(npixtot))
!   endif
!enddo   
!
!filename=trim(adjustl(pathout))//"fskymap.fits"
!call write_bintab(fskymap,npixtotout,1,mapheader,nlheader,"!"//filename)
!
!print*, "Read in all the disk masks"
!
!end subroutine read_apodized_diskmasks
!!########################################################################
!

end module gen_mask

