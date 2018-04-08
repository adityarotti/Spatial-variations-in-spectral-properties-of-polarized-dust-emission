module master
use healpix_modules
use healpix_types
use global_param


implicit none
integer(i4b),parameter :: nlheader=80
character(len=80) :: cloutheader(1:nlheader)


contains
!########################################################################
subroutine est_true_cl(clin,mllp,nbin,clout)

implicit none
integer(i4b), intent(in) :: nbin
real*8, intent(in) :: mllp(nbin,nbin)
real*8, intent(in) :: clin(nbin)
real*8, intent(out) :: clout(nbin)
real*8 :: LU(nbin,nbin)

real*8 :: wcl(nbin,1)
integer(i4b) :: info,ipiv(nbin)

info=0; ipiv=0 ; LU=mllp
call dgetrf(nbin, nbin, LU, nbin, ipiv, info)
wcl(:,1)=clin(:) ; !print*, "Converted 0"
call dgetrs("N",nbin,1, LU, nbin, ipiv, wcl, nbin, info) ; print*,"Hello"
clout(:)=wcl(:,1)

!call write_minimal_header(cloutheader,"CL")
!call write_asctab(wcl,lmax,1,cloutheader,nlheader,"!"//fileout)

end subroutine est_true_cl
!########################################################################

!########################################################################
subroutine calc_kernel(wl,lmin,lmax,masklmax,mllp)
implicit none

! Computes the coupling matrix and returns the matrix in LU factored form.

integer*4 :: i, j, k, ier
integer*4, intent(in) :: lmin,lmax,masklmax
integer*8, parameter:: ndim=5000
real*8 :: l, l1, k1, k1min, k1max,tempvar, wig3j(ndim)
real*8, parameter :: pi=3.1415d0
real*8, intent(in) :: wl(0:masklmax)
real*8, intent(out) :: mllp(lmax-lmin+1,lmax-lmin+1)

mllp=0.d0

do i=1,lmax-lmin+1
   l=float(lmin+i-1)
   do j=1,lmax-lmin+1
      l1=float(lmin+j-1)
      call drc3jj(l,l1,0.d0,0.d0,k1min,k1max,wig3j,ndim,ier)
      do k=1,int(min(k1max,masklmax*1.d0)-k1min)+1
         k1=k1min+float(k-1)
         tempvar=wl(int(k1))*(wig3j(k)**2.d0)*(2.d0*l1+1.d0)*(2.d0*k1+1.d0)
         mllp(i,j)=mllp(i,j)+tempvar/(4.d0*pi)
      enddo
   enddo
   !write(10,11) (mllp(i,j),j=0,lmax)
enddo
close(10)

!11 format(1000(x,e15.8))
end subroutine calc_kernel
!########################################################################

end module master
