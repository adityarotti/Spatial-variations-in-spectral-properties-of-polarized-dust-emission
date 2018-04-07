module master

implicit none

contains

!########################################################################
subroutine est_true_cl(cl,wl,lmax,masklmax,info)

implicit none
integer*8, intent(in) :: lmax,masklmax
real*8, intent(in) :: wl(0:masklmax)
real*8 :: ipiv(lmax+1)
real*8, intent(inout) :: cl(0:lmax)
real*8 :: wcl(0:lmax,1)
real*8 :: mllp(0:lmax,0:lmax)
integer*8, intent(out) :: info

wcl(:,1)=cl(:)
call calc_kernel(wl,lmax,masklmax,mllp,info)
!call dgetrs("N",lmax+1,1, mllp, lmax+1, ipiv, wcl, lmax+1, info)
cl(:)=wcl(:,1)

end subroutine est_true_cl
!########################################################################

!########################################################################
subroutine calc_kernel(wl,lmax,masklmax,mllp,info)
implicit none

! Computes the coupling matrix and returns the matrix in LU factored form.

integer*8, intent(out) :: info
integer*8 :: i, j, k, ier,lmax,masklmax
integer*8, parameter:: ndim=2000
real*8 :: l, l1, k1, k1min, k1max,tempvar, wig3j(ndim),ipiv(lmax+1)
real*8, parameter :: pi=3.1415d0
real*8, intent(in) :: wl(0:masklmax)
real*8, intent(out) :: mllp(0:lmax,0:lmax)

mllp=0.d0
do i=0,lmax
   l=float(i)
   do j=0,lmax
      l1=float(j)
      call drc3jj(l,l1,0.d0,0.d0,k1min,k1max,wig3j,ndim,ier)
      do k=1,int(min(k1max,masklmax*1.d0)-k1min)+1
         k1=k1min+float(k-1)
         tempvar=wl(int(k1))*(wig3j(k)**2.d0)*(2.d0*l1+1.d0)*(2.d0*k1+1.d0)
         mllp(i,j)=mllp(i,j)+tempvar/(4.d0*pi)
      enddo
   enddo
   !write(10,11) (mllp(i,j),j=0,lmax)
enddo

call dgetrf(lmax+1, lmax+1, mllp, lmax+1, ipiv, info)

! Computes the inverse which is not necessary.
!call dgetri(lmax+1, , lmax+1, ipiv, work, lmax+1, info)
end subroutine calc_kernel
!########################################################################

end module master
