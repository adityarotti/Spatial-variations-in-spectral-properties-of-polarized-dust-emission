subroutine divide(Lmin,Lmax,lmax1,nt,nmax,Li,Lf)
implicit none
integer l,l1,l2,m,r,nmax,nt,j
integer af,sum1(1:nt),Lmax,Lmin,lmax1,l2min,l2max
integer b,temp,temp1,a,counter
integer Li(nt),Lf(nt)
real*8 sum

sum=0
do l=lmin,lmax
        do m=0,l
                do l1=0,lmax1
                 l2min=l1
             If (Abs(l1-l).ge.l1) l2min=Abs(l1-l)
                 l2max=lmax1
             If ((l+l1).lt.lmax1) l2max=(l+l1)
                         do l2=l2min,l2max
!                                do r=1,2*l1+1   
                                        sum=sum+1
!                                enddo
                         enddo
                 enddo
        enddo          
!	print*,l,sum
enddo
sum=sum/nt

l=Lmin-1
j=1
sum1=0
        do while ((j.le.nt).and.(l.lt.Lmax))
                Li(j)=l+1
        do while ((sum1(j).lt.INT(sum)).and.(l.lt.Lmax))
                l=l+1
                do m=0,l
                 do l1=0,lmax1
                         l2min=l1
                             If (Abs(l1-l).ge.l1) l2min=Abs(l1-l)
                         l2max=lmax1
                             If ((l+l1).lt.lmax1) l2max=(l+l1)
                         do l2=l2min,l2max
!                                do r=1,2*l1+1
                                        sum1(j)=sum1(j)+1
!                                enddo
                         enddo
                 enddo
                enddo         
        enddo
                Lf(j)=l
                l=Lf(j)
                j=j+1
        enddo
nmax=j-1
               
     return
     END
   
     
     
