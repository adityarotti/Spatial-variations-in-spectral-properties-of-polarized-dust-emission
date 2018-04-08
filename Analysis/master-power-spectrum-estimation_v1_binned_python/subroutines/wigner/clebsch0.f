	subroutine clebsch0(NDIM,l,l1,m,m1,l2min,l2max,cleb)

	implicit none

	integer ier,i,m2,NDIM
        real*8 l,l1,l2,m,m1,l2min,l2max ! Wigner 3j variables
        real*8 cleb(NDIM),THRCOF(NDIM)

        CALL DRC3JJ(l, l1, m, m1, l2min, l2max, THRCOF, NDIM, ier)

        do i=1,int(l2max-l2min)+1
            l2=l2min+float(i-1)
            cleb(i)= ((-1.d0)**(l1-l2))*(sqrt(2.d0*l+1.d0))*THRCOF(i)
        enddo

        return
        end 
