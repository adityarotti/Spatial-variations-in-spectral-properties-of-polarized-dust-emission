	subroutine clebsch(NDIM,l,l1,l2,m,m1min,m1max,cleb)

	implicit none

	integer ier,i,m2,NDIM
        real*8 l,l1,l2,m,m1min,m1max ! Wigner 3j variables
        real*8 cleb(NDIM),THRCOF(NDIM)

        CALL DRC3JM (l, l1, l2, -m, m1min, m1max, THRCOF, NDIM, ier)

        do i=1,int(m1max-m1min)+1
        m2=int(m1min)+i-1
        if (mod(int(l1-l2+m),2).eq.0) then
            cleb(i)= (sqrt(2.d0*l+1.d0))*THRCOF(i)
        else
            cleb(i)= (-1.d0)*(sqrt(2.d0*l+1.d0))*THRCOF(i)
        endif
	enddo

        return
        end 
