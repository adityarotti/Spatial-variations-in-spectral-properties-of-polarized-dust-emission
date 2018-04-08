use global_param
use shared_data
use process_mask
use map_operations

implicit none

read(*,*) paramfile
call read_param(paramfile)
print*, "Done reading paramfile"

! INITIALIZATION
print*, "INITIALIZATION"
print*, "=========================================="
call allocate_data() ; print*, " --> Allocated data arrays"
call allocate_mask_arrays() ; print*, " --> Allocated mask arrays"
call return_corrections() ; print*, " --> Imported corrections"
print*, "=========================================="
print*, " "

call read_data_mask() ; print*, " --> Read in data mask"

do simnum=nstart,nstart+nrlz-1
    ! Initialization
    mapin=0.d0
    mapout=0.d0
    residual=0.d0

    write(s1,"(i4.1)") simnum

!   Reading in the data maps
    ! Filename "IN"
    if (swSUFFIX) then
       datafname=trim(adjustl(finprefix))//trim(adjustl(s1))//trim(adjustl(finsuffix))
    else
       datafname=trim(adjustl(finprefix))
    endif
    call return_data(datafname) ; print*, " --> Returned read-in data"
!   ------------------------------------------


!   Setting the filename for the file to be written out.
    ! Filename "OUT"
    if (swSUFFIX) then
       datafname=trim(adjustl(foutprefix))//trim(adjustl(s1))//trim(adjustl(foutsuffix))   
    else
       datafname=trim(adjustl(foutprefix))
    endif
    datafname=trim(adjustl(pathout))//datafname
!   ------------------------------------------
    
!   Fullsky E/B decomposition.
    if (swDOFS) then
       print*, "=========================================="
       print*, " ==> PERFORMING THE FULLSKY ANALYSIS"
       print*, "==========================================" 
    
       swMASK=.False.
    
    !  call convert_TQU2TEB() ; print*, " --> Converted TQU to TEB maps"
    !  This only works for the full sky, where there is no residue subtraction.
    !  For it to work with masked skies another module needs to be written which removes
    !  the prefactor even from the residue term. This is simple to do but not incorporated yet.
       call convert_TQU2TEB_tilde() ; print*, " --> Converted TQU to TEB maps without prefactor"
    
       call write_data(datafname) ; print*, " --> Writing out the TEB maps"
       print*, "=========================================="
       print*, " "
    else
   !   Masked E/B decomposition.
       print*, "=========================================="
       print*, " ==> PERFORMING THE MASKED ANALYSIS"
       print*, "=========================================="
   !   Masked sky E/B decomposition (cosine apodized band mask)
       swMASK=.True.
       call convert_TQU2TEB() ; print*, " --> Converted TQU to TEB maps"
       call calc_residual()
       call write_data(datafname) ; print*, " --> Writing out reconstructed TEB maps"
       print*, "=========================================="
       print*, " "
    endif
!   ------------------------------------------
enddo

print*, "DEALLOCATING ARRAYS"
print*, "=========================================="
call deallocate_data() ; print*, " --> Deallocated the data arrays"
call deallocate_mask_arrays() ; print*, " --> Deallocated the mask arrays"
print*, "=========================================="
print*, " "

print*, " "
print*, ">>>>> Done HURRAY !!!! <<<<< "
print*, "Check the results for this run in ==>",trim(adjustl(pathout))
print*, " "

stop
end
