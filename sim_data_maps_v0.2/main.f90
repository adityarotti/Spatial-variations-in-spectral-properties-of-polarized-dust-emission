use global_param
use generate_noise_map
use simulate_data_map

implicit none


read(*,*) paramfile
call read_param(paramfile)

call allocate_data()

call return_eff_bl() ! The pixel window correction is incorportated in the beam.
call read_TQU_cov_calc_LTM()

! Initializing for random number generator for CMB and noise.
! You have to ensure that year1 and year2 map pairs have the same CMB sky, 
! Hence the seeds for the random number initialization for CMB production is held constant,
! while the random seeds for noise generation are taken in as user inputs.
call rand_init(rng_handle_noise, seed1init+nstart, seed2init+nstart)
call rand_init(rng_handle_cmb, 12134+nstart, 23256+nstart)
call rand_init(rng_handle_frg, 31434+nstart, 76556+nstart)


! NOTE THAT YOU HAVE TO BE CAREFULLY ENSURE THAT YEAR 1 and YEAR 2 NOISE REALIZATIONS
! ARE GENERATED WITH DIFFERENT SEEDS.
do simnum=nstart,nrlz+nstart-1


!  INITIALIZATION
   TQU_map=0.d0 

   if (swNOISE) call return_noise_map()
   call gen_cmb_frg_add_noise_smooth_map()
   
enddo

if (swRECCOV) then
reccov=reccov/float(nrlz) ; call write_reconstructed_cov()
endif

call deallocate_data()

stop 
end


