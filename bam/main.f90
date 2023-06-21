	!> @mainpage
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@copyright 2018
	!>@brief
	!>Bulk Aerosol Module (BAM): 
	!>Solves for droplet activation using Abdul-Razzak et al. (1998) or 
	!> Fountoukis and Nenes (2004)
	!> compile using the Makefile and then run using: <br>
	!> ./main.exe  namelist.in
	!> <br><br>



	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>main programme reads namelist, allocates arrays and calls ctmm_activation to do
	!> calculation
	
    program main
	    use numerics_type
        use random, only : random_normal
        use bam, only : n_mode, n_sv, giant_flag, method_flag, sv_flag, &
        	n_aer1, d_aer1, sig_aer1, molw_core1, density_core1, nu_core1, org_content1, &
        	molw_org1, density_org1, delta_h_vap1, nu_org1, log_c_star1, p_test, t_test, &
    		w_test, act_frac1, &
    		a_eq_7, b_eq_7, &
    		r, mean_w, sigma_w, rs, seed, l, n_rand, rand_dist, &
    		allocate_arrays, ctmm_activation, initialise_arrays, &
    		read_in_bam_namelist, smax1, dcrit2
        implicit none

        character (len=200) :: nmlfile = ' '
        
		! for random number:
		integer(i4b) :: i






        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! read in namelists	and allocate arrays								   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call getarg(1,nmlfile)
        call read_in_bam_namelist(nmlfile)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! initialise arrays								                       !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		call initialise_arrays(n_mode,n_sv,p_test,t_test,w_test, &
					n_aer1,d_aer1,sig_aer1, molw_org1,density_core1)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! call activation code												   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(.not. rand_dist) then
        	print *,'act_frac1, n_act'
			call ctmm_activation(n_mode,n_sv,sv_flag, &
						n_aer1, d_aer1,sig_aer1,molw_core1, &
						density_core1, nu_core1, &
						org_content1,molw_org1, density_org1, delta_h_vap1, nu_org1,  &
						log_c_star1, &
						w_test, t_test,p_test, a_eq_7, b_eq_7, &
						act_frac1,smax1,dcrit2)
			print *,w_test,act_frac1,  sum(act_frac1*n_aer1)
		else
			! random number generator
			call random_seed(size=l)
			allocate(seed(1:l))
			seed(:)=2
			call random_seed(put=seed)
			
			

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! call activation code for n_rand samples							   !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			print *,'w,act_frac1, n_act'
			do i=1,n_rand
				r=random_normal() ! from the Netlib
				r=r*sigma_w
				r=r+mean_w
				r=abs(r) ! only positive values
				w_test=r
				
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! initialise arrays								                       !
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				call initialise_arrays(n_mode,n_sv,p_test,t_test,w_test, &
							n_aer1,d_aer1,sig_aer1, molw_org1,density_core1)
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				
				call ctmm_activation(n_mode,n_sv,sv_flag, &
							n_aer1, d_aer1,sig_aer1,molw_core1, &
							density_core1, nu_core1, &
							org_content1,molw_org1, density_org1, delta_h_vap1, nu_org1,  &
							log_c_star1, &
							w_test, t_test,p_test, a_eq_7, b_eq_7, &
							act_frac1,smax1,dcrit2)
				print *,w_test,act_frac1,  sum(act_frac1*n_aer1)
			enddo
			deallocate(seed)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


		
    end program main
