	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>code to allocate arrays, and call activation 
	module bam
	use numerics_type
	implicit none
		real(wp), parameter :: grav=9.8_wp, lv=2.5e6_wp,cp=1005._wp,molw_air=29.e-3_wp,     &
									   r_air=287._wp,r_vap=461._wp,r_gas=8.314_wp,molw_vap=18.e-3_wp, &
									   eps=r_air/r_vap,kappa=r_air/cp,rhow=1000._wp,sigma=72.e-3_wp
									   ! private
		real(wp) :: rhinit,tinit,pinit,w, ndrop_test, &
							mass_dummy,density_dummy,n_dummy,sig_dummy,d_dummy, &
							tcb, pcb,xmin,a,smax,smax1, &
							alpha_sup, sigma_sup, g, chi, sd_dummy, s,c0   
		! size n_mode
		real(wp), allocatable, dimension(:) :: n_aer1, sig_aer1, d_aer1, &
										n_aer, sig_aer, d_aer, d_aer_new, sgi, &
										density_final,mass_initial, & 
									   mass_final,sd,b,sm,eta,f1,f2 
							  
		! size n_mode
		real(wp), allocatable, dimension(:) :: density_core, & 
								  molw_core, & 
								  nu_core,act_frac,act_frac1, act_frac2 , &
								  density_core1, & 
								  molw_core1, & 
								  nu_core1, dcrit2
		! size n_sv
		real(wp), allocatable, dimension(:) :: molw_org, r_org, log_c_star, cstar, &
												org_content, org_content1, &
											  density_org,nu_org,mass_org_condensed, &
											  delta_h_vap, epsilon1, c_ions, &
											  molw_org1, log_c_star1, &
											  density_org1,nu_org1,&
											  delta_h_vap1
											  
		real(wp), dimension(6) :: c1	! private
		logical(lgt) :: check			! private
		integer(i4b) :: n_mode_s		! private
		real(wp) :: p_test, t_test, w_test, a_eq_7, b_eq_7 ! public
		
		integer(i4b) :: n_mode, n_sv, method_flag, giant_flag, sv_flag 
		! 1=abdul-razzak, ghan; 2=fountoukis and nenes; 3=fountoukis and nenes with quad
	
		! for random number:
		real(wp) :: r, mean_w, sigma_w
		real(wp), dimension(10,10) :: rs
		integer(i4b), allocatable, dimension(:) :: seed
		integer(i4b) :: l, n_rand
		logical :: rand_dist=.false.

	private 
	public :: ctmm_activation, allocate_arrays, initialise_arrays, &
		read_in_bam_namelist, find_d_and_s_crits, &
		n_mode, n_sv, method_flag, giant_flag, sv_flag, &
		p_test, t_test, w_test, a_eq_7, b_eq_7, n_aer1, d_aer1, sig_aer1, &
		org_content1, molw_org1, log_c_star1, density_org1, nu_org1, delta_h_vap1, &
		molw_core1, density_core1, nu_core1, act_frac1, &
		r, mean_w, sigma_w, rs, seed, l, n_rand, rand_dist, &
		smax1, dcrit2
				
	contains
	
	
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>root-find to solve for smax and dcrit needed so that drop number is achieved
	!>parameterisation developed at university of manchester
	!>@param[in] p,t
	!>@param[inout] ndrop drop number concentration
	!>@param[inout] scrit, dscrit
	subroutine find_d_and_s_crits(p,t,ndrop,w,smax,dcrit)
	    use numerics_type
	    use numerics, only : zeroin
	    implicit none
	    real(wp), intent(in) :: p,t
	    real(wp), intent(inout) :: ndrop, smax, w
	    real(wp), dimension(n_mode), intent(inout) :: dcrit
	    
	    
	    ! adjust ndrop to 99% of total aerosol number:
	    ndrop=min(ndrop,0.99_wp*sum(n_aer1))
        ndrop_test=ndrop
        
        ! for use inside the root-finder
        pinit=p
        tinit=t
        ! find the updraft speed required to activate these aerosol particles
        w=zeroin(1.e-40_wp,100._wp,find_wcbase,1.e-20_wp)
        ! call again so that dcrit1 and smax are properly set
        call initialise_arrays(n_mode,n_sv,pinit,tinit,w, &
                    n_aer1,d_aer1,sig_aer1, molw_org1,density_core1)

        call ctmm_activation(n_mode,n_sv,sv_flag, &
                    n_aer1, d_aer1,sig_aer1,molw_core1, &
                    density_core1, nu_core1, &
                    org_content1,molw_org1, density_org1, delta_h_vap1, nu_org1,  &
                    log_c_star1, &
                    w, tinit,pinit, a_eq_7, b_eq_7, &
                    act_frac1,smax,dcrit)
        
	    
	end subroutine find_d_and_s_crits
	
	
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>root-finding helper function
	!>@param[in] w
	function find_wcbase(w)
	    use numerics_type
	    implicit none
	    real(wp), intent(in) :: w
	    real(wp) :: find_wcbase, smax
	
        call initialise_arrays(n_mode,n_sv,pinit,tinit,w, &
                    n_aer1,d_aer1,sig_aer1, molw_org1,density_core1)

        call ctmm_activation(n_mode,n_sv,sv_flag, &
                    n_aer1, d_aer1,sig_aer1,molw_core1, &
                    density_core1, nu_core1, &
                    org_content1,molw_org1, density_org1, delta_h_vap1, nu_org1,  &
                    log_c_star1, &
                    w, tinit,pinit, a_eq_7, b_eq_7, &
                    act_frac1,smax,dcrit2)

	    find_wcbase=ndrop_test-sum(act_frac1*n_aer1)
	    
	end function find_wcbase
	
	
	
	
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculate the activated fraction of a lognormally distributed
	!>aerosol distribution including condensation of semi-volatile organics.
	!>code follows a paper by Connolly, Topping, Malavelle and Mcfiggans (in ACP 2014)
	!>parameterisation developed at university of manchester
	!>@param[in] n_modes1, n_sv1 : number of aerosol and volatility modes
	!>@param[in] sv_flag: flag for cocondensation
	!>@param[inout] n_aer, d_aer, sig_aer, molw_core, density_core, nu_core
	!>@param[in] org_content1: amount of organic in volatility bins
	!>@param[in] molw_org1: molecular weight in volatility bins
	!>@param[in] density_org1: density of organic in volatility bins
	!>@param[in] delta_h_vap1: enthalpy change in volatility bins
	!>@param[in] nu_org1: van hoff factor in volatility bins
	!>@param[in] log_c_star1: volatility bins
	!>@param[in] w1, t1, p1, a, b: vertical wind, temperature, pressure + params in ARG
	!>@param[inout] act_frac1: activated fraction in each mode
	!>@param[inout] smax1: maximum supersaturation
	!>@param[inout] dcrit1: critical diameters in each mode
	subroutine ctmm_activation(n_modes1,n_sv1,sv_flag, n_aer1,d_aer1,sig_aer1,molw_core1, &
							   density_core1, nu_core1, org_content1, &
							   molw_org1, density_org1, delta_h_vap1, nu_org1,  &
                               log_c_star1, &
                               w1, t1,p1,a_arg,b_arg, &
							   act_frac1,smax1,dcrit1)

	    use numerics_type
		use numerics, only : zeroin, fmin
		implicit none 
		      real(wp), dimension(:), intent(in) :: n_aer1
			  real(wp), dimension(:), intent(inout) :: d_aer1, sig_aer1, molw_core1, &
													density_core1, nu_core1
			  real(wp), dimension(:), intent(in) :: org_content1  , molw_org1, &
			  							density_org1, delta_h_vap1, nu_org1, log_c_star1                               
			  real(wp), intent(in) :: w1,t1,p1, a_arg, b_arg
			  integer, intent(in) :: n_modes1, n_sv1, sv_flag
			  real(wp), dimension(:), intent(inout) :: act_frac1, dcrit1
			  real(wp), intent(inout) :: smax1

		integer(i4b):: i
		
		n_mode_s=n_modes1
		n_aer=n_aer1
		d_aer=d_aer1
		sig_aer=sig_aer1
		molw_core=molw_core1
		density_core=density_core1
		nu_core=nu_core1
		

		org_content=org_content1*1e-9_wp/(p1/r_air/t1) ! kg/kg
		molw_org=molw_org1
		density_org=density_org1
		delta_h_vap=delta_h_vap1
		nu_org=nu_org1
		log_c_star=log_c_star1

		
		if(sv_flag.eq.1) then

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! Find how much semi-volatile is condensed
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			call solve_semivolatiles(n_modes1,n_sv1, &
					org_content, log_c_star, delta_h_vap, nu_org, molw_org, &
					mass_initial, nu_core, molw_core,rhinit, t1, &
					mass_org_condensed)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! distribute mass in proportion to number - this is wrong - Crooks has a 
			! better method
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			mass_final=mass_initial+sum(mass_org_condensed)* &
									n_aer/sum(n_aer(1:n_modes1)) 
									! final mass after co-condensation
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! calculate the new density - this is wrong - Crooks has a 
			! better method
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			density_final=mass_initial/density_core+sum(mass_org_condensed/density_org) * &
													n_aer/sum(n_aer(1:n_modes1))
			density_final=mass_final/density_final
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! calculate the arithmetic standard deviations                               !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			sd=exp(log(d_aer)+0.5_wp*sig_aer)*sqrt(exp(sig_aer**2)-1._wp)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! note that for multiple modes, assume sd remains constant for all modes and 
			! shift each mode by the same amount in diameter space, such that the total 
			! mass constraint is satisfied. (see Connolly et al, 2014 acp))
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! now calculate the new median diameter                                      !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			do i=1,n_modes1
				density_dummy=density_final(i)
				mass_dummy=mass_final(i)
				n_dummy=n_aer(i)
				sd_dummy=sd(i)
				d_aer_new(i)=fmin(d_aer(i),2000.e-9_wp, mass_integrate,1.e-30_wp)
				xmin = mass_integrate(d_aer_new(i))
			enddo
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! now calculate the new geometric standard deviation                         !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			d_aer=d_aer_new
			do i=1,n_modes1
				d_dummy=d_aer_new(i)
				sd_dummy=sd(i)
				sig_aer(i)=zeroin(1e-9_wp,2e0_wp, find_sig_aer,1.e-6_wp)
			enddo
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		else if(sv_flag.eq.0) then
			mass_org_condensed=0._wp
			mass_final=mass_initial
			density_final=density_core
		endif

		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! now calculate the activated fraction                                           !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if(method_flag.eq.1) then
			a=2._wp*sigma*molw_vap/(rhow*r_gas*tcb)             ! eq 5: abdul-razzak, ghan
			! above, 2 should be 4 for fountoukis and nenes
  
			b=density_final/ &
			  ( (molw_core*mass_initial/nu_core+ & 
			  sum(molw_org*mass_org_condensed/nu_org)* &
			    mass_initial/sum(mass_initial(1:n_modes1))) / &
			  mass_final)/(rhow/molw_vap)                       ! eq 6: abdul-razzak, ghan
			sm=2._wp/sqrt(b)*(a/(3._wp*d_aer/2._wp))**1.5_wp    ! eq 8: abdul-razzak, ghan 
																! or 9: of 2000 paper

			alpha_sup=grav*molw_vap*lv/(cp*r_gas*tcb**2)- &
					  grav*molw_air/(r_gas*tcb)                ! eq 11: abdul-razzak, ghan

			sigma_sup=r_gas*tcb/(svp(tcb)*molw_vap)+ &
					  molw_vap*lv**2/(cp*pcb*molw_air*tcb)     ! eq 12: abdul-razzak, ghan

			g=rhow*r_gas*tcb/(svp(tcb)*dd(tcb,pcb)*molw_vap) + &
			  lv*rhow/(ka(tcb)*tcb)* &
			  (lv*molw_vap/(r_gas*tcb)-1._wp)                     
			g=1._wp/g                                          ! eq 16: abdul-razzak, ghan

			eta=(alpha_sup*w/g)**1.5_wp/ &
				(2._wp*pi*rhow*sigma_sup*n_aer)                ! eq 22: abdul-razzak, ghan
															   ! or 11: of 2000 paper

			chi=(2._wp/3._wp)*(alpha_sup*w/g)**0.5_wp*a       ! eq 23: abdul-razzak, ghan
															   ! or 10: of 2000 paper
 
			! f1=1.5_wp*exp(2.25_wp*sig_aer**2)                ! eq 28: abdul-razzak, ghan
			f1=a_arg*exp(b_arg*sig_aer**2)                  ! or 7 : of 2000 paper 
															   ! a=0.5, b=2.5
										
		    f2=1._wp+0.25_wp*sig_aer                          ! eq 29: abdul-razzak, ghan
															   ! or 8 : of 2000 paper
 
			! act_frac=0.5_wp*erfc(log(f1*(chi/eta)**1.5_wp &
			!          +f2*(sm**2/(eta+3_wp*chi))**0.75_wp) / &
			!          (3._wp*sqrt(2_wp)*sig_aer))             ! eq 30: abdul-razzak, ghan

			smax=sum(1._wp/sm(1:n_modes1)**2* &
			      (f1(1:n_modes1)*(chi/eta(1:n_modes1))**1.5_wp+ &
				  f2(1:n_modes1)*(sm(1:n_modes1)**2._wp / &
				  (eta(1:n_modes1)+3._wp*chi))**0.75_wp ))**0.5_wp
			smax=1._wp/smax					               ! eq 6: of 2000 paper
													
			! smax=(f1*(chi/eta)**1.5_wp+ &
			!       f2*(sm**2/eta)**0.75)**0.5
			! smax=sm/smax                                     ! eq 31: abdul-razzak, ghan

			act_frac1=1._wp/sum(n_aer(1:n_modes1))* &
			          sum(n_aer(1:n_modes1)*5.e-1_wp*(1._wp- &
			          erf(2._wp*log(sm(1:n_modes1)/smax)/ &
			          (3._wp*sqrt(2._wp)*sig_aer(1:n_modes1)) )))   ! eq 13: of 2000 paper

			act_frac2=1._wp/(n_aer)*(n_aer*5.e-1_wp*(1._wp- &
		     erf(2._wp*log(sm/smax)/(3._wp*sqrt(2._wp)*sig_aer) ))) ! eq 13: of 2000 paper

            dcrit2=2._wp*a/3._wp*(2./smax/sqrt(b))**(2._wp/3._wp)
		     
		else if(method_flag.ge.2) then
			a=4._wp*sigma*molw_vap/(rhow*r_gas*tcb)             ! eq 5: abdul-razzak, ghan
			! above, 2 for abdul-razzak,4 for fountoukis and nenes
  
			b=density_final/ &
			  ( (molw_core*mass_initial/nu_core+ & 
			  sum(molw_org*mass_org_condensed/nu_org)* &
			  	mass_initial/sum(mass_initial(1:n_modes1))) / &
			  	mass_final)/(rhow/molw_vap)                   ! eq 6: abdul-razzak, ghan

			sgi=sqrt(4._wp*a**3/27._wp/(b*d_aer**3))           ! eq 17 fountoukis and nenes

			smax=max(zeroin(1.e-20_wp,100.e-2_wp, fountoukis_nenes,1.e-20_wp),1.e-20_wp)
			
			!act_frac=brent(10d-2,1.d-10,1.d-20,fountoukis_nenes,1.d-30,smax)
			!smax=max(smax,1.d-20)

			!act_frac=sum(0.5_wp*&
			!   erfc(2_wp*log(sgi/smax)/(3_wp*sqrt(2_wp)*sig_aer))) ! eq 8 and 9 f+n 
			act_frac1=1._wp/sum(n_aer(1:n_modes1))*&
			 sum(n_aer(1:n_modes1)*5.e-1_wp*(1._wp- &
			 erf(2._wp*log(sgi(1:n_modes1)/smax)/ &
			 (3._wp*sqrt(2._wp)*sig_aer(1:n_modes1)) )))   ! eq 13: of 2000 paper

			act_frac2=1._wp/(n_aer)*(n_aer*5.e-1_wp*(1._wp- &
			 erf(2._wp*log(sgi/smax)/(3._wp*sqrt(2._wp)*sig_aer) )))! eq 13: of 2000 paper

            dcrit2=( (4._wp*a**3) / (smax1**2 *(27._wp*b)) )**(1._wp/3._wp)
		 end if
 
		 act_frac1=act_frac2
		 smax1=smax
		 dcrit1(1:n_modes1)=dcrit2(1:n_modes1)
		 
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
	end subroutine ctmm_activation
	
	
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! diffusivity of water vapour in air										   !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the diffusivity of water vapour in air
	!>@param[in] t: temperature
	!>@param[in] p: pressure
	!>@return dd: diffusivity of water vapour in air
	function dd(t,p)
	    use numerics_type
		implicit none
		real(wp), intent(in) :: t, p
		real(wp) :: dd
		dd=2.11e-5_wp*(t/273.15_wp)**1.94*(101325._wp/p)
	end function dd

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! conductivity of water vapour												   !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the thermal conductivity of air
	!>@param[in] t: temperature
	!>@return ka: thermal conductivity of air
	function ka(t)
	    use numerics_type
		implicit none
		real(wp), intent(in) :: t
		real(wp) :: ka
		ka=(5.69_wp+0.017_wp*(t-273.15_wp))*1.e-3_wp*4.187_wp
	end function ka

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! dry potential temperature 												   !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the thermal conductivity of air
	!> Note uses tinit, pinit, rhinit from module
	!>@param[in] p: pressure (Pa)
	!>@return dry_potential temperature (K)
	function dry_potential(p)
	    use numerics_type
		implicit none
		real(wp), intent(in) :: p
		real(wp) :: dry_potential
		real(wp) :: total_water1, total_water2, tcalc 
		total_water1=rhinit*eps*svp(tinit)/(pinit-svp(tinit))

		!print *,'kappa',svp(tinit),eps
		tcalc=tinit*(p/pinit)**kappa

		total_water2=eps*svp(tcalc)/(p-svp(tcalc))

		dry_potential=total_water2-total_water1

	end function dry_potential

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! saturation vapour pressure over liquid                                       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the saturation vapour pressure over liquid water according to buck fit
	!>@param[in] t: temperature
	!>@return svp_liq: saturation vapour pressure over liquid water
	function svp(t)
	    use numerics_type
		implicit none
		real(wp), intent(in) :: t
		real(wp) :: svp
		svp = 100._wp*6.1121_wp* &
			  exp((18.678_wp - (t-273.15_wp)/ 234.5_wp)* &
			  (t-273.15_wp)/(257.14_wp + (t-273.15_wp)))
	end function svp


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 3rd moment - for integration												   !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the third moment of a lognormal
	!>@param[in] d: diameter (m)
	!>@return ln3: third moment in a size interval
	function ln3(d)
	    use numerics_type
		implicit none
		real(wp), dimension(:), intent(in) :: d
		real(wp), dimension(size(d)) :: ln3

		! add all modes together
		ln3=pi*d**2/(sqrt(twopi*sig_dummy**2)*6.)* &
		exp(-log(d/d_dummy)**2/(2.*sig_dummy**2))* &
		density_dummy*n_dummy
	end function ln3

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! mass in a lognormal       												   !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the mass across one lognormal distribution
	!>@param[in] d1: mode diameter (m)
	!>@return mass_integrate: total mass in the distribution
	function mass_integrate(d1)
	    use numerics_type
		use numerics, only : zeroin
		implicit none
		real(wp), intent(in) :: d1
		real(wp) :: mass_integrate

		integer(i4b):: i
		d_dummy=d1  ! guess at d_aer, used to calculate the new standard deviation
		sig_dummy=zeroin(1.e-9_wp,2._wp,find_sig_aer,1.e-6_wp)
		!  mass_integrate=abs(qromb(ln3,0.d-10,3.d-6)-mass_dummy)
		! moment generating function
		! http://www.mlahanas.de/math/lognormal.htm
		mass_integrate=n_dummy*exp(3._wp*log(d_dummy) + &
					3._wp**2*sig_dummy**2/2._wp) &
				   *density_dummy*pi/6._wp
		mass_integrate=abs(mass_integrate-mass_dummy)

	end function mass_integrate

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! find sigma for the size distribution										   !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the thermal conductivity of air
	!>@param[in] sig_aer_new: geometric standard deviation
	!>@return find_sig_aer: arithmetic standard deviation
	function find_sig_aer(sig_aer_new)
	    use numerics_type
		real(wp), intent(in) :: sig_aer_new
		real(wp) :: find_sig_aer
		real(wp) :: sd1
		sd1=exp(log(d_dummy)+0.5_wp*sig_aer_new)*sqrt(exp(sig_aer_new**2)-1._wp)
		find_sig_aer=sd1-sd_dummy
	end function find_sig_aer

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! fountoukis and nenes integrals											   !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the integral in the fountoukis and nenes method
	!>@param[in] smax1: max supersaturation
	!>@return fountoukis_nenes: integral for activation
	function fountoukis_nenes(smax1)
	    use numerics_type
		use numerics, only : romb
		implicit none
		real(wp), intent(in) :: smax1
		real(wp) :: fountoukis_nenes

		real(wp) :: integral, integral1,smax2,spart,upart,umax,i1,i2,discriminant, &
				  deq,dw3
		integer(i4b):: i

		smax2=max(smax1,1e-20_wp)
		alpha_sup=grav*molw_vap*lv/(cp*r_gas*tcb**2)- &
		grav*molw_air/(r_gas*tcb)                ! eq 11 arg
		sigma_sup=r_gas*tcb/(svp(tcb)*molw_vap)+ &
				molw_vap*lv**2/(cp*pcb*molw_air*tcb)  ! eq 11?
		g=rhow*r_gas*tcb/(svp(tcb)*dd(tcb,pcb)*molw_vap) + &
		lv*rhow/(ka(tcb)*tcb)*(lv*molw_vap/(r_gas*tcb)-1_wp)  ! eq 12
		g=4_wp/g
		! do the integral 
		! f-n (2005 use approximate form for integral, here we use quadrature)
		c1(1)=2._wp*a/3._wp
		c1(2)=g/alpha_sup/w
		c1(6)=smax2


		!  discriminant=smax2**4_wp-16._wp*a**2_wp*alpha_sup*w/(9_wp*g)
		discriminant=1._wp-16._wp*a**2._wp*alpha_sup*w/(9._wp*g)/smax2**4._wp

		if(discriminant.ge.0._wp) then
		 spart=smax2*(0.5_wp*(1._wp+ &
		   (1._wp-16_wp*a**2*alpha_sup*w/(9._wp*g*smax2**4 ))**0.5_wp) &
		   )**0.5_wp
		else
		 spart=smax2*min(2.e7_wp*a/3_wp*smax2**(-0.3824_wp),1._wp) 
		endif

		integral=0._wp
		integral1=0._wp
		do i=1,n_mode_s
		  c1(3)=2_wp*n_aer(i)/(3._wp*sqrt(2._wp*pi)*sig_aer(i))
		  c1(4)=sgi(i)
		  c1(5)=2._wp*sig_aer(i)**2
		  ! note for multiple modes, have to calculate the integral below for each mode
		  if(method_flag.eq.3) then
			 integral1=integral1+romb(integral3_fn,0._wp,smax2)
		  endif

		! approximate method
		  upart=2._wp*log(sgi(i)/spart)/(3._wp*sqrt(2._wp)*sig_aer(i))
		  umax =2._wp*log(sgi(i)/smax2)/(3._wp*sqrt(2._wp)*sig_aer(i))

		  i1=n_aer(i)/2._wp*sqrt(g/alpha_sup/w)*smax2* &
			 (erfc(upart)-5e-1_wp*(sgi(i)/smax2)**2_wp*exp(4.5_wp*sig_aer(i)**2)* &
			 erfc(upart+3_wp*sig_aer(i)/sqrt(2._wp)))
		  i2=a*n_aer(i)/(3_wp*sgi(i))*exp(9_wp/8_wp*sig_aer(i)**2_wp)* &
			 (erf(upart-3_wp*sig_aer(i)/(2._wp*sqrt(2._wp))) - &
			 erf(umax-3._wp*sig_aer(i)/(2._wp*sqrt(2._wp) )) )
		!      print *,'i1, i2',i1,i2, upart, umax

		  if(method_flag.eq.2) then
			 integral1=integral1+(i1+i2)
		  endif   

		  ! for the giant ccn - barahona et al (2010)
		  if(giant_flag.eq.1) then
			dw3  = (sqrt(2._wp)*log(sgi(i)/spart)/3_wp/sig_aer(i))- &
				   (1.5_wp*sig_aer(i)/sqrt(2._wp))
			deq= a*2_wp/sgi(i)/3_wp/sqrt(3._wp)      
			dw3=n_aer(i)*deq*exp(9_wp/8_wp*sig_aer(i)**2._wp)*smax2* &
			   (erfc(dw3))*((g*alpha_sup*w)**0.5_wp)   
		!         
			dw3=dw3/(2._wp*g*smax2)*(g/alpha_sup/w)**0.5_wp
			integral1=integral1 +dw3
		  endif
		enddo
		!  print *,'integral',integral,integral1
		! cost function - eq 10 fountoukis and nenes

		fountoukis_nenes=(2._wp*alpha_sup*w/(pi*sigma_sup*rhow)- &
					   g*smax2*integral1)

	end function fountoukis_nenes

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! integrand 3 in FN          												   !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the integrand in the FN method (for the quadrature case)
	!>@param[in] s: s parameter
	!>@return integral3_fn: integrand in FN (for use with quadrature)
	function integral3_fn(s)
	    use numerics_type
		implicit none
		real(wp), dimension(:), intent(in) :: s
		real(wp), dimension(size(s)) :: integral3_fn

		integral3_fn=((c1(1)/(s+1.e-50_wp))**2+c1(2)*(c1(6)**2-s**2))**0.5_wp* &
				   c1(3)/(s+1.e-50_wp)*exp(-log((c1(4)/(s+1.e-50_wp))**(2._wp/3._wp))**2/c1(5)) ! eq 14 f-n
	end function integral3_fn
		
		
	
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculate the condensed organics (in Crooks, GMD 2016)
	!>parameterisation developed at university of manchester
	!>@param[in] n_modes1: number of aerosol modes
	!>@param[in] n_sv1: number of volatility bins
	!>@param[in] org_content1: amount of organic in volatility bins
	!>@param[in] nu_org1: van hoff factor in volatility bins
	!>@param[in] molw_org1: molecular weight in volatility bins
	!>@param[in] mass_core: mass in aerosol modes
	!>@param[in] nu_core1: van hoff factor in core
	!>@param[in] molw_core1: molecular weight of modes (core only)
	!>@param[in] s1, t1: rh and temperature.
	!>@param[inout] mass_org_condensed1: condensed organics
	subroutine solve_semivolatiles(n_modes1,n_sv1, &
					org_content1, log_c_star1, delta_h_vap1, &
					nu_org1, molw_org1, &
					mass_core1, nu_core1, molw_core1,s1, t1, &
					mass_org_condensed1)
	    use numerics_type
		use numerics, only : zeroin
		implicit none
		integer(i4b), intent(in) :: n_modes1, n_sv1
		real(wp), dimension(n_sv1), intent(in) :: org_content1, nu_org1, molw_org1, &
							log_c_star1, delta_h_vap1
		real(wp), dimension(n_modes1), intent(in) :: mass_core1, nu_core1, molw_core1
		real(wp), intent(in) :: s1, t1
		real(wp), dimension(n_sv1), intent(inout) :: mass_org_condensed1
		
		
		real(wp) :: ct
		real(wp), dimension(n_sv1) :: c_c
		
		
		! set variables in module (for passing to optimizer)
		nu_org=nu_org1
		molw_org=molw_org1
		nu_core=nu_core1
		molw_core=molw_core1
		s=s1
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		cstar = 10._wp**log_c_star1* (298.15_wp/t1) * &
					exp(-delta_h_vap1*1.e3_wp/r_gas *(1._wp/t1-1._wp/298.15_wp))/1.e9_wp
										 ! c* needs to be adjusted by delta_h_vap / t
		c_ions=org_content1*nu_org1/molw_org1 ! c - all ions
		c0=sum(mass_core1*nu_core1/molw_core1)  ! number of "core" ions
		ct=1._wp/(1._wp-s)*(sum(c_ions)+c0)  ! equation 5 from Crooks et al. (2016, GMD)
										! basically saturation ratio is mole fraction
		! ct is the total concentration of all ions in the condensed phase
		! solve iteratively to find CT in matt's paper
		ct=abs(zeroin(ct,1._wp/(1._wp-s)*c0,partition01,1.e-8_wp))
		!xmin=brent(1.e-15_wp,1.e-8_wp,ct*5._wp, &
		!					partition01,1.e-10_wp,ct)
		
		epsilon1=(1._wp+cstar/ct)**(-1) ! partitioning coefficients
		c_c=c_ions*epsilon1   ! condensed

		mass_org_condensed1=c_c/nu_org1*molw_org1
		
	end subroutine solve_semivolatiles
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! partition01          										        		   !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the partitioning coefficient according to Crooks et al (2016 GMD)
	!>@param[in] ct: total ions, including water
	!>@return guess, minus calculated - for root-finder
	function partition01(ct)
	    use numerics_type
		implicit none
		real(wp), intent(in) :: ct
		real(wp), dimension(size(epsilon1)) :: c_c
		real(wp) :: partition01
		
		real(wp) :: ct1, ct2
		ct2=abs(ct)
		epsilon1 = (1._wp+cstar/(ct2))**(-1)
		c_c=(c_ions)*epsilon1
		ct1=1._wp/(1._wp-s)*(sum(c_c*nu_org)+c0)
		
		partition01=(ct1-ct2)
	
	end function partition01
	
	
		
		
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>allocate arrays for activation code
	!>@param[in] n_modes: number of aerosol modes
	!>@param[in] n_sv: number of organic / volatility modes
	!>@param[inout] n_aer1: number in modes
	!>@param[inout] d_aer1: diameter in modes
	!>@param[inout] sig_aer1: geo std in modes
	!>@param[inout] molw_core1:molw in core
	!>@param[inout] density_core1: solute density
	!>@param[inout] nu_core1: van hoff factor
	!>@param[inout] org_content1: organic content in vol bins
	!>@param[inout] molw_org1: molw in volatility bins
	!>@param[inout] density_org1: density in volatility bins
	!>@param[inout] delta_h_vap1: enthalpy in volatility bins
	!>@param[inout] nu_org1: van hoff factor in volatility bins
	!>@param[inout] log_c_star1: log_c_star in volatility bins
	!>@param[inout] act_frac1: activated fraction in modes
	!>@param[inout] dcrit1: critical diameter in modes
	subroutine allocate_arrays(n_mode,n_sv,n_aer1,d_aer1,sig_aer1, &
			molw_core1,density_core1,nu_core1,org_content1, &
			molw_org1, density_org1,delta_h_vap1,nu_org1,log_c_star1, act_frac1, dcrit1)
	    use numerics_type
		implicit none
		integer(i4b), intent(in) :: n_mode, n_sv
		real(wp), dimension(:), allocatable, intent(inout) :: n_aer1,d_aer1,sig_aer1, &
							molw_core1, density_core1, nu_core1, org_content1, &
							molw_org1, density_org1, delta_h_vap1, nu_org1, log_c_star1, &
							act_frac1, dcrit1
		
		integer(i4b) :: AllocateStatus
		allocate( n_aer(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( d_aer(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( sig_aer(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( n_aer1(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( d_aer1(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( sig_aer1(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( d_aer_new(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( sgi(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( density_final(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( mass_initial(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( mass_final(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( sd(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( b(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( sm(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( eta(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( f1(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( f2(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( density_core(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( density_core1(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( molw_core(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( molw_core1(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( nu_core(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( nu_core1(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( act_frac(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( act_frac1(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( act_frac2(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( dcrit1(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		
		allocate( molw_org(1:n_sv), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( molw_org1(1:n_sv), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( r_org(1:n_sv), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( log_c_star(1:n_sv), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( log_c_star1(1:n_sv), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( cstar(1:n_sv), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( c_ions(1:n_sv), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( epsilon1(1:n_sv), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( org_content(1:n_sv), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( org_content1(1:n_sv), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( density_org(1:n_sv), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( density_org1(1:n_sv), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( nu_org(1:n_sv), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( nu_org1(1:n_sv), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( mass_org_condensed(1:n_sv), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( delta_h_vap(1:n_sv), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( delta_h_vap1(1:n_sv), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		
	
	end subroutine allocate_arrays
	
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>initialise arrays for activation code
	!>@param[in] n_modes: number of aerosol modes
	!>@param[in] n_sv: number of volatility bins
	!>@param[in] p1: pressure (Pa)
	!>@param[in] t1: temperature (K)
	!>@param[in] w1: vertical wind (m/s)
	!>@param[in] n_aer1: number concentration in modes
	!>@param[in] d_aer1: diameter in modes
	!>@param[in] sig_aer1: geometric standard deviation in modes
	!>@param[in] molw_org1: molecular weight in volatility bins
	!>@param[in] density_core1: density in modes
	subroutine initialise_arrays(n_modes,n_sv,p1,t1,w1,n_aer1, &
								d_aer1,sig_aer1, molw_org1,density_core1)
	    use numerics_type
		implicit none
		integer(i4b), intent(in) :: n_modes, n_sv
		real(wp), intent(in) :: p1,t1,w1
		real(wp), dimension(n_modes), intent(in) :: n_aer1,d_aer1,sig_aer1, density_core1
		real(wp), dimension(n_sv), intent(in) :: molw_org1
		
		integer(i4b) :: i
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! define initial conditions                                                      !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		rhinit=0.999_wp                                ! as a fraction - assume we are
													   ! assume we are at cb.
		pinit=p1                                       ! pascals
		tinit=t1                                       ! kelvin
		w    =w1                                       ! m s-1
		n_aer=n_aer1
		d_aer=d_aer1
		sig_aer=sig_aer1
		molw_org=molw_org1
		density_core=density_core1
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		
		
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! define the organic properties                                                  !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!log_c_star=(/(i,i=-0,3)/)                      ! watch out for type?
		!nu_org=1._wp                                    ! disociation factor
		!molw_org=200e-3_wp                                ! kg per mol
		r_org=r_gas/molw_org
		!density_org=1500._wp                            ! kg m-3
		!delta_h_vap=150._wp                             ! enthalpy phase change (kj mol-1)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		
		

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! find the t & p at cloud base                                                   !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		pcb=pinit !zbrent(dry_potential,10000._wp,pinit,1.d-8)
		tcb=tinit !tinit*(pcb/pinit)**kappa
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! initial mass in ith distribution                                               !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		do i=1,n_modes
			density_dummy=density_core(i)                  ! the density of the core dist.
			n_dummy=n_aer(i)
			sd_dummy=sig_aer(i)
			d_dummy=d_aer(i)
			! initial mass in ith distribution
			! moment generating function
			! http://www.mlahanas.de/math/lognormal.htm
			mass_initial(i)=n_dummy*exp(3._wp*log(d_dummy) + &
							3._wp**2_wp*sd_dummy**2/2._wp) &
						   *density_dummy*pi/6._wp
		enddo
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	end subroutine initialise_arrays
	
	
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>read in the data from the namelists for the BAM module
	!>@param[in] nmlfile
	subroutine read_in_bam_namelist(nmlfile)
		implicit none
        character (len=200), intent(in) :: nmlfile
	    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! namelists                                                            !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! set up flags, etc
        namelist /bulk_aerosol_setup/ n_mode, n_sv, sv_flag, &
        			method_flag, giant_flag, &
        			a_eq_7, b_eq_7      
        ! run parameters / input arrays
        namelist /bulk_aerosol_spec/ n_aer1, d_aer1, sig_aer1, &
        					molw_core1, density_core1, &
        					nu_core1, org_content1, molw_org1, density_org1, &
        					delta_h_vap1, &
        					nu_org1, log_c_star1, p_test, t_test, w_test, &
        					rand_dist, n_rand, mean_w,sigma_w
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! read in namelists	and allocate arrays								   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        open(8,file=nmlfile,status='old', recl=80, delim='apostrophe')
        read(8,nml=bulk_aerosol_setup)
        ! allocate memory / init
		call allocate_arrays(n_mode,n_sv,n_aer1,d_aer1,sig_aer1, &
			molw_core1,density_core1,nu_core1,org_content1, &
			molw_org1, density_org1,delta_h_vap1,nu_org1,log_c_star1, &
			act_frac1,dcrit2)
        
        read(8,nml=bulk_aerosol_spec)
        close(8)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	end subroutine read_in_bam_namelist
	

	
	end module bam	

