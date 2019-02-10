	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>radiation solver routines 
    module radiation
    use nrtype
    implicit none


		!>@brief
		!>main model prognostic variables for radiation
        type radg
        	integer(i4b) :: year, mon, day, hour, sec
        	integer(i4b) :: doy, diy
        	integer(i4b) :: nprocv
        	real(sp) :: albedo, emiss, lat, lon
        	integer(i4b) :: ns, nl, ntot
            real(sp), dimension(:), allocatable ::	b_s_g, &
            	lambda, lambda_low, lambda_high, delta_lambda ! cross section for Rayleigh
            real(sp), dimension(:), allocatable ::	ext_s_g ! extinction for Rayleigh
            real(sp), dimension(:,:,:,:), allocatable :: flux_u, flux_d
            real(sp), dimension(:,:,:), allocatable :: rad_power
            real(sp), dimension(:), allocatable :: sflux_l	
            real(sp), dimension(:), allocatable :: mvrecv	
            real(sp), dimension(:), allocatable :: a,b,c,r,u
            integer(i4b) :: tdstart,tdend
            							 
        end type radg



		real(sp), parameter :: r_gas=8.314_sp, ma=29e-3_sp, ra=r_gas / ma, &
								t_stp=288._sp, p_stp=101300._sp, navog=6.022e23_sp, &
								sflux=1370._sp, mu1=1._sp/sqrt(3._sp)
		!>@brief
		!>variables for namelist radiation input
        type namelist_rad
        	integer(i4b) :: start_year, start_mon,start_day, &
        					start_hour, start_min, start_sec
        	real(sp) :: lat_ref,lon_ref
        	integer(i4b) :: quad_flag
        	character (len=200) :: albedo='Urban'
        	character (len=200) :: emissivity='Urban'
            integer(i4b) :: ns,nl
            real(sp), dimension(20) :: lambda_read_s
            real(sp), dimension(20) :: lambda_read_l
            real(sp) :: lambda_s_low, lambda_s_high, lambda_l_low, lambda_l_high
        end type namelist_rad

		integer(i4b), parameter :: nemiss=12, nalb=14

		!>@brief
		!>variables for radiation
        type rad_surf_vars
        	! emissivity
        	real(sp), dimension(nemiss) :: emiss=[1.0,0.99,0.82,0.25,0.1,0.96, &
        									0.94,0.925,0.88,0.96,0.80,0.96]
        	character(len=18), dimension(nemiss) :: emiss_des = &
        	   ["liquid water       ", &
        		"Fresh snow         ",&
        		"Old snow           ",&
        		"Liquid water clouds",&
        		"Cirrus clouds      ", &
        		"Ice                ", &
        		"Soil               ",&
        		"Grass              ",&
        		"Desert             ",&
        		"Forest             ",&
        		"Concrete           ",&
        		"Urban              "]
        	
        	! albedo
        	real(sp), dimension(nalb) :: albedo= &
        		[0.3,0.05,0.85,0.55,0.6,0.45,0.32,0.125,0.21,0.3,0.225,0.125,0.215,0.185]
        	character(len=20), dimension(nalb) :: albedo_des = &
        		["Earth and atmosphere", &
        		 "Liquid water        ", &
        		 "Fresh snow          ", &
        		 "Old snow            ", &
        		 "Thick clouds        ", &
        		 "Thin clouds         ", &
        		 "Sea ice             ", &
        		 "Soil                ", &
        		 "Grass               ", &
        		 "Desert              ", &
        		 "Forest              ", &
        		 "Asphalt             ", &
        		 "Concrete            ", &
        		 "Urban               "]
        end type rad_surf_vars
        
        
        type(rad_surf_vars) :: rsv

		! declare a namelist rad type
		type(namelist_rad) :: nm2
		! declare a radiation grid type
		type(radg) :: radg1

		integer(i4b) :: diy_s, doy_s
		real(sp), parameter :: h_planck=6.6256e-34_sp, &
							c_light=2.9979e8_sp, &	
							k_boltz=1.38e-23_sp, temp_sun=5778._sp, &
						sigma_sb=2._sp*pi**5*k_boltz**4/(15._sp*h_planck**3*c_light**2), &
						fac_planck=2._sp*k_boltz**4/(h_planck**3*c_light**2)*pi, &
			N_a_0 = (p_stp / (t_stp*ra)) / ma * navog ! number of air molecules per m^3 at stp


							
		real(sp), parameter :: diy_exact = 365.242199_sp
		
		real(sp), dimension(:), allocatable :: flux_down, flux_up
		real(sp), dimension(12) :: di_mn=[31, 28, 31,30,31,30,31,31,30,31,30,31], &
								   di_ml=[31, 29, 31,30,31,30,31,31,30,31,30,31]

		private
		public :: e_photon, plancks_law, real_refractive_air, solve_fluxes, solar_zenith, &
					nm2, radg1,allocate_and_set_radiation

		contains
		
		
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!>@author
		!>Paul J. Connolly, The University of Manchester
		!>@brief
		!>initialise variables for radiation
		!>@param[in] lambda
		!>@param[inout] b_s_g_lambda: the scattering cross section
		subroutine allocate_and_set_radiation(start_year, start_mon, &
		  				 				start_day, start_hour,start_min, start_sec, &
		  									ip,jp,kp,ns,nl,ns_r,nl_r, &
											ntot, &
		  									albedo_r,emiss_r, alb, emiss, &
		  									lat_ref,lon_ref, lat, lon, &
											lambda_read_s,lambda_read_l, &
											lambda_s_low, lambda_s_high, &
											lambda_l_low, lambda_l_high, &
											lambda,lambda_low,lambda_high, delta_lambda,&
											sflux_l, b_s_g, &
											ext_s_g, flux_u, flux_d, rad_power, &
											l_h, r_h, rhoan, thetan, &
											nprocv,mvrecv, &
											tdstart,tdend,a,b,c,r,u, &
											coords,dims, id, comm3d)
			use mpi
			implicit none
			character (len=200), intent(in) :: albedo_r, emiss_r
			real(sp), intent(inout) :: alb, emiss
			real(sp), intent(in) :: lat_ref, lon_ref
			real(sp), intent(inout) :: lat, lon
			integer(i4b), intent(in) :: start_year, start_mon, start_day, start_hour, &
									start_min, start_sec
			integer(i4b), intent(in) :: ns_r,nl_r, l_h, r_h
			integer(i4b), intent(inout) :: ip, jp, kp, ns,nl, ntot
			real(sp), intent(in), dimension(:) :: lambda_read_s, lambda_read_l
			real(sp), intent(in) :: lambda_s_low, lambda_s_high, &
									lambda_l_low, lambda_l_high
			real(sp), intent(inout), dimension(:), allocatable :: lambda, &
						lambda_low, lambda_high, delta_lambda, b_s_g, ext_s_g, sflux_l
			real(sp), intent(in), dimension(1-l_h:kp+r_h) :: rhoan, thetan
			real(sp), dimension(:,:,:,:), allocatable, intent(inout) :: flux_d,flux_u
			real(sp), dimension(:,:,:), allocatable, intent(inout) :: rad_power
			integer(i4b), intent(inout) :: nprocv,tdstart,tdend
			real(sp), dimension(:), allocatable, intent(inout) :: mvrecv,a,b,c,r,u
			integer(i4b), dimension(3), intent(inout) :: coords
			integer(i4b), dimension(3), intent(in) :: dims
			integer(i4b), intent(in) :: id, comm3d
			
		
			real(sp) :: sum_upper, sum_lower, x_upper, x_lower, fac
			integer(i4b) :: AllocateStatus, i, j,doy, error
			logical :: leap
			
! if the pe is not being used in the cartesian topology, do not use here
			if(id>=dims(1)*dims(2)*dims(3)) return 
		
				
		    ! find out how many procs in this communicator
		    call MPI_Comm_size(comm3d, nprocv,error)
            
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! array bounds for tridag solver                                             !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            tdend=kp*2
            if(coords(3)==0) then 
                ! bottom level: 1 is for lower BC
                tdstart=2
                tdend=tdend+1
            else 
                tdstart=1
            endif
            
            if(coords(3)==(dims(3)-1)) then
                tdend=tdend+1
            endif
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            

			! some calcs
			! set the lat / lon
			lat=lat_ref
			lon=lon_ref
			! find the surface emissivity
			do i=1,nemiss
				if(emiss_r == trim(rsv%emiss_des(i))) then
					exit
				endif
			enddo
			emiss=rsv%emiss(i)
			
			! find the surface albedo
			do i=1,nalb
				if(albedo_r == trim(rsv%albedo_des(i))) then
					exit
				endif
			enddo
			alb=rsv%albedo(i)
			
			
			
			! calculate days in year:
			call leap_year_calc(start_year,diy_s,leap)
			call doy_calc(start_year,start_mon,start_day,diy_s,doy_s, leap)
			
			! number of wave bands for radiation:
			nl=nl_r
			ns=ns_r
			ntot=ns+nl
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! allocate memory
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			allocate( a(1:tdend), STAT = AllocateStatus)
			if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
			allocate( b(1:tdend), STAT = AllocateStatus)
			if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
			allocate( c(1:tdend), STAT = AllocateStatus)
			if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
			allocate( r(1:tdend), STAT = AllocateStatus)
			if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
			allocate( u(1:tdend), STAT = AllocateStatus)
			if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
			allocate( mvrecv(1:nprocv), STAT = AllocateStatus)
			if (AllocateStatus /= 0) STOP "*** Not enough memory ***"

			allocate( lambda(1:ntot), STAT = AllocateStatus)
			if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
			allocate( lambda_low(1:ntot), STAT = AllocateStatus)
			if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
			allocate( lambda_high(1:ntot), STAT = AllocateStatus)
			if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
			allocate( delta_lambda(1:ntot), STAT = AllocateStatus)
			if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
			allocate( sflux_l(1:ntot), STAT = AllocateStatus)
			if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
			allocate( b_s_g(1:ntot), STAT = AllocateStatus)
			if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
			allocate( ext_s_g(1:kp), STAT = AllocateStatus)
			if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
			allocate( flux_d(1-r_h:kp+r_h,1-r_h:jp+r_h,1-r_h:ip+r_h,1:ntot), &
									STAT = AllocateStatus)
			if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
			allocate( flux_u(1-r_h:kp+r_h,1-r_h:jp+r_h,1-r_h:ip+r_h,1:ntot), &
									STAT = AllocateStatus)
			if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
			allocate( rad_power(1-r_h:kp+r_h,1-r_h:jp+r_h,1-r_h:ip+r_h), &
									STAT = AllocateStatus)
			if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
			

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! bands that we are solving for radiative transfer in:						 !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! short wave:
			if (ns > 0) then
				lambda(1:ns)=lambda_read_s(1:ns)
				lambda_low(1)=lambda_s_low
				lambda_low(2:ns)=(lambda(2:ns)+lambda(1:ns-1)) / 2._sp
				lambda_high(ns)=lambda_s_high
				lambda_high(1:ns-1)=(lambda(2:ns)+lambda(1:ns-1)) / 2._sp
				delta_lambda(1:ns) = lambda_high(1:ns)-lambda_low(1:ns)
			endif
			! long wave:
			if (nl > 0) then
				lambda(ns+1:ntot)=lambda_read_l(1:nl)
				lambda_low(ns+1)=lambda_l_low
				lambda_low(ns+2:ntot)=(lambda(ns+2:ntot)+lambda(ns+1:ntot-1)) / 2._sp
				lambda_high(ntot)=lambda_l_high
				lambda_high(ns+1:ntot-1)=(lambda(ns+2:ntot)+lambda(ns+1:ntot-1)) / 2._sp
				delta_lambda(ns+1:ntot) = lambda_high(ns+1:ntot)-lambda_low(ns+1:ntot)
			endif
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



			
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! calculate the planck function in these bins, and scale to TOA				 !
			! see http://www.spectralcalc.com/blackbody/inband_radiance.html             !
			! note, the integral result is the same no matter form of BB function        !
			! but the value of x is what is in the exponential                           !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			fac=fac_planck / sigma_sb * sflux
			do i=1,ntot
				x_upper=h_planck*c_light/(k_boltz*temp_sun*lambda_high(i))
				x_lower=h_planck*c_light/(k_boltz*temp_sun*lambda_low(i))
				sum_upper=0._sp
				do j=1,min(512 , 2+nint(20._sp/x_upper)) ! criteria on webpage
					sum_upper=sum_upper+(x_upper**3./real(j,sp) + &
							3._sp*x_upper**2./real(j,sp)**2 + &
							6._sp*x_upper/real(j,sp)**3 + &
							6._sp/real(j,sp)**4) * &
								exp(-real(j,sp)*x_upper)
				enddo
				
				sum_lower=0._sp
				do j=1,min(512 , 2+nint(20._sp/x_lower)) ! criteria on webpage
					sum_lower=sum_lower+(x_lower**3./real(j,sp) + &
							3._sp*x_lower**2./real(j,sp)**2 + &
							6._sp*x_lower/real(j,sp)**3 + &
							6._sp/real(j,sp)**4) * &
								exp(-real(j,sp)*x_lower)
				enddo
				sflux_l(i) = fac*(sum_upper-sum_lower)
			enddo
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


			! calculate scattering cross section due to rayleigh:
			do i=1,ntot
				call scattering_cs_rayleigh(lambda(i),b_s_g(i))
			enddo
	
			! now:
			! 3. asymmetry parameters, for aerosol, and hydrometeors - read
			! 4. scattering cross sections for aerosol, absorption cross sections?
			! 5. Tyndall, Mie, Geometric
			! 6. Optical depth

		end subroutine allocate_and_set_radiation
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!>@author
		!>Paul J. Connolly, The University of Manchester
		!>@brief
		!>calculate energy in a photon of wavelength lambda
		!>@param[in] lambda
		function e_photon(lambda)
			implicit none
			real(sp), intent(in) :: lambda
			real(sp) :: e_photon
		
			e_photon=h_planck*c_light / lambda
	
		end function e_photon
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!>@author
		!>Paul J. Connolly, The University of Manchester
		!>@brief
		!>calculate scattering cross section due to Rayleigh scattering
		!>@param[in] lambda
		!>@param[inout] b_s_g_lambda: the scattering cross section
		subroutine scattering_cs_rayleigh(lambda,b_s_g_lambda)
			implicit none
			real(sp), intent(in) :: lambda
			real(sp), intent(inout) :: b_s_g_lambda
			real(sp) :: n_a_l, delta, f_delta
		
			delta = 0.0279_sp ! Young (1980)
			! anisotropic correction term (in equation 9.36):
			f_delta = (6._sp+3._sp*delta) / (6._sp-7._sp*delta)
			
			
			call real_refractive_air(lambda,n_a_l)

			! scattering cross section due to Rayleigh scattering
			b_s_g_lambda=8._sp*pi**3._sp*(n_a_l**2._sp -1._sp)**2._sp &
						/ (3._sp*lambda**4._sp * N_a_0**2._sp) &
						* f_delta
		end subroutine scattering_cs_rayleigh
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!>@author
		!>Paul J. Connolly, The University of Manchester
		!>@brief
		!>radiant intensity, or radiance from Planck's law
		!>@param[in] lambda, T
		!>@param[inout] b_lambda_t
		subroutine plancks_law(lambda,t,b_lambda_t)
			implicit none
			real(sp), intent(in) :: lambda,t
			real(sp), intent(inout) :: b_lambda_t
		
			b_lambda_t = 2._sp*h_planck*c_light**2 / (lambda**5 * &
						(exp(h_planck*c_light/(lambda*k_boltz*T))-1._sp) )
	
		end subroutine plancks_law
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!>@author
		!>Paul J. Connolly, The University of Manchester
		!>@brief
		!>real refractive index in air vs wavelength
		!>Edlen (1966)
		!>@param[in] lambda
		!>@param[inout] na
		subroutine real_refractive_air(lambda,na)
			implicit none
			real(sp), intent(in) :: lambda
			real(sp), intent(inout) :: na
			
			real(sp) :: lambda_mu
			
			lambda_mu = lambda*1.e6_sp
		
			! Edlen, 1966:
! 			na=1._sp+1e-8_sp*( 8342.13_sp+(2406030._sp/(130._sp-lambda_mu**(-2))) + &
! 						(15997._sp/(38.9_sp-lambda_mu**(-2)))	)
			! Peck and Reeder, 1972:
			na=1._sp+1.e-8_sp*( 8060.51_sp+(2480990._sp/(132.274_sp-lambda_mu**(-2))) + &
						(17455.7_sp/(39.32957_sp-lambda_mu**(-2)))	)
	
		end subroutine real_refractive_air
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!>@author
		!>Paul J. Connolly, The University of Manchester
		!>@brief
		!>Calculate if current year is leap year or not
		!>See https://support.microsoft.com/en-gb/help/214019/
		!            method-to-determine-whether-a-year-is-a-leap-year
		!>@param[in] year
		!>@param[inout] diy, leap
		subroutine leap_year_calc(year,diy,leap)
			implicit none
			integer(i4b), intent(in) :: year
			integer(i4b), intent(inout) :: diy
			logical, intent(inout) :: leap
			
			integer(i4b) :: i
			
			if(modulo(year,4) == 0) then ! step 1:
				if(modulo(year,100) == 0) then ! step 2:
					if(modulo(year,400) == 0) then ! step 3:
						! step 4:
						leap=.true. 
					else
						! step 5:
						leap=.false.
					endif
				else
					! step 4:
					leap=.true.
				endif
			else ! step 5:
				leap=.false. 
			endif
					
			if(leap) then
				diy=366
			else
				diy=365
			endif
			
		end subroutine leap_year_calc	
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!>@author
		!>Paul J. Connolly, The University of Manchester
		!>@brief
		!>Calculate if current year is leap year or not
		!>See https://support.microsoft.com/en-gb/help/214019/
		!            method-to-determine-whether-a-year-is-a-leap-year
		!>@param[in] year
		!>@param[inout] diy, leap
		subroutine doy_calc(year,mon,day,diy,doy, leap)
			implicit none
			integer(i4b), intent(in) :: year,mon,day, diy
			integer(i4b), intent(inout) :: doy
			logical, intent(in) :: leap
			
			integer(i4b) :: i
			
			! calculates the day of year and whether it is a leap year / day in year
			doy=0
			if (leap ) then
				do i=1,mon-1
					doy=doy+di_ml(i)
				enddo
			else
				do i=1,mon-1
					doy=doy+di_mn(i)
				enddo
			endif
			doy=doy+day
			
			
			
		end subroutine doy_calc	
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!>@author
		!>Paul J. Connolly, The University of Manchester
		!>@brief
		!>Determine current year, mon, day, hour, min, sec
		!>@param[in] year
		!>@param[inout] diy, leap
		subroutine determine_current_time(start_year,start_hour,&
										start_min, start_sec, time, doy, tod,diy, &
										year)
			implicit none
			integer(i4b), intent(in) :: start_year, start_hour, &
										start_min, start_sec
			real(sp), intent(in) :: time
			real(sp), intent(inout) :: tod
			integer(i4b), intent(inout) :: doy,diy, year
			real(sp) :: days_to_add, doy_c
			integer(i4b) :: diy_c, i
			logical :: leap
			
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! doy                                                                        !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			days_to_add=time / 86400._sp
			doy_c=real(doy_s,sp)+real(start_hour,sp)/24._sp+ &
					real(start_min,sp)/3600._sp+real(start_sec,sp)/86400._sp+&
					days_to_add
			call leap_year_calc(start_year,diy_c,leap)

			i=1
			do while (doy_c.gt.real(diy_c,sp))
				doy_c=doy_c-real(diy_c,sp)
				call leap_year_calc(start_year+i,diy_c,leap)
				i=i+1
			end do
			year=start_year+i-1
			doy=floor(doy_c)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			tod = (doy_c-real(doy,sp))*86400._sp
			diy=diy_c
			
		end subroutine determine_current_time
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!										
										
										
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!>@author
		!>Paul J. Connolly, The University of Manchester
		!>@brief
		!>TOA irradiance
		!>See Jacobson, Fundamentals of Atmospheric Modelling, pp 325
		!>@param[in] day of year - doy
		!>@param[inout] frac - fraction of solar constant
		subroutine toa_solar_frac(doy,frac)
			implicit none
			real(sp), intent(in) :: doy
			real(sp), intent(inout) :: frac
			real(sp) :: th_j
			
			! equation 9.97 of Jacobson's book.
			th_j=2._sp*pi*(doy)/diy_exact ! doy=1 corresponds to 1st Jan
										  ! and not 22 dec.
										  ! note, diy is down as 365.25..., so not 100% 
										  ! accurate
			frac=1.00011_sp + 0.034221_sp*cos(th_j) + 0.00128_sp*sin(th_j) + &
					0.000719_sp*cos(2._sp*th_j)+0.000077_sp*sin(2._sp*th_j)
	
		end subroutine toa_solar_frac	
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!>@author
		!>Paul J. Connolly, The University of Manchester
		!>@brief
		!>solar zenith angle
		!>See Jacobson, Fundamentals of Atmospheric Modelling, pp 317
		!>an over simplification, but ok for non climate applications
		!>@param[inout] cosine of solar zenith angle
		!>@param[in] julian day of year - d_j 
		!            https://landweb.nascom.nasa.gov/browse/calendar.html
		!>@param[in] year - year
		!>@param[in] longitude, latitude
		!>@param[in] time - seconds in day (is adjusted to be past noon, in this func)
		subroutine solar_zenith(cos_theta_s,d_j,year, longitude, latitude, time)
			implicit none
			real(sp), intent(inout) :: cos_theta_s
			real(sp), intent(in) :: d_j, longitude, latitude
			real(sp), intent(in) :: year, time
			real(sp) :: th_j, d_l,n_jd, eps_ob, l_m, g_m, lambda_ec, time_local, ha_r, &
						delta_r

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! Obliquity of the ecliptic                                                  !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! equation 9.70 of Jacobson's book.
			if(year >= 2001) then
				d_l=int(year-2001._sp) / 4
			else
				d_l=int(year-2001._sp) / 4 - 1			
			endif		
			n_jd = 364.5_sp	+ (year-2001)*365._sp + d_l + d_j
			! equation 9.69 of Jacobson's book.
			! check this is correct - wikipedia different to Jacobson?:
			eps_ob = (23.439_sp-0.0000004_sp*n_jd)*pi/180._sp  
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! ecliptic longitude of the Sun                                              !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! equation 9.72 of Jacobson's book.
			l_m = (280.460_sp + 0.9856474_sp*n_jd)
			g_m = (357.528_sp + 0.9856003_sp*n_jd)*pi/180._sp
			! equation 9.71 of Jacobson
			lambda_ec = (l_m + 1.915_sp*sin(g_m) + 0.020_sp*sin(2._sp*g_m))*pi/180._sp
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! local hour angle                                              			 !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! longitude varies from -180 to 180
			time_local=time +longitude/180._sp*86400._sp
			! equation 9.73 of Jacobson's book.
			ha_r = 2._sp*pi*(time_local-86400._sp/2._sp)/86400._sp
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! solar declination angle                                              	     !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! equation 9.68 of Jacobson's book.
			delta_r=asin(sin(eps_ob)*sin(lambda_ec))
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			! equation 9.67 of Jacobson's book.
			cos_theta_s = sin(latitude*pi/180._sp)* &
						sin(delta_r)+cos(latitude*pi/180._sp)*cos(delta_r)*cos(ha_r)
	
		end subroutine solar_zenith
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!>@author
		!>Paul J. Connolly, The University of Manchester
		!>@brief
		!>calculate the in-band planck function at temperature t
		!>See http://www.spectralcalc.com/blackbody/inband_radiance.html 
		!>@param[inout] cosine of solar zenith angle
		!>@param[in] year - year
		!>@param[in] longitude, latitude
		!>@param[in] time - seconds in day (is adjusted to be past noon, in this func)
		subroutine inband_planck(nbands,lambda_low,lambda_high,temp,blt)
			implicit none
			integer(i4b), intent(in) :: nbands 
			real(sp), dimension(nbands), intent(inout) :: blt
			real(sp), dimension(nbands), intent(in) :: lambda_low, lambda_high
			real(sp), intent(in) :: temp
			
			integer(i4b) :: i,j
			real(sp) :: x_upper, x_lower, sum_upper, sum_lower, fac
			
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! calculate the planck function in these bins				 				 !
			! see http://www.spectralcalc.com/blackbody/inband_radiance.html             !
			! note, the integral result is the same no matter form of BB function        !
			! but the value of x is what is in the exponential                           !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			fac=fac_planck * temp**4
			do i=1,nbands
				x_upper=h_planck*c_light/(k_boltz*temp*lambda_high(i))
				x_lower=h_planck*c_light/(k_boltz*temp*lambda_low(i))
				sum_upper=0._sp
				do j=1,min(512 , 2+nint(20._sp/x_upper)) ! criteria on webpage
					sum_upper=sum_upper+(x_upper**3./real(j,sp) + &
							3._sp*x_upper**2./real(j,sp)**2 + &
							6._sp*x_upper/real(j,sp)**3 + &
							6._sp/real(j,sp)**4) * &
								exp(-real(j,sp)*x_upper)
				enddo
				
				sum_lower=0._sp
				do j=1,min(512 , 2+nint(20._sp/x_lower)) ! criteria on webpage
					sum_lower=sum_lower+(x_lower**3./real(j,sp) + &
							3._sp*x_lower**2./real(j,sp)**2 + &
							6._sp*x_lower/real(j,sp)**3 + &
							6._sp/real(j,sp)**4) * &
								exp(-real(j,sp)*x_lower)
				enddo
				blt(i) = fac*(sum_upper-sum_lower)
				
! 				call plancks_law(0.5_sp*(lambda_low(i)+lambda_high(i)),temp,blt(i))
! 				blt(i)=blt(i)*(lambda_high(i)-lambda_low(i))
			enddo
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		end subroutine inband_planck
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
			

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!>@author
		!>Paul J. Connolly, The University of Manchester
		!>@brief
		!> solve fluxes, using
		!>tridiagonal solver for two-stream method
		!> \f$ \begin{bmatrix} 
        !>        1  & -A_e & 0 & 0 & 0 & 0 & \cdots\\
        !>        A_k & B_k & C_k & 0 & 0 & 0 & \cdots\\
        !>        0 & E_k & F_k & G_k & 0 & 0  & \cdots\\
        !>        0 & 0 & A_{k+1} & B_{k+1} & C_{k+1} & 0  & \cdots\\
        !>        \vdots & \vdots & \vdots & \vdots & \vdots & \vdots & \vdots \\
        !>        0 & 0 & 0 & 0 &0 & 0 & 1   \\
        !>   \end{bmatrix}
        !>   \begin{bmatrix}
        !>   F^+_{k-{1/2}} \\F^-_{k-{1/2}} \\F^+_{k+{1/2}} \\
        !>   F^-_{k+{1/2}} \\\vdots \\F^-_{kp+{1/2}} \\
        !>   \end{bmatrix} = 
        !>   \begin{bmatrix}
        !>   A_e\mu_sF_s\exp\left(-\tau / \mu_s \right) +\epsilon \pi B_T\\
        !>   D_k \\H_k \\
        !>   D_{k+1} \\ \vdots \\\mu F_s \\
        !>   \end{bmatrix} 
        !>   \f$
        !>@param[in] comm3d,id, dims, coords: mpi stuff
		subroutine solve_fluxes(start_year,start_mon,start_day, &
								start_hour,start_min,start_sec, &
								lat, lon, &
								time, nbands,ns,nl,ip,jp,kp,r_h, &
								tdstart,tdend,a,b,c,r,u, &
								lambda_low, lambda_high, lambda,b_s_g, sflux_l,&
								rhoan, thetan, dz,dzn, albedo, emiss, &
								quad_flag, flux_u, flux_d, &
								nprocv,mvrecv, &
								coords,dims, id, comm3d)
			use nr, only : tridag
			use mpi
			use pts
			implicit none
			integer(i4b), intent(in) :: start_year, start_mon, start_day, start_hour, &
							start_min, start_sec
			real(sp), intent(in) :: time
			integer(i4b), intent(in) :: nbands, ns, nl, ip, jp, kp, r_h, nprocv
			real(sp), intent(in), dimension(nprocv) :: mvrecv
			real(sp), dimension(nbands), intent(in) :: &
					lambda_low, lambda_high, lambda, b_s_g, sflux_l
			real(sp), intent(in), dimension(1-r_h:kp+r_h) :: rhoan, thetan, dz,dzn
			real(sp), intent(inout), &
					dimension(1-r_h:kp+r_h, 1-r_h:jp+r_h, 1-r_h:ip+r_h,1:nbands) :: &
					flux_u, flux_d
			real(sp), intent(in) :: lat, lon, albedo, emiss
			integer(i4b), intent(in) :: quad_flag	

            integer(i4b), dimension(3), intent(in) :: coords
            integer(i4b), dimension(3), intent(in) :: dims
            integer(i4b), intent(in) :: id, comm3d
            
            integer(i4b), intent(in) :: tdstart,tdend
            real(sp), dimension(1:tdend), intent(inout) :: a,b,c,r,u

			
			! local variables:
			integer(i4b) :: doy,diy, year,error
			real(sp) :: tod, cos_theta_s, frac,msend
			real(sp), dimension(0:kp+1) :: sigma_ray, taun,dtau,dtaun
			real(sp), dimension(0:kp+1) :: tau,direct
			real(sp), dimension(kp+1,nbands) :: blt
			logical :: leap
			real(sp) :: omg_s=1._sp, ga=0._sp, gamma1, gamma2, gamma3
			
			integer(i4b) :: it, itermax=1, k, m,j,i,kend,kstart
			real(sp) :: Ac, Bc, Cc, Dc, alp, bet, gam, del, &
						Ak, Bk, Ck, Dk, Ek, Fk, Gk, Hk, Skp, Skm
			!real(sp), dimension(2*(kp+1)) :: b,r,u
			!real(sp), dimension(2*kp+1) :: a,c
			
			
			
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! determine parameters for sun position                                      !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			call determine_current_time(start_year,start_hour,&
						start_min, start_sec, time, doy, tod,diy, year)	
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! determine toa, solar zenith:                                               !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			call toa_solar_frac(real(doy,sp),frac)
			call solar_zenith(cos_theta_s,real(doy,sp),real(year,sp), lon, lat, tod)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! calculate planck function at all levels p-points:						     !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            do k=1,kp+1 
                call inband_planck(nbands,lambda_low,lambda_high,thetan(k),blt(k,:))
            enddo
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            if(dims(3)==(coords(3)+1)) then
                kend=kp
            else 
                kend=kp+1
            endif

			do m=1,nbands

				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! extinction:															 !
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				
				do k=1,kend 
					! number of molecules per cubic metre of air multiplied by 
					! collision cross-section:
					sigma_ray(k)=rhoan(k)*navog/ma*b_s_g(m)
				enddo
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
			
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! optical depth:														 !
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! sigma_ray is the extinction on p-points
				! tau is optical depth on w-points
				! mpi - add bottom taus to all and add cumulative
				tau(kp+1)=0._sp
				tau(kp)=0._sp ! for top level tau(kp) needs to be defined taking into 
				              ! account the atmosphere above domain
				do k=kend-1,0,-1 
					! cumulative extinction x dz
					tau(k)=tau(k+1)+sigma_ray(k+1)*dz(k) 
				enddo
				
				
				
			
			    ! for parallel solver tau and dtau need some parallel communication:
                ! http://mpitutorial.com/tutorials/mpi-scatter-gather-and-allgather/
                msend=tau(1)  ! message to send
                call MPI_allgather(msend,1,MPI_REAL8,mvrecv,1,MPI_REAL8,comm3d,error)
			    tau=tau+sum(mvrecv(  coords(3)+2:dims(3)  ))
			    
			    
			    
			    
			    ! change in optical depth on w-points
				do k=kp,0,-1 
				    dtau(k)=(tau(k+1)-tau(k) ) ! tau = 0 at top, so set to negative 
				                                ! to get derivative in correct direction
				enddo
                msend=dtau(1)  ! message to send
                call MPI_allgather(msend,1,MPI_REAL8,mvrecv,1,MPI_REAL8,comm3d,error)
                if((coords(3)+1).lt.(dims(3))) dtau(kp+1)=mvrecv(coords(3)+2)
				if(coords(3)==(dims(3)-1)) dtau(kp)=0._sp
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! store direct radiation                                                 !
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                direct=0._sp
                if(cos_theta_s > 0._sp) then
                    do k=0,kp
                        direct(k)=cos_theta_s*sflux_l(m)* &
                                exp(-0.5_sp*(tau(k+1)+tau(k))/cos_theta_s)
                    enddo			
                endif
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
                
                
                
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! diffuse phase function scaling on p-points:							 !
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				select case (quad_flag)
					case(0) ! quadrature
						gamma1=(1._sp-omg_s*(1._sp+ga)/2._sp)/mu1
						gamma2=omg_s*(1._sp-ga)/(2._sp*mu1)
						gamma3=(1._sp-3._sp*ga*mu1*cos_theta_s)/2._sp
					case(1) ! eddington
						gamma1=(7._sp-omg_s*(4._sp+3._sp*ga))/4._sp
						gamma2=-(1._sp-omg_s*(4._sp-3._sp*ga))/4._sp
						gamma3=(2._sp-3._sp*ga*cos_theta_s)/4._sp				
					case default
						print *,'not coded'
						stop
				end select
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			    
				
			
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! set up tridiagonal problem:											 !
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! dtau(0,kp+1) needs to be set
				! tau(0,kp+1) needs to be set				
				do k=1,kp
					! define coefficients for FDE:
					Ac=1._sp-0.5_sp*dtau(k)*gamma1
					Bc=-1._sp-0.5_sp*dtau(k)*gamma1
					Cc=0.5_sp*dtau(k)*gamma2
					Dc=Cc
					
					alp=-Cc
					bet=-Cc
					gam=-Bc
					del=-Ac
				
					! matrix coefficients:
					Ak=(Bc-bet*Cc/gam)
					Bk=(Dc-del*Cc/gam)
					Ck=(Ac-alp*Cc/gam)
					Skp=0._sp
					Skm=0._sp
					if(cos_theta_s > 0._sp) then
					    ! short-wave (see equations 9.121 Jacobson):
						Skp=-dtau(k)*gamma3*omg_s*sflux_l(m)* &
						    exp(-0.5_sp*(tau(k+1)+tau(k))/cos_theta_s)
						Skm=dtau(k)*(1._sp-gamma3)* &
								omg_s*sflux_l(m)*&
								exp(-0.5_sp*(tau(k+1)+tau(k))/cos_theta_s)
					endif
					! long wave (see equations 9.122 Jacobson):
					Skp=Skp-dtau(k)*2._sp*pi*(1._sp-omg_s)*blt(k+1,m)
					Skm=Skm+dtau(k)*2._sp*pi*(1._sp-omg_s)*blt(k+1,m)

					Dk=Skp-Cc/gam*Skm
					Ek=(del-Dc*bet/Bc)
					Fk=(alp-Ac*bet/Bc)
					Gk=(gam-Cc*bet/Bc)
					Hk=Skm-bet*Skp/Bc
						
						
						
					! lower diag:
					! a(1) needs to be set when not bottom
					if(coords(3)==0) then
                        a(2*k)     = Ak !2,4,6,8
                        a(2*k+1)   = Ek !3,5,7,9
                    else
                        a(2*k-1) = Ak !1,3,5,7
                        a(2*k)   = Ek !2,4,6,8
                    endif				


					! upper diag:
					! c(kpp) gets set to zero on top PE
					if(coords(3)==0) then
                        c(2*k) = Ck !2,4,6
                        c(2*k+1)   = Gk !3,5,7
                    else
                        c(2*k-1) = Ck !1,3,5
                        c(2*k)   = Gk !2,4,6                    
                    endif				
				
				
				
				
					! diagonal:
					if(coords(3)==0) then
                        b(2*k) = Bk !2,4,6
                        b(2*k+1)   = Fk !3,5,7
                    else
                        b(2*k-1) = Bk !1,3,5
                        b(2*k)   = Fk !2,4,6
                    endif                

					! rhs:
					if(coords(3)==0) then
                        r(2*k) = Dk !2,4,6
                        r(2*k+1)   = Hk !3,5,7
                    else
                        r(2*k-1) = Dk !1,3,5
                        r(2*k)   = Hk !2,4,6
                    endif
				enddo
				
				
				
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! boundary conditions
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				if(coords(3)==0) then
				    a(1)=0._sp ! this is specific to parallel
                    b(1)=1._sp ! diagonal, so fine
                    c(1)=-albedo ! ! c(1) has to be set in both serial and parallel
    				r(1)=0._sp
                    if(cos_theta_s > 0._sp) then
                        r(1)=albedo*cos_theta_s*sflux_l(m)*exp(-tau(0)/cos_theta_s)
                    endif
    				r(1)=r(1)+emiss*pi*blt(1,m)
				endif
				
                
                if(coords(3)==(dims(3)-1)) then
                    a(tdend)=0._sp ! no contribution from upward flux in last equation
    				b(tdend)=1._sp
    				r(tdend)=0._sp ! toa downward flux at night
    				c(tdend)=0._sp ! this is specific to parallel
                    if(cos_theta_s > 0._sp) then
                        !r(tdend)=cos_theta_s*sflux_l(m) ! toa downward flux - day
                                                        ! should be multiplied by tau(kp) 
                                                        ! for atmosphere above
                    endif
                endif                
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				
				
				
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! solve tridiagonal matrix:												 !
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				call parallel_tridag(a,b,c,r,u,tdend,coords,dims, id, comm3d)
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
                
                
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! assumes temperature is given by reference state:						 !
				! message passing done out of this routine                               !
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                kstart=1
                if(coords(3)==0) then
                    do i=1,ip	
                        do j=1,jp
                            ! flux_d(kp) recv
                            do k=1,kp
                                ! sets the surface fluxes too
                                flux_u(k-kstart,j,i,m)=u(2*(k)-1)
                                flux_d(k-kstart,j,i,m)=u(2*k)
                            enddo
                            flux_u(kp+1-kstart,j,i,m)=u(2*(kp+1)-1)
                            
                        enddo
                    enddo
                else if (coords(3)==(dims(3)-1)) then
                    do i=1,ip	
                        do j=1,jp
                            flux_d(0,j,i,m)=u(1) ! needs to be passed to below
                            do k=1,kp
                                flux_u(k,j,i,m)=u(2*(k)) ! 2,4,6
                                flux_d(k,j,i,m)=u(2*k+1) ! 3,5,7
                            enddo
                        enddo
                    enddo
                else
                    do i=1,ip	
                        do j=1,jp
                            flux_d(0,j,i,m)=u(1) ! needs to be passed to below
                            ! mpi_issend to one below
                            ! mpi_irecv from above flux_d(kp,j,i,m)=
                            do k=1,kp-1
                                flux_u(k,j,i,m)=u(2*(k))
                                flux_d(k,j,i,m)=u(2*k+1)
                            enddo
                            flux_u(kp,j,i,m)=u(2*kp)
                        enddo
                    enddo
                endif
                if((coords(3)==0).and.(coords(3)==(dims(3)-1))) then
                    do i=1,ip	
                        do j=1,jp
                            flux_d(kp+1-kstart,j,i,m)=u(2*(kp+1))                            
                        enddo
                    enddo
                endif
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if(cos_theta_s > 0._sp) then
                    do i=1,ip	
                        do j=1,jp
                            do k=0,kp
                                flux_d(k,j,i,m)=  flux_d(k,j,i,m)+direct(k)
                            enddo                        
                        enddo
                    enddo
                endif
				
			enddo
			
									
		end subroutine solve_fluxes
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	end module radiation
	