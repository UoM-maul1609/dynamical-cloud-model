	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>radiation solver routines 
    module radiation
    use nrtype
    use nr, only : locate, polint, qromb
    implicit none


		!>@brief
		!>main model prognostic variables for radiation
        type radg
        	integer(i4b) :: year, mon, day, hour, sec
        	integer(i4b) :: doy, diy
        	integer(i4b) :: nprocv
        	real(sp) :: albedo, emiss, lat, lon
        	integer(i4b) :: ns, nl, ntot
            real(sp), dimension(:), allocatable ::	b_s_g, & ! cross section for Rayleigh
            	lambda, lambda_low, lambda_high, delta_lambda, & 
            	nrwbin,niwbin ! bin averaged refractive indices
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
        	real(sp) :: asymmetry_water
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
		! number of entries in LUT for refrative index of liquid water
		integer(i4b), parameter :: n_rfrl=169
		real(sp), dimension(n_rfrl) :: lam_h2o=[0.200e-6_sp,0.225e-6_sp,0.250e-6_sp,&
		                    0.275e-6_sp,0.300e-6_sp,0.325e-6_sp,0.350e-6_sp,0.375e-6_sp,&
		                    0.400e-6_sp,0.425e-6_sp,0.450e-6_sp,0.475e-6_sp,0.500e-6_sp,&
		                    0.525e-6_sp,0.550e-6_sp,0.575e-6_sp,0.600e-6_sp,0.625e-6_sp,&
		                    0.650e-6_sp,0.675e-6_sp,0.700e-6_sp,0.725e-6_sp,0.750e-6_sp,&
		                    0.775e-6_sp,0.800e-6_sp,0.825e-6_sp,0.850e-6_sp,0.875e-6_sp,&
		                    0.900e-6_sp,0.925e-6_sp,0.950e-6_sp,0.975e-6_sp,1.000e-6_sp,&
		                    1.200e-6_sp,1.400e-6_sp,1.600e-6_sp,1.800e-6_sp,2.000e-6_sp,&
		                    2.200e-6_sp,2.400e-6_sp,2.600e-6_sp,2.650e-6_sp,2.700e-6_sp,&
		                    2.750e-6_sp,2.800e-6_sp,2.850e-6_sp,2.900e-6_sp,2.950e-6_sp,&
		                    3.000e-6_sp,3.050e-6_sp,3.100e-6_sp,3.150e-6_sp,3.200e-6_sp,&
		                    3.250e-6_sp,3.300e-6_sp,3.350e-6_sp,3.400e-6_sp,3.450e-6_sp,&
		                    3.500e-6_sp,3.600e-6_sp,3.700e-6_sp,3.800e-6_sp,3.900e-6_sp,&
		                    4.000e-6_sp,4.100e-6_sp,4.200e-6_sp,4.300e-6_sp,4.400e-6_sp,&
		                    4.500e-6_sp,4.600e-6_sp,4.700e-6_sp,4.800e-6_sp,4.900e-6_sp,&
		                    5.000e-6_sp,5.100e-6_sp,5.200e-6_sp,5.300e-6_sp,5.400e-6_sp,&
		                    5.500e-6_sp,5.600e-6_sp,5.700e-6_sp,5.800e-6_sp,5.900e-6_sp,&
		                    6.000e-6_sp,6.100e-6_sp,6.200e-6_sp,6.300e-6_sp,6.400e-6_sp,&
		                    6.500e-6_sp,6.600e-6_sp,6.700e-6_sp,6.800e-6_sp,6.900e-6_sp,&
		                    7.000e-6_sp,7.100e-6_sp,7.200e-6_sp,7.300e-6_sp,7.400e-6_sp,&
		                    7.500e-6_sp,7.600e-6_sp,7.700e-6_sp,7.800e-6_sp,7.900e-6_sp,&
		                    8.000e-6_sp,8.200e-6_sp,8.400e-6_sp,8.600e-6_sp,8.800e-6_sp,&
		                    9.000e-6_sp,9.200e-6_sp,9.400e-6_sp,9.600e-6_sp,9.800e-6_sp,&
		                    10.000e-6_sp,10.500e-6_sp,11.000e-6_sp,11.500e-6_sp,&
		                    12.000e-6_sp,12.500e-6_sp,13.000e-6_sp,13.500e-6_sp,&
		                    14.000e-6_sp,14.500e-6_sp,15.000e-6_sp,15.500e-6_sp,&
		                    16.000e-6_sp,16.500e-6_sp,17.000e-6_sp,17.500e-6_sp,&
		                    18.000e-6_sp,18.500e-6_sp,19.000e-6_sp,19.500e-6_sp,&
		                    20.000e-6_sp,21.000e-6_sp,22.000e-6_sp,23.000e-6_sp,&
		                    24.000e-6_sp,25.000e-6_sp,26.000e-6_sp,27.000e-6_sp,&
		                    28.000e-6_sp,29.000e-6_sp,30.000e-6_sp,32.000e-6_sp,&
		                    34.000e-6_sp,36.000e-6_sp,38.000e-6_sp,40.000e-6_sp,&
		                    42.000e-6_sp,44.000e-6_sp,46.000e-6_sp,48.000e-6_sp,&
		                    50.000e-6_sp,60.000e-6_sp,70.000e-6_sp,80.000e-6_sp,&
		                    90.000e-6_sp,100.000e-6_sp,110.000e-6_sp,120.000e-6_sp,&
		                    130.000e-6_sp,140.000e-6_sp,150.000e-6_sp,160.000e-6_sp,&
		                    170.000e-6_sp,180.000e-6_sp,190.000e-6_sp,200.000_sp]

		real(sp), dimension(n_rfrl) :: nr_h2o=[1.396_sp,1.373_sp,1.362_sp,1.354_sp,&
		    1.349_sp,1.346_sp,1.343_sp,1.341_sp,1.339_sp,1.338_sp,1.337_sp,1.336_sp,&
		    1.335_sp,1.334_sp,1.333_sp,1.333_sp,1.332_sp,1.332_sp,1.331_sp,1.331_sp,&
		    1.331_sp,1.33_sp,1.33_sp,1.33_sp,1.329_sp,1.329_sp,1.329_sp,1.328_sp,&
		    1.328_sp,1.328_sp,1.327_sp,1.327_sp,1.327_sp,1.324_sp,1.321_sp,1.317_sp,&
		    1.312_sp,1.306_sp,1.296_sp,1.279_sp,1.242_sp,1.219_sp,1.188_sp,1.157_sp,&
		    1.142_sp,1.149_sp,1.201_sp,1.292_sp,1.371_sp,1.426_sp,1.467_sp,1.483_sp,&
		    1.478_sp,1.467_sp,1.45_sp,1.432_sp,1.42_sp,1.41_sp,1.4_sp,1.385_sp,1.374_sp,&
		    1.364_sp,1.357_sp,1.351_sp,1.346_sp,1.342_sp,1.338_sp,1.334_sp,1.332_sp,&
		    1.33_sp,1.33_sp,1.33_sp,1.328_sp,1.325_sp,1.322_sp,1.317_sp,1.312_sp,&
		    1.305_sp,1.298_sp,1.289_sp,1.277_sp,1.262_sp,1.248_sp,1.265_sp,1.319_sp,&
		    1.363_sp,1.357_sp,1.347_sp,1.339_sp,1.334_sp,1.329_sp,1.324_sp,1.321_sp,&
		    1.317_sp,1.314_sp,1.312_sp,1.309_sp,1.307_sp,1.304_sp,1.302_sp,1.299_sp,&
		    1.297_sp,1.294_sp,1.291_sp,1.286_sp,1.281_sp,1.275_sp,1.269_sp,1.262_sp,&
		    1.255_sp,1.247_sp,1.239_sp,1.229_sp,1.218_sp,1.185_sp,1.153_sp,1.126_sp,&
		    1.111_sp,1.123_sp,1.146_sp,1.177_sp,1.21_sp,1.241_sp,1.27_sp,1.297_sp,&
		    1.325_sp,1.351_sp,1.376_sp,1.401_sp,1.423_sp,1.443_sp,1.461_sp,1.476_sp,&
		    1.48_sp,1.487_sp,1.5_sp,1.511_sp,1.521_sp,1.531_sp,1.539_sp,1.545_sp,&
		    1.549_sp,1.551_sp,1.551_sp,1.546_sp,1.536_sp,1.527_sp,1.522_sp,1.519_sp,&
		    1.522_sp,1.53_sp,1.541_sp,1.555_sp,1.587_sp,1.703_sp,1.821_sp,1.886_sp,&
		    1.924_sp,1.957_sp,1.966_sp,2.004_sp,2.036_sp,2.056_sp,2.069_sp,2.081_sp,&
		    2.094_sp,2.107_sp,2.119_sp,2.13_sp]

		real(sp), dimension(n_rfrl) :: ni_h2o=[1.10E-07_sp,4.90E-08_sp,3.35E-08_sp,&
		        2.35E-08_sp,1.60E-08_sp,1.08E-08_sp,6.50E-09_sp,3.50E-09_sp,1.86E-09_sp,&
		        1.30E-09_sp,1.02E-09_sp,9.35E-10_sp,1.00E-09_sp,1.32E-09_sp,1.96E-09_sp,&
		        3.60E-09_sp,1.09E-08_sp,1.39E-08_sp,1.64E-08_sp,2.23E-08_sp,3.35E-08_sp,&
		        9.15E-08_sp,1.56E-07_sp,1.48E-07_sp,1.25E-07_sp,1.82E-07_sp,2.93E-07_sp,&
		        3.91E-07_sp,4.86E-07_sp,1.06E-06_sp,2.93E-06_sp,3.48E-06_sp,2.89E-06_sp,&
		        9.89E-06_sp,1.38E-04_sp,8.55E-05_sp,1.15E-04_sp,1.10E-03_sp,2.89E-04_sp,&
		        9.56E-04_sp,3.17E-03_sp,6.70E-03_sp,0.019_sp,0.059_sp,0.115_sp,0.185_sp,&
		        0.268_sp,0.298_sp,0.272_sp,0.24_sp,0.192_sp,0.135_sp,0.0924_sp,0.061_sp,&
		        0.0368_sp,0.0261_sp,0.0195_sp,0.0132_sp,0.0094_sp,0.00515_sp,0.0036_sp,&
		        0.0034_sp,0.0038_sp,0.0046_sp,0.00562_sp,0.00688_sp,0.00845_sp,0.0103_sp,&
		        0.0134_sp,0.0147_sp,0.0157_sp,0.015_sp,0.0137_sp,0.0124_sp,0.0111_sp,&
		        0.0101_sp,0.0098_sp,0.0103_sp,0.0116_sp,0.0142_sp,0.0203_sp,0.033_sp,&
		        0.0622_sp,0.107_sp,0.131_sp,0.088_sp,0.057_sp,0.0449_sp,0.0392_sp,&
		        0.0356_sp,0.0337_sp,0.0327_sp,0.0322_sp,0.032_sp,0.032_sp,0.0321_sp,&
		        0.0322_sp,0.0324_sp,0.0326_sp,0.0328_sp,0.0331_sp,0.0335_sp,0.0339_sp,&
		        0.0343_sp,0.0351_sp,0.0361_sp,0.0372_sp,0.0385_sp,0.0399_sp,0.0415_sp,&
		        0.0433_sp,0.0454_sp,0.0479_sp,0.0508_sp,0.0662_sp,0.0968_sp,0.142_sp,&
		        0.199_sp,0.259_sp,0.305_sp,0.343_sp,0.37_sp,0.388_sp,0.402_sp,0.414_sp,&
		        0.422_sp,0.428_sp,0.429_sp,0.429_sp,0.426_sp,0.421_sp,0.414_sp,0.404_sp,&
		        0.393_sp,0.382_sp,0.373_sp,0.367_sp,0.361_sp,0.356_sp,0.35_sp,0.344_sp,&
		        0.338_sp,0.333_sp,0.328_sp,0.324_sp,0.329_sp,0.343_sp,0.361_sp,0.385_sp,&
		        0.409_sp,0.436_sp,0.462_sp,0.488_sp,0.514_sp,0.587_sp,0.576_sp,0.547_sp,&
		        0.536_sp,0.532_sp,0.531_sp,0.526_sp,0.514_sp,0.5_sp,0.495_sp,0.496_sp,&
		        0.497_sp,0.499_sp,0.501_sp,0.504_sp]

		private
		public :: e_photon, plancks_law, real_refractive_air, refractive_h2o, &
		            real_refractive_h2o, imag_refractive_h2o, &
		            solve_fluxes, solar_zenith, &
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
											nrwbin,niwbin, &
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
						lambda_low, lambda_high, delta_lambda, b_s_g, ext_s_g, sflux_l, &
						nrwbin,niwbin
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
			allocate( nrwbin(1:ntot), STAT = AllocateStatus)
			if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
			allocate( niwbin(1:ntot), STAT = AllocateStatus)
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
			do i=1,ntot
			    ! calculate the bin averaged real refractive index 
			    nrwbin(i)=qromb(real_refractive_h2o,lambda_low(i),lambda_high(i)) / &
			                (lambda_high(i)-lambda_low(i))
			    ! calculate the bin averaged imaginary refractive index 
			    niwbin(i)=qromb(imag_refractive_h2o,lambda_low(i),lambda_high(i)) / &
			                (lambda_high(i)-lambda_low(i))
			enddo

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
		!>real and imaginary refractive index of water vs wavelength
		!>just a lookup table from 
		!>https://en.wikipedia.org/wiki/Optical_properties_of_water_and_ice#Refractive_Index,_Real_and_Imaginary_Parts_for_Liquid_Water
		!>@param[in] lambda
		!>@param[inout] nr, ni
		subroutine refractive_h2o(lambda,nr,ni)
			implicit none
			real(sp), intent(in) :: lambda
			real(sp), intent(inout) :: nr, ni
			real(sp) :: var, dummy
			integer(i4b) :: iloc
			
            iloc=locate(lam_h2o(1:n_rfrl),lambda)
            iloc=min(n_rfrl-1,iloc)
            iloc=max(1,iloc)
            
            
            ! linear interp nr
            call polint(lam_h2o(iloc:iloc+1), nr_h2o(iloc:iloc+1), &
                        max(min(lambda,lam_h2o(n_rfrl)),lam_h2o(1)), nr,dummy)
            ! linear interp ni
            call polint(lam_h2o(iloc:iloc+1), ni_h2o(iloc:iloc+1), &
                        max(min(lambda,lam_h2o(n_rfrl)),lam_h2o(1)), ni,dummy)
	
		end subroutine refractive_h2o
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!>@author
		!>Paul J. Connolly, The University of Manchester
		!>@brief
		!>real refractive index of water vs wavelength
		!>just a lookup table from 
		!>https://en.wikipedia.org/wiki/Optical_properties_of_water_and_ice#Refractive_Index,_Real_and_Imaginary_Parts_for_Liquid_Water
		!>@param[in] lambda
		function real_refractive_h2o(lambda)
			implicit none
			real(sp), dimension(:), intent(in) :: lambda
			real(sp), dimension(size(lambda)) :: real_refractive_h2o
			real(sp) :: var, dummy
			integer(i4b) :: iloc, i
			
			do i=1,size(lambda)
                iloc=locate(lam_h2o(1:n_rfrl),lambda(i))
                iloc=min(n_rfrl-1,iloc)
                iloc=max(1,iloc)
            
                ! linear interp nr
                call polint(lam_h2o(iloc:iloc+1), nr_h2o(iloc:iloc+1), &
                            max(min(lambda(i),lam_h2o(n_rfrl)),lam_h2o(1)),&
                            real_refractive_h2o(i),dummy)
            enddo
            	
		end function real_refractive_h2o
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!>@author
		!>Paul J. Connolly, The University of Manchester
		!>@brief
		!>imaginary refractive index of water vs wavelength
		!>just a lookup table from 
		!>https://en.wikipedia.org/wiki/Optical_properties_of_water_and_ice#Refractive_Index,_Real_and_Imaginary_Parts_for_Liquid_Water
		!>@param[in] lambda
		function imag_refractive_h2o(lambda)
			implicit none
			real(sp), dimension(:), intent(in) :: lambda
			real(sp), dimension(size(lambda)) :: imag_refractive_h2o
		    real(sp) :: var, dummy
			integer(i4b) :: iloc, i
			
			do i=1,size(lambda)
                iloc=locate(lam_h2o(1:n_rfrl),lambda(i))
                iloc=min(n_rfrl-1,iloc)
                iloc=max(1,iloc)
            
                ! linear interp ni
                call polint(lam_h2o(iloc:iloc+1), ni_h2o(iloc:iloc+1), &
                            max(min(lambda(i),lam_h2o(n_rfrl)),lam_h2o(1)), &
                            imag_refractive_h2o(i),dummy)
            enddo	
		end function imag_refractive_h2o
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
		!>calculate the in-band planck function at temperature t
		!>See http://www.spectralcalc.com/blackbody/inband_radiance.html 
		!>Note that I fit an exponential curve to the approximation in inband_planck
		subroutine inband_planck_param(nbands,lambda_low,lambda_high,temp,blt)
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
				sum_upper=exp(-0.000000000633054_sp*x_upper**4 + &
				               0.000000950820139_sp*x_upper**3  &
				              -0.000489314838945_sp*x_upper**2  &
				              -0.890274259590316_sp*x_upper +  &
				              7.174105190082139_sp)
				
				sum_lower=exp(-0.000000000633054_sp*x_lower**4 + &
				               0.000000950820139_sp*x_lower**3  &
				              -0.000489314838945_sp*x_lower**2  &
				              -0.890274259590316_sp*x_lower +  &
				              7.174105190082139_sp)
				blt(i) = fac*(sum_upper-sum_lower)
				
			enddo
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		end subroutine inband_planck_param
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
			
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!>@author
		!>Paul J. Connolly, The University of Manchester
		!>@brief
		!> Use Mitchell (2000), JAS to calculate the extinction and absorption
        !> See equations 22 and 23 (although some changes needed)
        !>@param[in] ip,jp,kp,r_h, nbands,nrad
        !>@param[inout] sigma_s_clouds, sigma_a_clouds
        !>@param[in] ngs,mugs,lamgs,lambda,nrwbin,niwbin
        subroutine calculate_scattering_and_absorption(ip,jp,kp,r_h, &
                    nbands,nrad,sigma_s_clouds, &
                    sigma_a_clouds,ngs,mugs,lamgs,lambda,nrwbin,niwbin)
			implicit none
			integer(i4b), intent(in) :: nbands, nrad,ip, jp, kp, r_h
			real(sp), dimension(nbands), intent(in) :: lambda, nrwbin, niwbin
			real(sp), intent(in), &
					dimension(1-r_h:kp+r_h, 1-r_h:jp+r_h, 1-r_h:ip+r_h,1:nrad) :: &
					    ngs,lamgs, mugs
			real(sp), intent(inout), &
			        dimension(1-r_h:kp+r_h, 1-r_h:jp+r_h, 1-r_h:ip+r_h,1:nbands) :: &
			        sigma_s_clouds, sigma_a_clouds			
			! local variables:
			integer(i4b) :: i,j,k,m,n
			real(sp) :: gam_a_1,gam_a_2, gam_3,gam_4,gam_5,gam_6,gam_7, gam_m, &
			            g,e,e0,a1,a2,a3,a4,a5,a6,ra,rext, kmax, dmm, test
			real(sp), parameter :: mtun=0.5_sp
			complex(sp) :: q,ncom
			complex(sp), parameter :: imag1=complex(0._sp,1._sp)
			
			sigma_a_clouds=0._sp
			sigma_s_clouds=0._sp
			do m=1,nbands
			    g=8._sp*pi*niwbin(m)/(3._sp*lambda(m))
			    e0=0.25_sp+0.6_sp*(1._sp-exp(-8._sp*pi*niwbin(m)/3._sp))**2 ! eq 8
			    e=e0/lambda(m) ! equation 28
			    ra=0.7393_sp*nrwbin(m)-0.6069_sp ! equation 6
			    rext=ra/2._sp ! equation 11
			    a1=0.25_sp+0.25_sp*exp(-1167._sp*niwbin(m)) ! equation5
			    kmax=mtun/e0 ! equation 9
			    a2=ra/(kmax**mtun*exp(-mtun)*lambda(m)**mtun)
			    a3=rext/(kmax**mtun*exp(-mtun)*lambda(m)**mtun)
			    a4=0.06_sp*pi/lambda(m)
			    a5=(pi/lambda(m))**(-2._sp/3._sp)
			    a6=1._sp
			    ncom=complex(nrwbin(m),niwbin(m))
			    q=imag1*(ncom-complex(1._sp,0._sp))*complex(2._sp*pi/lambda(m),0._sp)
			    do n=1,2 ! cloud, rain, ice
                    do i=1,ip
                        do j=1,jp
                            do k=1-r_h,kp+r_h
                                if(isnan(lamgs(k,j,i,n))) cycle
                                if(isnan(mugs(k,j,i,n))) cycle
                                if(isnan(ngs(k,j,i,n))) cycle

                                ! babs here
                                gam_7=gamma(mtun+mugs(k,j,i,n)+1._sp)
                                gam_6=gam_7*(mtun+mugs(k,j,i,n)+1._sp) ! gamma m+mu+2
                                gam_5=gamma(mugs(k,j,i,n)+1._sp)
                                gam_4=gam_5*1._sp ! gamma mu+2
                                gam_a_1=gam_4*2._sp ! gamma mu+3
                                gam_m=gam_a_1*3._sp ! gamma mu+4
                                gam_a_2=gamma(3._sp+mtun+mugs(k,j,i,n))
                                gam_3=gamma(mugs(k,j,i,n)+7._sp/3._sp)
                            
                                ! equation 22 
                                test= &
                                    pi/4._sp*ngs(k,j,i,n) * &
                                        gam_a_1 / &
                                        (lamgs(k,j,i,n)**(3._sp+mugs(k,j,i,n))) - &
                                    pi/4._sp*ngs(k,j,i,n) * &
                                        gam_a_1 / &
                                        ((lamgs(k,j,i,n)+g)**(3._sp+mugs(k,j,i,n))) + &
                                    a1*pi/4._sp*ngs(k,j,i,n) * &
                                        gam_a_1 / &
                                        ((lamgs(k,j,i,n)+g)**(3._sp+mugs(k,j,i,n))) - &
                                    a1*pi/4._sp*ngs(k,j,i,n) * &
                                        gam_a_1 / &
                                        ((lamgs(k,j,i,n)+2._sp*g)**(3._sp+mugs(k,j,i,n))) + &
                                    a2*pi/4._sp*ngs(k,j,i,n) * &
                                        gam_a_2 / &
                                        ((lamgs(k,j,i,n)+e)**(3._sp+mugs(k,j,i,n)+mtun)) - &
                                    a2*pi/4._sp*ngs(k,j,i,n) * &
                                        gam_a_2 / &
                                        ((lamgs(k,j,i,n)+e+g)**(3._sp+mugs(k,j,i,n)+mtun))

                                if(.not.isnan(test)) &
                                    sigma_a_clouds(k,j,i,m)=sigma_a_clouds(k,j,i,m)+test
                                ! bext here
                            
                                ! equation 23 - there is a change of sign of 
                                ! first term in re function
                                test= &
                                    pi*ngs(k,j,i,n)*gam_a_1 / &
                                        (2._sp*lamgs(k,j,i,n)**(3._sp+mugs(k,j,i,n))) + &
                                    a3*pi*ngs(k,j,i,n)*gam_a_2 / &
                                        (2._sp*(lamgs(k,j,i,n)+&
                                        e)**(3._sp+mugs(k,j,i,n)+mtun)) + &
                                    a6*pi/4._sp*a5*ngs(k,j,i,n)*gam_3*&
                                        (lamgs(k,j,i,n)**(-(mugs(k,j,i,n)+7._sp/3._sp))- &
                                        (lamgs(k,j,i,n)+a4)**(-(mugs(k,j,i,n)+7._sp/3._sp)))+&
                                    pi*ngs(k,j,i,n)*real( &
                                     complex(gam_4,0._sp) / &
                                     (q*(complex(lamgs(k,j,i,n),0._sp)+q)**(&
                                     complex(mugs(k,j,i,n)+2._sp,0._sp))) +&
                                     complex(gam_5,0._sp) / &
                                        (q*q)*((complex(lamgs(k,j,i,n),0._sp)+q)**(&
                                        -(complex(mugs(k,j,i,n)+1._sp,0._sp))) - &
                                        complex(lamgs(k,j,i,n),0._sp)**(&
                                        -complex(mugs(k,j,i,n)+1._sp,0._sp))),sp)+&
                                    a3*pi*ngs(k,j,i,n)*real( &
                                     complex(gam_6,0._sp) / &
                                     (q*(complex(lamgs(k,j,i,n)+e,0._sp)+q)**&
                                     complex(mugs(k,j,i,n)+mtun+2._sp,0._sp)) +&
                                     complex(gam_7,0._sp) / &
                                     (q*q)* &
                                        ((complex(lamgs(k,j,i,n)+e,0._sp)+q)**(&
                                        -complex(mugs(k,j,i,n)+mtun+1._sp,0._sp)) - &
                                        complex(lamgs(k,j,i,n)+e,0._sp)**(&
                                        -complex(mugs(k,j,i,n)+1._sp,0._sp))),sp)

                                if(.not.isnan(test)) &
                                    sigma_s_clouds(k,j,i,m)=sigma_s_clouds(k,j,i,m)+test

                                ! Geometric approximation
!                                 sigma_s_clouds(k,j,i,m)=sigma_s_clouds(k,j,i,m)+&
!                                     pi/4._sp* &
!                                     2._sp*ngs(k,j,i,n)*gamma(mugs(k,j,i,n)+3._sp) / &
!                                     lamgs(k,j,i,n)**(mugs(k,j,i,n)+3._sp)
                                
                            enddo
                        enddo
                    enddo
                enddo
            enddo
            ! we need scattering, not extinction
            sigma_s_clouds=sigma_s_clouds-sigma_a_clouds
        end subroutine calculate_scattering_and_absorption
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
								lambda_low, lambda_high, lambda,nrwbin, niwbin, &
								b_s_g, sflux_l,&
								rhoan, trefn, dz,dzn, albedo, emiss, &
								quad_flag, th, flux_u, flux_d, &
								asymmetry_water, nrad, ngs,lamgs,mugs, &
								cloud_flag, &
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
					lambda_low, lambda_high, lambda, nrwbin, niwbin, b_s_g, sflux_l
			real(sp), intent(in), dimension(1-r_h:kp+r_h) :: rhoan, trefn, dz,dzn
			real(sp), intent(in), &
					dimension(1-r_h:kp+r_h, 1-r_h:jp+r_h, 1-r_h:ip+r_h) :: th
			real(sp), intent(inout), &
					dimension(1-r_h:kp+r_h, 1-r_h:jp+r_h, 1-r_h:ip+r_h,1:nbands) :: &
					flux_u, flux_d
			real(sp), intent(in), &
					dimension(1-r_h:kp+r_h, 1-r_h:jp+r_h, 1-r_h:ip+r_h,1:nrad) :: &
					    ngs,lamgs, mugs
			real(sp), intent(in) :: lat, lon, albedo, emiss, asymmetry_water
			integer(i4b), intent(in) :: quad_flag, nrad

            integer(i4b), dimension(3), intent(in) :: coords
            integer(i4b), dimension(3), intent(in) :: dims
            integer(i4b), intent(in) :: id, comm3d
            
            integer(i4b), intent(in) :: tdstart,tdend
            real(sp), dimension(1:tdend), intent(inout) :: a,b,c,r,u
            logical, intent(in) :: cloud_flag

			
			! local variables:
			integer(i4b) :: doy,diy, year,error
			real(sp) :: tod, cos_theta_s, frac,msend
			real(sp), dimension(0:kp+1) :: sigma_ray, sigma_tot,taun,dtau,dtaun
			real(sp), dimension(0:kp+1) :: tau,direct
			real(sp), dimension(kp+1,nbands) :: blt
			real(sp), dimension(1-r_h:kp+r_h, 1-r_h:jp+r_h, 1-r_h:ip+r_h,1:nbands) :: &
			        sigma_s_clouds, sigma_a_clouds
			real(sp), dimension(1:kp) :: gamma1, gamma2, gamma3, omg_s
			logical :: leap
			real(sp) :: ga=0._sp
			
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



            kend=kp+1
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Calculate the absorption and extinction for clouds
            if(cloud_flag) then
                call calculate_scattering_and_absorption(ip,jp,kp,r_h, &
                    nbands,nrad,sigma_s_clouds, &
                    sigma_a_clouds,ngs,mugs,lamgs,lambda,nrwbin,niwbin)
            endif
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


			do m=1,nbands
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! extinction:							    					 !
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                do k=1,kend 
                    ! number of molecules per cubic metre of air multiplied by 
                    ! collision cross-section:
                    sigma_ray(k)=rhoan(k)*navog/ma*b_s_g(m)
                enddo
                do i=1,ip
                    do j=1,jp
                
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        ! calculate planck function at all levels p-points:				 !
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        do k=1,kp+1 
                            call inband_planck(1,lambda_low(m),lambda_high(m), &
                                                    trefn(k)+th(k,j,i),blt(k,m))
                            !if(k==1.and.i==1.and.j==1) print *,blt(k,m)
                        enddo
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        ! extinction:							    					 !
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        do k=1,kend 
                            sigma_tot(k)=sigma_ray(k)
                            if(cloud_flag) then
                                sigma_tot(k)=sigma_tot(k)+sigma_s_clouds(k,j,i,m) + &
                                    sigma_a_clouds(k,j,i,m)
                            endif
                        enddo
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        ! optical depth:					    						 !
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        ! sigma_ray is the extinction on p-points
                        ! tau is optical depth on w-points
                        ! mpi - add bottom taus to all and add cumulative
!                         dtau(kp+1)=0._sp
                        tau(kp+1)=0._sp
                        tau(kp)=0._sp ! for top level tau(kp) needs to be defined taking into 
                                      ! account the atmosphere above domain
                        if(dims(3)==(coords(3)+1)) then
                            tau(kp+1)=1.e-2_sp
!                             dtau(kp+1)=-tau(kp+1)
                        endif
                        do k=kend-1,0,-1 
                            ! cumulative extinction x dz
                            tau(k)=tau(k+1)+sigma_tot(k+1)*dz(k)
                        enddo
                
            
                        ! for parallel solver tau and dtau need some parallel communication:
                        ! http://mpitutorial.com/tutorials/mpi-scatter-gather-and-allgather/
                        msend=tau(1)  ! message to send
                        call MPI_allgather(msend,1,MPI_REAL8,mvrecv,1,MPI_REAL8,&
                            comm3d,error)
                        tau=tau+sum(mvrecv(  coords(3)+2:dims(3)  ))
                
                
                
                
                        ! change in optical depth on w-points
!                         do k=kp,0,-1 
!                             dtau(k)=(tau(k+1)-tau(k) ) 
!                                                 ! tau = 0 at top, so set to negative 
!                                                 ! to get derivative in correct direction
!                         enddo
!                         msend=dtau(1)  ! message to send
!                         call MPI_allgather(msend,1,MPI_REAL8,mvrecv,1,MPI_REAL8,&
!                             comm3d,error)
!                         if((coords(3)+1).lt.(dims(3))) dtau(kp+1)=mvrecv(coords(3)+2)
                        !if(coords(3)==(dims(3)-1)) dtau(kp)=0._sp
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


                
                
                
                
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        ! diffuse phase function scaling on p-points:					 !
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        select case (quad_flag)
                            case(0) ! quadrature
                                do k=1,kp
                                    ! equation 9.103, Jacobson - single scattering albedo
                                    omg_s(k)=(sigma_ray(k)+sigma_s_clouds(k,j,i,m)) / &
                                        (sigma_ray(k)+sigma_s_clouds(k,j,i,m)+&
                                        sigma_a_clouds(k,j,i,m))
                                    ! equation 9.115, Jacobson - effective asymmetry param
                                    ga=(sigma_s_clouds(k,j,i,m)*asymmetry_water) / &
                                        (sigma_ray(k)+sigma_s_clouds(k,j,i,m))
                                    
                                    gamma1(k)=(1._sp-omg_s(k)*(1._sp+ga)/2._sp)/mu1
                                    gamma2(k)=omg_s(k)*(1._sp-ga)/(2._sp*mu1)
                                    gamma3(k)=(1._sp-3._sp*ga*mu1*cos_theta_s)/2._sp
                                enddo
                            case(1) ! eddington
                                do k=1,kp
                                    ! equation 9.103, Jacobson - single scattering albedo
                                    omg_s(k)=(sigma_ray(k)+sigma_s_clouds(k,j,i,m)) / &
                                        (sigma_ray(k)+sigma_s_clouds(k,j,i,m)+&
                                        sigma_a_clouds(k,j,i,m))
                                    ! equation 9.115, Jacobson - effective asymmetry param
                                    ga=(sigma_s_clouds(k,j,i,m)*asymmetry_water) / &
                                        (sigma_ray(k)+sigma_s_clouds(k,j,i,m))

                                    gamma1(k)=(7._sp-omg_s(k)*(4._sp+3._sp*ga))/4._sp
                                    gamma2(k)=-(1._sp-omg_s(k)*(4._sp-3._sp*ga))/4._sp
                                    gamma3(k)=(2._sp-3._sp*ga*cos_theta_s)/4._sp				
                                enddo
                            case default
                                print *,'not coded'
                                stop
                        end select
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
                
            
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        ! set up tridiagonal problem:									 !
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        ! dtau(0,kp+1) needs to be set
                        ! tau(0,kp+1) needs to be set				
                        do k=1,kp
                            ! define coefficients for FDE:
                            Ac=1._sp+0.5_sp*dz(k-1)*gamma1(k)*sigma_tot(k)
                            Bc=-1._sp+0.5_sp*dz(k-1)*gamma1(k)*sigma_tot(k)
                            Cc=-0.5_sp*dz(k-1)*gamma2(k)*sigma_tot(k)
                            Dc=Cc
                    
                            alp=-Cc
                            bet=alp
                            gam=-Ac
                            del=-Bc
                
                            ! matrix coefficients:
                            Ak=(Bc-bet*Cc/del)
                            Bk=(Dc-gam*Cc/del)
                            Ck=(Ac-alp*Cc/del)
                            Skp=0._sp
                            Skm=0._sp
                            if(cos_theta_s > 0._sp) then
                                ! short-wave (see equations 9.121 Jacobson):
                                Skp=dz(k-1)*sigma_tot(k)*gamma3(k)*omg_s(k)*sflux_l(m)* &
                                    exp(-0.5_sp*(tau(k+1)+tau(k))/cos_theta_s)
                                Skm=-dz(k-1)*sigma_tot(k)*(1._sp-gamma3(k))* &
                                        omg_s(k)*sflux_l(m)*&
                                        exp(-0.5_sp*(tau(k+1)+tau(k))/cos_theta_s)
                            endif
                            ! long wave (see equations 9.122 Jacobson):
                            Skp=Skp+&
                                dz(k-1)*sigma_tot(k)*2._sp*pi*(1._sp-omg_s(k))*blt(k,m)
                            Skm=Skm-& 
                                dz(k-1)*sigma_tot(k)*2._sp*pi*(1._sp-omg_s(k))*blt(k,m)

                            Dk=Skp-Cc/del*Skm
                            Ek=(gam-Dc*bet/Bc)
                            Fk=(alp-Ac*bet/Bc)
                            Gk=(del-Cc*bet/Bc)
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
                        
                
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        ! boundary conditions
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        if(coords(3)==0) then
                            a(1)=0._sp ! this is specific to parallel
                            b(1)=-1._sp ! diagonal, so fine
                            c(1)=albedo ! ! c(1) has to be set in both serial and parallel
                            r(1)=0._sp
                            if(cos_theta_s > 0._sp) then
                                r(1)=-albedo*cos_theta_s*sflux_l(m)* &
                                    exp(-tau(0)/cos_theta_s)
                            endif
                            r(1)=r(1)-emiss*pi*blt(1,m) 
                                                ! this needs to be for the surface, i.e. 
                                                ! blt at tref(0)
                        endif
                
                
                        if(coords(3)==(dims(3)-1)) then
                            a(tdend)=0._sp ! no contribution from upward flux in last equation
                            b(tdend)=-1._sp ! was 1
                            r(tdend)=0._sp ! toa downward flux at night
                            c(tdend)=0._sp ! this is specific to parallel
                            if(cos_theta_s > 0._sp) then
                                r(tdend)=-cos_theta_s*sflux_l(m) 
                                                        ! toa downward flux - day
                                                        ! should be multiplied by tau(kp) 
                                                        ! for atmosphere above
                            endif
                        endif                
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        ! solve tridiagonal matrix:										 !
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        call parallel_tridag(a,b,c,r,u,tdend,coords,dims, id, comm3d)
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        

                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        ! assumes temperature is given by reference state:				 !
                        ! message passing done out of this routine                       !
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        kstart=1
                        if(coords(3)==0) then ! bottom
                            ! flux_d(kp) recv
                            do k=1,kp
                                ! sets the surface fluxes too
                                flux_u(k-kstart,j,i,m)=u(2*(k)-1)
                                flux_d(k-kstart,j,i,m)=u(2*k)
                            enddo
                            flux_u(kp+1-kstart,j,i,m)=u(2*(kp+1)-1)
                        else if (coords(3)==(dims(3)-1)) then
                            flux_d(0,j,i,m)=u(1) ! needs to be passed to below
                            do k=1,kp
                                flux_u(k,j,i,m)=u(2*(k)) ! 2,4,6
                                flux_d(k,j,i,m)=u(2*k+1) ! 3,5,7
                            enddo
                        else
                            flux_d(0,j,i,m)=u(1) ! needs to be passed to below
                            ! mpi_issend to one below
                            ! mpi_irecv from above flux_d(kp,j,i,m)=
                            do k=1,kp-1
                                flux_u(k,j,i,m)=u(2*(k))
                                flux_d(k,j,i,m)=u(2*k+1)
                            enddo
                            flux_u(kp,j,i,m)=u(2*kp)
                        endif
                        ! bottom and top
                        if((coords(3)==0).and.(coords(3)==(dims(3)-1))) then
                            flux_d(kp+1-kstart,j,i,m)=u(2*(kp))                            
                        endif
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    enddo
                enddo
			enddo
			
									
		end subroutine solve_fluxes
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	end module radiation
	