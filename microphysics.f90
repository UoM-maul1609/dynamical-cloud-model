	!> @mainpage
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@copyright 2018
	!>@brief
	!>Prognostic Aerosol Microphysics Module (PAMM):
	!>A warm bulk cloud microphysics module, with prognostic aerosol
	!>  for use with other models.
	!> compile using the Makefile. Requires linking to other wrapper model for execution.

	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>microphysics code for the different cloud models
    module p_micro_module
    use numerics_type
    use numerics, only : find_pos, poly_int, erfinv

    use bam, only : n_mode, n_sv, giant_flag, method_flag, sv_flag, &
        	n_aer1, d_aer1, sig_aer1, molw_core1, density_core1, nu_core1, org_content1, &
        	molw_org1, density_org1, delta_h_vap1, nu_org1, log_c_star1, p_test, t_test, &
    		w_test, act_frac1, smax1, dcrit2, &
    		a_eq_7, b_eq_7, &
    		ctmm_activation, initialise_arrays, read_in_bam_namelist
    use numerics, only : quad2d_qgaus, invgammainc, romb
    
    private
    public :: p_microphysics_3d, &
            p_microphysics_2d, p_microphysics_1d, read_in_pamm_bam_namelist, &
            p_initialise_aerosol_3d,p_initialise_aerosol, p_initialise_aerosol_1d, &
            calculate_gamma_params
            
    real(wp) :: mrthresh, mrupper, miupper, f_mode2, lambda0r, lambda0i, n0r, n0i, &
            pthreshr, pthreshi
    real(wp), parameter :: phi_mode2=0.35_wp, probthresh=0.9999_wp, grav=9.81_wp
    
    ! Chen and Lamb (1994) Gamma variable fit (scaled and centred logarithm)
    integer(i4b), parameter :: n_cl=18
    real(wp), dimension(n_cl), parameter :: gam_cl=[-0.072328469664620_wp, &
        -0.324623262465577_wp, 0.363138099937540_wp, 3.323089908344732_wp, &
        0.874844989423720_wp, &
        -13.554426432462339_wp, -9.810322482346461_wp, 27.846739088352344_wp, &
        26.480447842355410_wp,&
         -29.890199206698309_wp, -32.327548996894521_wp, 15.827423311652167_wp, &
         18.466605783503052_wp, -4.158566361058538_wp, -5.039533848938808_wp, &
         1.477272813054374_wp, 1.038600921563425_wp, -0.457007828432810_wp]
    real(wp), dimension(2), parameter :: gam_mu_cl=[260.163817050062335_wp, &
                                                8.274747821396463_wp]
                                                
    ! Vardiman (1978) fits to figure 6 - note delta M in units on g cm s-1
    real(wp), dimension(3), parameter :: vard01=[0.000495314304309_wp, &
                                                 0.281199363154805_wp, &
                                                 3.380130133900658_wp], &
                                         vard02=[0.304288838581395_wp, &
                                                 4.452491028368538_wp, &
                                                 17.511640705855431_wp], &
                                         vard03=[1.549508244781713_wp, &
                                                 21.756014605694737_wp, &
                                                 77.539493556502251_wp], &
                                         vard04=[0.924318964759507_wp, &
                                                 15.774108106443462_wp, &
                                                 68.805308506959534_wp], &
                                         vard05=[0.162609020092783_wp, &
                                                 3.031949785103254_wp, &
                                                 15.296369750198556_wp]
                                          
    
    ! physical constants
    real(wp), parameter :: rhow=1000._wp, rhoi=920._wp,lv=2.5e6_wp,ls=2.8e6_wp,lf=ls-lv, &
    					   cp=1005._wp, cw=4187._wp, cice=2093._wp, r=8.314_wp, &
    						mw=18e-3_wp, ma=29e-3_wp, ra=r/ma,rv=r/mw, eps1=ra/rv, &
    						ttr=273.15_wp, t_hom=ttr-36._wp, joules_in_an_erg=1.0e-7_wp, &
    						joules_in_a_cal=4.187_wp, &
    						gamma_liq=0.072_wp, DEcrit=0.2_wp, &
    						oneoversix=1._wp/6._wp, dtt=10.e-6_wp, &
        					oneoverthree=1._wp/3._wp, oneovernine=1._wp/9._wp, &
        					oneoverpi=1._wp/pi, phi_phillips=3.5e-3_wp    

    ! mass-diameter and size spectra relations
    real(wp), parameter :: cr=523.6_wp, cc=523.6_wp, &
                        cs=52.36_wp, cg=261.8_wp, ci=104._wp, &
    					dr=3_wp, dc=3_wp, ds=3._wp, dg=3._wp, di=3._wp, &
    					alpha_r=2.5_wp, alpha_c=0.0_wp, & ! note, alpha_c is for a "mass" - number distribution
    					alpha_s=2.5_wp, alpha_g=2.5_wp, alpha_i=0._wp
    					
	! terminal fall-speed relations
	real(wp), parameter :: a_r=362._wp, a_c=362._wp, &
	                        a_s=4.84_wp, a_g=253._wp, a_i=71.34_wp, &
							b_r=0.65_wp, b_c=0.65_wp, &
							b_s=0.25_wp, b_g=0.734_wp, b_i=0.6635_wp, &
							f_r=0._wp, f_c=0._wp, f_s=0._wp, f_g=0._wp, f_i=0._wp
							
	! autoconversion
	real(wp), parameter :: aw0=1e-3_wp, dwa=20e-6_wp, nl=2.4e8_wp, &
						lw0=rhow*pi/6._wp*nl*dwa**3, &
						tsaut=60._wp, dimax=0.3e-3_wp, di2s=0.33e-3_wp, &
						lambda_imin=(1._wp+di+alpha_i)/dimax, &
						tsbreak=60._wp, lambda_s_break=1000._wp
	
    ! microphysical values:
    real(wp), parameter :: hm_rate=3.5e8_wp, nar=1.1e15_wp, nbr=0._wp, &
    						rho0=1.2_wp, bbigg=100._wp, abigg=0.66_wp
    real(wp) :: mi0=1.e-14_wp
    

	! coalescence efficiencies
	real(wp), parameter :: erw=1._wp, erg=1._wp, ers=1._wp, eri=1._wp, esw=1._wp, &
						egw=1._wp, eiw=1._wp, egs_wet=1._wp, egi_wet=1._wp
	
	! variables used in various process rates:
	real(wp) :: gam1r,gam2r,gam3r,gam1c, gam2c, gam3c, &
	            gam1i,gam2i, gam1s, gam2s,gam1g,gam2g, &
	            gam3ai,gam3bi,gam4ai,gam4bi,gam5ai,gam5bi, &
	            gam6ai, gam6bi, fall_q_i_hw,fall_n_i_hw, &
				fall_q_r, fall_q_c, fall_q_s, fall_q_g, fall_n_r, fall_n_s, fall_n_g, &
				fall_q_i, fall_n_i, fall_n_c, &
				phi_r, mass_iacr,num_iacr, mass_sacw_i, mass_iacw, &
				mass_raci1,mass_raci2,mass_raci3, &
				mass_racg1,mass_racg2,mass_racg3, &
				mass_iacr1,mass_iacr2,mass_iacr3, &
				mass_sacg1,mass_sacg2,mass_sacg3, &
				mass_gacr1,mass_gacr2,mass_gacr3, &
				mass_gacs1,mass_gacs2,mass_gacs3, &
				num_raci1,num_raci2,num_raci3, &
				num_iacr1,num_iacr2,num_iacr3, num_racg1,num_racg2,num_racg3, &
				num_sacg1, num_sacg2,num_sacg3, &
				mass_gacw, mass_gaci, &
				nu_r1,nu_r2,nu_i1, nu_i2, nu_s1, nu_s2, nu_g1, nu_g2, &
				mass_imm, num_imm, q0sat, &
				chi_rain, chi_cloud, chi_ice, chi_snow, chi_graupel, &
				chi_rain1, chi_cloud1, chi_ice1, chi_snow1, chi_graupel1, &
				chi_num_ice, chi_num_ice1, &
				gam1cr,gam2cr ! for radiation
				
	! some work space to transfer data between functions
	real(wp), dimension(5) :: phillips_br_workspace
				
	! Seifert and Beheng autoconversion
	real(wp) :: kc, kr, xstar
				
	real(wp), dimension(3) :: c=[1._wp,2._wp,1._wp]
	integer(i4b) :: k
	real(wp) :: isnow, iice, iice2, f1,f2,a,b, qsmall=1.e-30_wp
	
	! to send to integrator
	real(wp) :: a_hw_new, pre_hw_new, ci_new, t_send, lam_freeze, n0_freeze
	
	integer :: n_modes_prof, n_levels_s
	real(wp), allocatable, dimension(:,:) :: n_read, sig_read, d_read
	real(wp), allocatable, dimension(:) :: z_read
	real(wp) :: small_number
    contains
    
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>convert an integer to a string
	!>@param[in] i: input integer
	!>@param[out] res: the string
    function itoa(i) result(res)
      character(:),allocatable :: res
      integer,intent(in) :: i
      character(range(i)+2) :: tmp
      write(tmp,'(i0)') i
      res = trim(tmp)
    end function


	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculate the kth integral moment of a lognormal distribution
	!>see:
	!>https://en.wikipedia.org/wiki/Log-normal_distribution#Characteristic_function_and_moment_generating_function
	!>@param[in] k: moment to calculate
	!>@param[in] n,sig,d: parameters of the lognormal distribution
	!>@param[out] mom: the integral moment to calculate
    function ln_mom(k,n,sig,d) 
        implicit none
        integer(i4b), intent(in) :: k
        real(wp), intent(in) :: n,sig,d
        real(wp) :: ln_mom
        
        ln_mom=n*exp(real(k,wp)*log(d)+real(k,wp)**2*sig**2/2._wp)
    end function ln_mom
    
    

	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculate the partial kth integral moment of a lognormal distribution
	!>see:
	!>https://en.wikipedia.org/wiki/Log-normal_distribution#Partial_expectation
	!> https://math.stackexchange.com/questions/2055782/partial-expectations-of-lognormal-distributions
	!>@param[in] k: moment to calculate
	!>@param[in] n,sig,d: parameters of the lognormal distribution
	!>@param[out] mom: the integral moment to calculate
    function ln_part_mom(k,a,n,sig,d) 
        implicit none
        integer(i4b), intent(in) :: k
        real(wp), intent(in) :: a,n,sig,d
        real(wp) :: ln_part_mom
        
        real(wp) :: x1, phi1
        
        x1=(log(a)-log(d)-sig**2*real(k,wp))/(sig*sqrt(2._wp))
        phi1=0.5_wp*(1._wp+erf(x1))

        ln_part_mom= &
            max(n*exp(log(d)*real(k,wp)+sig**2*real(k,wp)**2*0.5_wp)* (1._wp-phi1), 0._wp)
        
    end function ln_part_mom
    
    
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculate the parameters of a lognormal distribution
	!>@param[in] n,s,m, rho
	!>@param[inout] sig_aer,d_aer
	subroutine ln_params_from_integral_moms(n,s,m,rho,sig_aer,d_aer)
	    implicit none
	    real(wp), intent(in) :: n,s,m, rho
	    real(wp), intent(inout) :: sig_aer, d_aer
	
	    if(m .gt. 0._wp) then
            ! this was derived by calculating moments of the distribution
            ! and solving to find dm and sig
            sig_aer=log( (36._wp*m**2*n*pi/rho**2)**(1._wp/3._wp) / s )
            sig_aer=sqrt(sig_aer)
        
            d_aer=(log(6._wp*m/(n*pi*rho))-4.5_wp*sig_aer**2) / 3._wp
            d_aer=exp(d_aer)	
        else
            sig_aer=0.3_wp
            d_aer=60.e-9_wp            
        endif	
	end subroutine ln_params_from_integral_moms
	

	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculate the parameters and properties of a lognormal distribution
	!>@param[in] n_mode, n,s,m
	!>@param[inout] rho, molw, nu, n_aer, sig_aer,d_aer,n_mix,s_mix,m_mix
	subroutine ln_params_and_props_from_integral_moms(n_mode,&
	                            n,s,m,n_aer,rho,molw,nu, &
	                            sig_aer,d_aer,n_mix,s_mix,m_mix)
	    implicit none
	    integer(i4b), intent(in) :: n_mode
	    real(wp), intent(in) :: n
	    real(wp), dimension(n_mode-1), intent(in) :: s,m
	    real(wp), dimension(n_mode), intent(inout) :: rho,molw,nu
	    real(wp), intent(inout) :: n_aer,sig_aer, d_aer,n_mix,s_mix,m_mix

        
	    n_aer=n
	    n_mix=n
	    s_mix=sum(s)
	    m_mix=sum(m)
        if(m_mix .gt. 0._wp) then
            ! conserve volume of particle:	
            rho(n_mode) = sum(m) / sum(m/rho(1:n_mode-1))
            ! conserve total number of moles in particle:	
            molw(n_mode) = sum(m) / sum(m/molw(1:n_mode-1))
            ! conserve total number of moles of ions in particle:	
            nu(n_mode) = molw(n_mode) * sum(m*nu(1:n_mode-1)/molw(1:n_mode-1)) / sum(m) 

            ! this was derived by calculating moments of the distribution
            ! and solving to find dm and sig
            sig_aer=log( (36._wp*m_mix**2*n*pi/rho(n_mode)**2)**(1._wp/3._wp) / s_mix )
            if(sig_aer.le.0._wp) then
                sig_aer=0.3_wp
                d_aer=60.e-9_wp
                n_aer=0._wp
            else
                sig_aer=sqrt(sig_aer)
        
                d_aer=(log(6._wp*m_mix/(n*pi*rho(n_mode)))-4.5_wp*sig_aer**2) / 3._wp
                d_aer=exp(d_aer)	
            endif

        else
            rho(n_mode)=rho(n_mode-1)    
            molw(n_mode)=molw(n_mode-1)    
            nu(n_mode)=nu(n_mode-1)    
            
            sig_aer=0.3_wp
            d_aer=60.e-9_wp
        endif
        
	end subroutine ln_params_and_props_from_integral_moms
	


	


	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>read in the data from the namelists for the BAM module
	!> and set variables for microphysics
	!>@param[in] nmlfile, aero_nmlfile
	!>@param[in] aero_prof_flag
	!>@param[in] ice_flag
	!>@param[inout] q_name, q_type, c_s, c_e
	!>@param[inout] nq,ncat, nprec, iqv, iqc, inc, iqr,inr,iqi,ini,iai
	!>$param[inout] n_modeg, cat_am,cat_c, cat_r
	!>@param[inout] cat_i
	subroutine read_in_pamm_bam_namelist(nmlfile, aero_nmlfile, &
	            aero_prof_flag, &
	            ice_flag, &
                q_name,q_type,c_s,c_e,nq,ncat,nprec,n_modeg, &
                iqv,iqc,inc, iqr,inr, iqi,ini, iai, cat_am,cat_c, cat_r, cat_i)
		use bam, only : read_in_bam_namelist, n_mode
		implicit none
        logical, intent(in) :: aero_prof_flag, ice_flag
        character (len=200), intent(in) :: nmlfile
        character (len=200), intent(in) :: aero_nmlfile
        integer(i4b), intent(inout) :: nq, ncat, nprec, iqv, iqc, inc, iqr,inr, &
                                         iqi, ini, cat_am,&
                                        cat_c, cat_r, cat_i, iai
        integer(i4b), intent(inout) :: n_modeg
        integer(i4b), intent(inout), dimension(:), allocatable :: q_type, c_s, c_e
        character(len=20), dimension(:), allocatable :: q_name
        
        integer(i4b) :: i
        ! define namelists for aerosol profile
        namelist /aerosol_profile/ n_modes_prof, n_levels_s
        namelist /aerosol_profile_data/ n_read,sig_read,d_read,z_read
        
        ! read in namelist
        call read_in_bam_namelist(nmlfile)
        
        if(aero_prof_flag) then
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! read in aerosol profile num modes									   !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            open(8,file=aero_nmlfile,status='old', recl=80, delim='apostrophe')
            read(8,nml=aerosol_profile)
            if(n_modes_prof .gt. n_mode) then
                !n_modes_prof=n_mode
            else
                n_mode=n_modes_prof
            endif
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! allocate and read aerosol profile data 							   !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            allocate(n_read(n_modes_prof,n_levels_s))
            allocate(sig_read(n_modes_prof,n_levels_s))
            allocate(d_read(n_modes_prof,n_levels_s))
            allocate(z_read(n_levels_s))
            read(8,nml=aerosol_profile_data)
            close(8)
            if(n_modes_prof .gt. n_mode) then
                n_modes_prof=n_mode
            endif
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        endif        
        
        n_modeg=n_mode  ! number of modes used in activation (>=1)
                        ! if n_mode==1 there there should be no mixed mode
        ncat=3+n_mode   ! number of categories that are advected separately

        
        nq=6+ &              ! vapour, qc,qr,nc,nr,mixed-mode number
            (n_mode-1)*3 + & ! aerosol
            (n_mode-1)*3 + & ! mixed-mode aerosol
            (n_mode-1)*3 + & ! aerosol in cloud water
            (n_mode-1)*3     ! aerosol in rain water

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! if you would like to calculate ice microphysics                                !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(ice_flag) then
            ncat=ncat+1             ! add the ice category
            nq=nq+(n_mode-1)*3 + &  ! aerosol in ice water
                6                   ! qi, ni, shape, density, number of mon, rime mass
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        allocate(q_name(nq))
        allocate(q_type(ncat))
        allocate(c_s(ncat))
        allocate(c_e(ncat))

        q_type(1)=0 ! vapour
        c_s(1)=1
        c_e(1)=1
        do i=1,n_mode-1
            q_type(i+1)=3   ! aerosol
            c_s(i+1)=(i-1)*3+2
            c_e(i+1)=(i)*3+2-1 
        enddo
        i=n_mode ! last mode is mixed-mode - 3x(n_mode-1)+1 for total number
        q_type(i+1)=3   ! aerosol
        c_s(i+1)=(i-1)*3+2
        c_e(i+1)=2+(n_mode-1)*6


        q_type(2+n_mode:2+n_mode)=1 
            ! cloud water - 3*(n_mode-1)+2 (for cloud number and mass)
        c_s(2+n_mode:2+n_mode)=(n_mode-1)*6+3
        c_e(2+n_mode:2+n_mode)=(n_mode-1)*6+3+1+3*(n_mode-1)

        q_type(3+n_mode:3+n_mode)=1 ! rain water        
        c_s(3+n_mode:3+n_mode)=(n_mode-1)*6+3+2+3*(n_mode-1)
        c_e(3+n_mode:3+n_mode)=(n_mode-1)*6+3+2+3*(n_mode-1)+1+3*(n_mode-1)
        
        
        ! name the categories
        q_name(1)="qv"
        ! next externally mixed aerosol particles
        do i=1,n_mode-1
            q_name((i-1)*3+2)="an_" // itoa(i)
            q_name((i-1)*3+3)="as_" // itoa(i)
            q_name((i-1)*3+4)="am_" // itoa(i)
        enddo
        
        ! internally mixed aerosol particles (total number, then n,sa,m for each)
        q_name((n_mode-1)*3+2)="an_m_t"
        do i=1,n_mode-1
            q_name((n_mode-1)*3+3+(i-1)*3) = "an_m_" // itoa(i)
            q_name((n_mode-1)*3+4+(i-1)*3) = "as_m_" // itoa(i)
            q_name((n_mode-1)*3+5+(i-1)*3) = "am_m_" // itoa(i)
        enddo

        ! cloud water
        q_name((n_mode-1)*6+3) = "nc"
        q_name((n_mode-1)*6+4) = "qc"
        ! aerosol particles in cloud water
        do i=1,n_mode-1
            q_name((n_mode-1)*6+5+(i-1)*3)="cn_" // itoa(i)
            q_name((n_mode-1)*6+6+(i-1)*3)="cs_" // itoa(i)
            q_name((n_mode-1)*6+7+(i-1)*3)="cm_" // itoa(i)
        enddo
        ! rain water
        q_name((n_mode-1)*6+3*(n_mode-1)+5) = "nr"
        q_name((n_mode-1)*6+3*(n_mode-1)+6) = "qr"
        ! aerosol particles in rain water
        do i=1,n_mode-1
            q_name((n_mode-1)*6+3*(n_mode-1)+7+(i-1)*3)="rn_" // itoa(i)
            q_name((n_mode-1)*6+3*(n_mode-1)+8+(i-1)*3)="rs_" // itoa(i)
            q_name((n_mode-1)*6+3*(n_mode-1)+9+(i-1)*3)="rm_" // itoa(i)
        enddo
        
        inc=(n_mode-1)*6+3
        iqc=(n_mode-1)*6+4       
        
        if(ice_flag) then
            nprec=2
        else
            nprec=1
        endif
        
        iqv=1
        
        cat_am=(n_mode-1)+2
        cat_c=cat_am+1
        cat_r=cat_c+1
        inr=(n_mode-1)*6+3*(n_mode-1)+5
        iqr=(n_mode-1)*6+3*(n_mode-1)+6
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! if you would like to calculate ice microphysics                                !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(ice_flag) then
            q_type(3+n_mode+1:3+n_mode+1)=1 ! ice water 
            c_s(3+n_mode+1:3+n_mode+1)=(n_mode-1)*6+3+2+3*(n_mode-1)+1+3*(n_mode-1)+1
            c_e(3+n_mode+1:3+n_mode+1)=(n_mode-1)*6+3+2+3*(n_mode-1)+1+3*(n_mode-1)+ &
                        (n_mode-1)*3 + 6 

            ! ice water
            q_name( 7+(n_mode-1)*12) = "ni"
            q_name( 8+(n_mode-1)*12) = "qi"
            q_name( 9+(n_mode-1)*12) = "phi"
            q_name(10+(n_mode-1)*12) = "vol"
            q_name(11+(n_mode-1)*12) = "nmon"
            q_name(12+(n_mode-1)*12) = "rmass"
            ! aerosol particles in ice water
            do i=1,n_mode-1
                q_name(13+(n_mode-1)*12+(i-1)*3)="in_" // itoa(i)
                q_name(14+(n_mode-1)*12+(i-1)*3)="is_" // itoa(i)
                q_name(15+(n_mode-1)*12+(i-1)*3)="im_" // itoa(i)
            enddo
            ini=7+(n_mode-1)*12
            iqi=8+(n_mode-1)*12
            iai=13+(n_mode-1)*12
            cat_i=cat_r+1
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
	end subroutine read_in_pamm_bam_namelist


    
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>initialises aerosol profile
	!>@param[in] aero_prof_flag
	!>@param[in] nq, ncat: number of q variables, categories
	!>@param[in] c_s, c_e: start and end pointers for categories
	!>@param[in] inc: pointer to drop number category
	!>@param[in] kp, o_halo: number of i, k and halo points
	!>@param[in] z,rho, p, t: grid values
	!>@param[inout] q: q_variables
    subroutine p_initialise_aerosol_1d(aero_prof_flag,nq,ncat,c_s,c_e, &
                inc, kp,o_halo, z,rho,p,t,q)
                
        use bam, only : n_mode, n_sv, n_aer1, d_aer1, sig_aer1, density_core1, &
                    nu_core1, molw_core1, org_content1, molw_org1, density_org1, &
                    delta_h_vap1, nu_org1, log_c_star1, a_eq_7, b_eq_7, &
                    initialise_arrays, ctmm_activation,find_d_and_s_crits
        implicit none
        ! arguments:
        logical :: aero_prof_flag
        integer(i4b), intent(in) :: nq, ncat, inc, kp, o_halo
        integer(i4b), dimension(ncat), intent(in) :: c_s, c_e
        real(wp), dimension(-o_halo+1:kp+o_halo), intent(in) :: z
        real(wp), dimension(-o_halo+1:kp+o_halo), intent(in) :: &
    					rho, p, t
        real(wp), dimension(-o_halo+1:kp+o_halo,nq), &
            intent(inout) :: q

        ! local variables
        integer(i4b) :: i, k, AllocateStatus, iloc
        real(wp) :: w, smax, phi, xx, kmom, var, dummy
        real(wp), dimension(:), allocatable :: act_frac1 , dcrit
         
        
        allocate(act_frac1(1:n_mode))
        if(AllocateStatus /= 0) STOP "*** Not enough memory ***"
        allocate(dcrit(1:n_mode))
        if(AllocateStatus /= 0) STOP "*** Not enough memory ***"
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! initialise prognostic aerosol profiles:                                        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! use linear interpolation to put sounding on grid:
		do k=1,kp
		
		
		    if(aero_prof_flag) then
                iloc=find_pos(z_read(1:n_levels_s),z(k))
                iloc=min(n_levels_s-1,iloc)
                iloc=max(1,iloc)
                do i=1,n_mode
                    ! linear interp n_aer
                    call poly_int(z_read(iloc:iloc+1), n_read(i,iloc:iloc+1), &
                                min(z(k),z_read(n_levels_s)), var,dummy)
                    n_aer1(i)=var
                    ! linear interp sig_aer
                    call poly_int(z_read(iloc:iloc+1), sig_read(i,iloc:iloc+1), &
                                min(z(k),z_read(n_levels_s)), var,dummy)
                    sig_aer1(i)=var
                    ! linear interp d_aer
                    call poly_int(z_read(iloc:iloc+1), d_read(i,iloc:iloc+1), &
                                min(z(k),z_read(n_levels_s)), var,dummy)
                    d_aer1(i)=var
                enddo
            endif 
            
                        
            do i=1,n_mode-1 ! only fill external mixtures
                ! zeroth moment:
                q(k,(i-1)*3+2)=n_aer1(i) 
                ! surface area: 2nd moment x pi:
                q(k,(i-1)*3+3)= pi* ln_mom(2,n_aer1(i),sig_aer1(i),d_aer1(i))
                ! mass: 3rd moment x pi/6*rho:
                q(k,(i-1)*3+4)= pi/6._wp*density_core1(i)* &
                    ln_mom(3,n_aer1(i),sig_aer1(i),d_aer1(i))

            enddo
            
            
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! 1. find the critical diameter of each aerosol mode, and                    !
            ! 2. perform integration to set the aerosol n,s,m, in cloud water            !                         
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! initialise aerosol in cloud water    
            !  
            if(q(k,(n_mode-1)*6+3) .gt. 0._wp) then 

                
                call find_d_and_s_crits(p(k),t(k),q(k,(n_mode-1)*6+3),w,smax,dcrit)
                ! dcrit is set now
                ! partial moments of a lognormal distribution:
                ! see:
                ! https://math.stackexchange.com/questions/2055782/partial_expectations_of_lognormal_distributions
                do i=1,n_mode-1 ! only for external mixtures
                    ! number
                    q(k,(n_mode-1)*6+5+(i-1)*3)= &
                        ln_part_mom(0,dcrit(i),n_aer1(i),sig_aer1(i),d_aer1(i))
                    q(k,(i-1)*3+2)=q(k,(i-1)*3+2)-q(k,(n_mode-1)*6+5+(i-1)*3)
                    
                    
                    ! surface area
                    q(k,(n_mode-1)*6+6+(i-1)*3)= pi* &
                        ln_part_mom(2,dcrit(i),n_aer1(i),sig_aer1(i),d_aer1(i))
                    q(k,(i-1)*3+3)=q(k,(i-1)*3+3)-q(k,(n_mode-1)*6+6+(i-1)*3)
                    
                    
                    ! mass
                    q(k,(n_mode-1)*6+7+(i-1)*3)= pi/6._wp*density_core1(i)* &
                        ln_part_mom(3,dcrit(i),n_aer1(i),sig_aer1(i),d_aer1(i))
                    q(k,(i-1)*3+4)=q(k,(i-1)*3+4)-q(k,(n_mode-1)*6+7+(i-1)*3)
                enddo
                
                
            else
                smax=0._wp
                dcrit=1000._wp
            endif
        enddo
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



        
        



        

        deallocate(act_frac1)
        deallocate(dcrit)
                
    end subroutine p_initialise_aerosol_1d
    
    


	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>initialises aerosol profile
	!>@param[in] aero_prof_flag
	!>@param[in] nq, ncat: number of q variables, categories
	!>@param[in] c_s, c_e: start and end pointers for categories
	!>@param[in] inc: pointer to drop number category
	!>@param[in] ip, kp, o_halo: number of i, k and halo points
	!>@param[in] x,z,rho, p, t: grid values
	!>@param[inout] q, q_old: q_variables
    subroutine p_initialise_aerosol(aero_prof_flag, nq,ncat,c_s,c_e, &
                inc, ip,kp,o_halo, x,z,rho,p,t,q,q_old)
                
        use bam, only : n_mode, n_sv, n_aer1, d_aer1, sig_aer1, density_core1, &
                    nu_core1, molw_core1, org_content1, molw_org1, density_org1, &
                    delta_h_vap1, nu_org1, log_c_star1, a_eq_7, b_eq_7, &
                    initialise_arrays, ctmm_activation,find_d_and_s_crits
        implicit none
        ! arguments:
        logical :: aero_prof_flag
        integer(i4b), intent(in) :: nq, ncat, inc, ip, kp, o_halo
        integer(i4b), dimension(ncat), intent(in) :: c_s, c_e
        real(wp), dimension(-o_halo+1:ip+o_halo), intent(in) :: x
        real(wp), dimension(-o_halo+1:kp+o_halo), intent(in) :: z
        real(wp), dimension(-o_halo+1:kp+o_halo,-o_halo+1:ip+o_halo), intent(in) :: &
    					rho, p, t
        real(wp), dimension(-o_halo+1:kp+o_halo,-o_halo+1:ip+o_halo,nq), &
            intent(inout) :: q, q_old

        ! local variables
        integer(i4b) :: i, k, AllocateStatus, iloc
        real(wp) :: w, smax, phi, xx, kmom, var, dummy
        real(wp), dimension(:), allocatable :: act_frac1 , dcrit
         
        
        allocate(act_frac1(1:n_mode))
        if(AllocateStatus /= 0) STOP "*** Not enough memory ***"
        allocate(dcrit(1:n_mode))
        if(AllocateStatus /= 0) STOP "*** Not enough memory ***"
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! initialise prognostic aerosol profiles:                                        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! use linear interpolation to put sounding on grid:
		do k=1,kp
		
		
		    if(aero_prof_flag) then
                iloc=find_pos(z_read(1:n_levels_s),z(k))
                iloc=min(n_levels_s-1,iloc)
                iloc=max(1,iloc)
                do i=1,n_mode
                    ! linear interp n_aer
                    call poly_int(z_read(iloc:iloc+1), n_read(i,iloc:iloc+1), &
                                min(z(k),z_read(n_levels_s)), var,dummy)
                    n_aer1(i)=var
                    ! linear interp sig_aer
                    call poly_int(z_read(iloc:iloc+1), sig_read(i,iloc:iloc+1), &
                                min(z(k),z_read(n_levels_s)), var,dummy)
                    sig_aer1(i)=var
                    ! linear interp d_aer
                    call poly_int(z_read(iloc:iloc+1), d_read(i,iloc:iloc+1), &
                                min(z(k),z_read(n_levels_s)), var,dummy)
                    d_aer1(i)=var
                enddo
            endif 
            
                        
            do i=1,n_mode-1 ! only fill external mixtures
                ! zeroth moment:
                q(k,:,(i-1)*3+2)=n_aer1(i) 
                ! surface area: 2nd moment x pi:
                q(k,:,(i-1)*3+3)= pi* ln_mom(2,n_aer1(i),sig_aer1(i),d_aer1(i))
                ! mass: 3rd moment x pi/6*rho:
                q(k,:,(i-1)*3+4)= pi/6._wp*density_core1(i)* &
                    ln_mom(3,n_aer1(i),sig_aer1(i),d_aer1(i))

            enddo
            
            
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! 1. find the critical diameter of each aerosol mode, and                    !
            ! 2. perform integration to set the aerosol n,s,m, in cloud water            !                         
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! initialise aerosol in cloud water    
            !  
            if(q(k,1,(n_mode-1)*6+3) .gt. 0._wp) then 

                
                call find_d_and_s_crits(p(k,1),t(k,1),q(k,1,(n_mode-1)*6+3),w,smax,dcrit)
                q(k,:,(n_mode-1)*6+3)=q(k,1,(n_mode-1)*6+3)
                ! dcrit is set now
                ! partial moments of a lognormal distribution:
                ! see:
                ! https://math.stackexchange.com/questions/2055782/partial_expectations_of_lognormal_distributions
                do i=1,n_mode-1 ! only fill external mixtures
                    ! number
                     ! number
                    q(k,:,(n_mode-1)*6+5+(i-1)*3)= &
                        ln_part_mom(0,dcrit(i),n_aer1(i),sig_aer1(i),d_aer1(i))
                    q(k,:,(i-1)*3+2)=q(k,:,(i-1)*3+2)-q(k,:,(n_mode-1)*6+5+(i-1)*3)
                    
                    
                    ! surface area
                    q(k,:,(n_mode-1)*6+6+(i-1)*3)= pi* &
                        ln_part_mom(2,dcrit(i),n_aer1(i),sig_aer1(i),d_aer1(i))
                    q(k,:,(i-1)*3+3)=q(k,:,(i-1)*3+3)-q(k,:,(n_mode-1)*6+6+(i-1)*3)
                    
                    
                    ! mass
                    q(k,:,(n_mode-1)*6+7+(i-1)*3)= pi/6._wp*density_core1(i)* &
                        ln_part_mom(3,dcrit(i),n_aer1(i),sig_aer1(i),d_aer1(i))
                    q(k,:,(i-1)*3+4)=q(k,:,(i-1)*3+4)-q(k,:,(n_mode-1)*6+7+(i-1)*3)
                enddo
                
                
            else
                smax=0._wp
                dcrit=1000._wp
            endif
        enddo
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





        
        q_old=q ! previous value same

        deallocate(act_frac1)
        deallocate(dcrit)
                
    end subroutine p_initialise_aerosol
    
    

	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>initialises aerosol profile - 3d version
	!>@param[in] aero_prof_flag
	!>@param[in] nq, ncat: number of q variables, categories
	!>@param[in] c_s, c_e: start and end pointers for categories
	!>@param[in] inc: pointer to drop number category
	!>@param[in] ip, jp,kp, o_halo: number of i, k and halo points
	!>@param[in] x,y,z,rho, p, t: grid values
	!>@param[inout] q: q_variables
    subroutine p_initialise_aerosol_3d(aero_prof_flag, nq,ncat,c_s,c_e, &
                inc, ip,jp,kp,o_halo, x,y,z,rho,p,t,q)
                
        use bam, only : n_mode, n_sv, n_aer1, d_aer1, sig_aer1, density_core1, &
                    nu_core1, molw_core1, org_content1, molw_org1, density_org1, &
                    delta_h_vap1, nu_org1, log_c_star1, a_eq_7, b_eq_7, &
                    initialise_arrays, ctmm_activation,find_d_and_s_crits
        implicit none
        ! arguments:
        logical :: aero_prof_flag
        integer(i4b), intent(in) :: nq, ncat, inc, ip, jp, kp, o_halo
        integer(i4b), dimension(ncat), intent(in) :: c_s, c_e
        real(wp), dimension(-o_halo+1:ip+o_halo), intent(in) :: x
        real(wp), dimension(-o_halo+1:jp+o_halo), intent(in) :: y
        real(wp), dimension(-o_halo+1:kp+o_halo), intent(in) :: z
        real(wp), dimension(-o_halo+1:kp+o_halo), intent(in) :: rho, p, t
        real(wp), dimension(-o_halo+1:kp+o_halo,-o_halo+1:jp+o_halo,-o_halo+1:ip+o_halo,nq), &
            intent(inout) :: q

        ! local variables
        integer(i4b) :: i, k, AllocateStatus, iloc
        real(wp) :: w, smax, phi, xx, kmom, var, dummy
        real(wp), dimension(:), allocatable :: act_frac1 , dcrit
         
        
        allocate(act_frac1(1:n_mode))
        if(AllocateStatus /= 0) STOP "*** Not enough memory ***"
        allocate(dcrit(1:n_mode))
        if(AllocateStatus /= 0) STOP "*** Not enough memory ***"
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! initialise prognostic aerosol profiles:                                        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! use linear interpolation to put sounding on grid:
		do k=1,kp
		
		
		    if(aero_prof_flag) then
                iloc=find_pos(z_read(1:n_levels_s),z(k))
                iloc=min(n_levels_s-1,iloc)
                iloc=max(1,iloc)
                do i=1,n_mode
                    ! linear interp n_aer
                    call poly_int(z_read(iloc:iloc+1), n_read(i,iloc:iloc+1), &
                                min(z(k),z_read(n_levels_s)), var,dummy)
                    n_aer1(i)=var
                    ! linear interp sig_aer
                    call poly_int(z_read(iloc:iloc+1), sig_read(i,iloc:iloc+1), &
                                min(z(k),z_read(n_levels_s)), var,dummy)
                    sig_aer1(i)=var
                    ! linear interp d_aer
                    call poly_int(z_read(iloc:iloc+1), d_read(i,iloc:iloc+1), &
                                min(z(k),z_read(n_levels_s)), var,dummy)
                    d_aer1(i)=var
                enddo
            endif 
            
                        
            do i=1,n_mode-1 ! only fill external mixtures
                ! zeroth moment:
                q(k,:,:,(i-1)*3+2)=n_aer1(i) 
                ! surface area: 2nd moment x pi:
                q(k,:,:,(i-1)*3+3)= pi* ln_mom(2,n_aer1(i),sig_aer1(i),d_aer1(i))
                ! mass: 3rd moment x pi/6*rho:
                q(k,:,:,(i-1)*3+4)= pi/6._wp*density_core1(i)* &
                    ln_mom(3,n_aer1(i),sig_aer1(i),d_aer1(i))

            enddo
            
            
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! 1. find the critical diameter of each aerosol mode, and                    !
            ! 2. perform integration to set the aerosol n,s,m, in cloud water            !                         
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! initialise aerosol in cloud water    
            !  
            if(q(k,1,1,(n_mode-1)*6+3) .gt. 0._wp) then 

                
                call find_d_and_s_crits(p(k),t(k),&
                    q(k,1,1,(n_mode-1)*6+3),w,smax,dcrit)
                
                q(k,:,:,(n_mode-1)*6+3)=q(k,1,1,(n_mode-1)*6+3)
                ! dcrit is set now
                ! partial moments of a lognormal distribution:
                ! see:
                ! https://math.stackexchange.com/questions/2055782/partial_expectations_of_lognormal_distributions
                do i=1,n_mode-1 ! only fill external mixtures
                    ! number
                     ! number
                    q(k,:,:,(n_mode-1)*6+5+(i-1)*3)= &
                        ln_part_mom(0,dcrit(i),n_aer1(i),sig_aer1(i),d_aer1(i))
                    q(k,:,:,(i-1)*3+2)=q(k,:,:,(i-1)*3+2)-q(k,:,:,(n_mode-1)*6+5+(i-1)*3)
                    
                    
                    ! surface area
                    q(k,:,:,(n_mode-1)*6+6+(i-1)*3)= pi* &
                        ln_part_mom(2,dcrit(i),n_aer1(i),sig_aer1(i),d_aer1(i))
                    q(k,:,:,(i-1)*3+3)=q(k,:,:,(i-1)*3+3)-q(k,:,:,(n_mode-1)*6+6+(i-1)*3)
                    
                    
                    ! mass
                    q(k,:,:,(n_mode-1)*6+7+(i-1)*3)= pi/6._wp*density_core1(i)* &
                        ln_part_mom(3,dcrit(i),n_aer1(i),sig_aer1(i),d_aer1(i))
                    q(k,:,:,(i-1)*3+4)=q(k,:,:,(i-1)*3+4)-q(k,:,:,(n_mode-1)*6+7+(i-1)*3)
                enddo
                
                
            else
                smax=0._wp
                dcrit=1000._wp
            endif
        enddo
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





        

        deallocate(act_frac1)
        deallocate(dcrit)
                
    end subroutine p_initialise_aerosol_3d
    
    






	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>initialises variables for use with the microphysics
    subroutine initialise_microphysics_vars
    use hypergeo, only : hygfx
    implicit none

    small_number=epsilon(small_number)
    
	! used to calculate intercept and slopes
	gam1r=gamma(1._wp+alpha_r)
	gam2r=gamma(1._wp+alpha_r+dr)
	gam3r=gamma(5._wp+alpha_r)
	gam1c=gamma(1._wp+alpha_c)
	gam2c=gamma(1._wp+alpha_c+1._wp) ! note the 1, instead of dc - drop distribution
	                                ! is a mass distribution
	gam3c=gamma(4._wp/dc+alpha_c+1._wp) ! note different to rain 
	                                    !because cloud is a mass distribution
	gam1i=gamma(1._wp+alpha_i)
	gam2i=gamma(1._wp+alpha_i+di)
	gam1s=gamma(1._wp+alpha_s)
	gam2s=gamma(1._wp+alpha_s+ds)
	gam1g=gamma(1._wp+alpha_g)
	gam2g=gamma(1._wp+alpha_g+dg)
	
	! Raphael's "Magic" function representation of fall-speeds
	! Derived from Heymsfield and Westbrook (2010)
	gam3ai=1._wp/gamma(alpha_i+di/2._wp)
	gam3bi=4._wp/gamma(alpha_i+di)
	gam4ai=1._wp/gamma(alpha_i+di+di/2._wp)
	gam4bi=4._wp/gamma(alpha_i+2._wp*di)
	gam5ai=1._wp/gamma(alpha_i+2._wp+0.5_wp*di)
	gam5bi=4._wp/gamma(alpha_i+2._wp+di)
	gam6ai=1._wp/gamma(alpha_i+di/2._wp+2._wp)
	gam6bi=4._wp/gamma(alpha_i+di+2._wp)
	fall_q_i_hw=1._wp/ gamma(1._wp+alpha_i+di)
	fall_n_i_hw=1._wp/ gamma(1._wp+alpha_i)
	

    ! mass weighted fall for r, c, s, g, i
    fall_q_r=a_r*gamma(1._wp+alpha_r+dr+b_r) / gamma(1._wp+alpha_r+dr)
    fall_q_c=a_c*gamma(1._wp+alpha_c+1._wp+b_c) / gamma(1._wp+alpha_c+1._wp)
    fall_q_s=a_s*gamma(1._wp+alpha_s+ds+b_s) / gamma(1._wp+alpha_s+ds)
    fall_q_g=a_g*gamma(1._wp+alpha_g+dg+b_g) / gamma(1._wp+alpha_g+dg)
    fall_q_i=a_i*gamma(1._wp+alpha_i+di+b_i) / gamma(1._wp+alpha_i+di)

    ! number weighted fall for r, c, s, g
    fall_n_r=a_r*gamma(1._wp+alpha_r+b_r) / gamma(1._wp+alpha_r)
    fall_n_c=a_c*gamma(1._wp+alpha_c+b_c) / gamma(1._wp+alpha_c)
    fall_n_s=a_s*gamma(1._wp+alpha_s+b_s) / gamma(1._wp+alpha_s)
    fall_n_g=a_g*gamma(1._wp+alpha_g+b_g) / gamma(1._wp+alpha_g)
    fall_n_i=a_i*gamma(1._wp+alpha_i+b_i) / gamma(1._wp+alpha_i)
    
    ! sweep out of rain
    phi_r=pi*a_r*gamma(3._wp+b_r+alpha_r) / 4._wp
    
    ! ice accreting rain
    mass_iacr=pi*eri*a_r*cr*gamma(3._wp+b_r+dr+alpha_r)/4._wp
    num_iacr =pi*eri*a_r*gamma(3._wp+b_r+alpha_r)/4._wp

	! collection of cloud by snow and ice
	mass_sacw_i=pi*a_s*gamma(3._wp+b_s+alpha_s)/4._wp
	mass_iacw=pi*a_i*gamma(3._wp+b_i+alpha_i)/4._wp
	
	! collisions between precipitating particles of different species
	! rain-snow
	mass_raci1=gamma(1._wp+alpha_r)*gamma(3._wp+alpha_i+di)
	mass_raci2=2._wp*gamma(2._wp+alpha_r)*gamma(2._wp+alpha_i+di)
	mass_raci3=gamma(3._wp+alpha_r)*gamma(1._wp+alpha_i+di)
	
	num_raci1=gamma(1._wp+alpha_r)*gamma(3._wp+alpha_i)
	num_raci2=2._wp*gamma(2._wp+alpha_r)*gamma(2._wp+alpha_i)
	num_raci3=gamma(3._wp+alpha_r)*gamma(1._wp+alpha_i)
    ! rain-graupel
	mass_racg1=gamma(1._wp+alpha_r)*gamma(3._wp+alpha_g+dg)
	mass_racg2=2._wp*gamma(2._wp+alpha_r)*gamma(2._wp+alpha_g+dg)
	mass_racg3=gamma(3._wp+alpha_r)*gamma(1._wp+alpha_g+dg)

	num_racg1=gamma(1._wp+alpha_r)*gamma(3._wp+alpha_g)
	num_racg2=2._wp*gamma(2._wp+alpha_r)*gamma(2._wp+alpha_g)
	num_racg3=gamma(3._wp+alpha_r)*gamma(1._wp+alpha_g)
    ! snow-rain
	mass_iacr1=gamma(1._wp+alpha_i)*gamma(3._wp+alpha_r+dr)
	mass_iacr2=2._wp*gamma(2._wp+alpha_i)*gamma(2._wp+alpha_r+dr)
	mass_iacr3=gamma(3._wp+alpha_i)*gamma(1._wp+alpha_r+dr)

	num_iacr1=gamma(1._wp+alpha_i)*gamma(3._wp+alpha_r)
	num_iacr2=2._wp*gamma(2._wp+alpha_i)*gamma(2._wp+alpha_r)
	num_iacr3=gamma(3._wp+alpha_i)*gamma(1._wp+alpha_r)


    ! snow-graupel
	mass_sacg1=gamma(1._wp+alpha_s)*gamma(3._wp+alpha_g+dg)
	mass_sacg2=2._wp*gamma(2._wp+alpha_s)*gamma(2._wp+alpha_g+dg)
	mass_sacg3=gamma(3._wp+alpha_s)*gamma(1._wp+alpha_g+dg)

	num_sacg1=gamma(1._wp+alpha_s)*gamma(3._wp+alpha_g)
	num_sacg2=2._wp*gamma(2._wp+alpha_s)*gamma(2._wp+alpha_g)
	num_sacg3=gamma(3._wp+alpha_s)*gamma(1._wp+alpha_g)
    ! graupel-rain
	mass_gacr1=gamma(1._wp+alpha_g)*gamma(3._wp+alpha_r+dr)
	mass_gacr2=2._wp*gamma(2._wp+alpha_g)*gamma(2._wp+alpha_r+dr)
	mass_gacr3=gamma(3._wp+alpha_g)*gamma(1._wp+alpha_r+dr)
    ! graupel-snow
	mass_gacs1=gamma(1._wp+alpha_g)*gamma(3._wp+alpha_s+ds)
	mass_gacs2=2._wp*gamma(2._wp+alpha_g)*gamma(2._wp+alpha_s+ds)
	mass_gacs3=gamma(3._wp+alpha_g)*gamma(1._wp+alpha_s+ds)
	
	
	! accretion and riming by graupel
	mass_gacw=pi*egw*a_g*gamma(3._wp+b_g+alpha_g)/4._wp
	mass_gaci=pi*a_g*gamma(3._wp+b_g+alpha_g)/4._wp
    
    ! gauss hypergeometric equations aggregation of ice with ice (and snow with snow)
    ! See Ferrier (1994, JAS part 1, equation B.21)
    ! snow:
    a=1._wp
    b=4._wp+2._wp*alpha_s+b_s
    isnow=0._wp
    do k=1,3
	    call hygfx(a, b, real(k,wp)+alpha_s+1.0_wp,0.5_wp,f1)
	    call hygfx(a, b, real(k,wp)+alpha_s+b_s+1.0_wp, 0.5_wp,f2)
	    isnow=isnow+c(k)*(f1/(real(k,wp)+alpha_s)-f2/(real(k,wp)+alpha_s+b_s))
	enddo
	isnow=a_s*pi*gamma(b)/(2._wp**(6._wp+2._wp*alpha_s+b_s)) * isnow
	
    ! ice:
    a=1._wp
    b=4._wp+2._wp*alpha_i+b_i
    iice=0._wp
    do k=1,3
	    call hygfx(a, b, real(k,wp)+alpha_i+1.0_wp,0.5_wp,f1)
	    call hygfx(a, b, real(k,wp)+alpha_i+b_i+1.0_wp, 0.5_wp,f2)
	    iice=iice+c(k)*(f1/(real(k,wp)+alpha_i)-f2/(real(k,wp)+alpha_i+b_i))
	enddo
	iice=a_i*pi*gamma(b)/(2._wp**(6._wp+2._wp*alpha_i+b_i)) * iice
	
    ! ice for Heymsfield and Westbrook (2010):
    a=1._wp
    b=3._wp+2._wp*alpha_i+di
    iice2=0._wp
    do k=1,3
	    call hygfx(a, b, real(k,wp)+alpha_i+1.0_wp,0.5_wp,f1)
	    call hygfx(a, b, real(k,wp)+alpha_i+di, 0.5_wp,f2)
	    iice2=iice2+c(k)*(f1/(real(k,wp)+alpha_i)-f2/(real(k,wp)+alpha_i+di-1._wp))
	enddo
	iice2=pi*gamma(b)/(2._wp**(5._wp+2._wp*alpha_i+di)) * iice2
		
	! ventilation
	! rain:
	nu_r1=0.78_wp*gamma(2._wp+alpha_r)
	nu_r2=gamma(0.5_wp*b_r+alpha_r+2.5)
	! ice:
	nu_i1=0.78_wp*gamma(2._wp+alpha_i)
	nu_i2=0.31_wp*gamma(0.5_wp*b_i+alpha_i+2.5)
	! snow:
	nu_s1=0.78_wp*gamma(2._wp+alpha_s)
	nu_s2=0.31_wp*gamma(0.5_wp*b_s+alpha_s+2.5)
	! graupel:
	nu_g1=0.78_wp*gamma(2._wp+alpha_g)
	nu_g2=0.31_wp*gamma(0.5_wp*b_g+alpha_g+2.5)
	
	! immersion freezing by bigg
	mass_imm=gamma(4._wp+dr+alpha_r)*pi*cr*bbigg/6._wp
	num_imm=gamma(4._wp+alpha_r)*pi*bbigg/6._wp
	
	
	! precipitation
	chi_rain=gamma(1._wp+alpha_r+b_r+dr)
	chi_cloud=gamma(1._wp+alpha_c+b_c+1._wp)
	chi_ice=gamma(1._wp+alpha_i+b_i+di)
	chi_num_ice=gamma(1._wp+alpha_i+b_i)
	chi_snow=gamma(1._wp+alpha_s+b_s+ds)
	chi_graupel=gamma(1._wp+alpha_g+b_g+dg)
	
	chi_rain1=gamma(1._wp+alpha_r+dr)
	chi_cloud1=gamma(1._wp+alpha_c+1._wp)
	chi_ice1=gamma(1._wp+alpha_i+di)
	chi_num_ice1=gamma(1._wp+alpha_i)
	chi_snow1=gamma(1._wp+alpha_s+ds)
	chi_graupel1=gamma(1._wp+alpha_g+dg)
	

    ! Seifert and Beheng autoconversion:
    kc=9.44e9_wp ! m3 kg-2 s-1
    kr=5.78e0_wp ! m3 kg-2 s-1
    xstar=2.6e-10_wp ! kg
    
    
    ! mode 2 multiplication - specify the limit of integration for gamma distribution
    pthreshr=invgammainc(probthresh,alpha_r+1.0_wp)
    pthreshi=invgammainc(probthresh,alpha_i+1.0_wp)
    
    
    ! for radiation - converting number-mass to number-diameter
    gam1cr=gamma(alpha_c+1._wp+1._wp/dc)
    gam2cr=gamma(alpha_c+1._wp+2._wp/dc)
    end subroutine initialise_microphysics_vars
    
    
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculate the size distribution parameters
	!>@param[in] nq: number of q-fields
	!>@param[in] ncat: number of categories
	!>@param[in] n_mode: number of aerosol modes
	!>@param[in] cst,cen: indices of categories
	!>@param[in] inc, iqc: index of cloud number, index of cloud mass
	!>@param[in] inr, iqr: index of rain number, index of rain mass
	!>@param[in] ini, iqi,iai: index of ice number, index of ice mass, and ice aerosol
	!>@param[in] cat_am,cat_c, cat_r, cat_i: category index for cloud and rain and ice
	!>@param[in] ip,jp: number of horizontal levels
	!>@param[in] kp: number of vertical levels
	!>@param[in] o_halo: extra points for advection
	!>@param[inout] q: q-variables 
	!>@param[inout] nrad,ngs,lamgs,mugs- needed for radiation
	!>@param[in] rho, rhon: density 
    subroutine calculate_gamma_params(nq,ncat,n_mode,cst,cen,inc,iqc, inr,iqr,ini,iqi,iai, &
                    cat_am,cat_c, cat_r, cat_i,&
                    ip,jp,kp,l_h,r_h,q,&
                    nrad,ngs,lamgs,mugs, &
                    rhoan,ice_flag)
    implicit none
    ! arguments:
    integer(i4b), intent(in) :: nq, ncat, n_mode, ip,jp,kp, inc, iqc, inr,iqr,&
        ini,iqi,iai, &
        cat_am,&
        cat_c, cat_r,cat_i,l_h,r_h
    integer(i4b), dimension(ncat), intent(in) :: cst,cen
    real(wp), dimension(-l_h+1:kp+r_h,-l_h+1:jp+r_h,-l_h+1:ip+r_h,nq), intent(inout) :: q

    real(wp), dimension(-l_h+1:kp+r_h), intent(in) :: rhoan
    logical, intent(in) :: ice_flag
    
    integer(i4b), intent(in) :: nrad
	real(wp), intent(inout), dimension(-l_h+1:kp+r_h,-l_h+1:jp+r_h,-l_h+1:ip+r_h,nrad) :: ngs,lamgs,mugs

	! locals
	integer(i4b) :: i,j,k
	real(wp) :: p
	real(wp), dimension(-l_h+1:kp+r_h,-l_h+1:jp+r_h,-l_h+1:ip+r_h,4) :: moms

    ! RAIN
    do i=1-l_h,ip+r_h
        do j=1-l_h,jp+r_h
            do k=1-l_h,kp+r_h
                ! rain n0, lambda
                lamgs(k,j,i,2)=(max(q(k,j,i,cst(cat_r)),1._wp)*cr*gam2r / &
                        (max(q(k,j,i,cst(cat_r)+1),1.e-10_wp)*gam1r))**(1._wp/dr)
            enddo
        enddo
    enddo
    do i=1-l_h,ip+r_h
        do j=1-l_h,jp+r_h
            do k=1-l_h,kp+r_h
                ! rain n0, lambda
                ngs(k,j,i,2)=rhoan(k)* &
                    max(q(k,j,i,cst(cat_r)),0._wp)*&
                    lamgs(k,j,i,2)**(1._wp+alpha_r) / gam1r
            enddo
        enddo
    enddo
    do i=1-l_h,ip+r_h
        do j=1-l_h,jp+r_h
            do k=1-l_h,kp+r_h
                ! rain n0, lambda
                mugs(k,j,i,2)=alpha_r
            enddo
        enddo
    enddo
               
               
    ! CLOUD - note these params are for a mass-distribution. They need to be converted
    ! to be applicable for diameter distribution: Calculate the 1st and 2nd moments 
    ! (diameter)
    ! which enables calculation of shape parameter by taking ratios
    do i=1-l_h,ip+r_h
        do j=1-l_h,jp+r_h
            do k=1-l_h,kp+r_h
                ! cloud n0, lambda    
                lamgs(k,j,i,1)=(max(q(k,j,i,inc),1._wp)*gam2c / &
                    (max(q(k,j,i,iqc),1.e-10_wp)*gam1c))**(1._wp/1._wp)
            enddo
        enddo
    enddo
                    
    do i=1-l_h,ip+r_h
        do j=1-l_h,jp+r_h
            do k=1-l_h,kp+r_h
                ! cloud n0, lambda    
                ngs(k,j,i,1)=rhoan(k)*max(q(k,j,i,inc),0._wp)*&
                    lamgs(k,j,i,1)**(1._wp+alpha_c) / gam1c
            enddo
        enddo
    enddo
    ! cloud moments
    ! see https://journals.ametsoc.org/view/journals/atsc/68/7/2011jas3645.1.xml
    ! equation 17 - modified gamma distribution
    moms=0._wp
    do i=1-l_h,ip+r_h
        do j=1-l_h,jp+r_h
            do k=1-l_h,kp+r_h
                ! zeroth - number
                if (q(k,j,i,inc).lt.1._wp) cycle
                if (q(k,j,i,iqc).lt.1.e-20_wp) cycle
                moms(k,j,i,1)=rhoan(k)*q(k,j,i,inc)
                ! first 
                moms(k,j,i,2)=ngs(k,j,i,1)*gam1cr / &
                    (cc**(1._wp/dc)*lamgs(k,j,i,1)**(alpha_c+1._wp+1._wp/dc))
                ! second  
                moms(k,j,i,3)=ngs(k,j,i,1)*gam2cr / &
                    (cc**(2._wp/dc)*lamgs(k,j,i,1)**(alpha_c+1._wp+2._wp/dc))
                ! third 
                moms(k,j,i,4)=rhoan(k)*q(k,j,i,iqc)/cc
            enddo
        enddo
    enddo
    ! now, convert these moments to parameters
    do i=1-l_h,ip+r_h
        do j=1-l_h,jp+r_h
            do k=1-l_h,kp+r_h
                ! calculate the ratio p M1/M0*M2/M3
                p=moms(k,j,i,2)/moms(k,j,i,1)* &
                    moms(k,j,i,3)/moms(k,j,i,4)
                ! p=(alpha + 1) / (alpha+3)
                ! mu    
                mugs(k,j,i,1)=(1._wp-3._wp*p) / (p-1._wp)
                ! now calculate ratio of m0 : m1
                p=moms(k,j,i,1) / moms(k,j,i,2)
                ! lambda from this ratio
                lamgs(k,j,i,1)=p*(mugs(k,j,i,1)+1._wp)
                ! n0
                ngs(k,j,i,1)=lamgs(k,j,i,1)**(mugs(k,j,i,1)+1._wp)*moms(k,j,i,1) / &
                    gamma(mugs(k,j,i,1)+1._wp)
            enddo
        enddo
    enddo


    ! ICE - note, this should be updated to take into account variable density
    if(ice_flag) then
        do i=1-l_h,ip+r_h
            do j=1-l_h,jp+r_h
                do k=1-l_h,kp+r_h
                    ! ice n0, lambda
                    lamgs(k,j,i,3)=(max(q(k,j,i,ini),1._wp)*ci*gam2i / &
                        (max(q(k,j,i,iqi),1.e-10_wp)*gam1i))**(1._wp/di)
                enddo
            enddo
        enddo
        do i=1-l_h,ip+r_h
            do j=1-l_h,jp+r_h
                do k=1-l_h,kp+r_h
                    ! ice n0, lambda
                    ngs(k,j,i,3)=rhoan(k)*max(q(k,j,i,ini),0._wp)*&
                        lamgs(k,j,i,3)**(1._wp+alpha_i) / gam1i
                enddo
            enddo
        enddo
        do i=1-l_h,ip+r_h
            do j=1-l_h,jp+r_h
                do k=1-l_h,kp+r_h
                    ! ice n0, lambda
                    mugs(k,j,i,3)=alpha_i
                enddo
            enddo
        enddo
    endif
    end subroutine calculate_gamma_params    
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
#if MPI_PAMM == 0
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>solves one time-step of the microphysics - calls p_microphysics_1d
	!>@param[in] nq: number of q-fields
	!>@param[in] ncat: number of categories
	!>@param[in] n_mode: number of aerosol modes
	!>@param[in] cst,cen: indices of categories
	!>@param[in] inc, iqc: index of cloud number, index of cloud mass
	!>@param[in] inr, iqr: index of rain number, index of rain mass
	!>@param[in] ini, iqi,iai: index of ice number, index of ice mass, and ice aerosol
	!>@param[in] cat_am,cat_c, cat_r, cat_i: category index for cloud and rain and ice
	!>@param[in] nprec
	!>@param[in] ip,jp: number of horizontal levels
	!>@param[in] kp: number of vertical levels
	!>@param[in] dt: time-step
	!>@param[in] dz: dz, dzn
	!>@param[in] o_halo: extra points for advection
	!>@param[inout] q: q-variables 
	!>@param[inout] nrad,ngs,lamgs,mugs- needed for radiation
	!>@param[inout] precip: precip in rain, snow, graupel, ice cats - diagnostic
	!>@param[inout] theta: theta 
	!>@param[inout] p: pressure
	!>@param[in] z: vertical levels 
	!>@param[inout] thetan: potential temperature 
	!>@param[inout] rho, rhon: density 
	!>@param[in] w: vertical wind 
	!>@param[inout] micro_init: boolean to initialise microphysics 
	!>@param[in] hm_flag, wr_flag: switch hm-process and warm rain on and off
	!>@param[in] mass_ice: mass of a single ice crystal (override)
	!>@param[in] ice_flag: ice microphysics
	!>@param[in] theta_flag: whether to alter theta
	!>@param[in] j_stochastic, ice_nuc_flag, mode2_ice_flag, 
	!>@param[in] coll_breakup_flag1,heyms_west, lawson, recycle, calc_params
    subroutine p_microphysics_3d(nq,ncat,n_mode,cst,cen,inc,iqc, inr,iqr,ini,iqi,iai, &
                    cat_am,cat_c, cat_r, cat_i,&
                    nprec, &
                    ip,jp,kp,l_h,r_h,dt,dz,dzn,q, &
                    nrad,ngs,lamgs,mugs, &
                    precip,th,prefn, z,thetan,rhoa,rhoan,w, &
    				micro_init,hm_flag, wr_flag, mass_ice, ice_flag, theta_flag, &
    				j_stochastic,ice_nuc_flag,mode1_ice_flag,mode2_ice_flag, &
    				coll_breakup_flag1, heyms_west, lawson, recycle, calc_params)
#else
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>solves one time-step of the microphysics - calls p_microphysics_1d
	!>@param[in] nq: number of q-fields
	!>@param[in] ncat: number of categories
	!>@param[in] n_mode: number of aerosol modes
	!>@param[in] cst,cen: indices of categories
	!>@param[in] inc, iqc: index of cloud number, index of cloud mass
	!>@param[in] inr, iqr: index of rain number, index of rain mass
	!>@param[in] ini, iqi,iai: index of ice number, index of ice mass, and ice aerosol
	!>@param[in] cat_am,cat_c, cat_r, cat_i: category index for cloud and rain and ice
	!>@param[in] nprec
	!>@param[in] ip,jp: number of horizontal levels
	!>@param[in] kp: number of vertical levels
	!>@param[in] dt: time-step
	!>@param[in] dz: dz, dzn
	!>@param[in] o_halo: extra points for advection
	!>@param[inout] q: q-variables 
	!>@param[inout] nrad,ngs,lamgs,mugs- needed for radiation
	!>@param[inout] precip: precip in rain, snow, graupel, ice cats - diagnostic
	!>@param[inout] theta: theta 
	!>@param[inout] p: pressure
	!>@param[in] z: vertical levels 
	!>@param[inout] thetan: potential temperature 
	!>@param[inout] rho, rhon: density 
	!>@param[in] w: vertical wind 
	!>@param[inout] micro_init: boolean to initialise microphysics 
	!>@param[in] hm_flag, wr_flag: switch hm-process and warm rain on and off
	!>@param[in] mass_ice: mass of a single ice crystal (override)
	!>@param[in] ice_flag: ice microphysics
	!>@param[in] theta_flag: whether to alter theta
	!>@param[in] j_stochastic, ice_nuc_flag, mode1_ice_flag, mode2_ice_flag,
	!>@param[in] coll_breakup_flag1, heyms_west, lawson, recycle, calc_parms
	!>@param[in] comm,comm_vert,id,dims,coords: MPI variables
    subroutine p_microphysics_3d(nq,ncat,n_mode,cst,cen,inc,iqc, inr,iqr,ini,iqi,iai, &
                    cat_am,cat_c, cat_r, cat_i,&
                    nprec, &
                    ip,jp,kp,l_h,r_h,dt,dz,dzn,q,precip,&
                    nrad,ngs,lamgs,mugs, &
                    th,prefn, z,thetan,rhoa,rhoan,w, &
    				micro_init,hm_flag, wr_flag, mass_ice, ice_flag, theta_flag, &
    				j_stochastic,ice_nuc_flag, mode1_ice_flag, &
    				mode2_ice_flag, coll_breakup_flag1, &
    				heyms_west, lawson, recycle, calc_params, &
    				comm,comm_vert,id,dims,coords)
    use mpi
	use advection_s_3d, only : mpdata_vec_vert_3d, mpdata_vert_3d
	use mpi_module
#endif
    implicit none
    ! arguments:
    integer(i4b), intent(in) :: nq, ncat, n_mode, ip,jp,kp, inc, iqc, inr,iqr,&
        ini,iqi,iai, &
        cat_am,&
        cat_c, cat_r,cat_i,l_h,r_h, nprec
    integer(i4b), dimension(ncat), intent(in) :: cst,cen
    real(wp), intent(in) :: dt
    real(wp), dimension(-l_h+1:kp+r_h,-l_h+1:jp+r_h,-l_h+1:ip+r_h,nq), intent(inout) :: q
    real(wp), dimension(1:kp,1-l_h:jp+r_h,1-l_h:ip+r_h,nprec), intent(inout) :: precip
    real(wp), dimension(-l_h+1:kp+r_h,-l_h+1:jp+r_h,-l_h+1:ip+r_h), intent(inout) :: &
    					th
    real(wp), dimension(-l_h+1:kp+r_h), intent(in) :: z, dz, dzn, rhoa,rhoan, thetan, &
        prefn
    real(wp), dimension(-l_h+1:kp+r_h,-l_h+1:jp+r_h,-l_h+1:ip+r_h), intent(in) :: w
    logical, intent(in) :: ice_flag, hm_flag, wr_flag, theta_flag, calc_params, &
        heyms_west, lawson, recycle
    integer(i4b), intent(in) :: ice_nuc_flag, mode1_ice_flag, mode2_ice_flag, &
                                coll_breakup_flag1
    logical , intent(inout) :: micro_init
    real(wp), intent(in) :: mass_ice, j_stochastic
    
    integer(i4b), intent(in) :: nrad
	real(wp), intent(inout), &
	    dimension(-l_h+1:kp+r_h,-l_h+1:jp+r_h,-l_h+1:ip+r_h,nrad) :: ngs,lamgs,mugs

	! locals
	integer(i4b) :: i,j,k,n, error,n1
#if MPI_PAMM == 1
	real(wp), dimension(-l_h+1:kp+r_h,-l_h+1:jp+r_h,-l_h+1:ip+r_h) :: & 
	                vqr, vqc, vqi, vni
	integer(i4b), dimension(3) :: n_step,n_step_o,n_step_g
    integer(i4b), intent(in) :: id, comm,comm_vert
    integer(i4b), dimension(3), intent(in) :: coords, dims
    real(wp), dimension(nq) :: lbc,ubc
    logical, dimension(3) :: adv_lg, adv_l=[.false.,.false.,.false.], &
                    adv_l_o=[.false.,.false.,.false.], adv_lgg

    n_step=1
    n_step_o=1
    lbc=0._wp
    ubc=0._wp
#endif
	
	do i=1,ip
	    do j=1,jp
#if MPI_PAMM == 0 
    		call p_microphysics_1d(nq,ncat,n_mode,cst,cen,inc,iqc, inr,iqr, ini,iqi,iai,&
		                cat_am,cat_c, cat_r, cat_i,nprec,&
		                kp,l_h,dt,dz,dzn,q(:,j,i,:),precip(:,j,i,:),th(:,j,i),&
		                    prefn, &
							z(:),thetan,rhoa(:),rhoan(:),w(:,j,i), &
    						micro_init,hm_flag, mass_ice, ice_flag, &
    						wr_flag,.true.,theta_flag, &
    						j_stochastic,ice_nuc_flag,mode1_ice_flag,mode2_ice_flag, &
    						coll_breakup_flag1, heyms_west, lawson, recycle)
#else

    		call p_microphysics_1d(nq,ncat,n_mode,cst,cen,inc,iqc,inr,iqr, ini,iqi,iai, &
		                cat_am,cat_c, cat_r, cat_i,nprec, &
		                kp,l_h,dt,dz,dzn,q(:,j,i,:),precip(:,j,i,:),th(:,j,i),&
		                    prefn, &
							z(:),thetan,rhoa(:),rhoan(:),w(:,j,i), &
							vqc(:,j,i),vqr(:,j,i),vqi(:,j,i),vni(:,j,i),n_step, adv_l, &
							coords, &
    						micro_init,hm_flag, mass_ice, ice_flag, &
    						wr_flag,.true.,theta_flag, &
    						j_stochastic,ice_nuc_flag,mode1_ice_flag, mode2_ice_flag, &
    						coll_breakup_flag1, heyms_west, lawson, recycle)
    		n_step_o=max(n_step,n_step_o)

    		adv_l_o=adv_l_o .or. adv_l ! if there has been a true at any point, 
    		                           ! set adv_l_o to true on this PE
#endif
    	enddo
	enddo
	
	! collective communication - no longer really need this because
	! vertical levels are not parallelised
	! (but you may want to keep the option of parallelising over the vertical)
	! normally the communicator should be comm_vert
	! but there is problem in the code  when exchange_full is called
	! because not all processors might be there (due to if statements)
#if MPI_PAMM == 1
	call mpi_allreduce(adv_l_o(1:3),adv_lg(1:3),3,MPI_LOGICAL,MPI_LOR, comm_vert,error)
	call mpi_allreduce(n_step_o(1:3),n_step_g(1:3),3,MPI_INTEGER,MPI_MAX, comm_vert,error)
	call mpi_allreduce(adv_l_o(1:3),adv_lgg(1:3),3,MPI_LOGICAL,MPI_LOR, comm,error)
#endif
    


#if MPI_PAMM == 1
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! full exchange needed                                                               !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call exchange_full(comm, id, kp, jp, ip, r_h,r_h,r_h,r_h,r_h,r_h, &
                            th,0._wp,0._wp,dims,coords)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! full exchange needed                                                               !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do n1=1,nq
        call exchange_full(comm, id, kp, jp, ip, r_h,r_h,r_h,r_h,r_h,r_h, &
                                q(:,:,:,n1),lbc(n1),ubc(n1),dims,coords)
    enddo
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Vert communications
    if(adv_lg(1)) then
        call exchange_along_z(comm_vert, id, kp, jp, ip, r_h,r_h,r_h,r_h,r_h,r_h, &
                                vqc(:,:,:),0._wp,0._wp,dims,coords)
        do n=1,n_step_g(1)
            call mpdata_vec_vert_3d(dt/real(n_step_g(1),wp),dz,dzn,&
                    rhoa,rhoan, &
                    ip,jp,kp,cen(cat_c)-cst(cat_c)+1,l_h,r_h,&
                    vqc,q(:,:,:,cst(cat_c):cen(cat_c)),&
                    lbc(cst(cat_c):cen(cat_c)),ubc(cst(cat_c):cen(cat_c)), &
                    1,.false., 2,comm_vert, id, &
                    dims,coords)
        enddo
    endif       
    if(adv_lg(2)) then
        call exchange_along_z(comm_vert, id, kp, jp, ip, r_h,r_h,r_h,r_h,r_h,r_h, &
                                vqr(:,:,:),0._wp,0._wp,dims,coords)
        do n=1,n_step_g(2)

            call mpdata_vec_vert_3d(dt/real(n_step_g(2),wp),dz,dzn,&
                    rhoa,rhoan, &
                    ip,jp,kp,cen(cat_r)-cst(cat_r)+1,l_h,r_h,&
                    vqr,q(:,:,:,cst(cat_r):cen(cat_r)),&
                    lbc(cst(cat_r):cen(cat_r)),ubc(cst(cat_r):cen(cat_r)), &
                    1,.false., 2,comm_vert, id, &
                    dims,coords) 

        enddo
    endif       
    if(adv_lg(3)) then
        call exchange_along_z(comm_vert, id, kp, jp, ip, r_h,r_h,r_h,r_h,r_h,r_h, &
                                vqi(:,:,:),0._wp,0._wp,dims,coords)
        do n=1,n_step_g(3)

            call mpdata_vec_vert_3d(dt/real(n_step_g(3),wp),dz,dzn,&
                    rhoa,rhoan, &
                    ip,jp,kp,cen(cat_i)-cst(cat_i)+1,l_h,r_h,&
                    vqi,q(:,:,:,cst(cat_i):cen(cat_i)),&
                    lbc(cst(cat_i):cen(cat_i)),ubc(cst(cat_i):cen(cat_i)), &
                    1,.false., 2,comm_vert, id, &
                    dims,coords) 

        enddo
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! comm - communications
    if(adv_lgg(1)) then
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! full exchange needed                                                           !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do n1=cst(cat_c),cen(cat_c)
            call exchange_full(comm, id, kp, jp, ip, r_h,r_h,r_h,r_h,r_h,r_h, &
                                    q(:,:,:,n1),lbc(n1),ubc(n1),dims,coords)
        enddo
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    endif
           
    if(adv_lgg(2)) then
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! full exchange needed                                                           !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do n1=cst(cat_r),cen(cat_r)
            call exchange_full(comm, id, kp, jp, ip, r_h,r_h,r_h,r_h,r_h,r_h, &
                                    q(:,:,:,n1),lbc(n1),ubc(n1),dims,coords)
        enddo
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    endif
           
    if(adv_lgg(3)) then
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! full exchange needed                                                           !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do n1=cst(cat_i),cen(cat_i)
            call exchange_full(comm, id, kp, jp, ip, r_h,r_h,r_h,r_h,r_h,r_h, &
                                    q(:,:,:,n1),lbc(n1),ubc(n1),dims,coords)
        enddo
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    endif
           
#endif

    if (calc_params) then
        call calculate_gamma_params(nq,ncat,n_mode,cst,cen,inc,iqc, inr,iqr,ini,iqi,iai, &
                    cat_am,cat_c, cat_r, cat_i,&
                    ip,jp,kp,l_h,r_h,q,&
                    nrad,ngs,lamgs,mugs, &
                    rhoan,ice_flag)
    endif
    
	end subroutine p_microphysics_3d
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	
	
    
    
    
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>solves one time-step of the microphysics - calls p_microphysics_1d
	!>@param[in] nq: number of q-fields
	!>@param[in] ncat: number of categories
	!>@param[in] n_mode: number of aerosol modes
	!>@param[in] cst,cen: indices of categories
	!>@param[in] inc, iqc: index of cloud number, index of cloud mass
	!>@param[in] inr, iqr: index of rain number, index of rain mass
	!>@param[in] ini, iqi,iai: index of ice number, index of ice mass, and ice aerosol
	!>@param[in] cat_am,cat_c, cat_r,cat_i: category index for cloud and rain and ice
	!>@param[in] nprec
	!>@param[in] ip: number of horizontal levels
	!>@param[in] kp: number of vertical levels
	!>@param[in] dt: time-step
	!>@param[in] dz: dz, dzn
	!>@param[in] o_halo: extra points for advection
	!>@param[inout] q: q-variables 
	!>@param[inout] precip: precip in rain, snow, graupel, ice cats - diagnostic
	!>@param[inout] theta: theta 
	!>@param[inout] p: pressure
	!>@param[in] z: vertical levels 
	!>@param[in] theta_ref: reference potential temperature 
	!>@param[inout] rho, rhon: density 
	!>@param[in] w: vertical wind 
	!>@param[inout] micro_init: boolean to initialise microphysics 
	!>@param[in] hm_flag: switch hm-process on and off
	!>@param[in] mass_ice: mass of a single ice crystal (override)
	!>@param[in] ice_flag: ice on /off
	!>@param[in] wr_flag: warm rain on /off
	!>@param[in] rm_flag: riming on /off
	!>@param[in] theta_flag: whether to alter theta
	!>@param[in] j_stochastic, ice_nuc_flag, mode1_ice_flag, mode2_ice_flag, 
	!>@param[in] coll_breakup_flag1, heyms_west, lawson, recycle
    subroutine p_microphysics_2d(nq,ncat,n_mode,cst,cen,inc,iqc,inr,iqr,ini,iqi,iai, &
                    cat_am,cat_c, cat_r, cat_i,nprec, &
                    ip,kp,o_halo,dt,dz,dzn,q,precip,theta,p, z,theta_ref,rho,rhon,w, &
    						micro_init,hm_flag, mass_ice, ice_flag, &
    						wr_flag, rm_flag, theta_flag, &
    						j_stochastic,ice_nuc_flag,mode1_ice_flag,mode2_ice_flag, &
    						coll_breakup_flag1, heyms_west, lawson, recycle)
    implicit none
    ! arguments:
    integer(i4b), intent(in) :: nq, ncat, n_mode, ip,kp, o_halo, inc, iqc, inr,iqr, &
        ini,iqi,iai, &
        cat_am,cat_c, cat_r, cat_i,nprec
    integer(i4b), dimension(ncat), intent(in) :: cst,cen
    real(wp), intent(in) :: dt
    real(wp), dimension(-o_halo+1:kp+o_halo,-o_halo+1:ip+o_halo,nq), intent(inout) :: q
    real(wp), dimension(1:kp,1:ip,1:nprec), intent(inout) :: precip
    real(wp), dimension(-o_halo+1:kp+o_halo,-o_halo+1:ip+o_halo), intent(inout) :: &
    					theta, p,rho
    real(wp), dimension(-o_halo+1:kp+o_halo), intent(in) :: z, dz, dzn, &
                    rhon, theta_ref
    real(wp), dimension(-o_halo+1:kp+o_halo,-o_halo+1:ip+o_halo), intent(in) :: w
    logical, intent(in) :: hm_flag, ice_flag, wr_flag, rm_flag, theta_flag, heyms_west, &
                        lawson, recycle
    integer(i4b), intent(in) :: ice_nuc_flag, mode1_ice_flag, &
                                mode2_ice_flag, coll_breakup_flag1
    logical , intent(inout) :: micro_init
    real(wp), intent(in) :: mass_ice, j_stochastic


	! locals
	integer(i4b) :: i
#if MPI_PAMM == 1
    real(wp), dimension(-o_halo+1:kp+o_halo,-o_halo+1:ip+o_halo) :: vqc,vqr, vqi, vni
    integer(i4b), dimension(3) :: n_step, n_step_o
    logical, dimension(3) :: adv_l=[.false.,.false.,.false.]
    integer(i4b), dimension(3) :: coords
    
    n_step_o=1
#endif
	
	
	
	do i=1,ip
#if MPI_PAMM == 0 
		call p_microphysics_1d(nq,ncat,n_mode,cst,cen,inc,iqc,inr,iqr,ini,iqi,iai, &
		                cat_am,cat_c, cat_r, cat_i,nprec,&
		                kp,o_halo,dt,dz,dzn,q(:,i,:),precip(:,i,:),theta(:,i),p(:,i), &
							z(:),theta_ref,rho(:,i),rhon(:),w(:,i), &
    						micro_init,hm_flag, mass_ice, ice_flag, &
    						wr_flag,rm_flag,theta_flag, &
    						j_stochastic,ice_nuc_flag,mode1_ice_flag,mode2_ice_flag, &
    						coll_breakup_flag1, heyms_west,lawson, recycle)
#else
		call p_microphysics_1d(nq,ncat,n_mode,cst,cen,inc,iqc,inr,iqr, ini,iqi,iai,&
		                cat_am,cat_c, cat_r, cat_i,nprec, &
		                kp,o_halo,dt,dz,dzn,q(:,i,:),precip(:,i,:),theta(:,i),p(:,i), &
							z(:),theta_ref,rho(:,i),rhon(:),w(:,i), &
							vqc(:,i),vqr(:,i),vqi(:,i), vni(:,i), n_step, adv_l, coords,&
    						micro_init,hm_flag, mass_ice, ice_flag, &
    						wr_flag,rm_flag,theta_flag, &
    						j_stochastic,ice_nuc_flag,mode1_ice_flag,mode2_ice_flag, &
    						coll_breakup_flag1, heyms_west,lawson, recycle)
    	n_step_o=max(n_step_o,n_step)
#endif	
	enddo


	end subroutine p_microphysics_2d
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	
	
	
#if MPI_PAMM == 0	    
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>solves one time-step of the microphysics
	!>@param[in] nq: number of q-fields
	!>@param[in] ncat: number of categories
	!>@param[in] n_mode: number of aerosol modes
	!>@param[in] cst,cen: indices of categories
	!>@param[in] inc, iqc: index of cloud number, index of cloud mass
	!>@param[in] inr, iqr: index of rain number, index of rain mass
	!>@param[in] ini, iqi,iai: index of ice number, index of ice mass, and ice aerosol
	!>@param[in] cat_am,cat_c, cat_r, cat_i: category index for cloud and rain and ice
	!>@param[in] nprec
	!>@param[in] kp: number of vertical levels
	!>@param[in] dt: time-step
	!>@param[in] dz: dz, dzn
	!>@param[in] o_halo: extra points for advection
	!>@param[inout] q: q-variables 
	!>@param[inout] precip: precip in rain, snow, graupel, ice cats - diagnostic
	!>@param[inout] th: theta perturbation
	!>@param[inout] p: pressure
	!>@param[inout] z: vertical levels 
	!>@param[inout] theta: potential temperature 
	!>@param[in] rhoa: density 
	!>@param[in] rhon: density
	!>@param[in] u: vertical wind 
	!>@param[inout] micro_init: boolean to initialise microphysics 
	!>@param[in] hm_flag: switch hm-process on and off
	!>@param[in] mass_ice: mass of a single ice crystal (override)
	!>@param[in] ice_flag: ice microphysics
	!>@param[in] wr_flag: warm rain
	!>@param[in] rm_flag: riming flag
	!>@param[in] theta_flag: whether to alter theta
	!>@param[in] j_stochastic, ice_nuc_flag,mode1_ice_flag,mode2_ice_flag, 
	!             coll_breakup_flag1, heyms_west
	!>@param[in] lawson, recycle
    subroutine p_microphysics_1d(nq,ncat,n_mode,cst,cen, inc, iqc, inr,iqr, ini,iqi,iai,&
                            cat_am,cat_c, cat_r, cat_i,nprec, &
                            kp,o_halo,dt,dz,dzn,q,precip,th,p, z,theta,rhoa,rhon,u, &
    						micro_init,hm_flag, mass_ice,ice_flag, &
    						wr_flag, rm_flag, theta_flag, &
    						j_stochastic,ice_nuc_flag,mode1_ice_flag,mode2_ice_flag, &
    						coll_breakup_flag1, heyms_west, lawson, recycle)
#else
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>solves one time-step of the microphysics
	!>@param[in] nq: number of q-fields
	!>@param[in] ncat: number of categories
	!>@param[in] n_mode: number of aerosol modes
	!>@param[in] cst,cen: indices of categories
	!>@param[in] inc, iqc: index of cloud number, index of cloud mass
	!>@param[in] inr, iqr: index of rain number, index of rain mass
	!>@param[in] ini, iqi,iai: index of ice number, index of ice mass, and ice aerosol
	!>@param[in] cat_am,cat_c, cat_r, cat_i: category index for cloud and rain and ice
	!>@param[in] nprec
	!>@param[in] kp: number of vertical levels
	!>@param[in] dt: time-step
	!>@param[in] dz: dz, dzn
	!>@param[in] o_halo: extra points for advection
	!>@param[inout] q: q-variables 
	!>@param[inout] precip: precip in rain, snow, graupel, ice cats - diagnostic
	!>@param[inout] th: theta perturbation
	!>@param[inout] p: pressure
	!>@param[inout] z: vertical levels 
	!>@param[inout] theta: potential temperature 
	!>@param[in] rho: density 
	!>@param[in] rhon: density
	!>@param[in] u: vertical wind 
	!>@param[inout] vqr,vqc, n_step, adv_l
	!>@param[in] coords
	!>@param[inout] micro_init: boolean to initialise microphysics 
	!>@param[in] hm_flag: switch hm-process on and off
	!>@param[in] mass_ice: mass of a single ice crystal (override)
	!>@param[in] ice_flag: ice microphysics
	!>@param[in] wr_flag: warm rain
	!>@param[in] rm_flag: riming flag
	!>@param[in] theta_flag: whether to alter theta
	!>@param[in] j_stochastic, ice_nuc_flag,mode1_ice_flag,mode2_ice_flag, 
	!           coll_breakup_flag1, heyms_west
	!>@param[in] lawson, recycle
    subroutine p_microphysics_1d(nq,ncat,n_mode,cst,cen, inc, iqc, inr,iqr,ini,iqi,iai,&
                            cat_am,cat_c, cat_r, cat_i,nprec,&
                            kp,o_halo,dt,dz,dzn,q,precip,th,p, z,theta,rhoa,rhon,u, &
                            vqc,vqr,vqi,vni,n_step, adv_l, coords,&
    						micro_init,hm_flag, mass_ice,ice_flag, &
    						wr_flag, rm_flag, theta_flag, &
    						j_stochastic,ice_nuc_flag,mode1_ice_flag,mode2_ice_flag, &
    						coll_breakup_flag1, heyms_west, lawson, recycle)
#endif

	use advection_1d
	use numerics, only : dfsid1
	use advection_s_1d, only : mpdata_vec_1d
    implicit none
    ! arguments:
    integer(i4b), intent(in) :: nq, ncat,n_mode, kp, o_halo, inc, iqc,inr,iqr,&
                             ini,iqi,iai, &
                            cat_am,cat_c, cat_r, cat_i,nprec
    integer(i4b), dimension(ncat), intent(in) :: cst,cen
    real(wp), intent(in) :: dt
    real(wp), dimension(-o_halo+1:kp+o_halo,nq), intent(inout) :: q
    real(wp), dimension(1:kp,nprec), intent(inout) :: precip
    real(wp), dimension(-o_halo+1:kp+o_halo), intent(inout) :: th
    real(wp), dimension(-o_halo+1:kp+o_halo), intent(in) :: dz, z, dzn, rhoa, &
                                                    rhon, theta,p
    real(wp), dimension(-o_halo+1:kp+o_halo), intent(in) :: u
    logical, intent(in) :: hm_flag, ice_flag, wr_flag, rm_flag, theta_flag, heyms_west, &
                            lawson, recycle
    integer(i4b), intent(in) :: ice_nuc_flag, mode1_ice_flag, mode2_ice_flag, &
                             coll_breakup_flag1
    logical , intent(inout) :: micro_init
    real(wp), intent(in) :: mass_ice, j_stochastic
    ! locals:
    integer(i4b) :: k,k1,iter, i
#if MPI_PAMM == 1
    integer(i4b), dimension(3), intent(inout) :: n_step
	logical, intent(inout), dimension(3) :: adv_l
    integer(i4b), dimension(3), intent(in) :: coords
#else
    integer(i4b), dimension(3) :: n_step
	logical, dimension(3) :: adv_l
#endif
    real(wp) :: temp, qtot,qaut, a, b, ab_ice, ab_liq, ice_dep,snow_dep,graup_dep, &
    			nu_ice, nu_snow, nu_graup, diff1, ktherm1, tc, nu_vis, sc, nu_rain, rain_evap, &
    			sb_aut, sb_acr, sb_cwaut, sb_cwacr, sb_raut, sb_rsel, sb_cwsel
    
    real(wp), dimension(-o_halo+1:kp+o_halo) :: rho
    real(wp), dimension(1-o_halo:kp+o_halo) :: smr, smr_i
    
    real(wp), dimension(1-o_halo:kp+o_halo) :: &
    			pgaci,  & ! accretion of cloud ice by graupel
				pgacr, & ! riming of graupel by rain
				pgacs,rgacs, & ! accretion of snow by graupel
				pgacw, & ! riming of graupel by liquid cloud
				pgaut, & ! autoconversion from snow to graupel due to heavy riming
				pgdep, & ! deposition of vapour onto graupel
				pgfr,rgfr,  & ! freezing of rain to form graupel
				pgmlt, & ! melting of graupel to form rain
				pgshd, & ! shedding of liquid from wet graupel to rain
				pgsub, & ! sublimation of graupel or evaporation from wet graupel
				riaci, & ! aggregation of ice crystals
				piacr_g,riacr_g, & ! riming of cloud ice by large rain drops to form graupel
				piacr_s,riacr_s, & ! riming of cloud ice by small rain drops to form snow
				piacw, & ! riming of ice cloud ice by liquid cloud 
				riacw, &
				picnt, & ! nucleation of ice crystals by contact freezing
				pidep, & ! deposition of water vapour onto cloud ice
				piprm, & ! primary nucleation of ice crystals by INPs
				pifrw, & ! nucleation of ice crystals by homogeneous freezing of cloud
				rihal, & ! production of ice crystals by hm process
				pimlt, & ! cloud ice melting to form rain
				rimlt, &
				pisub, & ! sublimation of cloud ice
				risub, &
				praci_g, & ! accretion of cloud ice by rain to form graupel
				praci_s, & ! accretion of cloud ice by rain to form snow
				pracs, & ! accretion of snow by rain to form graupel
				pracw, & ! accretion of liquid cloud by rain
				praut, & ! autoconversion from liquid cloud to rain (coalescence)
				prevp, & ! evaporation of rain
				rrevp, & ! evaporation of rain
				psacr,rsacr, & ! accretion of rain by snow to form graupel
				psaci, & ! accretion of cloud ice by snow
				rsacs, & ! aggregation of snowflakes
				psacw, & ! accretion of liquid cloud by snow
				psaut,rsaut, & ! conversion of large ice crystals to snow
				rsbrk, & ! break-up of large snowflakes
				psdep, & ! deposition of vapour onto snow
				psmlt, & ! melting of snow to form rain
				pssub, & ! sublimation of snow
				rcwaut, & ! reduction in cloud number due to autoconversion
				rcwacr, & ! cloud water number accreted onto rain
				rraut, &  ! increase in rain number due to autoconversion
    			rrsel, &     ! rain self accretion - number
    			rcwsel, &     ! cloud water self accretion - number
    			praci, & ! rain accreting ice - mass
    			rraci, &  ! rain accreting ice - number
    			piacr, & ! ice accreting rain - mass
    			riacr, & ! ice accreting rain - number
    			nfrag_m1c, & ! fragment rate due to mode 1 during collisions
    			nfrag_m2, & ! fragment rate due to mode 2 during collisions
    			nfrag_ii, &  ! fragment rate due to ice-ice collisions
    			nfrag_nucc, &
    			massc_nucc, &
    			nfrag_nucr, &
    			massr_nucr, &
    			frac_r, &
    			frac_x, &
    			vimlt, &
    			mono1, &
    			phi11
    				    
    real(wp) :: pgwet ! amount of liquid that graupel can freeze without shedding
    								

    real(wp), dimension(1-o_halo:kp+o_halo) :: n_r, lam_r, n_i, lam_i, n_s, &
                                                 lam_s, n_g, lam_g, lam_c, n_c, &
                                                 lam_i_star
    real(wp), dimension(1-o_halo:kp+o_halo) :: rho_fac
	real(wp), dimension(1-o_halo:kp+o_halo) :: vnr, vnc
#if MPI_PAMM == 0
	real(wp), dimension(1-o_halo:kp+o_halo) :: vqr, vqs, vqg, vqi, vns, vng, vni, &
	                                        vqc
#else
	real(wp), intent(inout), dimension(1-o_halo:kp+o_halo) :: vqr, vqc, vqi, vni
#endif
	real(wp), dimension(1-o_halo:kp+o_halo) :: t, a_hw,a_hw1,pre_hw
	! coalescence efficiencies
	real(wp), dimension(1-o_halo:kp+o_halo) :: egi_dry, egs_dry, esi, eii, ess
	real(wp) :: qold,des_dt,dqs_dt,err,cond,temp1, dummy1,dummy2, dummy3,&
	            dummy4, &
	            n_mix,s_mix,m_mix, nin_c, din_c,nin_r,din_r, n_tot, s_tot, m_tot
	
	real(wp), dimension(1-o_halo:kp+o_halo) :: gamma_t,dep_density, qold1
	real(wp) :: phi,vol, nfrag=0._wp
	
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! initialise some variables that do not depend on prognostics                        !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (micro_init) then 
    	call initialise_microphysics_vars
    	mi0=mass_ice
    	micro_init=.false.
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! zero arrays                                                                        !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    vqr=0._wp
    vqc=0._wp
	pgfr=0._wp
	rgfr=0._wp
	riaci=0._wp
	piacw=0._wp
	riacw=0._wp ! new
	pidep=0._wp
	pifrw=0._wp
	rihal=0._wp
	pimlt=0._wp
	rimlt=0._wp ! new
	pisub=0._wp
	risub=0._wp ! new
	pracw=0._wp
	praut=0._wp
	prevp=0._wp
	rrevp=0._wp
	rcwaut=0._wp
    rcwacr=0._wp
    rraut=0._wp
    rrsel=0._wp
    rcwsel=0._wp
    praci=0._wp
    rraci=0._wp
    piacr=0._wp
    riacr=0._wp
    nfrag_m1c=0._wp
    nfrag_m2=0._wp
    nfrag_ii=0._wp
    nfrag_nucc=0._wp
    massc_nucc=0._wp
    nfrag_nucr=0._wp
    massr_nucr=0._wp
    frac_r=0._wp ! new
    frac_x=0._wp ! new
    vimlt=0._wp ! new
    mono1=0._wp ! new
    phi11=0._wp ! new
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! some commonly used variables that depend on prognostics                            !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    t=(theta+th)*(p/1.e5_wp)**(ra/cp) ! temperature
    rho=rhon !p / (ra*t) ! air density    
    rho_fac=(rho0/rho(:))**0.5_wp
    ! rain n0, lambda
    lam_r=(max(q(:,cst(cat_r)),1._wp)*cr*gam2r / &
            (max(q(:,cst(cat_r)+1),1.e-10_wp)*gam1r))**(1._wp/dr)
    n_r=rho(:)*max(q(:,cst(cat_r)),1._wp)*lam_r**(1._wp+alpha_r) / gam1r
    ! cloud n0, lambda    
    lam_c=(max(q(:,inc),1._wp)*gam2c / (max(q(:,iqc),1.e-10_wp)*gam1c))**(1._wp/1._wp)
    n_c=rho(:)*max(q(:,inc),0._wp)*lam_c**(1._wp+alpha_c) / gam1c
    
    if(ice_flag) then
                                            
        ! ice n0, lambda
        lam_i=(max(q(:,ini),1._wp)*ci*gam2i / (max(q(:,iqi),1.e-10_wp)*gam1i))**(1._wp/di)
        n_i=rho(:)*max(q(:,ini),0._wp)*lam_i**(1._wp+alpha_i) / gam1i

        if(heyms_west) then
            ! HEYMSFIELD AND WESTBROOK (2010) FALL-SPEEDS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! this calculates the prefactor of the terminal velocity relation
            pre_hw=heymsfield_and_westbrook_fall_parameters2(t,rhoa)
            !
        
    !         a_hw = heymsfield_and_westbrook_fall_parameters1(vol,mass,nmon,&
    !                                         phi,rim,t(1:kp),rhoa)

            ! this calculates a, which is needed to implement the Magic Function.
            ! assumes Ar=1
            a_hw = heymsfield_and_westbrook_fall_parameters1(q(:,ini),q(:,iqi+2),&
                q(:,iqi),1._wp, 1._wp,q(:,iqi+4),t,rho)
            a_hw1 = a_hw
            a_hw = 1._wp/(a_hw**(2._wp/di))
                                        
            ! lambda star - Magic function.
            lam_i_star = lam_i*a_hw
        
            ! calculate fall-speeds
            vni(:)=0.7_wp*max(fall_n_i_hw*pre_hw*lam_i* &
                (1._wp/(gam3ai*lam_i_star**(di/2._wp)+gam3bi*lam_i_star**(di) )),0._wp)
            vqi(:)=0.7_wp*max(fall_q_i_hw*pre_hw*lam_i* &
                (1._wp/(gam4ai*lam_i_star**(di/2._wp)+gam4bi*lam_i_star**(di) )),0._wp)
                
            ! precipitation - this one is actually number flux
            ! the Magic Function - standard, but with n0
            precip(1:kp,2)=0.7_wp*n_i(1:kp)*pre_hw(1:kp)*lam_i(1:kp)**(-alpha_i)/rho(1:kp)* &
             (1._wp/(gam3ai*lam_i_star(1:kp)**(di/2._wp)+gam3bi*lam_i_star(1:kp)**(di) ))
        
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        else  
            ! ice
            vqi(:)=min(max(fall_q_i*rho_fac * lam_i**(1._wp+alpha_i+di) / &
                (lam_i+f_i)**(1._wp+alpha_i+di+b_i), 0._wp), 10._wp)
            vni(:)=max(fall_n_i*rho_fac * lam_i**(1._wp+alpha_i) / &
                (lam_i+f_i)**(1._wp+alpha_i+b_i), 0._wp)

            ! precipitation
            precip(1:kp,2)=n_i(1:kp)*(a_i*chi_num_ice/(lam_i(1:kp)**(alpha_i+b_i+1._wp)) - &
                            u(1:kp)*chi_num_ice1/(lam_i(1:kp)**(alpha_i+1._wp))) &
                            /rho(1:kp)
                        
        endif
        
        
    endif
    
    ! precipitation
	precip(1:kp,1)=cr*n_r(1:kp)*(a_r*chi_rain/(lam_r(1:kp)**(alpha_r+b_r+dr+1._wp)) - &
					u(1:kp)*chi_rain1/(lam_r(1:kp)**(alpha_r+dr+1._wp))) &
					/rho(1:kp) *3600._wp
    
    ! fall speeds
    ! rain
    vqr(:)=min(max(fall_q_r*rho_fac * lam_r**(1._wp+alpha_r+dr) / &
    	(lam_r+f_r)**(1._wp+alpha_r+dr+b_r), 0._wp),10._wp)
    
    vnr(:)=max(fall_n_r*rho_fac * lam_r**(1._wp+alpha_r) / &
    	(lam_r+f_r)**(1._wp+alpha_r+b_r), 0._wp)
    
    ! cloud
    vqc(:)=min(max(fall_q_c*rho_fac * lam_c**(1._wp+alpha_c+1._wp) / &
    	(lam_c+f_c)**(1._wp+alpha_c+1._wp+b_c), 1.e-3_wp), 10._wp)
    
    vnc(:)=max(fall_n_c*rho_fac * lam_c**(1._wp+alpha_c) / &
    	(lam_c+f_c)**(1._wp+alpha_c+b_c), 1.e-3_wp)
    ! coalescence efficiencies
    egi_dry(:)=0.2_wp*exp(0.08*(t(:)-ttr))
    egs_dry(:)=0.2_wp*exp(0.08*(t(:)-ttr))
    esi(:)=0.2_wp*exp(0.08*(t(:)-ttr))
    eii(:)=0.2_wp*exp(0.08*(t(:)-ttr))
    ess(:)=0.2_wp*exp(0.08*(t(:)-ttr))
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
     
    
   
    
    ! loop over all levels
    do k=-o_halo+1,kp+o_halo   
 		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! condensation of liquid water                                                   !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! smr at 0c
		!q0sat=eps1*svp_liq(ttr)/(p(k)-svp_liq(ttr))
		q0sat=eps1*svp_liq(ttr)/(p(k)-svp_liq(ttr))
    	smr(k)=eps1*svp_liq(t(k))/(p(k)-svp_liq(t(k))) ! saturation mixing ratio
        
        des_dt=dfsid1(svp_liq,t(k),1.e0_wp,1.e-8_wp,err)
        dqs_dt=eps1*p(k)*des_dt/(p(k)-svp_liq(t(k)))**2
        qold=q(k,iqc)
        qtot=q(k,1)+q(k,iqc)

		
        q(k,iqc)=q(k,1)+q(k,iqc)-smr(k)
        if (theta_flag) q(k,iqc)=(q(k,iqc)+(lv/cp*qold)*dqs_dt) / (1._wp+lv/cp*dqs_dt)
        q(k,iqc)=max(q(k,iqc),0._wp)
        t(k)=t(k)
        if(theta_flag) t(k)=t(k)+lv/cp*(q(k,iqc)-qold)
		
		tc=t(k)-ttr
    	smr(k)=eps1*svp_liq(t(k))/(p(k)-svp_liq(t(k))) ! saturation mixing ratio
    	q0sat=smr(k)	
    	smr_i(k)=eps1*svp_ice(t(k))/(p(k)-svp_ice(t(k))) ! saturation mixing ratio - ice	
    	
    	cond=(q(k,iqc)-qold)
    	q(k,1)=q(k,1)-cond
    	qold1(k)=qold
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    enddo
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate gamma_t and dep_density for ice growth model                             !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(ice_flag) then
        do k=-o_halo+1,kp+o_halo   
            call chen_and_lamb_anc(t(k),q(k,1),smr_i(k),rhoa(k), &
                                    gamma_t(k), dep_density(k))    
        enddo
    endif 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
    
    ! loop over all levels
    do k=-o_halo+1,kp+o_halo         

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! activation of cloud drops                                                      !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if MPI_PAMM == 1
        if(coords(3)==0) then
            k1=max(k-1,1)
        else
            k1=k-1
        endif
#else
        k1=max(k-1,1)
#endif
	    if((q(k1,iqc) .lt. qsmall) .and. (q(k,iqc) .gt. qsmall).and.(u(k)>0._wp)) then
	    
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Calculate the lognormal parameters                                         !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! this loop calculates the lognormal parameters for the external mixtures
            do i=1,n_mode-1
                call ln_params_from_integral_moms(&
                    q(k,cst(i+1)),q(k,cst(i+1)+1),q(k,cst(i+1)+2), &
                    density_core1(i),sig_aer1(i),d_aer1(i))
                n_aer1(i)=q(k,cst(i+1))
            enddo        
            
            ! calculate ln params and relevant terms for mixed-mode, density, etc
            ! note that we assume that surface area does not change. In reality it 
            ! does depending on aerosol type
            call ln_params_and_props_from_integral_moms( &
                n_mode, &
                q(k,cst(cat_am)), & ! total number
                q(k,cst(cat_am)+2:cen(cat_am)-1:3), & ! surface area
                q(k,cst(cat_am)+3:cen(cat_am):3), & ! mass 
                n_aer1(n_mode),density_core1, &
                molw_core1,nu_core1, & 
                sig_aer1(n_mode),d_aer1(n_mode),n_mix,s_mix,m_mix)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
!                 n_aer1=max(n_aer1,0.1_wp)

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Bulk Aerosol Activation - number of drops
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! calculate aerosol PSD parameters
            p_test=p(k)
            t_test=t(k)
            w_test=max(u(k),0.001_wp)
            call initialise_arrays(n_mode,n_sv,p_test,t_test,w_test, &
                        max(n_aer1,0.1e6),d_aer1,sig_aer1, molw_org1,density_core1)

            call ctmm_activation(n_mode,n_sv,sv_flag, &
                        max(n_aer1,0.1e6), d_aer1,sig_aer1,molw_core1, &
                        density_core1, nu_core1, &
                        org_content1,molw_org1, density_org1, delta_h_vap1, nu_org1, &
                        log_c_star1, &
                        w_test, t_test,p_test, a_eq_7, b_eq_7, &
                        act_frac1,smax1,dcrit2)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            act_frac1=max(act_frac1,0._wp)
            temp1=sum(n_aer1*act_frac1)
            ! put in-cloud aerosol into aerosol - i.e. remove it first
            do i=1,n_mode-1
                q(k,cst(i+1))=  q(k,cst(i+1))+q(k,cst(cat_c)+(i-1)*3+2)
                q(k,cst(i+1)+1)=q(k,cst(i+1)+1)+q(k,cst(cat_c)+(i-1)*3+3)
                q(k,cst(i+1)+2)=q(k,cst(i+1)+2)+q(k,cst(cat_c)+(i-1)*3+4)
                q(k,cst(cat_c)+(i-1)*3+2)=0._wp     ! aerosol num
                q(k,cst(cat_c)+(i-1)*3+3)=0._wp     ! aerosol sa
                q(k,cst(cat_c)+(i-1)*3+4)=0._wp     ! aerosol mass
            enddo
            ! cloud droplet number
            q(k  ,inc)=temp1
            ! now we take activated aerosol away from aerosol field and add to in-cloud
            do i=1,n_mode-1
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! remove from aerosol particles:                                         !
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! number in aerosol modes
                dummy1=ln_part_mom(0,dcrit2(i),q(k,cst(i+1)), sig_aer1(i),d_aer1(i))
                !print *,dummy1, dcrit2(i), act_frac1*q(k,cst(i+1))
                q(k,cst(i+1))=q(k,cst(i+1))-dummy1 
                    
                ! surface area in aerosol modes
                dummy2=pi*ln_part_mom(2,dcrit2(i),q(k,cst(i+1)), sig_aer1(i),d_aer1(i))
                q(k,cst(i+1)+1)=q(k,cst(i+1)+1)- dummy2 
                    
                ! mass in aerosol modes
                dummy3=pi/6._wp*density_core1(i)* &
                    ln_part_mom(3,dcrit2(i),q(k,cst(i+1)), sig_aer1(i),d_aer1(i))
                q(k,cst(i+1)+2)=q(k,cst(i+1)+2)- dummy3
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! add to aerosol particles in cloud water                                !
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! number in aerosol modes
                ! qv, n_mode aerosol + 1
                q(k,cst(cat_c)+(i-1)*3+2)=q(k,cst(cat_c)+(i-1)*3+2)+dummy1 
                    
                ! surface area in aerosol modes
                q(k,cst(cat_c)+(i-1)*3+3)=q(k,cst(cat_c)+(i-1)*3+3)+dummy2 
                    
                ! mass in aerosol modes
                q(k,cst(cat_c)+(i-1)*3+4)=q(k,cst(cat_c)+(i-1)*3+4)+dummy3
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            enddo
            ! deplete aerosol from mixed mode
            ! this calculates the total depletion. For each component
            ! deplete, base on fraction of each component
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! remove from aerosol particles:                                             !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! number in aerosol modes
            dummy1=ln_part_mom(0,dcrit2(n_mode),n_mix, &
                            sig_aer1(n_mode),d_aer1(n_mode))
            ! surface area in aerosol modes
            dummy2=pi*ln_part_mom(2,dcrit2(n_mode),n_mix, &
                            sig_aer1(n_mode),d_aer1(n_mode))
            ! mass in aerosol modes
            dummy3=pi/6._wp*density_core1(n_mode)* &
                ln_part_mom(3,dcrit2(n_mode),n_mix, &
                            sig_aer1(n_mode),d_aer1(n_mode))
                            
            q(k,cst(cat_am))=q(k,cst(cat_am))*(1._wp-min(dummy1/n_mix,1._wp))
            do i=1,n_mode-1 ! deplete aerosol
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! add to aerosol particles in cloud water                                !
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! number in aerosol modes
                ! qv, n_mode aerosol + 1
                q(k,cst(cat_c)+(i-1)*3+2)=q(k,cst(cat_c)+(i-1)*3+2)+ &
                    max(dummy1/(n_mix)*q(k,cst(cat_am)+3*(i-1)+1),0._wp)
                
                ! surface area in aerosol modes
                q(k,cst(cat_c)+(i-1)*3+3)=q(k,cst(cat_c)+(i-1)*3+3)+ &
                    max(dummy2/(s_mix)*q(k,cst(cat_am)+3*(i-1)+2),0._wp)
                
                ! mass in aerosol modes
                q(k,cst(cat_c)+(i-1)*3+4)=q(k,cst(cat_c)+(i-1)*3+4)+ &
                    max(dummy3/(m_mix)*q(k,cst(cat_am)+3*(i-1)+3),0._wp)
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


                ! number - remove aerosol particles
                q(k,cst(cat_am)+3*(i-1)+1)=q(k,cst(cat_am)+3*(i-1)+1) * &
                        (1._wp-max(dummy1/(n_mix),0._wp))
                
                ! surface area
                q(k,cst(cat_am)+3*(i-1)+2)=q(k,cst(cat_am)+3*(i-1)+2)* &
                        (1._wp-max(dummy2/(s_mix),0._wp) )
                ! mass
                q(k,cst(cat_am)+3*(i-1)+3)=q(k,cst(cat_am)+3*(i-1)+3)* &
                        (1._wp-max(dummy3/(m_mix),0._wp) )


            enddo 
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
        endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! end the activation of cloud drops                                              !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        
        
        
        
        
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Calculate process rates                                                        !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		nin_c=0._wp
		nin_r=0._wp
		if(ice_flag.and.(t(k).lt.ttr)) then
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! 0. define properties                                                       !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            frac_r(k) = q(k,ini+5) / (q(k,iqi)+qsmall) ! fraction of rime to total
            frac_x(k) = (q(k,iqi)-q(k,ini+5)) / (q(k,iqi)+qsmall) ! frac of xtal to total
            phi11(k) = q(k,iqi+1) / (q(k,ini)+qsmall) ! phi
            mono1(k) = q(k,iqi+3) / (q(k,ini)+qsmall) ! number of monomers
		
		
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! 1. collisions between precipitating particles of different species         !
            ! praci, rraci, piacr, riacr, rates                                          !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call collisions_between_precipitating_particles(praci(k),rraci(k),piacr(k), &
                riacr(k), n_r(k), n_i(k), rho(k), vqr(k), vqi(k), vnr(k), vni(k), &
                lam_r(k), lam_i(k),q(k,iqi), q(k,iqi+2), q(k,iqi+4), q(k,iqr), &
                a_hw1(k), pre_hw(k), heyms_west)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
                
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! 2. mode 1 SIP during collisions                                            !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if((mode1_ice_flag.eq.1).and.(t(k).lt.268._wp)) then
                ! n_frag_m1c (collisional) is a delta
                call mode1_sip_collisional(lam_r(k), lam_i(k), t(k), n_r(k), n_i(k), &
                    q(k,iqi), q(k,iqi+2), q(k,iqi+4), q(k,iqr), &
                    a_hw1(k), pre_hw(k), dt, heyms_west, nfrag_m1c(k))
            endif
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
                            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! 3. mode 2 SIP                                                              !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if((mode2_ice_flag.eq.1).and.(t(k).lt.268._wp)) then
                ! n_frag_m2 is a delta
                call mode2_sip_collisional(lam_r(k), lam_i(k), t(k), n_r(k), n_i(k), &
                    q(k,iqi), q(k,iqi+2), q(k,iqi+4), q(k,iqr), &
                    a_hw1(k), pre_hw(k), dt, heyms_west, nfrag_m2(k))
            endif
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! 4. Ice-ice collisions SIP                                                  !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! n_frag_ii is a delta
            if (t(k).lt.268._wp) &
            call ice2_sip_collisional(lam_i(k), t(k), n_i(k), &
                q(k,ini), q(k,iqi), q(k,iqi+2), q(k,iqi+4), q(k,iqi+1), &
                a_hw1(k), pre_hw(k), dt, heyms_west, nfrag_ii(k),coll_breakup_flag1)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! 5. ice nucleation with drop-freezing SIP                                   !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! cloud water
            if(q(k,inc) > 0._wp) then 
                ! nin_c, nfrag_nucc, massc_nucc are deltas, not rates
                call ice_nucleation_and_mode1(n_mode, q(k,inc), &
                    q(k,ini), q(k,iqc), t(k), q(k,cst(cat_c)+3:cen(cat_c)-1:3), & ! surface area
                    q(k,cst(cat_c)+4:cen(cat_c):3), & ! Mass aerosol 
                    density_core1, molw_core1, nu_core1, nin_c, din_c, &
                    n_mix, s_mix,m_mix, &
                    n_aer1(n_mode),sig_aer1(n_mode),d_aer1(n_mode), &
                    j_stochastic, dt, &
                    ice_nuc_flag, nfrag_nucc(k), massc_nucc(k), lawson,.true.,&
                    mode1_ice_flag)
            endif
            ! rain water
            if(q(k,inr) > 0._wp) then
                ! nin_r, nfrag_nucr, massc_nucr are deltas, not rates
                call ice_nucleation_and_mode1(n_mode, q(k,inr), &
                    q(k,ini), q(k,iqr), t(k), q(k,cst(cat_r)+3:cen(cat_r)-1:3), & ! surface area
                    q(k,cst(cat_r)+4:cen(cat_r):3), & ! Mass aerosol 
                    density_core1, molw_core1, nu_core1, nin_r, din_r, &
                    n_mix, s_mix,m_mix, &
                    n_aer1(n_mode),sig_aer1(n_mode),d_aer1(n_mode), &
                    j_stochastic, dt, &
                    ice_nuc_flag, nfrag_nucr(k), massr_nucr(k), lawson,.false.,&
                    mode1_ice_flag)
            endif
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! 6. collection of cloud by ice - riming                                     !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		    if(((q(k,iqc).gt.qsmall) .and. (q(k,iqi).gt.qsmall)) .and. rm_flag) then
		        ! piacw is a rate
		        if(heyms_west) then
		        
                    lambda0i=lam_i(k)
                    ci_new=pi/6._wp*min(910._wp, &
                        q(k,iqi)/(q(k,iqi+2)+q(k,iqi+4)/500._wp))
                    mrthresh=ci_new*1.e-6_wp**di
                    miupper=ci_new*(pthreshi/lambda0i)**di

                    n0i=n_i(k)
                    a_hw_new=a_hw1(k)
                    pre_hw_new=pre_hw(k)
                    ! only call integral if mrupper gt mrthresh
                    if(miupper.gt.mrthresh) then
                        ! riming
                        piacw(k)=romb(integral_rime_hw,mrthresh,miupper)*q(k,iqc)
                        miupper=((q(k,iqc)*6._wp)/ &
                            ((q(k,inc)+qsmall)*pi*rhow))**oneoverthree
                        piacw(k)=piacw(k)*min(max(0.0_wp, &
                            1._wp/(50.e-6_wp-10.e-6_wp)*miupper-0.25_wp ),1.0_wp)
                    endif
		        
		        
!                     piacw(k)=max(0.25_wp*pi*n_i(k)*pre_hw(k)*&
!                         lam_i(k)**(-2._wp-alpha_i) * eiw * 0.2 * q(k,iqc) &
!                       /(gam6ai*lam_i_star(k)**(di*0.5_wp)+gam6bi*lam_i_star(k)**(di)),&
!                         0._wp)		        
		        else
                    piacw(k)=max(mass_iacw * n_i(k)* eiw *q(k,iqc)*rho_fac(k) / &
                            (lam_i(k)+f_i)**(3._wp+b_i+alpha_i),0._wp)
                endif                                
                piacw(k)=max(min(piacw(k),q(k,iqc)/dt),0._wp)
                piacw(k)=min(piacw(k),max((ttr-t(k))*cp/lf/dt,0._wp))
            endif
		    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! 7. h-m process                                                             !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! h-m process                                                                !
            if(hm_flag.and. &
                ((t(k).le.(ttr-2.0_wp)).and.(t(k).ge.(ttr-9.0_wp)))) then
                rihal(k)=max(hm_rate*piacw(k)*hm_func(t(k)),0._wp)
            endif
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! 8. deposition and sublimation                                              !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            phi=1._wp
		    if(q(k,iqi).gt.qsmall) then
		        if(heyms_west) then
                    nu_ice=2._wp*pi*n_i(k) / rho(k) * &
                            (nu_i1 / lam_i(k)**(2._wp+alpha_i) + &
                            sqrt(pre_hw(k))*0.31_wp*sc**(1._wp/3._wp)* &
                            (rho(k)/nu_vis)**0.5_wp* &
                            
                            lam_i(k)**(-2._wp-alpha_i)/ &
                            (lam_i_star(k)**(0.25_wp*di)*gam5ai+&
                            lam_i_star(k)**(0.5_wp*di)*gam5bi) )
                else
                    nu_ice=2._wp*pi*n_i(k) / rho(k) * &
                            (nu_i1 / lam_i(k)**(2._wp+alpha_i) + &
                            (a_i/nu_vis)**0.5_wp*sc**(1._wp/3._wp)* &
                            (rho(k)*rho0)**0.25_wp*nu_i2 / &
                            (lam_i(k)+0.5_wp*f_i)**(0.5_wp*b_i+alpha_i+2.5_wp))                
                endif
                ab_ice=ls**2 / (ktherm1*rv*t(k)**2) + 1._wp/(rho(k)*smr_i(k)*diff1)
                ! chen and lamb growth rates
                phi=min(max(q(k,iqi+1) / (q(k,ini+4)+qsmall),1.e-5_wp),100._wp)
                nu_ice=nu_ice*chen_and_lamb_cap_fac(phi)
        
                ! non chen and lamb bit        
                ice_dep=(q(k,1)/smr_i(k)-1._wp) / (rho(k)*ab_ice)*nu_ice
        
                if(q(k,1).gt.smr_i(k)) then
                    pisub(k)=0._wp
                    pidep(k)=min(max(ice_dep,0._wp),(q(k,1)-smr_i(k))/dt)
                else
                    pidep(k)=0._wp
                    pisub(k)=min(-min(ice_dep,0._wp),-(q(k,1)-smr_i(k))/dt)
                endif
            endif
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! 9. ice aggregation see Ferrier (1994)                                      !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if(lam_i(k).lt.1.e5_wp) then
                ! riaci is a rate
                if(heyms_west) then
                    ! collisions
!                     dummy1=max(a_hw1(k)*a_hw1(k)*pre_hw(k)*iice2*n_i(k)*n_i(k) / &
!                             lam_i(k)**(3._wp+2.*wp*alpha_i+di),0._wp)
                    call collisions_between_ice_particles(dummy1,n_i(k), rho(k), &
                         lam_i(k), q(k,iqi),q(k,iqi+2),q(k,iqi+4),a_hw1(k),pre_hw(k))
                            
                else
                    ! collisions
                    dummy1=max(iice*n_i(k)**2._wp*rho_fac(k) / &
                            lam_i(k)**(4._wp+2.*wp*alpha_i+b_i),0._wp)
                endif                                    
                ! aggregation rate
                riaci(k)=eii(k)*dummy1
            endif
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
        endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! 10. melting of ice                                                             !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if((t(k).gt.ttr).and.ice_flag) then
            ! pimlt is a rate
            pimlt(k)=q(k,iqi)/dt+ &
                (pidep(k)-pisub(k)+piacw(k)) ! ice melts instantaneously
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! 11. evaporation of rain                                                        !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! prevp is a rate
		diff1=dd(t(k),p(k))
		ktherm1=ka(t(k))
		ab_liq=lv**2 / (ktherm1*rv*t(k)**2) + 1._wp/(rho(k)*smr(k)*diff1)
		nu_vis=viscosity_air(t(k)) / rho(k)
		sc=nu_vis / diff1
		nu_rain=2._wp*pi*n_r(k) / rho(k) * &
				(nu_r1 / lam_r(k)**(2._wp+alpha_r) + &
				(a_r/nu_vis)**0.5_wp*sc**(1._wp/3._wp)* &
				(rho(k)*rho0)**0.25_wp*nu_r2 / &
				(lam_r(k)+0.5_wp*f_r)**(0.5_wp*b_r+alpha_r+2.5_wp))

		rain_evap=(q(k,1)/smr(k)-1._wp) / (rho(k)*ab_liq)*nu_rain
		if(q(k,1).gt.smr(k)) then
			prevp(k)=0._wp
		else
			prevp(k)=-min(rain_evap,0._wp)
		endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! 12. warm rain autoconversion based on Seifert and Beheng (2006)                !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if (wr_flag) then
            call seifert_beheng(sb_aut,sb_acr, sb_cwaut, sb_cwacr, sb_raut, &
                        sb_rsel, sb_cwsel, q(k,cst(cat_c)+1),q(k,cst(cat_c)),&
                        q(k,cst(cat_r)+1),q(k,cst(cat_r)),rho(k),dt)
            ! all rates
            praut(k)=sb_aut
            pracw(k)=sb_acr
            rcwaut(k)=sb_cwaut
            rcwacr(k)=sb_cwacr
            rraut(k)=sb_raut
            rrsel(k)=sb_rsel
            rcwsel(k)=sb_cwsel
        endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        if(ice_flag) then
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Scale process rates so that cannot get negative values                     !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call scale_microphysics(praci(k),rraci(k), piacr(k), riacr(k), &
                nin_c, nin_r, massc_nucc(k), massr_nucr(k), &
                piacw(k), riacw(k), rihal(k), pidep(k), pisub(k), risub(k), riaci(k), &
                pimlt(k), rimlt(k), prevp(k), rrevp(k), &
                praut(k), pracw(k), rcwacr(k), rraut(k), rrsel(k), &
                rcwaut(k),rcwsel(k), &
                q(k,1),q(k,inc),q(k,iqc),q(k,inr), q(k,iqr),q(k,ini),q(k,iqi),t(k),dt)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        else
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Scale process rates so that cannot get negative values                     !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call scale_microphysics_warm(rrevp(k), prevp(k), praut(k), &
                pracw(k), rcwacr(k), rraut(k), rrsel(k), rcwaut(k),rcwsel(k), &
                q(k,inc),q(k,iqc),q(k,inr), q(k,iqr),dt)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
        endif

    
    
    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! ice nucleation block                                                           !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(ice_flag.and.(t(k)<ttr)) then
        
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! ice nucleation from cloud water                                            !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if((q(k,inc) > 0._wp).and.(nin_c.gt.0._wp)) then            
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! remove aerosol from cloud water and add to ice                         !
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                
                call ln_params_and_props_from_integral_moms(n_mode, q(k,inc), &
                    q(k,cst(cat_c)+3:cst(cat_c)+(n_mode-2)*3+3:3), &
                    q(k,cst(cat_c)+4:cst(cat_c)+(n_mode-2)*3+4:3), & ! mass 
                    n_aer1(n_mode),density_core1, molw_core1,nu_core1, &
                    sig_aer1(n_mode),d_aer1(n_mode),n_mix,s_mix,m_mix)
                    
                call move_aerosol_larger_than_size(n_mode, &
                    din_c,n_mix, sig_aer1(n_mode), &
                    d_aer1(n_mode), density_core1(n_mode), &
                    q(k,cst(cat_c)+2:cst(cat_c)+(n_mode-2)*3+2:3), & ! number in cw mode
                    q(k,cst(cat_c)+3:cst(cat_c)+(n_mode-2)*3+3:3), & ! sa in cw mode
                    q(k,cst(cat_c)+4:cst(cat_c)+(n_mode-2)*3+4:3), & ! mass in cw mode
                    q(k,iai:iai+(n_mode-2)*3:3), & ! number in ice mode
                    q(k,iai+1:iai+(n_mode-2)*3+1:3), & !sa in ice mode
                    q(k,iai+2:iai+(n_mode-2)*3+2:3))  !mass in ice mode

                pifrw(k)=pifrw(k)+massc_nucc(k)/dt
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            endif 
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! end ice nucleation from cloud water                                        !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! ice nucleation from rain water                                             !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if((q(k,cst(cat_r)) > 0._wp).and.(nin_r.gt.0._wp)) then
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! remove aerosol from rain water and add to ice                          !
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                call ln_params_and_props_from_integral_moms(n_mode, q(k,inr), &
                    q(k,cst(cat_r)+3:cst(cat_r)+(n_mode-2)*3+3:3), &
                    q(k,cst(cat_r)+4:cst(cat_r)+(n_mode-2)*3+4:3), & ! mass 
                    n_aer1(n_mode),density_core1, molw_core1,nu_core1, &
                    sig_aer1(n_mode),d_aer1(n_mode),n_mix,s_mix,m_mix)
                    
                call move_aerosol_larger_than_size(n_mode, &
                    din_r,n_mix, sig_aer1(n_mode), &
                    d_aer1(n_mode), density_core1(n_mode), &
                    q(k,cst(cat_r)+2:cst(cat_r)+(n_mode-2)*3+2:3), & ! number in rw mode
                    q(k,cst(cat_r)+3:cst(cat_r)+(n_mode-2)*3+3:3), & ! sa in rw mode
                    q(k,cst(cat_r)+4:cst(cat_r)+(n_mode-2)*3+4:3), & ! mass in rw mode
                    q(k,iai:iai+(n_mode-2)*3:3), & ! number in ice mode
                    q(k,iai+1:iai+(n_mode-2)*3+1:3), & !sa in ice mode
                    q(k,iai+2:iai+(n_mode-2)*3+2:3))  !mass in ice mode

                pgfr(k)=pgfr(k)+massr_nucr(k)/dt
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                    
            endif 
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! end ice nucleation from rain water                                         !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! end ice nucleation block                                                       !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! put aerosol from freezing rain in ice, 
        ! should be praci (assuming piacr goes to ice anyway, so doesn't need doing)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if(ice_flag.and.(t(k).lt.ttr).and.(((q(k,  iqr)-pgfr(k)*dt).gt.qsmall))) then
		    ! move aerosol from rain to ice due to praci
			call move_aerosol_proportional( n_mode, &
			    q(k,cst(cat_r)+2:cst(cat_r)+2+(n_mode-2)*3:3), &
			    q(k,cst(cat_r)+3:cst(cat_r)+3+(n_mode-2)*3:3), &
			    q(k,cst(cat_r)+4:cst(cat_r)+4+(n_mode-2)*3:3), &
			    q(k,iai:iai+(n_mode-2)*3:3), &
			    q(k,iai+1:iai+1+(n_mode-2)*3:3), &
			    q(k,iai+2:iai+2+(n_mode-2)*3:3), &
			    praci(k)*dt,q(k,  iqr)-pgfr(k)*dt ,.true.)
		endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! move cloud aerosol into rain aerosol
        ! aerosol going into rain, by coll-coal
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if((q(k,inc)-nin_c).gt.qsmall) then
            call move_aerosol_proportional( n_mode, &
                q(k,cst(cat_c)+2:cst(cat_c)+2+(n_mode-2)*3:3), &
                q(k,cst(cat_c)+3:cst(cat_c)+3+(n_mode-2)*3:3), &
                q(k,cst(cat_c)+4:cst(cat_c)+4+(n_mode-2)*3:3), &
                q(k,cst(cat_r)+2:cst(cat_r)+2+(n_mode-2)*3:3), &
                q(k,cst(cat_r)+3:cst(cat_r)+3+(n_mode-2)*3:3), &
                q(k,cst(cat_r)+4:cst(cat_r)+4+(n_mode-2)*3:3), &
                -(rcwaut(k)+rcwacr(k)+rcwsel(k))*dt,q(k,  inc)-nin_c, .true. )
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        ! evaporation of rain - similar to melting of ice
        if(prevp(k) .gt. 0._wp) then
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! calculate the number conc. of rain drops evaporated
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            dummy2=rrevp(k)*dt
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            if(dummy2 .gt. qsmall) then
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! add evaporated rain particles to mixed-mode aerosol
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if(recycle) &
                    q(k,cst(cat_am))=q(k,cst(cat_am))+dummy2 ! total number of the mixed-mode


                call move_aerosol_proportional( n_mode, &
                    q(k,cst(cat_r)+2:cst(cat_r)+2+(n_mode-2)*3:3), &
                    q(k,cst(cat_r)+3:cst(cat_r)+3+(n_mode-2)*3:3), &
                    q(k,cst(cat_r)+4:cst(cat_r)+4+(n_mode-2)*3:3), &
                    q(k,cst(cat_am)+1:cst(cat_am)+1+(n_mode-2)*3:3), &
                    q(k,cst(cat_am)+2:cst(cat_am)+2+(n_mode-2)*3:3), &
                    q(k,cst(cat_am)+3:cst(cat_am)+3+(n_mode-2)*3:3), &
                    prevp(k)*dt,q(k,  iqr)-pgfr(k)*dt-praci(k)*dt , recycle )
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            endif
        endif



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Update variables                                                               !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! vapour mass
        q(k,1)=q(k,1)+(pisub(k)-pidep(k))*dt
        
        ! liquid mass - riming not done here
        dummy3=q(k,iqc)
        q(k,iqc)=q(k,iqc)-(praut(k)+pracw(k)+pifrw(k))*dt
       
        ! rain number - rraut, rrsel are negative
        q(k,cst(cat_r))=q(k,cst(cat_r))+(rraut(k)+rrsel(k)-pgfr(k)* &
            q(k,cst(cat_r))/(q(k,cst(cat_r)+1)+qsmall))*dt

        ! rain mass
        q(k,cst(cat_r)+1)=q(k,cst(cat_r)+1)+(praut(k)+pracw(k)-(pgfr(k)))*dt
        ! treat rain evaporation separately - adjust
        prevp(k)=min(prevp(k),q(k,cst(cat_r)+1)/dt) 

        ! liquid number
        dummy4=q(k,inc)
        q(k,inc)=q(k,inc)- &
            min(max(-(rcwaut(k)+rcwacr(k)+rcwsel(k))*dt/(q(k,inc)+qsmall),0._wp),1._wp)* &
                                q(k,inc)


        if(ice_flag) then
            ! rime mass divided by mass
            dummy1=q(k,iqi+4)/(q(k,iqi)+qsmall)

            ! ****ALTERING Q-VARIABLES / PROPERTIES****
            ! mode-1 SIP         
            ! increase ice crystal number
            q(k  ,ini)=q(k  ,ini)+nfrag_m1c(k)
            ! increase ice crystal shape factor
            q(k  ,iqi+1)=q(k  ,iqi+1)+nfrag_m1c(k)
            ! increase ice crystal monomers
            q(k  ,iqi+3)=q(k  ,iqi+3)+nfrag_m1c(k)

            ! mode-2 SIP        
            ! increase ice crystal number
            q(k  ,ini)=q(k  ,ini)+nfrag_m2(k)
            ! increase ice crystal shape factor
            q(k  ,iqi+1)=q(k  ,iqi+1)+nfrag_m2(k)
            ! increase ice crystal monomers
            q(k  ,iqi+3)=q(k  ,iqi+3)+nfrag_m2(k)
        
            ! ice-ice collisions
            ! increase ice crystal number
            q(k  ,ini)=q(k  ,ini)+nfrag_ii(k)
            ! increase ice crystal shape factor
            q(k  ,iqi+1)=q(k  ,iqi+1)+nfrag_ii(k)
            ! increase ice crystal monomers
            q(k  ,iqi+3)=q(k  ,iqi+3)+nfrag_ii(k)
            
            ! ****ALTERING Q-VARIABLES / PROPERTIES****
            ! increase ice crystal number
            q(k  ,ini)=q(k  ,ini)+nin_c+nfrag_nucc(k)
            ! increase ice crystal mass - added divided by number of cloud, 
                            ! multiplied by mass of cloud
            q(k  ,iqi)=q(k  ,iqi)+massc_nucc(k)
            ! increase ice crystal shape factor
            q(k  ,iqi+1)=q(k  ,iqi+1)+nin_c+nfrag_nucc(k)
            ! increase ice crystal volume factor
            q(k  ,iqi+2)=q(k  ,iqi+2)+massc_nucc(k)/rhoi
            ! increase ice crystal monomers
            q(k  ,iqi+3)=q(k  ,iqi+3)+nin_c+nfrag_nucc(k)                
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            ! deplete cloudnc
            ! ****ALTERING Q-VARIABLES / PROPERTIES****
            q(k,  inc)=q(k,inc)-nin_c
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

            ! ****ALTERING Q-VARIABLES / PROPERTIES****
            ! increase in ice crystal number
            q(k  ,ini)=q(k  ,ini)+nin_r+nfrag_nucr(k)
            ! increase in ice crystal mass
            q(k  ,iqi)=q(k  ,iqi)+massr_nucr(k)
            
            ! increase ice crystal shape factor
            q(k  ,iqi+1)=q(k  ,iqi+1)+nin_r+nfrag_nucr(k)
            ! increase ice crystal volume factor
            q(k  ,iqi+2)=q(k  ,iqi+2)+massr_nucr(k)/rhoi
            ! increase ice crystal monomers
            q(k  ,iqi+3)=q(k  ,iqi+3)+nin_r+nfrag_nucr(k)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! deplete rain
            ! ****ALTERING Q-VARIABLES / PROPERTIES****
            q(k,  inr)=q(k,inr)-nin_r
            
            ! ****ALTERING Q-VARIABLES / PROPERTIES****
			! iacr is a source of ice mass only - not number
            ! increase ice crystal mass
            q(k  ,iqi)  =q(k  ,iqi)+praci(k)*dt
            ! increase rime mass of ice
            q(k,  iqi+4)=q(k  ,iqi+4)+praci(k)*dt
            ! iacr is a sink of rain mass and number
            q(k,  iqr)  =q(k, iqr)-praci(k)*dt
            q(k,  inr)  =q(k, inr)-rraci(k)*dt


            ! rime mass - could do this in proportion i.e. rm/q*dm
            ! NOTE, iqi has been changed prior
            q(k,iqi+4)=q(k,iqi+4)+dummy1*(pidep(k)-pisub(k))*dt
            ! ice mass - total
            q(k,iqi)=q(k,iqi)+(pidep(k)-pisub(k))*dt

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! do the riming here                                                         !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if(rm_flag.and.(q(k,iqc).ge.qsmall)) then
                !dummy3 is initial cloud water mass, dummy4 is initial nc

                ! riming
                q(k,iqi)=q(k,iqi)+piacw(k)*dt
                q(k,iqi+4)=q(k,iqi+4)+piacw(k)*dt
                
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! add add aerosol in cw to ice during riming
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                call move_aerosol_proportional( n_mode, &
                    q(k,cst(cat_c)+2:cst(cat_c)+2+(n_mode-2)*3:3), &
                    q(k,cst(cat_c)+3:cst(cat_c)+3+(n_mode-2)*3:3), &
                    q(k,cst(cat_c)+4:cst(cat_c)+4+(n_mode-2)*3:3), &
                    q(k,iai:iai+(n_mode-2)*3:3), &
                    q(k,iai+1:iai+1+(n_mode-2)*3:3), &
                    q(k,iai+2:iai+2+(n_mode-2)*3:3), &
                    piacw(k)*dt,q(k,iqc) ,.true.)
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                
                q(k,inc)=q(k,inc)-riacw(k)*dt
                q(k,iqc)=q(k,iqc)-piacw(k)*dt
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! h-m process                                                            !
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if((t(k).le.(ttr-2.0_wp)).or.(t(k).ge.(ttr-9.0_wp))) then

                    if(hm_flag) then
                        rihal(k)=max(hm_rate*piacw(k)*hm_func(t(k)),0._wp)
                        dummy3=rihal(k)*dt
                        ! increase ice number
                        q(k,ini) = q(k,ini) + dummy3
                        ! increase phi
                        q(k,iqi+1) = q(k,iqi+1) + dummy3
                        ! increase monomers
                        q(k,iqi+3) = q(k,iqi+3) + dummy3
                    endif
                endif
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            endif
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



            
            ! add the aerosol in ice into the mixed-mode aerosol
            ! fraction that number reduces by
            q(k,ini)=q(k,ini)-(riaci(k))*dt
            if(q(k,iqi)<qsmall) then
                dummy2=q(k,ini)
                q(k,1)=q(k,1)+q(k,iqi)
                q(k,ini:ini+5)=0._wp ! all properties, except aerosol
                
            
                if(recycle) &
                    q(k,cst(cat_am))=q(k,cst(cat_am))+dummy2 ! total number of the mixed-mode
                    
                    
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! move all aerosol in ice to the mixed-mode
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                call move_aerosol_proportional( n_mode, &
                    q(k,iai:iai+(n_mode-2)*3:3), & ! number in ice mode
                    q(k,iai+1:iai+(n_mode-2)*3+1:3), & !sa in ice mode
                    q(k,iai+2:iai+(n_mode-2)*3+2:3), & ! mass in ice mode
                    q(k,cst(cat_am)+2:cst(cat_am)+2+(n_mode-2)*3:3), &
                    q(k,cst(cat_am)+3:cst(cat_am)+3+(n_mode-2)*3:3), &
                    q(k,cst(cat_am)+4:cst(cat_am)+4+(n_mode-2)*3:3), &
                    1._wp,1._wp ,recycle)
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    
            endif
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! deposition & sublimation onto ice                                          !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
		    if((q(k,iqi).gt.qsmall).and.(q(k,iqi+2).gt.0._wp)) then
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! chen and lamb                                                          !
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! current volume
                vol=q(k,iqi+2)
                ! dummy1 is the rime mass fraction
                ! dm of the crystal
                call chen_and_lamb_prop((1._wp-dummy1)*(pidep(k)-pisub(k))*dt,gamma_t(k), &
                    vol,phi, dep_density(k))
                ! this is the new volume of the crystals
                ! NB, IQI has changed prior
                vol=min(max(vol,(q(k,iqi)-q(k,iqi+4))/rhoi),(q(k,iqi)-q(k,iqi+4))/10._wp)
                q(k,iqi+2)=vol
                q(k,iqi+1)=phi*q(k,ini+4)
            endif
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! end deposition & sublimation onto ice                                      !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
            
        endif


        if(q(k,iqc) .lt. qsmall) then ! if evaporated
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! move all aerosol in cloud to the aerosol
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call move_aerosol_proportional( n_mode, &
                q(k,cst(cat_c)+2:cst(cat_c)+2+(n_mode-2)*3:3), &
                q(k,cst(cat_c)+3:cst(cat_c)+3+(n_mode-2)*3:3), &
                q(k,cst(cat_c)+4:cst(cat_c)+4+(n_mode-2)*3:3), &
                q(k,cst(2):cst(n_mode):3), &
                q(k,cst(2)+1:cst(n_mode)+1:3), &
                q(k,cst(2)+2:cst(n_mode)+2:3), &
                1._wp,1._wp ,recycle)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            q(k,1)=q(k,1)+q(k,iqc)
            q(k,inc:iqc) = 0.0_wp
        endif
        


        

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! adjust temperature 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        t(k)=t(k)-lv/cp*prevp(k)*dt+lf/cp*pifrw(k)*dt+lf/cp*pgfr(k)*dt+ &
            lf/cp*praci(k)*dt + &
            ls/cp*(pidep(k)-pisub(k))*dt + & !-lf/cp*(pimlt(k))*dt + &
            lf/cp*(piacw(k))*dt
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! melting ice!
        if(ice_flag) then
            if((t(k)+lf/cp*q(k,iqi) ).gt.ttr) then
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! calculate the number conc. of ice melted
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                pimlt(k)=q(k,iqi)/dt
                dummy2=rimlt(k)*(q(k,ini)/(qsmall+q(k,iqi)))*dt
                dummy2=min(dummy2,q(k,ini))
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                if(dummy2 .gt. qsmall) then
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! move aerosol in melting ice to the aerosol in rain
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    call move_aerosol_proportional( n_mode, &
                        q(k,iai:iai+(n_mode-2)*3:3), & ! number in ice mode
                        q(k,iai+1:iai+(n_mode-2)*3+1:3), & !sa in ice mode
                        q(k,iai+2:iai+(n_mode-2)*3+2:3), & ! mass in ice mode
                        q(k,cst(cat_r)+2:cst(cat_r)+2+(n_mode-2)*3:3), &
                        q(k,cst(cat_r)+3:cst(cat_r)+3+(n_mode-2)*3:3), &
                        q(k,cst(cat_r)+4:cst(cat_r)+4+(n_mode-2)*3:3), &
                        dummy2,q(k,ini) ,recycle)
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! ice properties
                    q(k,cst(cat_i):cst(cat_i)+5) = q(k,cst(cat_i):cst(cat_i)+5) * &
                        (1._wp - min(dummy2/(q(k,ini)+qsmall),1._wp ))
        
                    ! add the number of ice and mass to the rain
                    q(k,inr)=q(k,inr)+dummy2

                    ! mass already added
                    q(k,iqr)=q(k,iqr)+pimlt(k)*dt
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                endif
                t(k)=t(k)-lf/cp*pimlt(k)*dt   
            endif    
        endif

        q(k,cst(cat_r)+1)=q(k,cst(cat_r)+1)-prevp(k)*dt
    
        q(k,cst(cat_r))=q(k,cst(cat_r))-rrevp(k)*dt
    
        q(k,1)=q(k,1)+prevp(k)*dt
        
!         if(any(q(k,:)< 0._wp)) then 
!             print *, t(k),k,q(k,:)
!             stop
!         endif


             



        q(k,:)=max(q(k,:),0._wp)	

     
        if (theta_flag) th(k)=t(k)*(1.e5_wp/p(k))**(ra/cp)-theta(k)

    enddo
 	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 





    

    rho=1._wp ! fudge for advection conservation
 	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! advection rain 0th order Bott, a.k.a. upstream advection                           !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! cloud 
    if(sum(q(1:kp,cst(cat_c)+1)).gt.qsmall) then
        adv_l(1)=.true.
		where(isnan(vqc))
			vqc=0._wp
		end where
		vqc(kp+1:kp+o_halo)=vqc(kp+o_halo)
		n_step(1)=max(ceiling(maxval(vqc(-o_halo+1:kp+o_halo)*dt/dz*2._wp)),1)
		vqc(1-o_halo:kp+o_halo-1)=-vqc(-o_halo+2:kp+o_halo)
#if MPI_PAMM == 0
		do iter=1,n_step(1)
            call mpdata_vec_1d(dt/real(n_step(1),wp),dz,dzn,&
                            rho,rhon,kp,cen(cat_c)-cst(cat_c)+1,o_halo,o_halo,&
                            vqc(-o_halo+1:kp+o_halo),&
                            q(:,cst(cat_c):cen(cat_c)),1,.false.,0)		
        enddo
#endif	
	endif
    ! rain 
    !vqr=0._wp
    if(sum(q(1:kp,cst(cat_r)+1)).gt.qsmall) then
        adv_l(2)=.true.
		where(isnan(vqr))
			vqr=0._wp
		end where
		vqr(kp+1:kp+o_halo)=vqr(kp+o_halo)
		n_step(2)=max(ceiling(maxval(vqr(-o_halo+1:kp+o_halo)*dt/dz*2._wp)),1)
		vqr(1-o_halo:kp+o_halo-1)=-vqr(-o_halo+2:kp+o_halo)
#if MPI_PAMM == 0
		do iter=1,n_step(2)
            call mpdata_vec_1d(dt/real(n_step(2),wp),dz,dzn,&
                            rho,rhon,kp,cen(cat_r)-cst(cat_r)+1,o_halo,o_halo,&
                            vqr(-o_halo+1:kp+o_halo),&
                            q(:,cst(cat_r):cen(cat_r)),1,.false.,0)		
        enddo
#endif	
	endif
	if(ice_flag) then
        ! ice 
        if(sum(q(1:kp,cst(cat_i)+1)).gt.qsmall) then
            adv_l(3)=.true.
            where(isnan(vqi))
                vqi=0._wp
            end where
            vqi(kp+1:kp+o_halo)=vqi(kp+o_halo)
            n_step(3)=max(ceiling(maxval(vqi(-o_halo+1:kp+o_halo)*dt/dz*2._wp)),1)
            vqi(1-o_halo:kp+o_halo-1)=-vqi(-o_halo+2:kp+o_halo)
#if MPI_PAMM == 0
            do iter=1,n_step(3)
                call mpdata_vec_1d(dt/real(n_step(3),wp),dz,dzn,&
                                rho,rhon,kp,cen(cat_i)-cst(cat_i)+1,o_halo,o_halo,&
                                vqi(-o_halo+1:kp+o_halo),&
                                q(:,cst(cat_i):cen(cat_i)),1,.false.,0)		
            enddo
#endif	
        endif
        
    endif
 	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    

    end subroutine p_microphysics_1d
    
    
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates variables for Seifert and Beheng (2001) autoconversion scheme
	!>@param[inout] praut, pracw, rcwaut, rcwacr, rraut, rrsel, rwsel: variables for SB 
	!>@param[in] qc: qcloud
	!>@param[in] nc: ncloud
	!>@param[in] qr: qrain
	!>@param[in] nr: nrain
	!>@param[in] rho: density of air
	!>@param[in] dt: timestep
    subroutine seifert_beheng(praut,pracw, rcwaut, rcwacr, rraut, &
		                    rrsel, rcwsel, qc,nc,qr,nr,rho,dt)
	use advection_1d
    implicit none
    ! arguments:
    real(wp), intent(inout) :: praut,pracw, rcwaut, rcwacr, rraut, &
		                    rrsel, rcwsel
	real(wp), intent(in) :: qc,nc,qr,nr,rho, dt
	real(wp) :: lc, lr, nc1,nr1,xc_bar, phi_au, phi_ac, b_slope, tau, factor1, factor2, &
	            test
	
	
	
	    praut=0._wp
	    pracw=0._wp
	    rcwaut=0._wp
	    rcwacr=0._wp
	    rraut=0._wp
	    rrsel=0._wp
	    rcwsel=0._wp
	    
	    lc=qc*rho
	    lr=qr*rho
	    
	    nc1=max(nc*rho,lc/xstar)
	    nr1=max(nr*rho,1._wp)
	    
		b_slope=((nc1+1.e-20_wp)/(lc+1.e-20_wp))*(gam2c+1.e-20_wp)/(gam1c +1.e-20_wp)
		xc_bar=gam2c/(gam1c*b_slope+1.e-20_wp)
		
		tau=1._wp-lc/(lc+lr+1.e-20_wp)
		tau=max(tau,1.e-6_wp)
		phi_au=600._wp*tau**0.68*(1._wp-tau**0.68)**3
		phi_ac=(tau/(tau+5.e-4_wp))**4

		if (lc .gt. qsmall) then
            ! autoconversion: equation a1 in Seifert and Beheng (2001, atmos res)
            praut = kc/(20._wp*xstar)*(alpha_c+2._wp)*(alpha_c+4._wp)/(alpha_c+1._wp)**2 * &
                    (lc*xc_bar)**2*(1._wp+phi_au/(1._wp-tau+1.e-20_wp)**2)
            
            ! accretion: equation a2 in Seifert and Beheng (2001, atmos res)
            pracw=kr*lc*lr*phi_ac
        
            ! cloud number autoconversion: equation a5 in Seifert and Beheng (2001, atmos res)
            rcwaut=-2._wp/xstar*praut
            ! cloud num accretion: equation a6 in Seifert and Beheng (2001, atmos res)
            rcwacr=-1._wp/xc_bar*pracw
            ! rain num autoconversion: equation a7 in Seifert and Beheng (2001, atmos res)
            rraut=-1._wp/2._wp*rcwaut
            ! rain num self collection: equation a8 in Seifert and Beheng (2001, atmos res)
            rrsel=-kr*nr*lr
            ! cloud num self collection: equation a9 in Seifert and Beheng (2001, atmos res)
            rcwsel=-kr*(alpha_c+2._wp)/(alpha_c+1._wp)*lc**2-rcwaut

            factor1=min(lc/dt,praut+pracw)/(praut+pracw)
            factor2=min(nc1/dt,-(rcwaut+rcwacr+rcwsel))/(-(rcwaut+rcwacr+rcwsel))

            praut=praut*factor1
            pracw=pracw*factor1
            rcwaut=rcwaut*factor2
            rcwacr=rcwacr*factor2
            rcwsel=rcwsel*factor2

            
            praut=praut/rho
            pracw=pracw/rho
            rcwaut=rcwaut/rho
            rcwacr=rcwacr/rho
            rraut=rraut/rho   
            rrsel=rrsel/rho
            rcwsel=rcwsel/rho
        endif
        

    end subroutine seifert_beheng
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
      
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! saturation vapour pressure over liquid                                       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the saturation vapour pressure over liquid water according to buck fit
	!>@param[in] t: temperature
	!>@return svp_liq: saturation vapour pressure over liquid water
	function svp_liq(t)
		use numerics_type
		implicit none
		real(wp), intent(in) :: t
		real(wp) :: svp_liq
		svp_liq = 100._wp*6.1121_wp* &
			  exp((18.678_wp - (t-ttr)/ 234.5_wp)* &
			  (t-ttr)/(257.14_wp + (t-ttr)))
	end function svp_liq
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! saturation vapour pressure over ice										   !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the saturation vapour pressure over ice water according to buck fit
	!>@param[in] t: temperature
	!>@return svp_ice: saturation vapour pressure over ice water
	function svp_ice(t)
		use numerics_type
		implicit none
		real(wp), intent(in) :: t
		real(wp) :: svp_ice
		svp_ice = 100._wp*6.1115_wp* &
			  exp((23.036_wp - (t-ttr)/ 333.7_wp)* &
			  (t-ttr)/(279.82_wp + (t-ttr)))
	end function svp_ice
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
    
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Viscosity of air - Page 417 Pruppacher and Klett							   !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the viscosity of air vs temperature
	!>@param[in] t: temperature
	!>@return viscosity_air: viscosity of air
	elemental function viscosity_air(t)
		use numerics_type
		implicit none
		real(wp), intent(in) :: t
		real(wp) :: viscosity_air
		real(wp) :: tc

		tc = t-ttr
		tc = max(tc,-200._wp)

		if( tc.ge.0._wp) then
			viscosity_air = (1.718_wp+0.0049_wp*tc) * 1E-5_wp ! the 1d-5 converts from poise to si units
		else
			viscosity_air = (1.718_wp+0.0049_wp*tc-1.2e-5_wp*tc**2) * 1e-5_wp
		end if
	end function viscosity_air
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


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
	  real(wp) :: dd, t1
	  t1=max(t,200._wp)
	  dd=2.11e-5_wp*(t1/ttr)**1.94_wp*(101325_wp/p)
	end function dd
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
	  real(wp) :: ka, t1
	  t1=max(t,200._wp)
	  ka=(5.69_wp+0.017_wp*(t1-ttr))*1e-3_wp*joules_in_a_cal
	end function ka
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! hm factor                                                                    !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the factor to multiply the peak production rate by
	!>@param[in] t: temperature
	!>@return hm_func: factor to multiple hm by
	function hm_func(t)
		use numerics_type
		implicit none
		real(wp), intent(in) :: t
		real(wp) :: hm_func
		hm_func=(min(max((t-265.65) / 2.5_wp,0._wp),1._wp) + &
		        min(max((270.65-t) / 2.5_wp,0._wp),1._wp)) -1.0_wp
	end function hm_func
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! INP source function                                                          !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the number concentration of INPs using the DeMott et al. (2010)
	!> parameterisation (Predicting global atmospheric ice...
	!>                    https://doi.org/10.1073/pnas.0910818107)
	!>@param[in] t, naer05: temperature, number concentration of aerosols > 0.5 um
	!>@return demott_2010: number concentration of INPs
	function demott_2010(t,naer05)
		use numerics_type
		implicit none
		real(wp), intent(in) :: t,naer05
		real(wp) :: demott_2010
		real(wp) :: tc
		tc=ttr-t
		! equation 1 from
		! https://www.pnas.org/content/107/25/11217
		! number per std m^3
		demott_2010=min(0.0594_wp*(tc)**3.33_wp * &
		    (naer05/1.e6_wp)**(0.0264_wp*tc+0.0033_wp),naer05)
		
	end function demott_2010
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! ice nucleation from aerosol                                                  !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief calculate the number of active INPs and the threshold diameter for 
	!> activation
	!>@param[in] flag for type of nucleation 1=demott, 2=stochastic....
	!>@param[in] n_aer,sig_aer,d_aer,T,icen, nc, qc, j_stochastic, dt
	!>@param[inout] nin,din
    subroutine ice_nucleation_aerosol(nin,din, &
                n_aer, &    ! number
                sig_aer, &  ! sigma 
                d_aer, t ,icen, nc, qc, &
                j_stochastic, dt,ice_nuc_flag)     ! d
    use numerics_type
    implicit none
    integer(i4b), intent(in) :: ice_nuc_flag
    real(wp), intent(inout) :: nin, din
    real(wp), intent(in) :: n_aer,sig_aer,d_aer, t,icen, nc, qc, &
                        j_stochastic,dt

    real(wp) :: naer05, x, arg, dq

    if(t>= t_hom) then  
        if(ice_nuc_flag.eq.1) then
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! DeMott 2010 nucleation                         !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            naer05=ln_part_mom(0,0.5e-6_wp,n_aer,sig_aer,d_aer)
            ! source function
            nin=demott_2010(t,naer05)
            ! limit nucleation
            nin=max(nin-icen,0._wp)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        elseif(ice_nuc_flag.eq.2) then
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! basic stochastic nucleation                    !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            nin=j_stochastic*n_aer*dt
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        endif  
    else
        ! homogeneous nucleation
        dq = min(max((t_hom-t)*cp/lf*dt, 0._wp), qc)
        nin = dq/qc*nc
        
    endif  
    
    
    

    ! deplete aerosol up to this diameter - using erfinv
    ! limit the argument so that it is not equal to -1 or +1
    arg=max(min(((1._wp-nin/n_aer)*2._wp-1._wp),1._wp-small_number),-1._wp+small_number)

    ! re-calculate nin based on limited value of arg
    nin=(1._wp-(1._wp+arg)/2._wp)*n_aer
    
    ! inverse erf
    call erfinv(arg,x)
    

    ! but x is equal to log(d/dm)/(sig_aer*sqrt(2))
    din=exp(x*sig_aer*sqrt(2._wp)+log(d_aer))
    
    if((din<0.5e-6_wp) .and. ((t>t_hom).and.(ice_nuc_flag.eq.1))) din=0.5e-6_wp
    
    end subroutine ice_nucleation_aerosol
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! calculate new volume and phi                                                       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief calculates ice growth model of Chen and Lamb (1994) 
	!>@param[in] t,qv, qvsat,rhoa,dm,gamma_t,dep_density
	!>@param[inout] v,phi
    subroutine chen_and_lamb_prop(dm,gamma_t,v,phi, dep_density)
        use numerics_type
        implicit none
        real(wp), intent(in) :: dm, gamma_t,dep_density
        real(wp), intent(inout) :: v, phi

        real(wp) :: deltaV,v_old,rgamma_tp2,ln_vn_vo
        integer(i4b) :: i
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! increment to volume of crystals - see equation 41
        ! note that this will be per kg of air, rather than crystal but, since we are 
        ! taking the ratio to determine c and a-axes, it should not matter
        deltaV=dm/dep_density
        v_old=v
        v=v+deltaV
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! solving equations 43 and 43 over dV
        ! i.e. d (ln c /a) = (gam-1)/(gam+2) *d ln v
        !      1/phi * d phi = (gam-1)/(gam+2) / v * dv
        !      ln (phi2/phi1) = (gam-1)/(gam+2) * ln(v2/v1)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        rgamma_tp2=1._wp/(gamma_t+2._wp)
        ln_vn_vo=log(v/v_old)
        phi=phi*exp((gamma_t-1._wp)*rgamma_tp2*ln_vn_vo)       
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
    end subroutine chen_and_lamb_prop
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Chen and Lamb (1994) ancillary variables                                           !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief calculates ancillary variables for ice growth model of Chen and Lamb (1994) 
	!>@param[in] t,qv, qvsat,rhoa
	!>@param[inout] v,phi,gamma_t,dep_density
    subroutine chen_and_lamb_anc(t,qv,qvsat,rhoa,gamma_t, dep_density)
        use numerics_type
        implicit none
        real(wp), intent(in) :: t,qv,qvsat,rhoa
        real(wp), intent(inout) :: gamma_t,dep_density

        real(wp) :: delta_rho,t1
        integer(i4b) :: i
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! calculate the inherent growth ratio - this is from a 17th order polynomial
        gamma_t=0._wp
        t1=min(max(t,243.15),273.15) ! range of fit
        do i=1,n_cl
            gamma_t=gamma_t+((t1-gam_mu_cl(1))/gam_mu_cl(2))**(n_cl-i)*gam_cl(i)
        enddo
        gamma_t=10._wp**gamma_t
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! equation 42 from Chen and Lamb (1994, JAS: The Theoretical Basis for 
        !   Parameterisation of Ice Crystal Habits)
        delta_rho=(qv-qvsat)*rhoa*1000._wp ! g/m^3
        dep_density=rhoi*exp(-3._wp*max(delta_rho-0.05_wp,0._wp)/gamma_t)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end subroutine chen_and_lamb_anc
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Chen and Lamb (1994) capacitance factors                                           !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief calculates the ratio of capacitance to that of an equivalent spehre
	!  Chen and Lamb (1994) 
	!>@param[in] phi
	!>@return cap_fac
    function chen_and_lamb_cap_fac(phi)
        use numerics_type
        implicit none
        real(wp), intent(in) :: phi
        real(wp) :: chen_and_lamb_cap_fac
        real(wp) :: fac1,fac2,ecc

        
        ! factor to convert between R and a - derived from equating volume of sphere to 
        ! volume of spheroid and taking the ratio of a / r
        fac1=(1._wp/(phi))**(1/3)
        
        ! factor to convert between a and capacitance
        if(phi<0.99_wp) then
            ! see equation 39 of Chen and Lamb (1994)
            ecc=sqrt(1._wp-phi**2)
            fac2=ecc/asin(ecc)
        elseif(phi>1.01_wp) then
            ! see equation 40 of Chen and Lamb (1994)
            ecc=sqrt(1._wp-(1._wp/phi)**2)
            fac2=(ecc)/log((1._wp+ecc)*phi)*phi/fac1
        else
            fac2=1._wp
        endif
        
        ! total factor
        chen_and_lamb_cap_fac=fac2
        
    end function chen_and_lamb_cap_fac
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! riming integral over size distribution                                       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates integral over distribution for riming
    function integral_rime_hw(x)
        use numerics_type, only : wp, i4b
        implicit none
		real(wp), dimension(:), intent(in) :: x
		real(wp), dimension(size(x)) :: integral_rime_hw
		
		real(wp), dimension(size(x)) :: diami, vi
		
		diami=(x/ci_new)**(1.0_wp/di)
        vi=pre_hw_new*(diami**-1)*((1._wp+a_hw_new*diami**(0.5_wp*di))**0.5_wp-1._wp)**2
        
        integral_rime_hw=0.25_wp*pi*diami**(2+alpha_i)*n0i*exp(-lambda0i*diami)*vi* &
            (diami**(1.0_wp-di)) / (ci_new*di) 
            
            
    end function integral_rime_hw
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This evaluates the integrand                                                       !
    ! for collisions - number weighted                                                   !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function dintegral_coll_num_hw(x,y)
        use numerics_type, only : wp, i4b
        implicit none
        real(wp), intent(in) :: x
        real(wp), dimension(:), intent(in) :: y
        real(wp), dimension(size(y)) :: dintegral_coll_num_hw
        real(wp) :: diamr, mr, vr, n,nt,nb,mb,mt
        real(wp), dimension(size(y)) :: mi, diami, delv, vi
        integer(i4b) :: i


        mr=x
        mi=y
        diamr=(mr/cr)**(1.0_wp/dr)
        diami=(mi/ci_new)**(1.0_wp/di)
        ! fall-speeds
        vr=a_r*diamr**b_r
!         vi=a_i*diami**b_i
        vi=pre_hw_new*(diami**-1)*((1._wp+a_hw_new*diami**(0.5_wp*di))**0.5_wp-1._wp)**2
        delv=abs(vr-vi)
        !delv=max((vx+vy)/8.0,abs(vx-vy))
        ! last bit is to convert to integral over m
        dintegral_coll_num_hw=eri*pi/4.0_wp*(diamr+diami)**2* &
            delv*n0r*diamr**alpha_r* &
            exp(-lambda0r*diamr)*n0i*diami**alpha_i*exp(-lambda0i*diami)* &
            (diamr**(1.0_wp-dr)) / (cr*dr)*(diami**(1.0_wp-di)) / (ci_new*di)
        
    end function dintegral_coll_num_hw
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This evaluates the integrand                                                       !
    ! for collisions - number weighted                                                   !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function dintegral_coll_ice_num_hw(x,y)
        use numerics_type, only : wp, i4b
        implicit none
        real(wp), intent(in) :: x
        real(wp), dimension(:), intent(in) :: y
        real(wp), dimension(size(y)) :: dintegral_coll_ice_num_hw
        real(wp) :: diamr, mr, vr, n,nt,nb,mb,mt
        real(wp), dimension(size(y)) :: mi, diami, delv, vi
        integer(i4b) :: i


        mr=x
        mi=y
        diami=(mi/ci_new)**(1.0_wp/di)
        diamr=(mr/ci_new)**(1.0_wp/di)
        ! fall-speeds
        vr=pre_hw_new*(diamr**-1)*((1._wp+a_hw_new*diamr**(0.5_wp*di))**0.5_wp-1._wp)**2
        vi=pre_hw_new*(diami**-1)*((1._wp+a_hw_new*diami**(0.5_wp*di))**0.5_wp-1._wp)**2
        delv=abs(vr-vi)
        !delv=max((vx+vy)/8.0,abs(vx-vy))
        ! last bit is to convert to integral over m
        dintegral_coll_ice_num_hw=pi/4.0_wp*(diamr+diami)**2* &
            delv*n0i*diamr**alpha_i* &
            exp(-lambda0i*diamr)*n0i*diami**alpha_i*exp(-lambda0i*diami)* &
            (diamr**(1.0_wp-di)) / (ci_new*di)*(diami**(1.0_wp-di)) / (ci_new*di)
        
    end function dintegral_coll_ice_num_hw
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This evaluates the integrand                                                       !
    ! for collisions - mass weighted praci                                               !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function dintegral_coll_mass1_hw(x,y)
        use numerics_type, only : wp, i4b
        implicit none
        real(wp), intent(in) :: x
        real(wp), dimension(:), intent(in) :: y
        real(wp), dimension(size(y)) :: dintegral_coll_mass1_hw
        real(wp) :: diamr, mr, vr, n,nt,nb,mb,mt
        real(wp), dimension(size(y)) :: mi, diami, delv, vi
        integer(i4b) :: i


        mr=x
        mi=y
        diamr=(mr/cr)**(1.0_wp/dr)
        diami=(mi/ci_new)**(1.0_wp/di)
        ! fall-speeds
        vr=a_r*diamr**b_r
!         vi=a_i*diami**b_i
        vi=pre_hw_new*(diami**-1)*((1._wp+a_hw_new*diami**(0.5_wp*di))**0.5_wp-1._wp)**2
        delv=abs(vr-vi)
        !delv=max((vx+vy)/8.0,abs(vx-vy))
        ! last bit is to convert to integral over m
        dintegral_coll_mass1_hw=ci_new*diami**di*eri*pi/4.0_wp*(diamr+diami)**2* &
            delv*n0r*diamr**alpha_r* &
            exp(-lambda0r*diamr)*n0i*diami**alpha_i*exp(-lambda0i*diami)* &
            (diamr**(1.0_wp-dr)) / (cr*dr)*(diami**(1.0_wp-di)) / (ci_new*di)
        
    end function dintegral_coll_mass1_hw
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This evaluates the integrand                                                       !
    ! for collisions - mass weighted piacr                                               !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function dintegral_coll_mass2_hw(x,y)
        use numerics_type, only : wp, i4b
        implicit none
        real(wp), intent(in) :: x
        real(wp), dimension(:), intent(in) :: y
        real(wp), dimension(size(y)) :: dintegral_coll_mass2_hw
        real(wp) :: diamr, mr, vr, n,nt,nb,mb,mt
        real(wp), dimension(size(y)) :: mi, diami, delv, vi
        integer(i4b) :: i


        mr=x
        mi=y
        diamr=(mr/cr)**(1.0_wp/dr)
        diami=(mi/ci_new)**(1.0_wp/di)
        ! fall-speeds
        vr=a_r*diamr**b_r
!         vi=a_i*diami**b_i
        vi=pre_hw_new*(diami**-1)*((1._wp+a_hw_new*diami**(0.5_wp*di))**0.5_wp-1._wp)**2
        delv=abs(vr-vi)
        !delv=max((vx+vy)/8.0,abs(vx-vy))
        ! last bit is to convert to integral over m
        dintegral_coll_mass2_hw=cr*diamr**dr*eri*pi/4.0_wp*(diamr+diami)**2* &
            delv*n0r*diamr**alpha_r* &
            exp(-lambda0r*diamr)*n0i*diami**alpha_i*exp(-lambda0i*diami)* &
            (diamr**(1.0_wp-dr)) / (cr*dr)*(diami**(1.0_wp-di)) / (ci_new*di)
        
    end function dintegral_coll_mass2_hw
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! mode 1 fragmentation                                                         !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the number of fragments and their mass
    subroutine calculate_mode1(min1,min2,t,n,nt,nb,mb,mt)
        use numerics_type
        implicit none
        real(wp), intent(in) :: min1,min2,t
        real(wp), intent(inout) :: n, nt, nb, mb, mt
        real(wp) :: tc, dthresh, x, beta1,log10zeta, log10nabla, t0, zetab, nablab, tb0, &
            sigma, omega, m,d, fac1
        
        if((min2>min1).or.(min1<=6.55e-11_wp)) then
            ! the ice is more massive than the drop or drop small, don't do it
            n=0._wp
            nt=0._wp
            nb=0._wp
            return
        endif
    
        d = (6._wp*min1/rhow)**(1._wp/3._wp)
        tc=t-ttr
        dthresh = min(d,1.6e-3)
        x = log10(dthresh*1000._wp)
        
        ! table 3, phillips et al.
        beta1 = 0.
        log10zeta = 2.4268_wp*x*x*x + 3.3274_wp*x*x + 2.0783_wp*x + 1.2927_wp
        log10nabla = 0.1242_wp*x*x*x - 0.2316_wp*x*x - 0.9874_wp*x - 0.0827_wp
        t0 = -1.3999_wp*x*x*x - 5.3285_wp*x*x - 3.9847_wp*x - 15.0332_wp
        
        ! table 4, phillips et al. 
        zetab = -0.4651_wp*x*x*x - 1.1072_wp*x*x - 0.4539_wp*x+0.5137_wp
        nablab = 28.5888*x*x*x + 49.8504_wp*x*x + 22.4873_wp*x + 8.0481_wp
        tb0 = 13.3588_wp*x*x*x + 15.7432_wp*x*x - 2.6545_wp*x - 18.4875_wp
        
        sigma = min(max((d-50.e-6_wp)/10.e-6_wp,0._wp), 1._wp)
        omega = min(max((-3._wp-tc)/3._wp,0._wp),1._wp)
        
        n = sigma*omega*(10._wp**log10zeta *(10**log10nabla)**2) / &
            ((tc-t0)**2+(10._wp*log10nabla)**2+beta1*tc)
        
        ! total number of fragments
        n=n*d/dthresh
        ! number of large fragments
        nb = min(sigma*omega*(zetab*nablab**2/((tc-tb0)**2+nablab**2)),n)
        ! number of small fragments
        nt = n-nb
        
        m=oneoversix*rhow*pi*d**3
        
        ! mass of large fragments
        mb=0.4_wp*m
        
        ! mass of small fragments
        mt=oneoversix*rhoi*pi*dtt**3
        
        fac1=min((mt*nt+mb*nb)/min1,1._wp)
        nt = nt *fac1
        nb = nb *fac1
        n=nt+nb
        
    end subroutine calculate_mode1
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! mode 1 fragmentation integral over size distribution                         !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the number of fragments and their mass
    function integral_m1(x)
        use numerics_type, only : wp, i4b
        implicit none
		real(wp), dimension(:), intent(in) :: x
		real(wp), dimension(size(x)) :: integral_m1
		
		real(wp), dimension(size(x)) :: nfrag
		real(wp) :: n,nt,nb,mb,mt
		integer(i4b) :: i
		
		do i=1,size(x)
		    call calculate_mode1(x(i),0._wp,t_send,n,nt,nb,mb,mt)
		    nfrag(i)=n
		enddo
		
    
        integral_m1=n0_freeze*exp(-lam_freeze*x)*x**alpha_r*nfrag
    end function integral_m1
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This evaluates the integrand                                                       !
    ! for mode 1 ice multiplication                                                      !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function dintegral_mode1(x,y)
        use numerics_type, only : wp, i4b
        implicit none
        real(wp), intent(in) :: x
        real(wp), dimension(:), intent(in) :: y
        real(wp), dimension(size(y)) :: dintegral_mode1
        real(wp) :: diamr, mr, vr, n,nt,nb,mb,mt
        real(wp), dimension(size(y)) :: mi, diami, delv, vi
        integer(i4b) :: i


        mr=x
        mi=y
        diamr=(mr/cr)**(1.0_wp/dr)
        diami=(mi/ci)**(1.0_wp/di)
        ! fall-speeds
        vr=a_r*diamr**b_r
        vi=a_i*diami**b_i
        delv=abs(vr-vi)
        !delv=max((vx+vy)/8.0,abs(vx-vy))
        ! last bit is to convert to integral over m
        dintegral_mode1=eri*pi/4.0_wp*(diamr+diami)**2* &
            delv*n0r*diamr**alpha_r* &
            exp(-lambda0r*diamr)*n0i*diami**alpha_i*exp(-lambda0i*diami)* &
            (diamr**(1.0_wp-dr)) / (cr*dr)*(diami**(1.0_wp-di)) / (ci*di)
        
        do i=1,size(y)
            call calculate_mode1(mr,mi(i),t_send,n,nt,nb,mb,mt)
            dintegral_mode1(i)=dintegral_mode1(i)*n
        enddo
    
    end function dintegral_mode1
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This evaluates the integrand                                                       !
    ! for mode 1 ice multiplication                                                      !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function dintegral_mode1_hw(x,y)
        use numerics_type, only : wp, i4b
        implicit none
        real(wp), intent(in) :: x
        real(wp), dimension(:), intent(in) :: y
        real(wp), dimension(size(y)) :: dintegral_mode1_hw
        real(wp) :: diamr, mr, vr, n,nt,nb,mb,mt
        real(wp), dimension(size(y)) :: mi, diami, delv, vi
        integer(i4b) :: i


        mr=x
        mi=y
        diamr=(mr/cr)**(1.0_wp/dr)
        diami=(mi/ci_new)**(1.0_wp/di)
        ! fall-speeds
        vr=a_r*diamr**b_r
!         vi=a_i*diami**b_i
        vi=pre_hw_new*(diami**-1)*((1._wp+a_hw_new*diami**(0.5_wp*di))**0.5_wp-1._wp)**2
        delv=abs(vr-vi)
        !delv=max((vx+vy)/8.0,abs(vx-vy))
        ! last bit is to convert to integral over m
        dintegral_mode1_hw=eri*pi/4.0_wp*(diamr+diami)**2* &
            delv*n0r*diamr**alpha_r* &
            exp(-lambda0r*diamr)*n0i*diami**alpha_i*exp(-lambda0i*diami)* &
            (diamr**(1.0_wp-dr)) / (cr*dr)*(diami**(1.0_wp-di)) / (ci_new*di)
        
        do i=1,size(y)
            call calculate_mode1(mr,mi(i),t_send,n,nt,nb,mb,mt)
            dintegral_mode1_hw(i)=dintegral_mode1_hw(i)*n
        enddo

    
    end function dintegral_mode1_hw
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This evaluates the integrand                                                       !
    ! for mode 2 ice multiplication                                                      !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function dintegral_mode2(x,y)
        use numerics_type, only : wp
        implicit none
        real(wp), intent(in) :: x
        real(wp), dimension(:), intent(in) :: y
        real(wp), dimension(size(y)) :: dintegral_mode2
        real(wp) :: diamr, mr, vr
        real(wp), dimension(size(y)) :: mi, diami, delv, vi, k0, de, nfrag, nfrag_freeze1, &
            nfrag_freeze2


        mr=x
        mi=y
        diamr=(mr/cr)**(1.0_wp/dr)
        diami=(mi/ci)**(1.0_wp/di)
        ! fall-speeds
        vr=a_r*diamr**b_r
        vi=a_i*diami**b_i
        delv=abs(vr-vi)
        !delv=max((vx+vy)/8.0,abs(vx-vy))
        ! last bit is to convert to integral over m
        dintegral_mode2=eri*pi/4.0_wp*(diamr+diami)**2* &
            delv*n0r*diamr**alpha_r* &
            exp(-lambda0r*diamr)*n0i*diami**alpha_i*exp(-lambda0i*diami)* &
            (diamr**(1.0_wp-dr)) / (cr*dr)*(diami**(1.0_wp-di)) / (ci*di)
        
        ! cke from equation 6
        k0=0.5_wp*(mr*mi/(mr+mi))*(vr-vi)**2
        ! de parameter
        de=k0/(gamma_liq*pi*diamr**2)
        ! number of fragments in spalsh
        nfrag=3.0_wp*max(de-decrit, 0.0_wp)
        ! number of fragments in splash that freeze due to mode 1
        nfrag_freeze1=nfrag*f_mode2
        ! number of fragments in splash that freeze due to mode 2
        nfrag_freeze2=nfrag*(1.0_wp-f_mode2)*phi_mode2
        
        dintegral_mode2=dintegral_mode2*nfrag_freeze2
    
    end function dintegral_mode2
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This evaluates the integrand                                                       !
    ! for mode 2 ice multiplication                                                      !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function dintegral_mode2_hw(x,y)
        use numerics_type, only : wp
        implicit none
        real(wp), intent(in) :: x
        real(wp), dimension(:), intent(in) :: y
        real(wp), dimension(size(y)) :: dintegral_mode2_hw
        real(wp) :: diamr, mr, vr
        real(wp), dimension(size(y)) :: mi, diami, delv, vi, k0, de, nfrag, nfrag_freeze1, &
            nfrag_freeze2


        mr=x
        mi=y
        diamr=(mr/cr)**(1.0_wp/dr)
        diami=(mi/ci_new)**(1.0_wp/di)
        ! fall-speeds
        vr=a_r*diamr**b_r
!         vi=a_i*diami**b_i
        vi=pre_hw_new*(diami**-1)*((1._wp+a_hw_new*diami**(0.5_wp*di))**0.5_wp-1._wp)**2
        delv=abs(vr-vi)
        !delv=max((vx+vy)/8.0,abs(vx-vy))
        ! last bit is to convert to integral over m
        dintegral_mode2_hw=eri*pi/4.0_wp*(diamr+diami)**2* &
            delv*n0r*diamr**alpha_r* &
            exp(-lambda0r*diamr)*n0i*diami**alpha_i*exp(-lambda0i*diami)* &
            (diamr**(1.0_wp-dr)) / (cr*dr)*(diami**(1.0_wp-di)) / (ci_new*di)
        
        ! cke from equation 6
        k0=0.5_wp*(mr*mi/(mr+mi))*(vr-vi)**2
        ! de parameter
        de=k0/(gamma_liq*pi*diamr**2)
        ! number of fragments in spalsh
        nfrag=3.0_wp*max(de-decrit, 0.0_wp)
        ! number of fragments in splash that freeze due to mode 1
        nfrag_freeze1=nfrag*f_mode2
        ! number of fragments in splash that freeze due to mode 2
        nfrag_freeze2=nfrag*(1.0_wp-f_mode2)*phi_mode2
        
        dintegral_mode2_hw=dintegral_mode2_hw*nfrag_freeze2
    
    end function dintegral_mode2_hw
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function limit1_coll(x)
        use numerics_type, only : wp
        implicit none
        real(wp), intent(in) :: x
        real(wp) :: limit1_coll
        limit1_coll=mrthresh
    end function limit1_coll
!
    function limit2_coll(x)
        use numerics_type, only : wp
        implicit none
        real(wp), intent(in) :: x
        real(wp) :: limit2_coll
        limit2_coll=mrupper
    end function limit2_coll


    function limit1_mode1(x)
        use numerics_type, only : wp
        implicit none
        real(wp), intent(in) :: x
        real(wp) :: limit1_mode1
        limit1_mode1=mrthresh
    end function limit1_mode1
!
    function limit2_mode1(x)
        use numerics_type, only : wp
        implicit none
        real(wp), intent(in) :: x
        real(wp) :: limit2_mode1
        limit2_mode1=x
    end function limit2_mode1


    function limit1_mode2(x)
        use numerics_type, only : wp
        implicit none
        real(wp), intent(in) :: x
        real(wp) :: limit1_mode2
        limit1_mode2=x
    end function limit1_mode2
!
    function limit1_collisional(x)
        use numerics_type, only : wp
        implicit none
        real(wp), intent(in) :: x
        real(wp) :: limit1_collisional
        limit1_collisional=1.e-6_wp
    end function limit1_collisional
!
    function limit2_mode2(x)
        use numerics_type, only : wp
        implicit none
        real(wp), intent(in) :: x
        real(wp) :: limit2_mode2
        limit2_mode2=miupper
    end function limit2_mode2

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This evaluates the integrand                                                       !
    ! for collisional breakup following vardiman 1978                                    !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function dintegral_collisional_breakup(x,y)
        use numerics_type, only : wp
        implicit none
        real(wp), intent(in) :: x
        real(wp), dimension(:), intent(in) :: y
        real(wp), dimension(size(y)) :: dintegral_collisional_breakup
        real(wp) :: diamr, mr, vr
        real(wp), dimension(size(y)) :: mi, diami, delv, vi, k0, de, nfrag, delm


        mr=x
        mi=y
        diamr=(mr/ci_new)**(1.0_wp/di)
        diami=(mi/ci_new)**(1.0_wp/di)
        ! fall-speeds
        vr=a_i*diamr**b_i
        vi=a_i*diami**b_i
        delv=abs(vr-vi)

        ! last bit is to convert to integral over m
        ! need to multiply by eii(k) - f_mode2
        dintegral_collisional_breakup=f_mode2*pi/4.0_wp*(diamr+diami)**2* &
            delv*n0i*diamr**alpha_i* &
            exp(-lambda0i*diamr)*n0i*diami**alpha_i*exp(-lambda0i*diami)* &
            (diamr**(1.0_wp-di)) / (ci_new*di)*(diami**(1.0_wp-di)) / (ci_new*di)
        
        ! calculate the change in momentum - equation 7
        ! assume a coefficient of restitution of 0.5?
        ! units are g cm s-1
        delm = 0.25_wp*pi*mr*mi/(mr+mi)*(1._wp+0.5_wp)*abs(vr-vi)*1.e5_wp 
        nfrag = vard02(1)*(log(delm)**2)+vard02(2)*log(delm)+vard02(3)
        
        dintegral_collisional_breakup = dintegral_collisional_breakup * nfrag
            
             
    end function dintegral_collisional_breakup
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This evaluates the integrand                                                       !
    ! for collisional breakup following vardiman 1978                                    !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function dintegral_collisional_breakup_hw(x,y)
        use numerics_type, only : wp
        implicit none
        real(wp), intent(in) :: x
        real(wp), dimension(:), intent(in) :: y
        real(wp), dimension(size(y)) :: dintegral_collisional_breakup_hw
        real(wp) :: diamr, mr, vr
        real(wp), dimension(size(y)) :: mi, diami, delv, vi, k0, de, nfrag, delm


        mr=x
        mi=y
        diamr=(mr/ci_new)**(1.0_wp/di)
        diami=(mi/ci_new)**(1.0_wp/di)
        ! fall-speeds
        vr=pre_hw_new*(diamr**-1)*((1._wp+a_hw_new*diamr**(0.5_wp*di))**0.5_wp-1._wp)**2
        vi=pre_hw_new*(diami**-1)*((1._wp+a_hw_new*diami**(0.5_wp*di))**0.5_wp-1._wp)**2
        vr=max(vr,0._wp)
        vi=max(vi,0._wp)
        delv=abs(vr-vi)

        ! last bit is to convert to integral over m
        ! need to multiply by eii(k) - f_mode2
        dintegral_collisional_breakup_hw=f_mode2*pi/4.0_wp*(diamr+diami)**2* &
            delv*n0i*diamr**alpha_i* &
            exp(-lambda0i*diamr)*n0i*diami**alpha_i*exp(-lambda0i*diami)* &
            (diamr**(1.0_wp-di)) / (ci_new*di)*(diami**(1.0_wp-di)) / (ci_new*di)
        
        ! calculate the change in momentum - equation 7
        ! assume a coefficient of restitution of 0.5?
        ! units are g cm s-1
        delm = 0.125_wp*pi*mr*mi/(mr+mi)*(1._wp+0.5_wp)*delv*1.e5_wp 
        delm = min(max(delm,exp(-vard02(2)/(2*vard02(1)))),1._wp)
        
        nfrag = vard02(1)*(log(delm)**2)+vard02(2)*log(delm)+vard02(3)
        
        dintegral_collisional_breakup_hw = dintegral_collisional_breakup_hw * nfrag
            
             
    end function dintegral_collisional_breakup_hw
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This evaluates the integrand                                                       !
    ! for collisional breakup following Phillips et al (2017)                            !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function dintegral_collisional_breakup2_hw(x,y)
        use numerics_type, only : wp
        implicit none
        real(wp), intent(in) :: x
        real(wp), dimension(:), intent(in) :: y

        real(wp), dimension(size(y)) :: dintegral_collisional_breakup2_hw
        real(wp) :: diamr, mr, vr
        real(wp), dimension(size(y)) :: mi, diami, delv, vi, k0, nfrag, delm,vol2, &
                    twicea2, dmax2
        real(wp) :: dmax1, twicea1, dsmall, frimes, frimel, &
                    rhois, phis, &
                    alpha, a0, t0, a, c, nmax, gamma, zeta, &
                    rhoi1,rhoi2,frime1,frime2, phi1, phi2, vol1, t=0._wp
        integer(i4b) :: i

		t=t_send
        mr=x
        mi=y
        diamr=(mr/ci_new)**(1.0_wp/di)
        diami=(mi/ci_new)**(1.0_wp/di)
        ! fall-speeds
        vr=pre_hw_new*(diamr**-1)*((1._wp+a_hw_new*diamr**(0.5_wp*di))**0.5_wp-1._wp)**2
        vi=pre_hw_new*(diami**-1)*((1._wp+a_hw_new*diami**(0.5_wp*di))**0.5_wp-1._wp)**2
        vr=max(vr,0._wp)
        vi=max(vi,0._wp)
        delv=abs(vr-vi)

        ! last bit is to convert to integral over m
        ! need to multiply by eii(k) - f_mode2
        dintegral_collisional_breakup2_hw=f_mode2*pi/4.0_wp*(diamr+diami)**2* &
            delv*n0i*diamr**alpha_i* &
            exp(-lambda0i*diamr)*n0i*diami**alpha_i*exp(-lambda0i*diami)* &
            (diamr**(1.0_wp-di)) / (ci_new*di)*(diami**(1.0_wp-di)) / (ci_new*di)
        
        
        ! now Phillips et al++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! calculate particle properties from moments
        ! 1=vol, 2=mass, 3=rime mass, 4=phi, 5=number
        frime1=phillips_br_workspace(3)/phillips_br_workspace(2)
        frime2=frime1
        phi1=phillips_br_workspace(4)/phillips_br_workspace(5)
        phi2=phi1
        rhoi1 = min(rhoi,(phillips_br_workspace(2)-phillips_br_workspace(3)) / &
                phillips_br_workspace(1))
        rhoi2 = rhoi1      
        vol1=mr/rhoi1
        vol2=mi/rhoi1
        
        
        
        ! calculate the max length of the ice crystals
        twicea1 = (6._wp*vol1 / (pi*phi1))**oneoverthree
        dmax1 = max(twicea1, dmax1*phi1  )
        
        twicea2=(6._wp*vol2 / (pi*phi2))**oneoverthree
        dmax2 = max(twicea2, dmax2*phi2  )
        
        ! calculate the max dimension of the particle assuming rime fills in like a sphere
        dmax1=max(dmax1,  &
         (6._wp*oneoverpi*(pi*twicea1**3*phi1*oneoversix)+ &
            mr*frime1/rhoi  )**oneoverthree)
        dmax2=max(dmax2,  &
         (6._wp*oneoverpi*(pi*twicea2**3*phi2*oneoversix)+&
            mi*frime2/rhoi  )**oneoverthree)
        
        
        
        do i=1,size(y)
            ! some swapping
            if(dmax1<dmax2(i)) then
                dsmall = dmax1
                frimes = frime1
                frimel = frime2
                rhois = rhoi1
                phis = phi1
            else
                dsmall = dmax2(i)
                frimes = frime2
                frimel = frime1
                rhois = rhoi2
                phis = phi2
            endif
        
        
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! table 1: phillips et al. (2017)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! type 1
            alpha=pi*dsmall**2
            if ((dsmall> 5.e-4_wp).and.(dsmall<5.e-3_wp).and.(frimes>=0.5_wp) &
                .and.(frimes<0.9_wp).and.(frimel>=0.5_wp)) then
            
                ! collisions of graupel size 500 micron to 5 mm with graupel and hail
                a0=3.78e4_wp*(1._wp+0.0079_wp/dsmall**1.5_wp)
                t0=-15._wp
                A=a0*oneoverthree+max(2.*a0*oneoverthree- &
                        a0*oneovernine*abs(t-ttr-t0),0._wp)
                c=6.30e6_wp*phi_phillips
                nmax=100._wp
                gamma=0.3_wp
                zeta=0.001_wp
        
            ! type 1
            elseif ((frimes>=0.9_wp).and.(frimel>=0.9_wp)) then
            
                ! collisions of hail and hail - no size constraint
                a0=4.35e5_wp
                t0=-15._wp
                A=a0*oneoverthree+max(2.*a0*oneoverthree- &
                    a0*oneovernine*abs(t-ttr-t0),0._wp)
                c=3.31e5_wp
                nmax=1000._wp
                gamma=0.54_wp
                zeta=1.e-6_wp
            
            ! types 2 or 3
            elseif ((dsmall> 5.e-4_wp).and.(dsmall<5.e-3_wp).and.(frimes<0.5_wp) &
                .and. (phis < 1._wp)) then
                ! collisions of ice /snow size 500 micron to 5 mm with any ice
            
                ! seems like columnar habits dont fragment?
                if(rhois < 400._wp) then
                    ! dendrites
                    A=1.41e6_wp*(1._wp+100._wp*frimes**2)* &
                        (1._wp+3.98e-5_wp/dsmall**1.5_wp)
                    c=3.09e6_wp*frimes
                    nmax=100._wp
                    gamma=0.50_wp - 0.25_wp*frimes
                    zeta=0.001_wp
                
                else
                    ! spatial planar
                    A=1.58e7_wp*(1._wp+100._wp*frimes**2)* &
                        (1._wp+1.33e-4_wp/dsmall**1.5_wp)
                    c=7.08e6_wp*frimes
                    nmax=100._wp
                    gamma=0.50_wp - 0.25_wp*frimes
                    zeta=0.001_wp
            
                endif
            else
                alpha=0._wp
                nmax=0._wp            
            endif
            ! CKE
            k0(i) = 0.5_wp*(mr*mi(i)/(mr+mi(i)))*(vr-vi(i))**2
            ! finally apply equation 13
            nfrag(i) = min(alpha*A*(1._wp-exp(-(C*K0(i)/(alpha*A))**gamma )), nmax)
            !---------------------------------------------------------------------------------
        enddo
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        dintegral_collisional_breakup2_hw = &
            dintegral_collisional_breakup2_hw * nfrag
                    
             
    end function dintegral_collisional_breakup2_hw
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate the a parameter in the heymsfield and westbrook fall-speed scheme        !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester, 2021
	!>@brief calculates the "constant term" in the reynolds number for Heymsfield and 
	! Westbrook (2010)
	!>@param[in] num,vol, mass, nmon, phi, rim, t, rhoa
	!>@param[out] a
    elemental function &
        heymsfield_and_westbrook_fall_parameters1(num,vol,mass,nmon,&
                                        phi,rim,t,rhoa) result(a)
        use numerics_type, only : wp
        implicit none
        real(wp), intent(in) :: num,vol, mass, nmon, phi, rim, t, rhoa
        real(wp) :: a
        
        
!         a=0.0625_wp/sqrt(0.175_wp) / &
!             (viscosity_air(t))*2._wp*sqrt(rhoa*grav/(6._wp)*(920._wp))
        a=0.0625_wp/(sqrt(0.175_wp) *viscosity_air(t))* &
            2._wp*sqrt(rhoa*grav*0.1667_wp * &
                min((mass)/(vol+rim/920._wp),920._wp))
        
    
    end function heymsfield_and_westbrook_fall_parameters1
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate the pre-factor in the heymsfield and westbrook fall-speed scheme         !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester, 2021
	!>@brief calculates the "pre-factor" in the reynolds number for Heymsfield and 
	! Westbrook (2010)
	!>@param[in] t, rhoa
	!>@param[out] a
    elemental function &
        heymsfield_and_westbrook_fall_parameters2(t,rhoa) result(pre)
        use numerics_type, only : wp
        implicit none
        real(wp), intent(in) :: t, rhoa
        real(wp) :: pre
        
        
        pre=16._wp*viscosity_air(t)/rhoa
        
    
    end function heymsfield_and_westbrook_fall_parameters2
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate the collisions between ice particles                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester, 2023
	!>@brief calculates the collisions between ice particles
	!>@param[in] n_i, rho, lam_i, qi,a_hw1,pre_hw,rime_m,vol_xtal
	!>@param[inout] riaci
    subroutine collisions_between_ice_particles(riaci, &
                n_i, rho, lam_i, &
                qi,vol_xtal,rime_m,a_hw1,pre_hw)
        implicit none
        real(wp), intent(in) :: n_i, rho, lam_i, &
            qi,a_hw1,pre_hw,rime_m,vol_xtal
        real(wp), intent(inout) :: riaci
                 
                
        ! actual integrals
        riaci=0._wp
        ci_new=pi/6._wp*min(910._wp, &
            qi/(vol_xtal+rime_m/500._wp))
        lambda0i=lam_i
        mrthresh=ci_new*1.e-6_wp**di
        mrupper=ci_new*(pthreshi/lambda0i)**di
        ! only call integral if mrupper gt mrthresh
        if((mrupper.gt.mrthresh).and.(qi.gt.qsmall)) then
            n0i=n_i

            a_hw_new=a_hw1
            pre_hw_new=pre_hw
            call quad2d_qgaus(dintegral_coll_ice_num_hw, &
                limit1_coll,limit2_coll,mrthresh,mrupper,riaci)
         endif
                
			
    end subroutine collisions_between_ice_particles
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
    
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate the collisions between precipitating particles of different species      !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester, 2023
	!>@brief calculates the collisions between precipitating particles of different species
	!>@param[in] nr,ni,rho, vqi, vqr, vni,vnr,lamr, lami
	!>@param[inout] praci, rraci, piacr, riacr
    subroutine collisions_between_precipitating_particles(praci,rraci,piacr, &
                riacr, n_r, n_i, rho, vqr, vqi,vnr, vni, lam_r, lam_i, &
                qi,vol_xtal,rime_m,qr,a_hw1,pre_hw,heyms_west)
        implicit none
        real(wp), intent(in) :: n_r, n_i, rho,vqr, vqi, vnr,vni, lam_r, lam_i, &
            qr,qi,a_hw1,pre_hw,rime_m,vol_xtal
        logical, intent(in) :: heyms_west
        real(wp), intent(inout) :: praci, rraci, piacr, riacr
        
        ! rain mass collected by ice
        praci=max(n_r*n_i*pi/(4._wp*rho)*eri*ci*max((vqi+vqr)/8._wp,abs(vqi-vqr)) * &
                ( &
                mass_raci1/(lam_r**(1._wp+alpha_r) *lam_i**(3._wp+alpha_i+di)) + &
                mass_raci2/(lam_r**(2._wp+alpha_r) *lam_i**(2._wp+alpha_i+di)) + &
                mass_raci3/(lam_r**(3._wp+alpha_r) *lam_i**(1._wp+alpha_i+di))  &
                ) , 0._wp)  
        ! rain number collected by ice
        rraci=max(n_i*n_r*pi/(4._wp*rho)*eri*max((vnr+vni)/8._wp,abs(vnr-vni)) * &
                ( &
                num_raci1/(lam_r**(1._wp+alpha_r) *lam_i**(3._wp+alpha_i)) + &
                num_raci2/(lam_r**(2._wp+alpha_r) *lam_i**(2._wp+alpha_i)) + &
                num_raci3/(lam_r**(3._wp+alpha_r) *lam_i**(1._wp+alpha_i))  &
                )    , 0._wp)    

        ! ice mass collected by rain
        piacr=max(n_i*n_r*pi/(4._wp*rho)*eri*cr*max((vqr+vqi)/8._wp,abs(vqr-vqi)) * &
                ( &
                mass_iacr1/(lam_i**(1._wp+alpha_i) *lam_r**(3._wp+alpha_r+dr)) + &
                mass_iacr2/(lam_i**(2._wp+alpha_i) *lam_r**(2._wp+alpha_r+dr)) + &
                mass_iacr3/(lam_i**(3._wp+alpha_i) *lam_r**(1._wp+alpha_r+dr))  &
                )    , 0._wp)    
        ! ice number collected by rain
        riacr=max(n_i*n_r*pi/(4._wp*rho)*eri*max((vnr+vni)/8._wp,abs(vnr-vni)) * &
                ( &
                num_iacr1/(lam_i**(1._wp+alpha_i) *lam_r**(3._wp+alpha_r)) + &
                num_iacr2/(lam_i**(2._wp+alpha_i) *lam_r**(2._wp+alpha_r)) + &
                num_iacr3/(lam_i**(3._wp+alpha_i) *lam_r**(1._wp+alpha_r))  &
                )    , 0._wp)  
                
                
 
 
            
                
                
        ! actual integrals
        riacr=0._wp
        piacr=0._wp
        praci=0._wp
        lambda0r=lam_r
        lambda0i=lam_i
        mrthresh=cr*1.e-6_wp**dr
        mrupper=cr*(pthreshr/lambda0r)**dr
        miupper=ci*(pthreshi/lambda0i)**di
        mrupper=min(mrupper,miupper)
        ! only call integral if mrupper gt mrthresh
        if((mrupper.gt.mrthresh).and.(qr.gt.qsmall) .and.(qi.gt.qsmall)) then
            n0r=n_r
            n0i=n_i

            if(heyms_west) then
                ci_new=pi/6._wp*min(910._wp, &
                    qi/(vol_xtal+rime_m/500._wp))
                a_hw_new=a_hw1
                pre_hw_new=pre_hw
                call quad2d_qgaus(dintegral_coll_num_hw, &
                    limit1_coll,limit2_coll,mrthresh,mrupper,riacr)
                call quad2d_qgaus(dintegral_coll_mass1_hw, &
                    limit1_coll,limit2_coll,mrthresh,mrupper,praci)
                call quad2d_qgaus(dintegral_coll_mass2_hw, &
                    limit1_coll,limit2_coll,mrthresh,mrupper,piacr)
!             else
!                 call quad2d_qgaus(dintegral_mode1, &
!                     limit1_mode1,limit2_mode1,mrthresh,mrupper,dummy3)                    
            endif
         endif
        rraci=riacr
                
			
    end subroutine collisions_between_precipitating_particles
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate the number of fragments due to mode 1 drop fragmentation when collisions !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester, 2023
	!>@brief calculates mode 1 due to collisions
	!>@param[in]lam_r, lam_i, t, n_r, n_i, qi, vol_xtal, rime_m, qr, a_hw1, pre_hm
	!>@param[inout] nfrag_m1c
    subroutine mode1_sip_collisional(lam_r, lam_i, t, n_r, n_i, &
                    qi, vol_xtal, rime_m, qr, a_hw1,pre_hw, dt, heyms_west,nfrag_m1c)
        implicit none
        logical :: heyms_west
        real(wp), intent(in) :: lam_r, lam_i, t, n_r, n_i, qi, vol_xtal, rime_m, qr, &
            a_hw1, pre_hw, dt
        real(wp), intent(inout) :: nfrag_m1c
        real(wp) :: dummy3
        
        ! calculate the number of fragments
        lambda0r=lam_r
        lambda0i=lam_i
        mrthresh=cr*1.e-6_wp**dr
        mrupper=cr*(pthreshr/lambda0r)**dr
        miupper=ci*(pthreshi/lambda0i)**di
        mrupper=min(mrupper,miupper)
        t_send=t
        nfrag_m1c=0._wp
        ! only call integral if mrupper gt mrthresh
        if((mrupper.gt.mrthresh).and.(qr.gt.qsmall) .and.(qi.gt.qsmall)) then

            n0r=n_r
            n0i=n_i

            if(heyms_west) then
                ci_new=pi/6._wp*min(910._wp, &
                    qi/(vol_xtal+rime_m/500._wp))
                a_hw_new=a_hw1
                pre_hw_new=pre_hw
                call quad2d_qgaus(dintegral_mode1_hw, &
                    limit1_mode1,limit2_mode1,mrthresh,mrupper,dummy3)
            else
                call quad2d_qgaus(dintegral_mode1, &
                    limit1_mode1,limit2_mode1,mrthresh,mrupper,dummy3)                    
            endif
            ! multiplication according to mode-1
            nfrag_m1c=dummy3*dt
         endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end subroutine mode1_sip_collisional
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate the number of fragments due to mode 2 drop fragmentation when collisions !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester, 2023
	!>@brief calculates mode 2 due to collisions
	!>@param[in]lam_r, lam_i, t, n_r, n_i, qi, vol_xtal, rime_m, qr, a_hw1, pre_hm
	!>@param[inout] nfrag_m1c
    subroutine mode2_sip_collisional(lam_r, lam_i, t, n_r, n_i, &
                    qi, vol_xtal, rime_m, qr, a_hw1,pre_hw, dt, heyms_west,nfrag_m2)
        implicit none
        logical :: heyms_west
        real(wp), intent(in) :: lam_r, lam_i, t, n_r, n_i, qi, vol_xtal, rime_m, qr, &
            a_hw1, pre_hw, dt
        real(wp), intent(inout) :: nfrag_m2
        real(wp) :: dummy3
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Mode 2 secondary ice                                                           !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! calculate the number of fragments
        lambda0r=lam_r
        lambda0i=lam_i
        mrthresh=cr*150.e-6_wp**dr
        mrupper=cr*(pthreshr/lambda0r)**dr
        miupper=ci*(pthreshi/lambda0i)**di
        mrupper=min(mrupper,miupper)
    
        ! only call integral if mrupper gt mrthresh
        if((mrupper.gt.mrthresh).and.(qr.gt.qsmall).and.(qi.gt.qsmall)) then

            f_mode2=min(-cw*(t-ttr)/lf,1.0_wp)
            n0r=n_r
            n0i=n_i

            if(heyms_west) then
                ci_new=pi/6._wp*min(910._wp, &
                    qi/(vol_xtal+rime_m/500._wp))
                a_hw_new=a_hw1
                pre_hw_new=pre_hw
                call quad2d_qgaus(dintegral_mode2_hw, &
                    limit1_mode2,limit2_mode2,mrthresh,mrupper,dummy3)
            else
                call quad2d_qgaus(dintegral_mode2, &
                    limit1_mode2,limit2_mode2,mrthresh,mrupper,dummy3)                    
            endif
            ! multiplication according to new lab results...
            nfrag_m2=dummy3*dt
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end subroutine mode2_sip_collisional
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate the number of fragments due to mode 2 drop fragmentation when collisions !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester, 2023
	!>@brief calculates mode 2 due to collisions
	!>@param[in]lam_i, t, n_i, ni, qi, vol_xtal, rime_m, phi1, a_hw1, pre_hm
	!>@param[inout] nfrag_m1c
    subroutine ice2_sip_collisional(lam_i, t, n_i, &
                ni, qi, vol_xtal, rime_m, phi1, a_hw1,pre_hw, dt, heyms_west,nfrag_ii, &
                coll_breakup_flag1)
        implicit none
        logical :: heyms_west
        real(wp), intent(in) :: lam_i, t, n_i, ni,qi, vol_xtal, rime_m, phi1, &
            a_hw1, pre_hw, dt
        integer(i4b), intent(in) :: coll_breakup_flag1
        real(wp), intent(inout) :: nfrag_ii
        real(wp) :: dummy3
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Collisional break-up of ice                                                    !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(coll_breakup_flag1==1) then
        
            ! calculate the number of fragments
            lambda0r=lam_i
            lambda0i=lam_i
            mrthresh=cr*1.e-6_wp**dr
            mrupper=ci*(pthreshi/lambda0r)**di
            miupper=ci*(pthreshi/lambda0i)**di
            mrupper=min(mrupper,miupper)
        
            ! only call integral if mrupper gt mrthresh
            if((mrupper.gt.mrthresh).and.(qi.gt.qsmall)) then

                f_mode2=1._wp !eii(k) ! only need to consider collision, not sticking
                n0r=n_i
                n0i=n_i
                
                if(heyms_west) then
                    ci_new=pi/6._wp*min(910._wp, &
                        qi/(vol_xtal+rime_m/500._wp))
                    a_hw_new=a_hw1
                    pre_hw_new=pre_hw
                    call quad2d_qgaus(dintegral_collisional_breakup_hw, &
                        limit1_collisional,limit2_mode2,mrthresh,mrupper,dummy3)
                else
                    call quad2d_qgaus(dintegral_collisional_breakup, &
                        limit1_collisional,limit2_mode2,mrthresh,mrupper,dummy3)                    
                endif
                ! multiplication according to Vardiman (1978)
                nfrag_ii=dummy3*dt
            endif
        elseif(coll_breakup_flag1==2) then
            ! calculate the number of fragments
            lambda0r=lam_i
            lambda0i=lam_i
            mrthresh=cr*1.e-6_wp**dr
            mrupper=ci*(pthreshi/lambda0r)**di
            miupper=ci*(pthreshi/lambda0i)**di
            mrupper=min(mrupper,miupper)
        	t_send=t
            ! only call integral if mrupper gt mrthresh
            if((mrupper.gt.mrthresh).and.(qi.gt.qsmall)) then

                ! total volume
                phillips_br_workspace(1)=max(vol_xtal,0._wp)
                ! total mass
                phillips_br_workspace(2)=max(qi,0._wp)
                ! total rime mass
                phillips_br_workspace(3)=max(rime_m,0._wp)
                ! total phi
                phillips_br_workspace(4)=max(phi1,0._wp)
                ! total number
                phillips_br_workspace(5)=max(ni,0._wp)


                f_mode2=1._wp !eii(k) ! only need to consider collision, not sticking
                n0r=max(n_i,0._wp)
                n0i=n0r
                
                if(ni>1._wp) then
                
                    ci_new=pi/6._wp*min(910._wp, &
                        qi/(vol_xtal+rime_m/500._wp))
                    a_hw_new=a_hw1
                    pre_hw_new=pre_hw
                    call quad2d_qgaus(dintegral_collisional_breakup2_hw, &
                        limit1_collisional,limit2_mode2,mrthresh,mrupper,dummy3)

                    ! multiplication according to Phillips et al (2017)
                    nfrag_ii=dummy3*dt
                else
                    nfrag_ii=0._wp
                endif
            endif
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end subroutine ice2_sip_collisional
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate the number of ice particles from primary nucleation and drop-frag        !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester, 2023
	!>@brief calculates ice nucleation and mode1
	!>@param[in]lawson, ice_nuc_flag, n_mode, nc, ni, qc, t, j_stochatic, dt, sa, ma,
	!           cloud_flag, mode1_ice_flag
	!>@param[inout] nin_c, din_c, nfrag_nucc, massc_nucc, n_mix, s_mix,
	!>       m_mix, rho, molw, nu, n_aer, sig_aer, d_aer
    subroutine ice_nucleation_and_mode1(n_mode, nc, ni, qc, &
        t, sa, ma, rho, molw, nu, &
        nin_c, din_c, n_mix,s_mix,m_mix, &
        n_aer,sig_aer,d_aer,j_stochastic, dt, ice_nuc_flag, nfrag_nucc, &
        massc_nucc, lawson, cloud_flag, mode1_ice_flag)
        implicit none
        logical, intent(in) :: lawson, cloud_flag
        integer(i4b), intent(in) :: ice_nuc_flag, mode1_ice_flag
        integer(i4b), intent(in) :: n_mode
        real(wp), intent(in) :: nc, ni, qc, t, j_stochastic, dt
        real(wp), intent(inout) :: nin_c, din_c, nfrag_nucc, massc_nucc, n_mix,s_mix,m_mix

	    real(wp), dimension(n_mode-1), intent(in) :: sa,ma
	    real(wp), dimension(n_mode), intent(inout) :: rho,molw,nu
	    real(wp), intent(inout) :: n_aer,sig_aer, d_aer
	    
        real(wp) :: lam_freeze, n0_freeze
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! ice nucleation via immersion                                                   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! set the log normal parameters for the aerosol in cloud-water
        ! to n_mix, s_mix, m_mix 
        ! (and save corresponding sig_aer1(n_mode) and d_aer1(n_mode))
        ! this is the aerosol distribution that we would have if we just took the 
        ! drops and evaporated them
        call ln_params_and_props_from_integral_moms(n_mode, nc, sa, ma, & ! mass 
            n_aer,rho, molw,nu, sig_aer,d_aer,n_mix,s_mix,m_mix)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Ice nucleation                                                                 !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! using the aerosol parameters work out how many nucleate ice 
        ! (nin_c) and down to which size (din_c)
        call ice_nucleation_aerosol(nin_c,din_c, n_mix, sig_aer, & 
            d_aer, t ,ni, nc, qc, j_stochastic,dt,ice_nuc_flag)    
        
        nin_c=min(nin_c,nc)
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! drop fragmentation                                                             !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        nfrag_nucc=0._wp
        massc_nucc=nin_c/(nc+qsmall)*qc ! mass of cloud water frozen
        if(massc_nucc.gt.qsmall) then
            if(cloud_flag) then
                lam_freeze=(nin_c/massc_nucc*gam2c/gam1c)
                n0_freeze = nin_c/gam1c*lam_freeze**(alpha_c+1)
                ! lawson et al
                if (lawson) then
                    nfrag_nucc = 2.5e13_wp*n0_freeze/(cc**(4._wp/dc))* &
                         gam3c/(lam_freeze**(4._wp/dc+1._wp+alpha_c))
                endif
            else
                ! rain water
                lam_freeze=(nin_c/massc_nucc*cr*gam2r/gam1r)**(1._wp/dr)
                n0_freeze = nin_c/gam1r*lam_freeze**(alpha_r+1)
                ! lawson et al
                if(lawson) then
                    nfrag_nucc = 2.5e13_wp*n0_freeze* &
                         gam3r/(lam_freeze**(5._wp+alpha_r))
                elseif(mode1_ice_flag.eq.1) then 
                    ! mode-1 fragmentation
                    mrthresh=cr*1.e-6_wp**dr
                    mrupper=cr*(pthreshr/lam_freeze)**dr
                    t_send=t
                    ! only call integral if mrupper gt mrthresh
                    if((mrupper.gt.mrthresh).and.(qc.gt.qsmall)) then
                        ! multiplication according to mode-1
                        nfrag_nucc=romb(integral_m1,mrthresh,mrupper)
                    endif
                endif
            endif
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end subroutine ice_nucleation_and_mode1
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Scale process rates so that cannot get negative values                             !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester, 2023
	!>@brief calculates ice nucleation and mode1
	!>@param[in] nc, qc, nr, qr, ni, qi , t,dt
	!>@param[inout] praci, rraci, piacr, riacr, nin_c, nin_r, &
	!>          massc_nucc, massr_nucr, 
    !>        piacw, rihal, pidep, pisub, riaci, pimlt, prevp, praut, pracw, &
    !>        rcwaut, rcwacr, rraut, rrsel, rcwsel
    subroutine scale_microphysics(praci,rraci, piacr, riacr, &
        nin_c, nin_r, massc_nucc, massr_nucr, piacw, riacw, rihal, pidep, pisub, risub, &
        riaci, &
        pimlt, rimlt, prevp, rrevp, praut, pracw, rcwacr, rraut, rrsel, rcwaut, &
        rcwsel, qv,nc,qc,nr, qr,ni,qi,t,dt)
        implicit none
        real(wp), intent(inout) :: praci,rraci, piacr, riacr, &
            nin_c, nin_r, massc_nucc, massr_nucr, piacw, riacw, &
            rihal, pidep, pisub, risub, riaci, pimlt, rimlt, prevp, rrevp, &
            praut, pracw, rcwacr, rraut, rrsel, rcwaut, rcwsel
        real(wp), intent(in) :: qv,nc, qc, nr, qr, ni, qi,t,dt


        real(wp) :: factor1, factor2, factor3

        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! set some rates to zero, or max                                                 !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(t.gt.ttr) then
            pimlt=qi/dt
            pisub=0._wp
            pidep=0._wp
            piacr=0._wp
            piacw=0._wp
            massc_nucc=0._wp
            massr_nucr=0._wp
            nin_c=0._wp
            nin_r=0._wp
            riaci=0._wp
            praci=0._wp
            rraci=0._wp
            piacr=0._wp
            riacr=0._wp
            rihal=0._wp
        endif
        pidep=min(qv/dt,pidep)
        riacw = piacw *nc/(qc+qsmall)
        rrevp = prevp * nr/(qr+qsmall)
        rimlt = pimlt *ni / (qi+qsmall)
        risub = pisub *ni / (qi+qsmall)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! sinks of q-cloud                                                               !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        factor1 = (praut+pracw+massc_nucc/dt+piacw)
        if(factor1 .gt. 0._wp) then
            factor2 = min(factor1,qc/dt) / factor1
            praut = praut*factor2
            pracw = pracw*factor2
            massc_nucc = massc_nucc*factor2
            piacw = piacw*factor2
        endif    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! sinks of n-cloud - note number rates are negative                              !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        factor1 = (-rcwaut-rcwacr-rcwsel+nin_c/dt+riacw)
        if(factor1 .gt. 0._wp) then
            factor2 = min(factor1,nc/dt) / factor1
            rcwaut = rcwaut*factor2
            rcwacr = rcwacr*factor2
            rcwsel = rcwsel*factor2
            nin_c = nin_c*factor2
            riacw = riacw*factor2
        endif    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! sinks of q-rain                                                                !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        factor1 = (massr_nucr/dt+prevp+praci)
        if(factor1 .gt. 0._wp) then
            factor2 = min(factor1,qr/dt) / factor1
            massr_nucr = massr_nucr*factor2
            prevp = prevp*factor2
            praci = praci*factor2
        endif    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! sinks of n-rain                                                                !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        factor1 = (-rraut-rrsel+nin_r/dt+rraci+rrevp)
        if(factor1 .gt. 0._wp) then
            factor2 = min(factor1,nr/dt) / factor1
            rraut = rraut*factor2
            rrsel = rrsel*factor2
            nin_r = nin_r*factor2
            rraci = rraci*factor2
            rrevp = rrevp*factor2
        endif    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! sinks of q-ice                                                                 !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        factor1 = (pisub+pimlt)
        if(factor1 .gt. 0._wp) then
            factor2 = min(factor1,qi/dt) / factor1
            pisub = pisub*factor2
            pimlt = pimlt*factor2
        endif    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! sinks of n-ice                                                                 !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        factor1 = (riaci+rimlt+risub)
        if(factor1 .gt. 0._wp) then
            factor2 = min(factor1,ni/dt) / factor1
            riaci = riaci*factor2
            rimlt = rimlt*factor2
            risub = risub*factor2
        endif    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 
               
    end subroutine scale_microphysics
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Scale process rates so that cannot get negative values                             !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester, 2023
	!>@brief calculates ice nucleation and mode1
	!>@param[in] nc, qc, nr, qr, dt
	!>@param[inout] prevp, praut, pracw, &
    !>        rcwaut, rcwacr, rraut, rrsel, rcwsel
    subroutine scale_microphysics_warm(prevp, rrevp, &
        praut, pracw, rcwaut, rcwacr, rraut, rrsel, &
        rcwsel, nc,qc,nr, qr,dt)
        implicit none
        real(wp), intent(inout) :: rrevp, prevp, praut, pracw, &
            rcwaut, rcwacr, rraut, rrsel, rcwsel
        real(wp), intent(in) :: nc, qc, nr, qr, dt


        real(wp) :: factor1, factor2, factor3
        
        rrevp = prevp * nr/(qr+qsmall)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! sinks of q-cloud                                                               !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        factor1 = (praut+pracw)
        if(factor1 .gt. 0._wp) then
            factor2 = min(factor1,qc/dt) / factor1
            praut = praut*factor2
            pracw = pracw*factor2
        endif    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! sinks of n-cloud - note number rates are negative                              !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        factor1 = (-rcwaut-rcwacr-rcwsel)
        if(factor1 .gt. 0._wp) then
            factor2 = min(factor1,nc/dt) / factor1
            rcwaut = rcwaut*factor2
            rcwacr = rcwacr*factor2
            rcwsel = rcwsel*factor2
        endif    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! sinks of q-rain                                                                !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        factor1 = (prevp)
        if(factor1 .gt. 0._wp) then
            factor2 = min(factor1,qr/dt) / factor1
            prevp = prevp*factor2
        endif    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! sinks of n-rain                                                                !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        factor1 = (-rraut-rrsel+rrevp)
        if(factor1 .gt. 0._wp) then
            factor2 = min(factor1,nr/dt) / factor1
            rraut = rraut*factor2
            rrsel = rrsel*factor2
            rrevp = rrevp*factor2
        endif    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





    end subroutine scale_microphysics_warm
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Move aerosol from one category to another based on a process                       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester, 2023
	!>@brief calculates ice nucleation and mode1
	!>@param[in] n_mode, dq, q
	!>@param[inout] n_aer_a, s_aer_a, m_aer_a, n_aer_b, s_aer_b, m_aer_b
    subroutine move_aerosol_proportional( n_mode,n_aer_a, s_aer_a, m_aer_a, &
        n_aer_b, s_aer_b, m_aer_b, dq,q ,recycle)
            implicit none
            integer(i4b), intent(in) :: n_mode
            real(wp), dimension(n_mode-1), intent(inout) :: n_aer_a, s_aer_a, m_aer_a, &
                n_aer_b, s_aer_b, m_aer_b
            real(wp), intent(in) :: dq, q
            logical, intent(in) :: recycle
            
            integer(i4b) :: i
            real(wp) :: dummy1, dummy2, dummy3, rat
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! remove aerosol from one and add to other due to process
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            rat = dq/q
            do i=1,n_mode-1
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! remove from a:                                                !
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! number in aerosol modes
                dummy1=n_aer_a(i)*rat
                n_aer_a(i)=n_aer_a(i)-dummy1 
                
                ! surface area in aerosol modes
                dummy2=s_aer_a(i)*rat
                s_aer_a(i)=s_aer_a(i)-dummy2

                ! mass in aerosol modes
                dummy3=m_aer_a(i)*rat
                m_aer_a(i)=m_aer_a(i)-dummy3
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                if(recycle) then
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! add to aerosol particles in b                                      !
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! number in aerosol modes
                    ! qv, n_mode aerosol + 1
                    n_aer_b(i)=n_aer_b(i)+dummy1 
                
                    ! surface area in aerosol modes
                    s_aer_b(i)=s_aer_b(i)+dummy2
                
                    ! mass in aerosol modes
                    m_aer_b(i)=m_aer_b(i)+dummy3
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                endif
            enddo
    end subroutine move_aerosol_proportional
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Move aerosol from one category to another based on size                            !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester, 2023
	!>@brief calculates ice nucleation and mode1
	!>@param[in] n_mode, din_c, n_mix, sig_aer1, d_aer1, density_core1
	!>@param[inout] n_aer_a, s_aer_a, m_aer_a, n_aer_b, s_aer_b, m_aer_b
    subroutine move_aerosol_larger_than_size(n_mode, &
                    din_c,n_mix, sig_aer1, d_aer1, density_core1, &
                    n_aer_a, s_aer_a, m_aer_a, &
                    n_aer_b, s_aer_b, m_aer_b)
                    
        implicit none
        integer(i4b), intent(in) :: n_mode
        real(wp), intent(in) :: din_c, n_mix, sig_aer1, d_aer1, density_core1
        real(wp), dimension(n_mode-1), intent(inout) :: &
            n_aer_a, s_aer_a, m_aer_a, n_aer_b, s_aer_b, m_aer_b
        
        integer(i4b) :: i
        real(wp) :: n_tot, s_tot, m_tot, dummy1, dummy2, dummy3
        
        n_tot=sum(n_aer_a) ! total number in cw
        s_tot=sum(s_aer_a)
        m_tot=sum(m_aer_a)
        do i=1,n_mode-1
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! remove from cloud water:                                           !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! number in aerosol modes
            dummy1=ln_part_mom(0,din_c,n_mix, sig_aer1,d_aer1) * n_aer_a(i)/n_tot
            n_aer_a(i)=n_aer_a(i)-dummy1
            ! surface area in aerosol modes
            dummy2=pi*ln_part_mom(2,din_c,n_mix, sig_aer1,d_aer1) * s_aer_a(i)/s_tot
            s_aer_a(i)=s_aer_a(i)- dummy2 
            
            ! mass in aerosol modes
            dummy3=pi/6._wp* &
                density_core1*ln_part_mom(3,din_c,n_mix, sig_aer1,d_aer1) * &
                m_aer_a(i)/m_tot
            m_aer_a(i)=m_aer_a(i)- dummy3
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! add to aerosol particles in ice water                              !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! number in aerosol modes
            ! qv, n_mode aerosol + 1
            n_aer_b(i)=n_aer_b(i)+dummy1 
            
            ! surface area in aerosol modes
            s_aer_b(i)=s_aer_b(i)+dummy2 
            
            ! mass in aerosol modes
            m_aer_b(i)=m_aer_b(i)+dummy3
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        enddo
    end subroutine move_aerosol_larger_than_size
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    end module p_micro_module
    
