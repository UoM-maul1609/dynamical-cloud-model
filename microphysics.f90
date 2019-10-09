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
    use nrtype
    use nr, only : locate, polint
    use erfinv_m, only : erfinv

    use bam, only : n_mode, n_sv, giant_flag, method_flag, sv_flag, &
        	n_aer1, d_aer1, sig_aer1, molw_core1, density_core1, nu_core1, org_content1, &
        	molw_org1, density_org1, delta_h_vap1, nu_org1, log_c_star1, p_test, t_test, &
    		w_test, act_frac1, smax1, dcrit2, &
    		a_eq_7, b_eq_7, &
    		ctmm_activation, initialise_arrays, read_in_bam_namelist
    
    private
    public :: p_microphysics_3d, &
            p_microphysics_2d, p_microphysics_1d, read_in_pamm_bam_namelist, &
            p_initialise_aerosol_3d,p_initialise_aerosol, p_initialise_aerosol_1d
            
    ! Chen and Lamb (1994) Gamma variable fit (scaled and centred logarithm)
    integer(i4b), parameter :: n_cl=18
    real(sp), dimension(n_cl), parameter :: gam_cl=[-0.072328469664620_sp, &
        -0.324623262465577_sp, 0.363138099937540_sp, 3.323089908344732_sp, &
        0.874844989423720_sp, &
        -13.554426432462339_sp, -9.810322482346461_sp, 27.846739088352344_sp, &
        26.480447842355410_sp,&
         -29.890199206698309_sp, -32.327548996894521_sp, 15.827423311652167_sp, &
         18.466605783503052_sp, -4.158566361058538_sp, -5.039533848938808_sp, &
         1.477272813054374_sp, 1.038600921563425_sp, -0.457007828432810_sp]
    real(sp), dimension(2), parameter :: gam_mu_cl=[260.163817050062335_sp, &
                                                8.274747821396463_sp]
    
    ! physical constants
    real(sp), parameter :: rhow=1000._sp, rhoi=920._sp,lv=2.5e6_sp,ls=2.8e6_sp,lf=ls-lv, &
    					   cp=1005._sp, cw=4187._sp, cice=2093._sp, r=8.314_sp, &
    						mw=18e-3_sp, ma=29e-3_sp, ra=r/ma,rv=r/mw, eps1=ra/rv, &
    						ttr=273.15_sp, joules_in_an_erg=1.0e-7_sp, &
    						joules_in_a_cal=4.187_sp
    						
    						
    ! mass-diameter and size spectra relations
    real(sp), parameter :: cr=523.6_sp, cc=523.6_sp, &
                        cs=52.36_sp, cg=261.8_sp, ci=104._sp, &
    					dr=3_sp, dc=3_sp, ds=3._sp, dg=3._sp, di=3._sp, &
    					alpha_r=2.5_sp, alpha_c=0.0_sp, & ! note, alpha_c is for a "mass" - number distribution
    					alpha_s=2.5_sp, alpha_g=2.5_sp, alpha_i=0._sp
    					
	! terminal fall-speed relations
	real(sp), parameter :: a_r=362._sp, a_c=362._sp, &
	                        a_s=4.84_sp, a_g=253._sp, a_i=71.34_sp, &
							b_r=0.65_sp, b_c=0.65_sp, &
							b_s=0.25_sp, b_g=0.734_sp, b_i=0.6635_sp, &
							f_r=0._sp, f_c=0._sp, f_s=0._sp, f_g=0._sp, f_i=0._sp
							
	! autoconversion
	real(sp), parameter :: aw0=1e-3_sp, dwa=20e-6_sp, nl=2.4e8_sp, &
						lw0=rhow*pi/6._sp*nl*dwa**3, &
						tsaut=60._sp, dimax=0.3e-3_sp, di2s=0.33e-3_sp, &
						lambda_imin=(1._sp+di+alpha_i)/dimax, &
						tsbreak=60._sp, lambda_s_break=1000._sp
	
    ! microphysical values:
    real(sp), parameter :: hm_rate=3.5e8_sp, nar=1.1e15_sp, nbr=0._sp, &
    						rho0=1.2_sp, bbigg=100._sp, abigg=0.66_sp
    real(sp) :: mi0=1.e-14_sp

	! coalescence efficiencies
	real(sp), parameter :: erw=1._sp, erg=1._sp, ers=1._sp, eri=1._sp, esw=1._sp, &
						egw=1._sp, eiw=1._sp, egs_wet=1._sp, egi_wet=1._sp
	
	! variables used in various process rates:
	real(sp) :: gam1r,gam2r,gam1c, gam2c, gam1i,gam2i, gam1s, gam2s,gam1g,gam2g, &
				fall_q_r, fall_q_c, fall_q_s, fall_q_g, fall_n_r, fall_n_s, fall_n_g, &
				fall_q_i, fall_n_i, fall_n_c, &
				phi_r, mass_iacr,num_iacr, mass_sacw_i, mass_iacw, &
				mass_racs1,mass_racs2,mass_racs3, &
				mass_racg1,mass_racg2,mass_racg3, &
				mass_sacr1,mass_sacr2,mass_sacr3, &
				mass_sacg1,mass_sacg2,mass_sacg3, &
				mass_gacr1,mass_gacr2,mass_gacr3, &
				mass_gacs1,mass_gacs2,mass_gacs3, &
				num_racs1,num_racs2,num_racs3, num_racg1,num_racg2,num_racg3, &
				num_sacg1, num_sacg2,num_sacg3, &
				mass_gacw, mass_gaci, &
				nu_r1,nu_r2,nu_i1, nu_i2, nu_s1, nu_s2, nu_g1, nu_g2, &
				mass_imm, num_imm, q0sat, &
				chi_rain, chi_cloud, chi_ice, chi_snow, chi_graupel, &
				chi_rain1, chi_cloud1, chi_ice1, chi_snow1, chi_graupel1, &
				chi_num_ice, chi_num_ice1
				
	! Seifert and Beheng autoconversion
	real(sp) :: kc, kr, xstar
				
	real(sp), dimension(3) :: c=[1._sp,2._sp,1._sp]
	integer(i4b) :: k
	real(sp) :: isnow, iice, f1,f2,a,b, qsmall=1.e-30_sp
	
	integer :: n_modes_prof, n_levels_s
	real(sp), allocatable, dimension(:,:) :: n_read, sig_read, d_read
	real(sp), allocatable, dimension(:) :: z_read
	real(sp) :: small_number
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
        real(sp), intent(in) :: n,sig,d
        real(sp) :: ln_mom
        
        ln_mom=n*exp(real(k,sp)*log(d)+real(k,sp)**2*sig**2/2._sp)
    end function ln_mom
    
    

	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculate the partial kth integral moment of a lognormal distribution
	!>see:
	!>https://en.wikipedia.org/wiki/Log-normal_distribution#Partial_expectation
	!>@param[in] k: moment to calculate
	!>@param[in] n,sig,d: parameters of the lognormal distribution
	!>@param[out] mom: the integral moment to calculate
    function ln_part_mom(k,a,n,sig,d) 
        implicit none
        integer(i4b), intent(in) :: k
        real(sp), intent(in) :: a,n,sig,d
        real(sp) :: ln_part_mom
        
        real(sp) :: x, phi
        
        x=(log(a)-log(d)-sig**2*real(k,sp))/sig
        phi=0.5_sp*(1._sp+erf(x/sqrt(2._sp)))
        ln_part_mom= &
            n*(exp(log(d)*real(k,sp)+sig**2*real(k,sp)**2/2._sp)* (1._sp-phi))
        
    end function ln_part_mom
    
    
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculate the parameters of a lognormal distribution
	!>@param[in] n,s,m, rho
	!>@param[inout] sig_aer,d_aer
	subroutine ln_params_from_integral_moms(n,s,m,rho,sig_aer,d_aer)
	    implicit none
	    real(sp), intent(in) :: n,s,m, rho
	    real(sp), intent(inout) :: sig_aer, d_aer
	
	    if(m .gt. 0._sp) then
            ! this was derived by calculating moments of the distribution
            ! and solving to find dm and sig
            sig_aer=log( (36._sp*m**2*n*pi/rho**2)**(1._sp/3._sp) / s )
            sig_aer=sqrt(sig_aer)
        
            d_aer=(log(6._sp*m/(n*pi*rho))-4.5_sp*sig_aer**2) / 3._sp
            d_aer=exp(d_aer)	
        else
            sig_aer=0.3_sp
            d_aer=60.e-9_sp            
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
	    real(sp), intent(in) :: n
	    real(sp), dimension(n_mode-1), intent(in) :: s,m
	    real(sp), dimension(n_mode), intent(inout) :: rho,molw,nu
	    real(sp), intent(inout) :: n_aer,sig_aer, d_aer,n_mix,s_mix,m_mix

        
	    n_aer=n
	    n_mix=n
	    s_mix=sum(s)
	    m_mix=sum(m)
        if(m_mix .gt. 0._sp) then
            ! conserve volume of particle:	
            rho(n_mode) = sum(m) / sum(m/rho(1:n_mode-1))
            ! conserve total number of moles in particle:	
            molw(n_mode) = sum(m) / sum(m/molw(1:n_mode-1))
            ! conserve total number of moles of ions in particle:	
            nu(n_mode) = molw(n_mode) * sum(m*nu(1:n_mode-1)/molw(1:n_mode-1)) / sum(m) 

            ! this was derived by calculating moments of the distribution
            ! and solving to find dm and sig
            sig_aer=log( (36._sp*m_mix**2*n*pi/rho(n_mode)**2)**(1._sp/3._sp) / s_mix )
            if(sig_aer.le.0._sp) then
                sig_aer=0.3_sp
                d_aer=60.e-9_sp
                n_aer=0._sp
            else
                sig_aer=sqrt(sig_aer)
        
                d_aer=(log(6._sp*m_mix/(n*pi*rho(n_mode)))-4.5_sp*sig_aer**2) / 3._sp
                d_aer=exp(d_aer)	
            endif

        else
            rho(n_mode)=rho(n_mode-1)    
            molw(n_mode)=molw(n_mode-1)    
            nu(n_mode)=nu(n_mode-1)    
            
            sig_aer=0.3_sp
            d_aer=60.e-9_sp
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
        real(sp), dimension(-o_halo+1:kp+o_halo), intent(in) :: z
        real(sp), dimension(-o_halo+1:kp+o_halo), intent(in) :: &
    					rho, p, t
        real(sp), dimension(-o_halo+1:kp+o_halo,nq), &
            intent(inout) :: q

        ! local variables
        integer(i4b) :: i, k, AllocateStatus, iloc
        real(sp) :: w, smax, phi, xx, kmom, var, dummy
        real(sp), dimension(:), allocatable :: act_frac1 , dcrit
         
        
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
                iloc=locate(z_read(1:n_levels_s),z(k))
                iloc=min(n_levels_s-1,iloc)
                iloc=max(1,iloc)
                do i=1,n_mode
                    ! linear interp n_aer
                    call polint(z_read(iloc:iloc+1), n_read(i,iloc:iloc+1), &
                                min(z(k),z_read(n_levels_s)), var,dummy)
                    n_aer1(i)=var
                    ! linear interp sig_aer
                    call polint(z_read(iloc:iloc+1), sig_read(i,iloc:iloc+1), &
                                min(z(k),z_read(n_levels_s)), var,dummy)
                    sig_aer1(i)=var
                    ! linear interp d_aer
                    call polint(z_read(iloc:iloc+1), d_read(i,iloc:iloc+1), &
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
                q(k,(i-1)*3+4)= pi/6._sp*density_core1(i)* &
                    ln_mom(3,n_aer1(i),sig_aer1(i),d_aer1(i))

            enddo
            
            
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! 1. find the critical diameter of each aerosol mode, and                    !
            ! 2. perform integration to set the aerosol n,s,m, in cloud water            !                         
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! initialise aerosol in cloud water    
            !  
            if(q(k,(n_mode-1)*6+3) .gt. 0._sp) then 

                
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
                    q(k,(n_mode-1)*6+7+(i-1)*3)= pi/6._sp*density_core1(i)* &
                        ln_part_mom(3,dcrit(i),n_aer1(i),sig_aer1(i),d_aer1(i))
                    q(k,(i-1)*3+4)=q(k,(i-1)*3+4)-q(k,(n_mode-1)*6+7+(i-1)*3)
                enddo
                
                
            else
                smax=0._sp
                dcrit=1000._sp
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
        real(sp), dimension(-o_halo+1:ip+o_halo), intent(in) :: x
        real(sp), dimension(-o_halo+1:kp+o_halo), intent(in) :: z
        real(sp), dimension(-o_halo+1:kp+o_halo,-o_halo+1:ip+o_halo), intent(in) :: &
    					rho, p, t
        real(sp), dimension(-o_halo+1:kp+o_halo,-o_halo+1:ip+o_halo,nq), &
            intent(inout) :: q, q_old

        ! local variables
        integer(i4b) :: i, k, AllocateStatus, iloc
        real(sp) :: w, smax, phi, xx, kmom, var, dummy
        real(sp), dimension(:), allocatable :: act_frac1 , dcrit
         
        
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
                iloc=locate(z_read(1:n_levels_s),z(k))
                iloc=min(n_levels_s-1,iloc)
                iloc=max(1,iloc)
                do i=1,n_mode
                    ! linear interp n_aer
                    call polint(z_read(iloc:iloc+1), n_read(i,iloc:iloc+1), &
                                min(z(k),z_read(n_levels_s)), var,dummy)
                    n_aer1(i)=var
                    ! linear interp sig_aer
                    call polint(z_read(iloc:iloc+1), sig_read(i,iloc:iloc+1), &
                                min(z(k),z_read(n_levels_s)), var,dummy)
                    sig_aer1(i)=var
                    ! linear interp d_aer
                    call polint(z_read(iloc:iloc+1), d_read(i,iloc:iloc+1), &
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
                q(k,:,(i-1)*3+4)= pi/6._sp*density_core1(i)* &
                    ln_mom(3,n_aer1(i),sig_aer1(i),d_aer1(i))

            enddo
            
            
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! 1. find the critical diameter of each aerosol mode, and                    !
            ! 2. perform integration to set the aerosol n,s,m, in cloud water            !                         
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! initialise aerosol in cloud water    
            !  
            if(q(k,1,(n_mode-1)*6+3) .gt. 0._sp) then 

                
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
                    q(k,:,(n_mode-1)*6+7+(i-1)*3)= pi/6._sp*density_core1(i)* &
                        ln_part_mom(3,dcrit(i),n_aer1(i),sig_aer1(i),d_aer1(i))
                    q(k,:,(i-1)*3+4)=q(k,:,(i-1)*3+4)-q(k,:,(n_mode-1)*6+7+(i-1)*3)
                enddo
                
                
            else
                smax=0._sp
                dcrit=1000._sp
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
        real(sp), dimension(-o_halo+1:ip+o_halo), intent(in) :: x
        real(sp), dimension(-o_halo+1:jp+o_halo), intent(in) :: y
        real(sp), dimension(-o_halo+1:kp+o_halo), intent(in) :: z
        real(sp), dimension(-o_halo+1:kp+o_halo), intent(in) :: rho, p, t
        real(sp), dimension(-o_halo+1:kp+o_halo,-o_halo+1:jp+o_halo,-o_halo+1:ip+o_halo,nq), &
            intent(inout) :: q

        ! local variables
        integer(i4b) :: i, k, AllocateStatus, iloc
        real(sp) :: w, smax, phi, xx, kmom, var, dummy
        real(sp), dimension(:), allocatable :: act_frac1 , dcrit
         
        
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
                iloc=locate(z_read(1:n_levels_s),z(k))
                iloc=min(n_levels_s-1,iloc)
                iloc=max(1,iloc)
                do i=1,n_mode
                    ! linear interp n_aer
                    call polint(z_read(iloc:iloc+1), n_read(i,iloc:iloc+1), &
                                min(z(k),z_read(n_levels_s)), var,dummy)
                    n_aer1(i)=var
                    ! linear interp sig_aer
                    call polint(z_read(iloc:iloc+1), sig_read(i,iloc:iloc+1), &
                                min(z(k),z_read(n_levels_s)), var,dummy)
                    sig_aer1(i)=var
                    ! linear interp d_aer
                    call polint(z_read(iloc:iloc+1), d_read(i,iloc:iloc+1), &
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
                q(k,:,:,(i-1)*3+4)= pi/6._sp*density_core1(i)* &
                    ln_mom(3,n_aer1(i),sig_aer1(i),d_aer1(i))

            enddo
            
            
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! 1. find the critical diameter of each aerosol mode, and                    !
            ! 2. perform integration to set the aerosol n,s,m, in cloud water            !                         
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! initialise aerosol in cloud water    
            !  
            if(q(k,1,1,(n_mode-1)*6+3) .gt. 0._sp) then 

                
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
                    q(k,:,:,(n_mode-1)*6+7+(i-1)*3)= pi/6._sp*density_core1(i)* &
                        ln_part_mom(3,dcrit(i),n_aer1(i),sig_aer1(i),d_aer1(i))
                    q(k,:,:,(i-1)*3+4)=q(k,:,:,(i-1)*3+4)-q(k,:,:,(n_mode-1)*6+7+(i-1)*3)
                enddo
                
                
            else
                smax=0._sp
                dcrit=1000._sp
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
	gam1r=gamma(1._sp+alpha_r)
	gam2r=gamma(1._sp+alpha_r+dr)
	gam1c=gamma(1._sp+alpha_c)
	gam2c=gamma(1._sp+alpha_c+1._sp) ! note the 1, instead of dc - drop distribution
	                                ! is a mass distribution
	gam1i=gamma(1._sp+alpha_i)
	gam2i=gamma(1._sp+alpha_i+di)
	gam1s=gamma(1._sp+alpha_s)
	gam2s=gamma(1._sp+alpha_s+ds)
	gam1g=gamma(1._sp+alpha_g)
	gam2g=gamma(1._sp+alpha_g+dg)

    ! mass weighted fall for r, c, s, g, i
    fall_q_r=a_r*gamma(1._sp+alpha_r+dr+b_r) / gamma(1._sp+alpha_r+dr)
    fall_q_c=a_c*gamma(1._sp+alpha_c+1._sp+b_c) / gamma(1._sp+alpha_c+1._sp)
    fall_q_s=a_s*gamma(1._sp+alpha_s+ds+b_s) / gamma(1._sp+alpha_s+ds)
    fall_q_g=a_g*gamma(1._sp+alpha_g+dg+b_g) / gamma(1._sp+alpha_g+dg)
    fall_q_i=a_i*gamma(1._sp+alpha_i+di+b_i) / gamma(1._sp+alpha_i+di)

    ! number weighted fall for r, c, s, g
    fall_n_r=a_r*gamma(1._sp+alpha_r+b_r) / gamma(1._sp+alpha_r)
    fall_n_c=a_c*gamma(1._sp+alpha_c+b_c) / gamma(1._sp+alpha_c)
    fall_n_s=a_s*gamma(1._sp+alpha_s+b_s) / gamma(1._sp+alpha_s)
    fall_n_g=a_g*gamma(1._sp+alpha_g+b_g) / gamma(1._sp+alpha_g)
    fall_n_i=a_i*gamma(1._sp+alpha_i+b_i) / gamma(1._sp+alpha_i)
    
    ! sweep out of rain
    phi_r=pi*a_r*gamma(3._sp+b_r+alpha_r) / 4._sp
    
    ! ice accreting rain
    mass_iacr=pi*eri*a_r*cr*gamma(3._sp+b_r+dr+alpha_r)/4._sp
    num_iacr =pi*eri*a_r*gamma(3._sp+b_r+alpha_r)/4._sp

	! collection of cloud by snow and ice
	mass_sacw_i=pi*a_s*gamma(3._sp+b_s+alpha_s)/4._sp
	mass_iacw=pi*a_i*gamma(3._sp+b_i+alpha_i)/4._sp
	
	! collisions between precipitating particles of different species
	! rain-snow
	mass_racs1=gamma(1._sp+alpha_r)*gamma(3._sp+alpha_s+ds)
	mass_racs2=2._sp*gamma(2._sp+alpha_r)*gamma(2._sp+alpha_s+ds)
	mass_racs3=gamma(3._sp+alpha_r)*gamma(1._sp+alpha_s+ds)
	
	num_racs1=gamma(1._sp+alpha_r)*gamma(3._sp+alpha_s)
	num_racs2=2._sp*gamma(2._sp+alpha_r)*gamma(2._sp+alpha_s)
	num_racs3=gamma(3._sp+alpha_r)*gamma(1._sp+alpha_s)
    ! rain-graupel
	mass_racg1=gamma(1._sp+alpha_r)*gamma(3._sp+alpha_g+dg)
	mass_racg2=2._sp*gamma(2._sp+alpha_r)*gamma(2._sp+alpha_g+dg)
	mass_racg3=gamma(3._sp+alpha_r)*gamma(1._sp+alpha_g+dg)

	num_racg1=gamma(1._sp+alpha_r)*gamma(3._sp+alpha_g)
	num_racg2=2._sp*gamma(2._sp+alpha_r)*gamma(2._sp+alpha_g)
	num_racg3=gamma(3._sp+alpha_r)*gamma(1._sp+alpha_g)
    ! snow-rain
	mass_sacr1=gamma(1._sp+alpha_s)*gamma(3._sp+alpha_r+dr)
	mass_sacr2=2._sp*gamma(2._sp+alpha_s)*gamma(2._sp+alpha_r+dr)
	mass_sacr3=gamma(3._sp+alpha_s)*gamma(1._sp+alpha_r+dr)
    ! snow-graupel
	mass_sacg1=gamma(1._sp+alpha_s)*gamma(3._sp+alpha_g+dg)
	mass_sacg2=2._sp*gamma(2._sp+alpha_s)*gamma(2._sp+alpha_g+dg)
	mass_sacg3=gamma(3._sp+alpha_s)*gamma(1._sp+alpha_g+dg)

	num_sacg1=gamma(1._sp+alpha_s)*gamma(3._sp+alpha_g)
	num_sacg2=2._sp*gamma(2._sp+alpha_s)*gamma(2._sp+alpha_g)
	num_sacg3=gamma(3._sp+alpha_s)*gamma(1._sp+alpha_g)
    ! graupel-rain
	mass_gacr1=gamma(1._sp+alpha_g)*gamma(3._sp+alpha_r+dr)
	mass_gacr2=2._sp*gamma(2._sp+alpha_g)*gamma(2._sp+alpha_r+dr)
	mass_gacr3=gamma(3._sp+alpha_g)*gamma(1._sp+alpha_r+dr)
    ! graupel-snow
	mass_gacs1=gamma(1._sp+alpha_g)*gamma(3._sp+alpha_s+ds)
	mass_gacs2=2._sp*gamma(2._sp+alpha_g)*gamma(2._sp+alpha_s+ds)
	mass_gacs3=gamma(3._sp+alpha_g)*gamma(1._sp+alpha_s+ds)
	
	
	! accretion and riming by graupel
	mass_gacw=pi*egw*a_g*gamma(3._sp+b_g+alpha_g)/4._sp
	mass_gaci=pi*a_g*gamma(3._sp+b_g+alpha_g)/4._sp
    
    ! gauss hypergeometric equations aggregation of ice with ice (and snow with snow)
    ! See Ferrier (1994, JAS part 1, equation B.21)
    ! snow:
    a=1._sp
    b=4._sp+2._sp*alpha_s+b_s
    isnow=0._sp
    do k=1,3
	    call hygfx(a, b, real(k,sp)+alpha_s+1.0_sp,0.5_sp,f1)
	    call hygfx(a, b, real(k,sp)+alpha_s+b_s+1.0_sp, 0.5_sp,f2)
	    isnow=isnow+c(k)*(f1/(real(k,sp)+alpha_s)-f2/(real(k,sp)+alpha_s+b_s))
	enddo
	isnow=a_s*pi*gamma(b)/(2._sp**(6._sp+2._sp*alpha_s+b_s)) * isnow
	
    ! ice:
    a=1._sp
    b=4._sp+2._sp*alpha_i+b_i
    iice=0._sp
    do k=1,3
	    call hygfx(a, b, real(k,sp)+alpha_i+1.0_sp,0.5_sp,f1)
	    call hygfx(a, b, real(k,sp)+alpha_i+b_i+1.0_sp, 0.5_sp,f2)
	    iice=iice+c(k)*(f1/(real(k,sp)+alpha_i)-f2/(real(k,sp)+alpha_i+b_i))
	enddo
	iice=a_i*pi*gamma(b)/(2._sp**(6._sp+2._sp*alpha_i+b_i)) * iice
	
	
	! ventilation
	! rain:
	nu_r1=0.78_sp*gamma(2._sp+alpha_r)
	nu_r2=gamma(0.5_sp*b_r+alpha_r+2.5)
	! ice:
	nu_i1=0.78_sp*gamma(2._sp+alpha_i)
	nu_i2=0.31_sp*gamma(0.5_sp*b_i+alpha_i+2.5)
	! snow:
	nu_s1=0.78_sp*gamma(2._sp+alpha_s)
	nu_s2=0.31_sp*gamma(0.5_sp*b_s+alpha_s+2.5)
	! graupel:
	nu_g1=0.78_sp*gamma(2._sp+alpha_g)
	nu_g2=0.31_sp*gamma(0.5_sp*b_g+alpha_g+2.5)
	
	! immersion freezing by bigg
	mass_imm=gamma(4._sp+dr+alpha_r)*pi*cr*bbigg/6._sp
	num_imm=gamma(4._sp+alpha_r)*pi*bbigg/6._sp
	
	
	! precipitation
	chi_rain=gamma(1._sp+alpha_r+b_r+dr)
	chi_cloud=gamma(1._sp+alpha_c+b_c+1._sp)
	chi_ice=gamma(1._sp+alpha_i+b_i+di)
	chi_num_ice=gamma(1._sp+alpha_i+b_i)
	chi_snow=gamma(1._sp+alpha_s+b_s+ds)
	chi_graupel=gamma(1._sp+alpha_g+b_g+dg)
	
	chi_rain1=gamma(1._sp+alpha_r+dr)
	chi_cloud1=gamma(1._sp+alpha_c+1._sp)
	chi_ice1=gamma(1._sp+alpha_i+di)
	chi_num_ice1=gamma(1._sp+alpha_i)
	chi_snow1=gamma(1._sp+alpha_s+ds)
	chi_graupel1=gamma(1._sp+alpha_g+dg)
	

    ! Seifert and Beheng autoconversion:
    kc=9.44e9_sp ! m3 kg-2 s-1
    kr=5.78e0_sp ! m3 kg-2 s-1
    xstar=2.6e-10_sp ! kg
    end subroutine initialise_microphysics_vars
    
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
	!>@param[inout] precip: precip in rain, snow, graupel, ice cats - diagnostic
	!>@param[inout] theta: theta 
	!>@param[inout] p: pressure
	!>@param[in] z: vertical levels 
	!>@param[inout] thetan: potential temperature 
	!>@param[inout] rho, rhon: density 
	!>@param[in] w: vertical wind 
	!>@param[inout] micro_init: boolean to initialise microphysics 
	!>@param[in] hm_flag: switch hm-process on and off
	!>@param[in] mass_ice: mass of a single ice crystal (override)
	!>@param[in] ice_flag: ice microphysics
	!>@param[in] theta_flag: whether to alter theta
    subroutine p_microphysics_3d(nq,ncat,n_mode,cst,cen,inc,iqc, inr,iqr,ini,iqi,iai, &
                    cat_am,cat_c, cat_r, cat_i,&
                    nprec, &
                    ip,jp,kp,l_h,r_h,dt,dz,dzn,q,precip,th,prefn, z,thetan,rhoa,rhoan,w, &
    				micro_init,hm_flag, mass_ice, ice_flag, theta_flag)
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
	!>@param[inout] precip: precip in rain, snow, graupel, ice cats - diagnostic
	!>@param[inout] theta: theta 
	!>@param[inout] p: pressure
	!>@param[in] z: vertical levels 
	!>@param[inout] thetan: potential temperature 
	!>@param[inout] rho, rhon: density 
	!>@param[in] w: vertical wind 
	!>@param[inout] micro_init: boolean to initialise microphysics 
	!>@param[in] hm_flag: switch hm-process on and off
	!>@param[in] mass_ice: mass of a single ice crystal (override)
	!>@param[in] ice_flag: ice microphysics
	!>@param[in] theta_flag: whether to alter theta
	!>@param[in] comm,comm_vert,id,dims,coords: MPI variables
    subroutine p_microphysics_3d(nq,ncat,n_mode,cst,cen,inc,iqc, inr,iqr,ini,iqi,iai, &
                    cat_am,cat_c, cat_r, cat_i,&
                    nprec, &
                    ip,jp,kp,l_h,r_h,dt,dz,dzn,q,precip,th,prefn, z,thetan,rhoa,rhoan,w, &
    				micro_init,hm_flag, mass_ice, ice_flag, theta_flag, &
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
    real(sp), intent(in) :: dt
    real(sp), dimension(-l_h+1:kp+r_h,-l_h+1:jp+r_h,-l_h+1:ip+r_h,nq), intent(inout) :: q
    real(sp), dimension(1:kp,1-l_h:jp+r_h,1-l_h:ip+r_h,nprec), intent(inout) :: precip
    real(sp), dimension(-l_h+1:kp+r_h,-l_h+1:jp+r_h,-l_h+1:ip+r_h), intent(inout) :: &
    					th
    real(sp), dimension(-l_h+1:kp+r_h), intent(in) :: z, dz, dzn, rhoa,rhoan, thetan, &
        prefn
    real(sp), dimension(-l_h+1:kp+r_h,-l_h+1:jp+r_h,-l_h+1:ip+r_h), intent(in) :: w
    logical, intent(in) :: ice_flag, hm_flag, theta_flag
    logical , intent(inout) :: micro_init
    real(sp), intent(in) :: mass_ice

	! locals
	integer(i4b) :: i,j,n, error,n1
#if MPI_PAMM == 1
	real(sp), dimension(-l_h+1:kp+r_h,-l_h+1:jp+r_h,-l_h+1:ip+r_h) :: & 
	                vqr, vqc, vqi
	integer(i4b), dimension(3) :: n_step,n_step_o,n_step_g
    integer(i4b), intent(in) :: id, comm,comm_vert
    integer(i4b), dimension(3), intent(in) :: coords, dims
    real(sp), dimension(nq) :: lbc,ubc
    logical, dimension(3) :: adv_lg, adv_l=[.false.,.false.,.false.], &
                    adv_l_o=[.false.,.false.,.false.]

    n_step=1
    n_step_o=1
    lbc=0._sp
    ubc=0._sp
#endif
	
	do i=1,ip
	    do j=1,jp
#if MPI_PAMM == 0 
    		call p_microphysics_1d(nq,ncat,n_mode,cst,cen,inc,iqc, inr,iqr, ini,iqi,iai,&
		                cat_am,cat_c, cat_r, cat_i,nprec,&
		                kp,l_h,dt,dz,dzn,q(:,j,i,:),precip(:,j,i,:),th(:,j,i),&
		                    prefn, &
							z(:),thetan,rhoa(:),rhoan(:),w(:,j,i), &
    						micro_init,hm_flag, mass_ice, ice_flag, theta_flag)
#else

    		call p_microphysics_1d(nq,ncat,n_mode,cst,cen,inc,iqc,inr,iqr, ini,iqi,iai, &
		                cat_am,cat_c, cat_r, cat_i,nprec, &
		                kp,l_h,dt,dz,dzn,q(:,j,i,:),precip(:,j,i,:),th(:,j,i),&
		                    prefn, &
							z(:),thetan,rhoa(:),rhoan(:),w(:,j,i), &
							vqc(:,j,i),vqr(:,j,i),vqi(:,j,i),n_step, adv_l, &
							coords, &
    						micro_init,hm_flag, mass_ice, ice_flag, theta_flag)
    		n_step_o=max(n_step,n_step_o)

    		adv_l_o=adv_l_o .or. adv_l ! if there has been a true at any point, 
    		                           ! set adv_l_o to true on this PE
#endif
    	enddo
	enddo
	
	! collective communication
#if MPI_PAMM == 1
	call mpi_allreduce(adv_l_o(1:3),adv_lg(1:3),3,MPI_LOGICAL,MPI_LOR, comm_vert,error)
	call mpi_allreduce(n_step_o(1:3),n_step_g(1:3),3,MPI_INTEGER,MPI_MAX, comm_vert,error)
#endif
    
    
    
    
    


#if MPI_PAMM == 1
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! full exchange needed                                                               !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call exchange_full(comm, id, kp, jp, ip, r_h,r_h,r_h,r_h,r_h,r_h, &
                            th,0._sp,0._sp,dims,coords)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! full exchange needed                                                               !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do n1=1,nq
        call exchange_full(comm, id, kp, jp, ip, r_h,r_h,r_h,r_h,r_h,r_h, &
                                q(:,:,:,n1),lbc(n1),ubc(n1),dims,coords)
    enddo
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    if(adv_lg(1)) then
        call exchange_along_z(comm_vert, id, kp, jp, ip, r_h,r_h,r_h,r_h,r_h,r_h, &
                                vqc(:,:,:),0._sp,0._sp,dims,coords)
        do n=1,n_step_g(1)
            call mpdata_vec_vert_3d(dt/real(n_step_g(1),sp),dz,dzn,&
                    rhoa,rhoan, &
                    ip,jp,kp,cen(cat_c)-cst(cat_c)+1,l_h,r_h,&
                    vqc,q(:,:,:,cst(cat_c):cen(cat_c)),&
                    lbc(cst(cat_c):cen(cat_c)),ubc(cst(cat_c):cen(cat_c)), &
                    1,.false., 2,comm_vert, id, &
                    dims,coords)
        enddo
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! full exchange needed                                                           !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do n1=cst(cat_c),cen(cat_c)
            call exchange_full(comm, id, kp, jp, ip, r_h,r_h,r_h,r_h,r_h,r_h, &
                                    q(:,:,:,n1),lbc(n1),ubc(n1),dims,coords)
        enddo
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    endif       
    if(adv_lg(2)) then
        call exchange_along_z(comm_vert, id, kp, jp, ip, r_h,r_h,r_h,r_h,r_h,r_h, &
                                vqr(:,:,:),0._sp,0._sp,dims,coords)
        do n=1,n_step_g(2)

            call mpdata_vec_vert_3d(dt/real(n_step_g(2),sp),dz,dzn,&
                    rhoa,rhoan, &
                    ip,jp,kp,cen(cat_r)-cst(cat_r)+1,l_h,r_h,&
                    vqr,q(:,:,:,cst(cat_r):cen(cat_r)),&
                    lbc(cst(cat_r):cen(cat_r)),ubc(cst(cat_r):cen(cat_r)), &
                    1,.false., 2,comm_vert, id, &
                    dims,coords) 

        enddo
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! full exchange needed                                                           !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do n1=cst(cat_r),cen(cat_r)
            call exchange_full(comm, id, kp, jp, ip, r_h,r_h,r_h,r_h,r_h,r_h, &
                                    q(:,:,:,n1),lbc(n1),ubc(n1),dims,coords)
        enddo
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    endif       
    if(adv_lg(3)) then
        call exchange_along_z(comm_vert, id, kp, jp, ip, r_h,r_h,r_h,r_h,r_h,r_h, &
                                vqi(:,:,:),0._sp,0._sp,dims,coords)
        do n=1,n_step_g(3)

            call mpdata_vec_vert_3d(dt/real(n_step_g(3),sp),dz,dzn,&
                    rhoa,rhoan, &
                    ip,jp,kp,cen(cat_i)-cst(cat_i)+1,l_h,r_h,&
                    vqi,q(:,:,:,cst(cat_i):cen(cat_i)),&
                    lbc(cst(cat_i):cen(cat_i)),ubc(cst(cat_i):cen(cat_i)), &
                    1,.false., 2,comm_vert, id, &
                    dims,coords) 

        enddo
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
	!>@param[in] theta_flag: whether to alter theta
    subroutine p_microphysics_2d(nq,ncat,n_mode,cst,cen,inc,iqc,inr,iqr,ini,iqi,iai, &
                    cat_am,cat_c, cat_r, cat_i,nprec, &
                    ip,kp,o_halo,dt,dz,dzn,q,precip,theta,p, z,theta_ref,rho,rhon,w, &
    						micro_init,hm_flag, mass_ice, theta_flag)
    implicit none
    ! arguments:
    integer(i4b), intent(in) :: nq, ncat, n_mode, ip,kp, o_halo, inc, iqc, inr,iqr, &
        ini,iqi,iai, &
        cat_am,cat_c, cat_r, cat_i,nprec
    integer(i4b), dimension(ncat), intent(in) :: cst,cen
    real(sp), intent(in) :: dt
    real(sp), dimension(-o_halo+1:kp+o_halo,-o_halo+1:ip+o_halo,nq), intent(inout) :: q
    real(sp), dimension(1:kp,1:ip,1:nprec), intent(inout) :: precip
    real(sp), dimension(-o_halo+1:kp+o_halo,-o_halo+1:ip+o_halo), intent(inout) :: &
    					theta, p,rho
    real(sp), dimension(-o_halo+1:kp+o_halo), intent(in) :: z, dz, dzn, &
                    rhon, theta_ref
    real(sp), dimension(-o_halo+1:kp+o_halo,-o_halo+1:ip+o_halo), intent(in) :: w
    logical, intent(in) :: hm_flag, theta_flag
    logical , intent(inout) :: micro_init
    real(sp), intent(in) :: mass_ice


	! locals
	integer(i4b) :: i
#if MPI_PAMM == 1
    real(sp), dimension(-o_halo+1:kp+o_halo,-o_halo+1:ip+o_halo) :: vqc,vqr, vqi
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
    						micro_init,hm_flag, mass_ice, .false., theta_flag)
#else
		call p_microphysics_1d(nq,ncat,n_mode,cst,cen,inc,iqc,inr,iqr, ini,iqi,iai,&
		                cat_am,cat_c, cat_r, cat_i,nprec, &
		                kp,o_halo,dt,dz,dzn,q(:,i,:),precip(:,i,:),theta(:,i),p(:,i), &
							z(:),theta_ref,rho(:,i),rhon(:),w(:,i), &
							vqc(:,i),vqr(:,i),vqi(:,i), n_step, adv_l, coords,&
    						micro_init,hm_flag, mass_ice, .false., theta_flag)
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
	!>@param[in] theta_flag: whether to alter theta
    subroutine p_microphysics_1d(nq,ncat,n_mode,cst,cen, inc, iqc, inr,iqr, ini,iqi,iai,&
                            cat_am,cat_c, cat_r, cat_i,nprec, &
                            kp,o_halo,dt,dz,dzn,q,precip,th,p, z,theta,rhoa,rhon,u, &
    						micro_init,hm_flag, mass_ice,ice_flag, theta_flag)
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
	!>@param[in] theta_flag: whether to alter theta
    subroutine p_microphysics_1d(nq,ncat,n_mode,cst,cen, inc, iqc, inr,iqr,ini,iqi,iai,&
                            cat_am,cat_c, cat_r, cat_i,nprec,&
                            kp,o_halo,dt,dz,dzn,q,precip,th,p, z,theta,rhoa,rhon,u, &
                            vqc,vqr,vqi,n_step, adv_l, coords,&
    						micro_init,hm_flag, mass_ice,ice_flag, theta_flag)
#endif

	use advection_1d
	use nr, only : dfridr
	use advection_s_1d, only : mpdata_vec_1d
    implicit none
    ! arguments:
    integer(i4b), intent(in) :: nq, ncat,n_mode, kp, o_halo, inc, iqc,inr,iqr,&
                             ini,iqi,iai, &
                            cat_am,cat_c, cat_r, cat_i,nprec
    integer(i4b), dimension(ncat), intent(in) :: cst,cen
    real(sp), intent(in) :: dt
    real(sp), dimension(-o_halo+1:kp+o_halo,nq), intent(inout) :: q
    real(sp), dimension(1:kp,nprec), intent(inout) :: precip
    real(sp), dimension(-o_halo+1:kp+o_halo), intent(inout) :: th
    real(sp), dimension(-o_halo+1:kp+o_halo), intent(in) :: dz, z, dzn, rhoa, &
                                                    rhon, theta,p
    real(sp), dimension(-o_halo+1:kp+o_halo), intent(in) :: u
    logical, intent(in) :: hm_flag, ice_flag, theta_flag
    logical , intent(inout) :: micro_init
    real(sp), intent(in) :: mass_ice
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
    real(sp) :: temp, qtot,qaut, a, b, ab_ice, ab_liq, ice_dep,snow_dep,graup_dep, &
    			nu_ice, nu_snow, nu_graup, diff1, ktherm1, tc, nu_vis, sc, nu_rain, rain_evap, &
    			sb_aut, sb_acr, sb_cwaut, sb_cwacr, sb_raut, sb_rsel, sb_cwsel
    
    real(sp), dimension(-o_halo+1:kp+o_halo) :: rho
    real(sp), dimension(kp) :: smr, smr_i
    
    real(sp), dimension(kp) :: &
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
				picnt, & ! nucleation of ice crystals by contact freezing
				pidep, & ! deposition of water vapour onto cloud ice
				piprm, & ! primary nucleation of ice crystals by INPs
				pifrw, & ! nucleation of ice crystals by homogeneous freezing of cloud
				pihal, & ! production of ice crystals by hm process
				pimlt, & ! cloud ice melting to form rain
				pisub, & ! sublimation of cloud ice
				praci_g, & ! accretion of cloud ice by rain to form graupel
				praci_s, & ! accretion of cloud ice by rain to form snow
				pracs, & ! accretion of snow by rain to form graupel
				pracw, & ! accretion of liquid cloud by rain
				praut, & ! autoconversion from liquid cloud to rain (coalescence)
				prevp, & ! evaporation of rain
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
    			rcwsel     ! cloud water self accretion - number
    				    
    real(sp) :: pgwet ! amount of liquid that graupel can freeze without shedding
    								

    real(sp), dimension(1-o_halo:kp+o_halo) :: n_r, lam_r, n_i, lam_i, n_s, &
                                                 lam_s, n_g, lam_g, lam_c, n_c
    real(sp), dimension(1-o_halo:kp+o_halo) :: rho_fac
	real(sp), dimension(1-o_halo:kp+o_halo) :: vnr, vnc
#if MPI_PAMM == 0
	real(sp), dimension(1-o_halo:kp+o_halo) :: vqr, vqs, vqg, vqi, vns, vng, vni, &
	                                        vqc
#else
	real(sp), intent(inout), dimension(1-o_halo:kp+o_halo) :: vqr, vqc, vqi
#endif
	real(sp), dimension(1-o_halo:kp+o_halo) :: t
	! coalescence efficiencies
	real(sp), dimension(kp) :: egi_dry, egs_dry, esi, eii, ess
	real(sp) :: qold,des_dt,dqs_dt,err,cond,temp1, dummy1,dummy2, dummy3,&
	            n_mix,s_mix,m_mix, nin_c, din_c,nin_r,din_r, n_tot, s_tot, m_tot
	
	real(sp), dimension(1-o_halo:kp+o_halo) :: gamma_t,dep_density
	real(sp) :: phi,vol
	
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
    vqr=0._sp
    vqc=0._sp
	pgaci=0._sp
	pgacr=0._sp
	pgacs=0._sp
	rgacs=0._sp
	pgacw=0._sp
	pgaut=0._sp
	pgdep=0._sp
	pgfr=0._sp
	rgfr=0._sp
	pgmlt=0._sp
	pgshd=0._sp
	pgsub=0._sp
	riaci=0._sp
	piacr_g=0._sp
	riacr_g=0._sp
	piacr_s=0._sp
	riacr_s=0._sp
	piacw=0._sp
	picnt=0._sp
	pidep=0._sp
	piprm=0._sp
	pifrw=0._sp
	pihal=0._sp
	pimlt=0._sp
	pisub=0._sp
	praci_g=0._sp
	praci_s=0._sp
	pracs=0._sp
	pracw=0._sp
	praut=0._sp
	prevp=0._sp
	psacr=0._sp
	rsacr=0._sp
	psaci=0._sp
	rsacs=0._sp
	psacw=0._sp
	psaut=0._sp
	rsaut=0._sp
	rsbrk=0._sp
	psdep=0._sp
	psmlt=0._sp
	pssub=0._sp
	rcwaut=0._sp
    rcwacr=0._sp
    rraut=0._sp
    rrsel=0._sp
    rcwsel=0._sp
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! some commonly used variables that depend on prognostics                            !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    t=(theta+th)*(p/1.e5_sp)**(ra/cp) ! temperature
    rho=rhoa !p / (ra*t) ! air density    
    rho_fac=(rho0/rho(:))**0.5_sp
    ! rain n0, lambda
    lam_r=(max(q(:,cst(cat_r)),1._sp)*cr*gam2r / &
            (max(q(:,cst(cat_r)+1),1.e-10_sp)*gam1r))**(1._sp/dr)
    n_r=rho(:)*max(q(:,cst(cat_r)),0._sp)*lam_r**(1._sp+alpha_r) / gam1r
    ! cloud n0, lambda    
    lam_c=(max(q(:,inc),1._sp)*cc*gam2c / (max(q(:,iqc),1.e-10_sp)*gam1c))**(1._sp/1._sp)
    n_c=rho(:)*max(q(:,inc),0._sp)*lam_c**(1._sp+alpha_c) / gam1c
    
    if(ice_flag) then
        ! ice n0, lambda
        lam_i=(max(q(:,ini),1._sp)*ci*gam2i / (max(q(:,iqi),1.e-10_sp)*gam1i))**(1._sp/di)
        n_i=rho(:)*max(q(:,ini),0._sp)*lam_i**(1._sp+alpha_i) / gam1i

        ! ice
        vqi(:)=min(max(fall_q_i*rho_fac * lam_i**(1._sp+alpha_i+di) / &
            (lam_i+f_i)**(1._sp+alpha_i+di+b_i), 0._sp), 10._sp)

        ! precipitation
        precip(1:kp,2)=n_i(1:kp)*(a_i*chi_num_ice/(lam_i(1:kp)**(alpha_i+b_i+1._sp)) - &
                        u(1:kp)*chi_num_ice1/(lam_i(1:kp)**(alpha_i+1._sp))) &
                        /rho(1:kp)
    endif
    
    ! precipitation
	precip(1:kp,1)=cr*n_r(1:kp)*(a_r*chi_rain/(lam_r(1:kp)**(alpha_r+b_r+dr+1._sp)) - &
					u(1:kp)*chi_rain1/(lam_r(1:kp)**(alpha_r+dr+1._sp))) &
					/rho(1:kp) *3600._sp
    
    ! fall speeds
    ! rain
    vqr(:)=min(max(fall_q_r*rho_fac * lam_r**(1._sp+alpha_r+dr) / &
    	(lam_r+f_r)**(1._sp+alpha_r+dr+b_r), 0._sp),10._sp)
    
    vnr(:)=max(fall_n_r*rho_fac * lam_r**(1._sp+alpha_r) / &
    	(lam_r+f_r)**(1._sp+alpha_r+b_r), 0._sp)
    
    ! cloud
    vqc(:)=min(max(fall_q_c*rho_fac * lam_c**(1._sp+alpha_c+1._sp) / &
    	(lam_c+f_c)**(1._sp+alpha_c+1._sp+b_c), 1.e-3_sp), 10._sp)
    
    vnc(:)=max(fall_n_c*rho_fac * lam_c**(1._sp+alpha_c) / &
    	(lam_c+f_c)**(1._sp+alpha_c+b_c), 1.e-3_sp)
    ! coalescence efficiencies
    egi_dry=0.2_sp*exp(0.08*(t(1:kp)-ttr))
    egs_dry=0.2_sp*exp(0.08*(t(1:kp)-ttr))
    esi=0.2_sp*exp(0.08*(t(1:kp)-ttr))
    eii=0.2_sp*exp(0.08*(t(1:kp)-ttr))
    ess=0.2_sp*exp(0.08*(t(1:kp)-ttr))
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
     
    
   
    
    ! loop over all levels
    do k=1,kp
 		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! condensation of liquid water                                                   !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! smr at 0c
		!q0sat=eps1*svp_liq(ttr)/(p(k)-svp_liq(ttr))
		q0sat=eps1*svp_liq(ttr)/(p(k)-svp_liq(ttr))
    	smr(k)=eps1*svp_liq(t(k))/(p(k)-svp_liq(t(k))) ! saturation mixing ratio
        
        des_dt=dfridr(svp_liq,t(k),1.e0_sp,err)
        dqs_dt=eps1*p(k)*des_dt/(p(k)-svp_liq(t(k)))**2
        qold=q(k,iqc)
        qtot=q(k,1)+q(k,iqc)

		
        q(k,iqc)=q(k,1)+q(k,iqc)-smr(k)
        if (theta_flag) q(k,iqc)=(q(k,iqc)+(lv/cp*qold)*dqs_dt) / (1._sp+lv/cp*dqs_dt)
        q(k,iqc)=max(q(k,iqc),0._sp)
        t(k)=t(k)
        if(theta_flag) t(k)=t(k)+lv/cp*(q(k,iqc)-qold)
		
		tc=t(k)-ttr
    	smr(k)=eps1*svp_liq(t(k))/(p(k)-svp_liq(t(k))) ! saturation mixing ratio
    	q0sat=smr(k)	
    	smr_i(k)=eps1*svp_ice(t(k))/(p(k)-svp_ice(t(k))) ! saturation mixing ratio - ice	
    	
    	cond=(q(k,iqc)-qold)
    	q(k,1)=q(k,1)-cond
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    enddo
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate gamma_t and dep_density for ice growth model                             !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(ice_flag) then
        do k=1,kp
            call chen_and_lamb_anc(t(k),q(k,1),smr_i(k),rhoa(k), &
                                    gamma_t(k), dep_density(k))    
        enddo
    endif 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
    
    ! loop over all levels
    do k=1,kp      
        ! inhomogeneous mixing -https://journals.ametsoc.org/doi/pdf/10.1175/2007JAS2374.1
 !        if(q(2,k)<qold) then
!             q(4,k)=q(4,k)*(q(2,k)/qold)**1._sp
!         endif
        
        
        
        

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
	    if((q(k,iqc) .gt. qsmall) .and. (q(k1,iqc) .le. qsmall)) then
	    
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
            
!                 n_aer1=max(n_aer1,0.1_sp)

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Bulk Aerosol Activation - number of drops
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! calculate aerosol PSD parameters
            p_test=p(k)
            t_test=t(k)
            w_test=max(u(k),0.001_sp)
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

            act_frac1=max(act_frac1,0._sp)
            temp1=sum(n_aer1*act_frac1)
            ! put in-cloud aerosol into aerosol - i.e. remove it first
            do i=1,n_mode-1
                q(k,cst(i+1))=  q(k,cst(i+1))+q(k,cst(cat_c)+(i-1)*3+2)
                q(k,cst(i+1)+1)=q(k,cst(i+1)+1)+q(k,cst(cat_c)+(i-1)*3+3)
                q(k,cst(i+1)+2)=q(k,cst(i+1)+2)+q(k,cst(cat_c)+(i-1)*3+4)
                q(k,cst(cat_c)+(i-1)*3+2)=0._sp
                q(k,cst(cat_c)+(i-1)*3+3)=0._sp
                q(k,cst(cat_c)+(i-1)*3+4)=0._sp
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
                q(k,cst(i+1))=q(k,cst(i+1))-dummy1 
                    
                ! surface area in aerosol modes
                dummy2=pi*ln_part_mom(2,dcrit2(i),q(k,cst(i+1)), sig_aer1(i),d_aer1(i))
                q(k,cst(i+1)+1)=q(k,cst(i+1)+1)- dummy2 
                    
                ! mass in aerosol modes
                dummy3=pi/6._sp*density_core1(i)* &
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
            dummy3=pi/6._sp*density_core1(n_mode)* &
                ln_part_mom(3,dcrit2(n_mode),n_mix, &
                            sig_aer1(n_mode),d_aer1(n_mode))
                            
            q(k,cst(cat_am))=q(k,cst(cat_am))*(1._sp-min(dummy1/n_mix,1._sp))
            do i=1,n_mode-1 ! deplete aerosol
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! add to aerosol particles in cloud water                                !
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! number in aerosol modes
                ! qv, n_mode aerosol + 1
                q(k,cst(cat_c)+(i-1)*3+2)=q(k,cst(cat_c)+(i-1)*3+2)+ &
                    max(dummy1/(n_mix)*q(k,cst(cat_am)+3*(i-1)+1),0._sp)
                
                ! surface area in aerosol modes
                q(k,cst(cat_c)+(i-1)*3+3)=q(k,cst(cat_c)+(i-1)*3+3)+ &
                    max(dummy2/(s_mix)*q(k,cst(cat_am)+3*(i-1)+2),0._sp)
                
                ! mass in aerosol modes
                q(k,cst(cat_c)+(i-1)*3+4)=q(k,cst(cat_c)+(i-1)*3+4)+ &
                    max(dummy3/(m_mix)*q(k,cst(cat_am)+3*(i-1)+3),0._sp)
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


                ! number - remove aerosol particles
                q(k,cst(cat_am)+3*(i-1)+1)=q(k,cst(cat_am)+3*(i-1)+1) * &
                        (1._sp-max(dummy1/(n_mix),0._sp))
                
                ! surface area
                q(k,cst(cat_am)+3*(i-1)+2)=q(k,cst(cat_am)+3*(i-1)+2)* &
                        (1._sp-max(dummy2/(s_mix),0._sp) )
                ! mass
                q(k,cst(cat_am)+3*(i-1)+3)=q(k,cst(cat_am)+3*(i-1)+3)* &
                        (1._sp-max(dummy3/(m_mix),0._sp) )


            enddo 
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
        endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! end the activation of cloud drops                                              !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        
        
        
        
        
        
    
    
    
    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! ice nucleation block                                                           !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(ice_flag.and.(t(k)<ttr)) then
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! ice nucleation from cloud water                                            !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if(q(k,cst(cat_c)) > 0._sp) then
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! ice nucleation via immersion                                           !
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! set the log normal parameters for the aerosol in cloud-water
                ! to n_mix, s_mix, m_mix 
                ! (and save corresponding sig_aer1(n_mode) and d_aer1(n_mode))
                ! this is the aerosol distribution that we would have if we just took the 
                ! drops and evaporated them
                call ln_params_and_props_from_integral_moms( &
                    n_mode, &
                    q(k,cst(cat_c)), & ! total number
                    q(k,cst(cat_c)+3:cen(cat_c)-1:3), & ! surface area
                    q(k,cst(cat_c)+4:cen(cat_c):3), & ! mass 
                    n_aer1(n_mode),density_core1, &
                    molw_core1,nu_core1, & 
                    sig_aer1(n_mode),d_aer1(n_mode),n_mix,s_mix,m_mix)

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! Ice nucleation                                                         !
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! using the aerosol parameters work out how many nucleate ice 
                ! (nin_c) and down to which size (din_c)
                call ice_nucleation_aerosol(nin_c,din_c, &
                    n_mix, &    ! number
                    sig_aer1(n_mode), &  ! sigma 
                    d_aer1(n_mode), t(k) ,q(k  ,ini))     ! d
             
                ! increase ice crystal number
                q(k  ,ini)=q(k  ,ini)+nin_c
                ! increase ice crystal mass
                dummy1=nin_c/q(k  ,inc)*q(k, iqc)
                q(k  ,iqi)=q(k  ,iqi)+dummy1
                
                ! increase ice crystal shape factor
                q(k  ,iqi+1)=q(k  ,iqi+1)+nin_c
                ! increase ice crystal volume factor
                q(k  ,iqi+2)=q(k  ,iqi+2)+dummy1/rhoi
                ! increase ice crystal monomers
                q(k  ,iqi+3)=q(k  ,iqi+3)+nin_c
                
                
                ! deplete cloudnc
                !q(k  ,iqc)=q(k  ,iqc)-dummy1
                q(k,  inc)=q(k,inc)-nin_c
                pifrw(k)=pifrw(k)+dummy1/dt
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! remove aerosol from cloud water and add to ice                         !
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                n_tot=sum(q(k,cst(cat_c)+2:cen(cat_c)-2:3)) ! total number in cw
                s_tot=sum(q(k,cst(cat_c)+3:cen(cat_c)-1:3))
                m_tot=sum(q(k,cst(cat_c)+4:cen(cat_c):3))
                do i=1,n_mode-1
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! remove from cloud water:                                           !
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! number in aerosol modes
                    dummy1=ln_part_mom(0,din_c,n_mix, sig_aer1(n_mode),d_aer1(n_mode)) * &
                        q(k,cst(cat_c)+(i-1)*3+2)/n_tot
                    q(k,cst(cat_c)+(i-1)*3+2)=q(k,cst(cat_c)+(i-1)*3+2)-dummy1
                    ! surface area in aerosol modes
                    dummy2=pi* &
                        ln_part_mom(2,din_c,n_mix, sig_aer1(n_mode),d_aer1(n_mode)) * &
                        q(k,cst(cat_c)+(i-1)*3+3)/s_tot
                    q(k,cst(cat_c)+(i-1)*3+3)=q(k,cst(cat_c)+(i-1)*3+3)- dummy2 
                    
                    ! mass in aerosol modes
                    dummy3=pi/6._sp*density_core1(n_mode)* &
                        ln_part_mom(3,din_c,n_mix, sig_aer1(n_mode),d_aer1(n_mode)) * &
                        q(k,cst(cat_c)+(i-1)*3+4)/m_tot
                    q(k,cst(cat_c)+(i-1)*3+4)=q(k,cst(cat_c)+(i-1)*3+4)- dummy3
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! add to aerosol particles in ice water                              !
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! number in aerosol modes
                    ! qv, n_mode aerosol + 1
                    q(k,iai+(i-1)*3)=q(k,iai+(i-1)*3)+dummy1 
                    
                    ! surface area in aerosol modes
                    q(k,iai+(i-1)*3+1)=q(k,iai+(i-1)*3+1)+dummy2 
                    
                    ! mass in aerosol modes
                    q(k,iai+(i-1)*3+2)=q(k,iai+(i-1)*3+2)+dummy3
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                enddo
            endif 
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! end ice nucleation from cloud water                                        !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! ice nucleation from rain water                                             !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if(q(k,cst(cat_r)) > 0._sp) then
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! ice nucleation via immersion                                           !
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! set the log normal parameters for the aerosol in rain-water
                ! to n_mix, s_mix, m_mix 
                ! (and save corresponding sig_aer1(n_mode) and d_aer1(n_mode))
                ! this is the aerosol distribution that we would have if we just took the 
                ! rain drops and evaporated them
                call ln_params_and_props_from_integral_moms( &
                    n_mode, &
                    q(k,cst(cat_r)), & ! total number
                    q(k,cst(cat_r)+3:cen(cat_r)-1:3), & ! surface area
                    q(k,cst(cat_r)+4:cen(cat_r):3), & ! mass 
                    n_aer1(n_mode),density_core1, &
                    molw_core1,nu_core1, & 
                    sig_aer1(n_mode),d_aer1(n_mode),n_mix,s_mix,m_mix)

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! Ice nucleation                                                         !
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! using the aerosol parameters work out how many nucleate ice 
                ! (nin_r) and down to which size (din_r)
                call ice_nucleation_aerosol(nin_r,din_r, &
                    n_mix, &    ! number
                    sig_aer1(n_mode), &  ! sigma 
                    d_aer1(n_mode), t(k), q(k  ,ini) )     ! d
            
                ! increase in ice crystal number
                q(k  ,ini)=q(k  ,ini)+nin_r
                ! increase in ice crystal mass
                dummy1=nin_r/q(k  ,inr)*q(k, iqr)
                q(k  ,iqi)=q(k  ,iqi)+dummy1
                
                ! increase ice crystal shape factor
                q(k  ,iqi+1)=q(k  ,iqi+1)+nin_r
                ! increase ice crystal volume factor
                q(k  ,iqi+2)=q(k  ,iqi+2)+dummy1/rhoi
                ! increase ice crystal monomers
                q(k  ,iqi+3)=q(k  ,iqi+3)+nin_r


                ! deplete rain
                !q(k  ,iqr)=q(k  ,iqr)-dummy1
                q(k,  inr)=q(k,inr)-nin_r
                pgfr(k)=pgfr(k)+dummy1/dt
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! remove aerosol from rain water and add to ice                          !
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                n_tot=sum(q(k,cst(cat_r)+2:cen(cat_r)-2:3))
                s_tot=sum(q(k,cst(cat_r)+3:cen(cat_r)-1:3))
                m_tot=sum(q(k,cst(cat_r)+4:cen(cat_r):3))
                do i=1,n_mode-1
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! remove from rain water:                                            !
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! number in aerosol modes
                    dummy1=ln_part_mom(0,din_r,n_mix, sig_aer1(n_mode),d_aer1(n_mode)) * &
                        q(k,cst(cat_r)+(i-1)*3+2)/n_tot
                    q(k,cst(cat_r)+(i-1)*3+2)=q(k,cst(cat_r)+(i-1)*3+2)-dummy1 
                    
                    ! surface area in aerosol modes
                    dummy2=pi* &
                        ln_part_mom(2,din_r,n_mix, sig_aer1(n_mode),d_aer1(n_mode)) * &
                        q(k,cst(cat_r)+(i-1)*3+3)/s_tot
                    q(k,cst(cat_r)+(i-1)*3+3)=q(k,cst(cat_r)+(i-1)*3+3)- dummy2 
                    
                    ! mass in aerosol modes
                    dummy3=pi/6._sp*density_core1(n_mode)* &
                        ln_part_mom(3,din_r,n_mix, sig_aer1(n_mode),d_aer1(n_mode)) *&
                        q(k,cst(cat_r)+(i-1)*3+4)/m_tot
                    q(k,cst(cat_r)+(i-1)*3+4)=q(k,cst(cat_r)+(i-1)*3+4)- dummy3
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! add to aerosol particles in ice water                              !
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! number in aerosol modes
                    ! qv, n_mode aerosol + 1
                    q(k,iai+(i-1)*3)=q(k,iai+(i-1)*3)+dummy1 
                    
                    ! surface area in aerosol modes
                    q(k,iai+(i-1)*3+1)=q(k,iai+(i-1)*3+1)+dummy2 
                    
                    ! mass in aerosol modes
                    q(k,iai+(i-1)*3+2)=q(k,iai+(i-1)*3+2)+dummy3
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                enddo
            endif 
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! end ice nucleation from rain water                                         !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! end ice nucleation block                                                       !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!









		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! evaporation of rain                                                            !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		diff1=dd(t(k),p(k))
		ktherm1=ka(t(k))
		ab_liq=lv**2 / (ktherm1*rv*t(k)**2) + 1._sp/(rho(k)*smr(k)*diff1)
		nu_vis=viscosity_air(t(k)) / rho(k)
		sc=nu_vis / diff1
		nu_rain=2._sp*pi*n_r(k) / rho(k) * &
				(nu_r1 / lam_r(k)**(2._sp+alpha_r) + &
				(a_r/nu_vis)**0.5_sp*sc**(1._sp/3._sp)* &
				(rho(k)*rho0)**0.25_sp*nu_r2 / &
				(lam_r(k)+0.5_sp*f_r)**(0.5_sp*b_r+alpha_r+2.5_sp))

	
		rain_evap=(q(k,1)/smr(k)-1._sp) / (rho(k)*ab_liq)*nu_rain
		if(q(k,1).gt.smr(k)) then
			prevp(k)=0._sp
		else
			prevp(k)=-min(rain_evap,0._sp)
		endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! end evaporation of rain                                                        !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!









		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! deposition & sublimation onto ice                                              !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
		if ((t(k).le.ttr) .and. ice_flag) then
		    if(q(k,iqi).gt.qsmall) then
                nu_ice=2._sp*pi*n_i(k) / rho(k) * &
                        (nu_i1 / lam_i(k)**(2._sp+alpha_i) + &
                        (a_i/nu_vis)**0.5_sp*sc**(1._sp/3._sp)* &
                        (rho(k)*rho0)**0.25_sp*nu_i2 / &
                        (lam_i(k)+0.5_sp*f_i)**(0.5_sp*b_i+alpha_i+2.5_sp))
            
                ab_ice=ls**2 / (ktherm1*rv*t(k)**2) + 1._sp/(rho(k)*smr_i(k)*diff1)

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! chen and lamb                                                          !
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                vol=q(k,iqi+2)
                if(vol.gt.0._sp) then
                    phi=min(max(q(k,iqi+1) / q(k,ini),1.e-5_sp),100._sp)
                    nu_ice=nu_ice*chen_and_lamb_cap_fac(phi)
            
                    ! non chen and lamb bit        
                    ice_dep=(q(k,1)/smr_i(k)-1._sp) / (rho(k)*ab_ice)*nu_ice
            
                    if(q(k,1).gt.smr_i(k)) then
                        pisub(k)=0._sp
                        pidep(k)=min(max(ice_dep,0._sp),q(k,1)-smr_i(k))
                    else
                        pidep(k)=0._sp
                        pisub(k)=min(-min(ice_dep,0._sp),q(k,iqi)/dt)
                    endif
                    !!!
                
                        
                    call chen_and_lamb_prop((pidep(k)-pisub(k))*dt,gamma_t(k), &
                        vol,phi, dep_density(k))
                    q(k,iqi+2)=vol
                    q(k,iqi+1)=phi*q(k,ini)
                endif
            endif
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! end deposition & sublimation onto ice                                          !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	




		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! warm rain autoconversion based on Seifert and Beheng (2006)                    !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		call seifert_beheng(sb_aut,sb_acr, sb_cwaut, sb_cwacr, sb_raut, &
		            sb_rsel, sb_cwsel, q(k,cst(cat_c)+1),q(k,cst(cat_c)),&
		            q(k,cst(cat_r)+1),q(k,cst(cat_r)),rho(k),dt)
		praut(k)=sb_aut
		pracw(k)=sb_acr
		rcwaut(k)=sb_cwaut
		rcwacr(k)=sb_cwacr
		rraut(k)=sb_raut
		rrsel(k)=sb_rsel
		rcwsel(k)=sb_cwsel
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! end rain auto-conversion                                                       !		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! collection of cloud by ice                                                     !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if ((t(k).le.ttr).and.ice_flag) then
		    if((q(k,iqc).gt.qsmall) .and. (q(k,iqi).gt.qsmall)) then
                piacw(k)=max(mass_iacw * n_i(k)* eiw *q(k,iqc)*rho_fac(k) / &
                        (lam_i(k)+f_i)**(3._sp+b_i+alpha_i),0._sp)
                
                dummy1=q(k,iqi) ! total mass of particle
                
                piacw(k)=max(min(piacw(k),q(k,iqc)/dt),0._sp)
                
            endif
		endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! ice aggregation see Ferrier (1994)                                             !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if(ice_flag.and.(t(k).lt.ttr)) then
		    ! collisions
		    dummy1=max(iice*n_i(k)**2._sp*rho_fac(k) / &
                    max(lam_i(k),1.e5_sp)**(4._sp+2.*sp*alpha_i+b_i),0._sp)
                        
            
            ! aggregation rate
            riaci(k)=eii(k)*dummy1
            
            !dummy2=dummy1*(1._sp-eii(K))*100._sp*dt

            !q(k,ini)=q(k,ini)+dummy2
            !q(k,iqi+3)=q(k,iqi+3)+dummy2 ! update the number of monomers
            
            
            
        endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! end ice aggregation                                                            !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! melting of ice                                                                 !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if((t(k).gt.ttr).and.ice_flag) then
			pimlt(k)=q(k,iqi)/dt+ &
			    (pidep(k)-pisub(k)+piacw(k)) ! ice melts instantaneously
        endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! end melting of ice                                                             !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    enddo
    
    
    
    
    ! update variables
    ! vapour mass
    q(1:kp,1)=q(1:kp,1)+(pgsub+pssub+pisub-(psdep+pidep+piprm+pgdep))*dt
    ! ice mass and number
    if(ice_flag) then 
        q(1:kp,iqi)=q(1:kp,iqi)+(pidep-pisub)*dt
        q(1:kp,iqi+4)=q(1:kp,iqi)+(-pisub)*dt
        q(1:kp,ini)=q(1:kp,ini)-(riaci)*dt
    endif

    		
    ! liquid mass - riming not done here
    q(1:kp,iqc)=q(1:kp,iqc)-((pgacw+praut+psacw+pracw+pihal+picnt+pifrw))*dt
    


    ! rain number
    q(1:kp,cst(cat_r))=q(1:kp,cst(cat_r))+ &
        (rraut+rrsel-pgfr*q(1:kp,cst(cat_r))/(q(1:kp,cst(cat_r)+1)+1.e-20_sp))*dt

    ! rain mass
    q(1:kp,cst(cat_r)+1)=q(1:kp,cst(cat_r)+1)+(pgmlt+praut+pgshd+pracw+psmlt+pimlt- &
    			(pgacr+pgfr+psacr+piacr_g+piacr_s))*dt
    ! treat rain evaporation separately - adjust
    prevp=min(prevp,q(1:kp,cst(cat_r)+1)/dt) 
    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! rain aerosol
    ! aerosol going into rain, by coll-coal
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do k=1,kp
        q(k,cst(cat_r)+2:cen(cat_r))=q(k,cst(cat_r)+2:cen(cat_r))+ &
           min(max(-(rcwaut(k)+rcwacr(k)+rcwsel(k))*dt / (q(k,cst(cat_c))+qsmall),0._sp),1._sp) * &
                q(k,cst(cat_c)+2:cen(cat_c))
    enddo
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
    ! liquid aerosol
    do k=1,kp
        q(k,cst(cat_c)+2:cen(cat_c))=q(k,cst(cat_c)+2:cen(cat_c))- &
         min(max(-(rcwaut(k)+rcwacr(k)+rcwsel(k))*dt / (q(k,cst(cat_c))+qsmall),0._sp),1._sp) * &
                q(k,cst(cat_c)+2:cen(cat_c))
    enddo
    ! liquid number
    q(1:kp,inc)=q(1:kp,inc)- &
        min(max(-(rcwaut+rcwacr+rcwsel)*dt/(q(1:kp,inc)+qsmall),0._sp),1._sp)* &
                            q(1:kp,inc)
        
        
    rho=1._sp ! fudge for advection conservation
    do k=-o_halo+1,kp+o_halo
        if(q(k,iqc) .lt. qsmall) then ! if evaporated
            do i=1,n_mode-1
                ! add aerosol in cloud water back to aerosol
                q(k,cst(i+1))      =q(k,cst(i+1))   +q(k,cst(cat_c)+(i-1)*3+2)
                q(k,cst(i+1)+1)    =q(k,cst(i+1)+1) +q(k,cst(cat_c)+(i-1)*3+3)
                q(k,cst(i+1)+2)    =q(k,cst(i+1)+2) +q(k,cst(cat_c)+(i-1)*3+4)
                q(k,cst(cat_c)+(i-1)*3+2)=0._sp
                q(k,cst(cat_c)+(i-1)*3+3)=0._sp
                q(k,cst(cat_c)+(i-1)*3+4)=0._sp
            enddo
            q(k,inc) = 0.0_sp
        endif
    enddo





    do k=1,kp
        if((prevp(k) .gt. 0._sp)) then
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! add up the total number of aerosol in rain - all modes
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            dummy1=0._sp
            do i=1,n_mode-1
                dummy1=dummy1+q(k,cst(cat_r)+(i-1)*3+2)
            enddo
            !dummy1=max(dummy1,q(k,cst(cat_r)),1._sp)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! calculate the number conc. of rain drops evaporated
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            dummy2=prevp(k)*(q(k,cst(cat_r))/(qsmall+q(k,cst(cat_r)+1)))*dt
            dummy2=min(dummy2,q(k,cst(cat_r)))
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





            if(dummy2 .lt. qsmall) cycle



        
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! add evaporated rain particles to mixed-mode aerosol
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            q(k,cst(cat_am))=q(k,cst(cat_am))+dummy2 ! total number of the mixed-mode

            do i=1,n_mode-1
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! add aerosol in evaporating rain water back to aerosol
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                ! this is number in aerosol, plus number in rain 
                !  (scaled by fraction in composition category)
                ! aer_in_rain * rain_num_evap / rain_num
                q(k,cst(cat_am)+(i-1)*3+1)      = &
                    q(k,cst(cat_am)+(i-1)*3+1)   + &
                    q(k,cst(cat_r)+(i-1)*3+2) * &
                    min(dummy2/(q(k,cst(cat_r))+qsmall),1._sp)
                    
                ! this is surface area going into aerosol
                q(k,cst(cat_am)+(i-1)*3+2)    = &
                    q(k,cst(cat_am)+(i-1)*3+2) +&
                    q(k,cst(cat_r)+(i-1)*3+3) * &
                    min(dummy2/(q(k,cst(cat_r))+qsmall),1._sp)
                ! this is mass going into aerosol
                q(k,cst(cat_am)+(i-1)*3+3)    = &
                    q(k,cst(cat_am)+(i-1)*3+3) +&
                    q(k,cst(cat_r)+(i-1)*3+4) * &
                    min(dummy2/(q(k,cst(cat_r))+qsmall),1._sp)
                
                ! aerosol in rain
                q(k,cst(cat_r)+(i-1)*3+2)=q(k,cst(cat_r)+(i-1)*3+2)* &
                    (1._sp - min(dummy2/(q(k,cst(cat_r))+qsmall),1._sp))
                    
                q(k,cst(cat_r)+(i-1)*3+3)=q(k,cst(cat_r)+(i-1)*3+3)* &
                    (1._sp - min(dummy2/(q(k,cst(cat_r))+qsmall),1._sp))
                    
                q(k,cst(cat_r)+(i-1)*3+4)=q(k,cst(cat_r)+(i-1)*3+4)* &
                    (1._sp - min(dummy2/(q(k,cst(cat_r))+qsmall),1._sp))
                    
            enddo
            ! rain number
            q(k,cst(cat_r)) = q(k,cst(cat_r)) * &
                (1._sp - min(dummy2/(q(k,cst(cat_r))+qsmall),1._sp ))
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        endif
        
    enddo
    

    ! metling ice!
    do k=1,kp
        if((pimlt(k) .gt. 0._sp)) then
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! add up the total number of aerosol in ice - all modes
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            dummy1=0._sp
            do i=1,n_mode-1
                dummy1=dummy1+q(k,iai+(i-1)*3)
            enddo
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! calculate the number conc. of ice melted
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            dummy2=pimlt(k)*(q(k,cst(cat_i))/(qsmall+q(k,cst(cat_i)+1)))*dt
            dummy2=min(dummy2,q(k,cst(cat_i)))
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





            if(dummy2 .lt. qsmall) cycle



        

            do i=1,n_mode-1
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! add aerosol in melting ice water back to aerosol
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                ! this is number in aerosol, plus number in ice 
                !  (scaled by fraction in composition category)
                ! aer_in_ice * ice_num_evap / ice_num
                q(k,cst(cat_r)+(i-1)*3+2)      = &
                    q(k,cst(cat_r)+(i-1)*3+2)   + &
                    q(k,iai+(i-1)*3) * &
                    min(dummy2/(q(k,cst(cat_i))+qsmall),1._sp)
                    
                ! this is surface area going into rain aerosol
                q(k,cst(cat_r)+(i-1)*3+3)    = &
                    q(k,cst(cat_r)+(i-1)*3+3) +&
                    q(k,iai+(i-1)*3+1) * &
                    min(dummy2/(q(k,cst(cat_i))+qsmall),1._sp)
                ! this is mass going into rain aerosol
                q(k,cst(cat_r)+(i-1)*3+4)    = &
                    q(k,cst(cat_r)+(i-1)*3+4) +&
                    q(k,iai+(i-1)*3+2) * &
                    min(dummy2/(q(k,cst(cat_i))+qsmall),1._sp)
                
                ! aerosol in ice
!                 q(k,iai+(i-1)*3)=q(k,iai+(i-1)*3)* &
!                     (1._sp - min(dummy2/(q(k,cst(cat_i))+qsmall),1._sp))
!                     
!                 q(k,iai+(i-1)*3+1)=q(k,iai+(i-1)*3+1)* &
!                     (1._sp - min(dummy2/(q(k,cst(cat_i))+qsmall),1._sp))
!                     
!                 q(k,iai+(i-1)*3+2)=q(k,iai+(i-1)*3+2)* &
!                     (1._sp - min(dummy2/(q(k,cst(cat_i))+qsmall),1._sp))
                    
            enddo
            ! add the number of ice and mass to the rain
            q(k,cst(cat_r))=q(k,cst(cat_r))+dummy2
            ! mass already added
!             q(k,cst(cat_r)+1)=q(k,cst(cat_r)+1)+pimlt(k)*dt
            ! ice properties
            q(k,cst(cat_i):cen(cat_i)) = q(k,cst(cat_i):cen(cat_i)) * &
                (1._sp - min(dummy2/(q(k,cst(cat_i))+qsmall),1._sp ))
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if(t(k) > ttr) q(k,cst(cat_i):cen(cat_i))=0._sp
        endif
        
    enddo

 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! trial and error shows that riming can sometimes cause problems if it warms above 
    ! ttr. So adjust temperature in two places
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    t(1:kp)=t(1:kp)-lv/cp*prevp*dt+lf/cp*pifrw*dt+lf/cp*pgfr*dt+ls/cp*(pidep-pisub)*dt- &
                lf/cp*(pimlt)*dt
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! do the riming here                                                                 !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(ice_flag) then
        do k=1,kp
            if((q(k,iqc).lt.qsmall).or.(t(k).gt.ttr)) cycle

            dummy1=q(k,iqc)
            
            ! limit riming so 0.1 times liquid water or to a temperature 
            ! that will not raise above the melting point
            piacw(k)=max(min(piacw(k),dummy1/dt*0.1_sp,(ttr-t(k))*cp/lf)/dt*0.5_sp,0._sp)
            ! riming
            q(k,iqi)=q(k,iqi)+piacw(k)*dt
            q(k,iqi+4)=q(k,iqi+4)+piacw(k)*dt

            ! and the aerosol
            ! increase ice aerosol
            q(k,iai:cen(cat_i))=q(k,iai:cen(cat_i))+ &
                q(k,cst(cat_c)+2:cen(cat_c))*piacw(k)*dt/dummy1
            ! reduce cloud props
            q(k,cst(cat_c):cen(cat_c))=q(k,cst(cat_c):cen(cat_c))* &
                max(1._sp- piacw(k)*dt/dummy1,0._sp)
        enddo
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! adjust temperature due to riming                                                   !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    t(1:kp)=t(1:kp)+lf/cp*(piacw)*dt
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    q(1:kp,cst(cat_r)+1)=q(1:kp,cst(cat_r)+1)-prevp*dt
    
    q(1:kp,1)=q(1:kp,1)+(prevp)*dt

    q=max(q,0._sp)	


     
    if (theta_flag) th=t*(1.e5_sp/p)**(ra/cp)-theta



 	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! advection rain 0th order Bott, a.k.a. upstream advection                           !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! cloud 
    if(sum(q(1:kp,cst(cat_c)+1)).gt.qsmall) then
        adv_l(1)=.true.
		where(isnan(vqc))
			vqc=0._sp
		end where
		vqc(kp+1:kp+o_halo)=vqc(kp+o_halo)
		n_step(1)=max(ceiling(maxval(vqc(-o_halo+1:kp+o_halo)*dt/dz*2._sp)),1)
		vqc(1-o_halo:kp+o_halo-1)=-vqc(-o_halo+2:kp+o_halo)
#if MPI_PAMM == 0
		do iter=1,n_step(1)
            call mpdata_vec_1d(dt/real(n_step(1),sp),dz,dzn,&
                            rho,rhon,kp,cen(cat_c)-cst(cat_c)+1,o_halo,o_halo,&
                            vqc(-o_halo+1:kp+o_halo),&
                            q(:,cst(cat_c):cen(cat_c)),1,.false.,0)		
        enddo
#endif	
	endif
    ! rain 
    if(sum(q(1:kp,cst(cat_r)+1)).gt.qsmall) then
        adv_l(2)=.true.
		where(isnan(vqr))
			vqr=0._sp
		end where
		vqr(kp+1:kp+o_halo)=vqr(kp+o_halo)
		n_step(2)=max(ceiling(maxval(vqr(-o_halo+1:kp+o_halo)*dt/dz*2._sp)),1)
		vqr(1-o_halo:kp+o_halo-1)=-vqr(-o_halo+2:kp+o_halo)
#if MPI_PAMM == 0
		do iter=1,n_step(2)
            call mpdata_vec_1d(dt/real(n_step(2),sp),dz,dzn,&
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
                vqi=0._sp
            end where
            vqi(kp+1:kp+o_halo)=vqi(kp+o_halo)
            n_step(3)=max(ceiling(maxval(vqi(-o_halo+1:kp+o_halo)*dt/dz*2._sp)),1)
            vqi(1-o_halo:kp+o_halo-1)=-vqi(-o_halo+2:kp+o_halo)
#if MPI_PAMM == 0
            do iter=1,n_step(3)
                call mpdata_vec_1d(dt/real(n_step(3),sp),dz,dzn,&
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
    real(sp), intent(inout) :: praut,pracw, rcwaut, rcwacr, rraut, &
		                    rrsel, rcwsel
	real(sp), intent(in) :: qc,nc,qr,nr,rho, dt
	real(sp) :: lc, lr, nc1,nr1,xc_bar, phi_au, phi_ac, b_slope, tau, factor1, factor2, &
	            test
	
	
	
	    praut=0._sp
	    pracw=0._sp
	    rcwaut=0._sp
	    rcwacr=0._sp
	    rraut=0._sp
	    rrsel=0._sp
	    rcwsel=0._sp
	    
	    lc=qc*rho
	    lr=qr*rho
	    
	    nc1=max(nc*rho,lc/xstar)
	    nr1=max(nr*rho,1._sp)
	    
		b_slope=((nc1+1.e-20_sp)/(lc+1.e-20_sp))*(gam2c+1.e-20_sp)/(gam1c +1.e-20_sp)
		xc_bar=gam2c/(gam1c*b_slope+1.e-20_sp)
		
		tau=1._sp-lc/(lc+lr+1.e-20_sp)
		tau=max(tau,1.e-6_sp)
		phi_au=600._sp*tau**0.68*(1._sp-tau**0.68)**3
		phi_ac=(tau/(tau+5.e-4_sp))**4

		if (lc .gt. qsmall) then
            ! autoconversion: equation a1 in Seifert and Beheng (2001, atmos res)
            praut = kc/(20._sp*xstar)*(alpha_c+2._sp)*(alpha_c+4._sp)/(alpha_c+1._sp)**2 * &
                    (lc*xc_bar)**2*(1._sp+phi_au/(1._sp-tau+1.e-20_sp)**2)
            
            ! accretion: equation a2 in Seifert and Beheng (2001, atmos res)
            pracw=kr*lc*lr*phi_ac
        
            ! cloud number autoconversion: equation a5 in Seifert and Beheng (2001, atmos res)
            rcwaut=-2._sp/xstar*praut
            ! cloud num accretion: equation a6 in Seifert and Beheng (2001, atmos res)
            rcwacr=-1._sp/xc_bar*pracw
            ! rain num autoconversion: equation a7 in Seifert and Beheng (2001, atmos res)
            rraut=-1._sp/2._sp*rcwaut
            ! rain num self collection: equation a8 in Seifert and Beheng (2001, atmos res)
            rrsel=-kr*nr*lr
            ! cloud num self collection: equation a9 in Seifert and Beheng (2001, atmos res)
            rcwsel=-kr*(alpha_c+2._sp)/(alpha_c+1._sp)*lc**2-rcwaut

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
		use nrtype
		implicit none
		real(sp), intent(in) :: t
		real(sp) :: svp_liq
		svp_liq = 100._sp*6.1121_sp* &
			  exp((18.678_sp - (t-ttr)/ 234.5_sp)* &
			  (t-ttr)/(257.14_sp + (t-ttr)))
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
		use nrtype
		implicit none
		real(sp), intent(in) :: t
		real(sp) :: svp_ice
		svp_ice = 100._sp*6.1115_sp* &
			  exp((23.036_sp - (t-ttr)/ 333.7_sp)* &
			  (t-ttr)/(279.82_sp + (t-ttr)))
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
	function viscosity_air(t)
		use nrtype
		implicit none
		real(sp), intent(in) :: t
		real(sp) :: viscosity_air
		real(sp) :: tc

		tc = t-ttr
		tc = max(tc,-200._sp)

		if( tc.ge.0._sp) then
			viscosity_air = (1.718_sp+0.0049_sp*tc) * 1E-5_sp ! the 1d-5 converts from poise to si units
		else
			viscosity_air = (1.718_sp+0.0049_sp*tc-1.2e-5_sp*tc**2) * 1e-5_sp
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
	  use nrtype
	  implicit none
	  real(sp), intent(in) :: t, p
	  real(sp) :: dd, t1
	  t1=max(t,200._sp)
	  dd=2.11e-5_sp*(t1/ttr)**1.94_sp*(101325_sp/p)
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
	  use nrtype
	  implicit none
	  real(sp), intent(in) :: t
	  real(sp) :: ka, t1
	  t1=max(t,200._sp)
	  ka=(5.69_sp+0.017_sp*(t1-ttr))*1e-3_sp*joules_in_a_cal
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
		use nrtype
		implicit none
		real(sp), intent(in) :: t
		real(sp) :: hm_func
		hm_func=(min(max((t-265.65) / 2.5_sp,0._sp),1._sp) + &
		        min(max((270.65-t) / 2.5_sp,0._sp),1._sp)) -1.0_sp
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
		use nrtype
		implicit none
		real(sp), intent(in) :: t,naer05
		real(sp) :: demott_2010
		real(sp) :: tc
		tc=ttr-t
		! equation 1 from
		! https://www.pnas.org/content/107/25/11217
		! number per std m^3
		demott_2010=min(0.0594_sp*(tc)**3.33_sp * &
		    (naer05/1.e6_sp)**(0.0264_sp*tc+0.0033_sp),naer05)
		
	end function demott_2010
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! ice nucleation from aerosol                                                  !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief calculate the number of active INPs and the threshold diameter for 
	!> activation
	!>@param[in] n_aer,sig_aer,d_aer,T,icen
	!>@param[inout] nin,din
    subroutine ice_nucleation_aerosol(nin,din, &
                n_aer, &    ! number
                sig_aer, &  ! sigma 
                d_aer, t ,icen)     ! d
    use nrtype
    implicit none
    real(sp), intent(inout) :: nin, din
    real(sp), intent(in) :: n_aer,sig_aer,d_aer, t,icen

    real(sp) :: naer05, x, arg


    naer05=ln_part_mom(0,0.5e-6_sp,n_aer,sig_aer,d_aer)
    ! source function
    nin=demott_2010(t,naer05)
    ! limit nucleation
    nin=max(nin-icen,0._sp)

    ! deplete aerosol up to this diameter - using erfinv
    ! limit the argument so that it is not equal to -1 or +1
    arg=max(min(((1._sp-nin/n_aer)*2._sp-1._sp),1._sp-small_number),-1._sp+small_number)

    ! re-calculate nin based on limited value of arg
    nin=(1._sp-(1._sp+arg)/2._sp)*n_aer
    
    ! inverse erf
    call erfinv(arg,x)
    

    ! but x is equal to log(d/dm)/(sig_aer*sqrt(2))
    din=exp(x*sig_aer*sqrt(2._sp)+log(d_aer))
    
    if(din<0.5e-6_sp) din=0.5e-6_sp
    
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
        use nrtype
        implicit none
        real(sp), intent(in) :: dm, gamma_t,dep_density
        real(sp), intent(inout) :: v, phi

        real(sp) :: deltaV,v_old,rgamma_tp2,ln_vn_vo
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
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        rgamma_tp2=1._sp/(gamma_t+2._sp)
        ln_vn_vo=log(v/v_old)
        phi=phi*exp((gamma_t-1._sp)*rgamma_tp2*ln_vn_vo)       
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
        use nrtype
        implicit none
        real(sp), intent(in) :: t,qv,qvsat,rhoa
        real(sp), intent(inout) :: gamma_t,dep_density

        real(sp) :: delta_rho,t1
        integer(i4b) :: i
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! calculate the inherent growth ratio - this is from a 17th order polynomial
        gamma_t=0._sp
        t1=min(max(t,243.15),273.15) ! range of fit
        do i=1,n_cl
            gamma_t=gamma_t+((t1-gam_mu_cl(1))/gam_mu_cl(2))**(n_cl-i)*gam_cl(i)
        enddo
        gamma_t=10._sp**gamma_t
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! equation 42 from Chen and Lamb (1994, JAS: The Theoretical Basis for 
        !   Parameterisation of Ice Crystal Habits)
        delta_rho=(qv-qvsat)*rhoa*1000._sp ! g/m^3
        dep_density=rhoi*exp(-3._sp*max(delta_rho-0.05_sp,0._sp)/gamma_t)
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
        use nrtype
        implicit none
        real(sp), intent(in) :: phi
        real(sp) :: chen_and_lamb_cap_fac
        real(sp) :: fac1,fac2,ecc

        
        ! factor to convert between R and a - derived from equating volume of sphere to 
        ! volume of spheroid and taking the ratio of a / r
        fac1=(1._sp/(phi))**(1/3)
        
        ! factor to convert between a and capacitance
        if(phi<0.99_sp) then
            ! see equation 39 of Chen and Lamb (1994)
            ecc=sqrt(1._sp-phi**2)
            fac2=ecc/asin(ecc)
        elseif(phi>1.01_sp) then
            ! see equation 40 of Chen and Lamb (1994)
            ecc=sqrt(1._sp-(1._sp/phi)**2)
            fac2=1._sp/phi/log((1._sp+ecc)*phi)
        else
            fac2=1._sp
        endif
        
        ! total factor
        chen_and_lamb_cap_fac=fac1*fac2
        
    end function chen_and_lamb_cap_fac
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




    end module p_micro_module
    
