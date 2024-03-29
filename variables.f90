	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>variables for the dynamical cloud model
    module variables
    use numerics_type
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>variables and types for the dynamical cloud model

    implicit none
    
	integer(i4b), parameter :: n_lev=1000, n_q=9
    
		!>@brief
		!>main model prognostic variables
        type grid
            ! variables for grid
            integer(i4b) :: ip, jp, kp, ntim, l_halo, r_halo, ipstart, jpstart, kpstart, &
                            nq,ncat, nprec, &
                            iqv, iqc, iqr, iqi, iqs, iqg, inc, inr, ini, ins, ing, &
                            cat_am, cat_c, cat_r, cat_i, iai
            integer(i4b), dimension(3) :: coords
            real(wp) :: f, re, g, dt, forcing_tau
            real(wp), dimension(:,:,:), allocatable :: u,v, w, rho, th, p, &
            										su, sv, sw, psrc,zu,zv,zw,tu,tv,tw, &
            										sth,strain,vism,vist, div

            real(wp), dimension(:,:,:,:), allocatable :: q, sq, viss, &
            										precip, ngs,lamgs,mugs
            real(wp), dimension(:), allocatable ::	dx, dy, dz, dxn,dyn,dzn, &
            										x, y, z, xn,yn,zn, theta, thetan, &
            										rhoa, rhoan, lamsq, lamsqn, &
            										ubar,vbar,wbar,thbar, &
            										dampfacn,dampfac, &
            										u_force,v_force, w_subs, &
            										pref,prefn,tref,trefn
            ! point to the start and end of a category
            integer(i4b), dimension(:), allocatable :: c_s, c_e
            character(len=20), dimension(:), allocatable :: q_name
            integer(i4b) :: n_mode, nrad=0
            real(wp), dimension(:,:), allocatable :: qbar
            real(wp) :: thbase, thtop
            
            integer(i4b), dimension(:), allocatable :: q_type
            logical :: micro_init=.true.
            										 
        end type grid



    											
				
	

		!>@brief
		!>variables for namelist input
        type namelist_input
            character (len=200) :: inputfile='input'
            character (len=200) :: outputfile='output'
            character (len=200) :: bam_nmlfile='input'
            character (len=200) :: aero_nmlfile='input'
            character (len=200) :: rad_nmlfile='input'
            character (len=200) :: lsm_nmlfile='input'
            logical :: add_random_height_noise, &
            			initially_geostrophic, &
            			viscous_dissipation, &
            			radiation, land_surface, &
            			dissipate_h, nudge, restart, &
            			monotone, moisture, &
            			damping_layer,forcing, aero_prof_flag,drop_num_init, theta_flag, &
            			hm_flag=.false.,wr_flag=.false., ice_flag=.false., &
            			adiabatic_prof=.false.,divergence, &
            			param_wind=.false., heyms_west=.false., lawson=.false., &
            			recycle=.true., param_theta=.true., bubble=.false.
            integer(i4b) :: nq,ip, jp, kp, subgrid_model, advection_scheme, kord, &
                        microphysics_flag,nprec, ice_nuc_flag=1, mode1_ice_flag=0, &
                        mode2_ice_flag=0, coll_breakup_flag1=0, land_surface_init=0
            real(wp) :: vis, &
            			runtime, dt, output_interval, &
            			rotation_period_hours, &
            			nudge_timescale, &
            			cvis,  &
            			dx, dy, dz, &
            			damping_thickness, damping_tau, forcing_tau, &
            			divergence_val,divergence_hgt, &
            			adiabatic_frac,t_cbase,t_ctop,rh_above,th_jump,th_grad, &
            			param_vmax, param_z, param_sigz, param_delz
            real(wp) :: psurf,tsurf,z0,z0th, ptol=1.e-8_wp, drop_num, &
                        j_stochastic=0.5e-9_wp
            integer(i4b) :: n_levels
        end type namelist_input


        real(wp), allocatable, dimension(:) :: theta_read, z_read, u_read, v_read
        real(wp), allocatable, dimension(:,:) :: q_read
        integer(i4b), dimension(:), allocatable :: q_type
        logical, dimension(:), allocatable :: q_init
        



        
		!>@brief
		!>variables for NetCDF file output
        type io
            ! variables for io
            integer(i4b) :: ncid, varid, x_dimid, y_dimid, z_dimid, &
                            dimids(2), a_dimid, xx_dimid, yy_dimid, &
                            zz_dimid, i_dimid, j_dimid, k_dimid, nq_dimid, nprec_dimid
            integer(i4b) :: icur=1
            logical :: new_file=.true.
        end type io




		! declare a namelist type
		type(namelist_input) :: nm1
        ! declare a grid type
        type(grid) :: grid1
        ! declare an io type
        type(io) :: io1


        
		
    end module variables




