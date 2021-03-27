	!> @mainpage
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@copyright 2017
	!>@brief
	!>Dynamical Cloud Model (DCM): 
	!>3-D cloud model on Cartesian grid
	!> <br> <b>Transport:</b> <br>
	!>\f$	\frac{\partial \psi}{\partial t}+\nabla \cdot \psi\bar{v}=0 
    !>\f$
	!> <br> <b>Momentum:</b> <br>
	!>\f$ \rho \frac{\partial u_i}{\partial t}+\rho \vec{v}\cdot 
	!> \nabla u_i =-\frac{\partial P}{\partial x_i}-\delta _{i,3}\rho g
    !>\f$ <br><br>
    !> This is a 3-D dynamical model based on the anelastic equations of motion.
    !> Scalar transport, advection of momentum, and pressure terms are all calculated.
    !> Subgrid terms are also calculated.
    !>
    !>
	!> <br><br>
	!> compile using the Makefile (note requires netcdf) and then run using: <br>
	!> mpiexec -n 4 ./main.exe namelist.in
	!> <br><br>
	!> (namelist used for initialisation).
	!> <br><br>



	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>main programme reads in information, allocates arrays, then calls the model driver

    program main
        use nrtype
        use variables
        use mpi
        use mpi_module
        use initialisation, only : allocate_and_set, allocate_nml_qs
        use drivers
        use p_micro_module, only : read_in_pamm_bam_namelist, p_initialise_aerosol, &
                p_initialise_aerosol_3d, calculate_gamma_params
        use radiation, only : allocate_and_set_radiation, nm2, radg1
        use pts
        
        implicit none
        real(sp) :: var
        
        character (len=200) :: nmlfile = ' '
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! namelist for run variables                                           !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        namelist /run_vars/ nm1
        namelist /rad_vars/ nm2
        namelist /sounding_vars/ nm1,q_read,theta_read,z_read
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! MPI initialisation                                                   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		call MPI_Init ( mp1%error )
		call MPI_Comm_rank ( MPI_COMM_WORLD, mp1%id, mp1%error )
		call MPI_Comm_size ( MPI_COMM_WORLD, mp1%rank, mp1%error )
		mp1%wtime = MPI_Wtime ( )	
		print *,'MPI running with ID: ',mp1%id,' and rank ',mp1%rank
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! read in namelist 1												   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call getarg(1,nmlfile)
        open(8,file=nmlfile,status='old', recl=80, delim='apostrophe')
        read(8,nml=run_vars)
        close(8)
		grid1%l_halo=1 
		grid1%r_halo=1 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! set initial q-variable properties                                    !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(nm1%moisture) then
            ! set up microphysics
            select case(nm1%microphysics_flag)
                case(3) ! pamm
                    ! read in the aerosol properties
                    call read_in_pamm_bam_namelist(nm1%bam_nmlfile,&
                        nm1%aero_nmlfile, &
                        nm1%aero_prof_flag, &
                        nm1%ice_flag, &
                        grid1%q_name,grid1%q_type,grid1%c_s,grid1%c_e,grid1%nq,&
                        grid1%ncat, &
                        grid1%nprec, grid1%n_mode, &
                        grid1%iqv, grid1%iqc, grid1%inc, grid1%iqr, grid1%inr, &
                        grid1%iqi,grid1%ini,grid1%iai,grid1%cat_am, &
                        grid1%cat_c, grid1%cat_r,grid1%cat_i)    

                    nm1%nq=grid1%nq
                    nm1%nprec=grid1%nprec
                    if (nm1%radiation) grid1%nrad=2
                    if (nm1%radiation.and.nm1%ice_flag) grid1%nrad=3
                case default
                    print *, 'error'
                    stop
            end select
            
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! allocate arrays for namelist q-variables							   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call allocate_nml_qs(nm1%nq,nm1%n_levels,q_type,q_init,&
            q_read, z_read,theta_read) 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! read in namelist 2												   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call getarg(1,nmlfile)
        open(8,file=nmlfile,status='old', recl=80, delim='apostrophe')
        read(8,nml=sounding_vars)
        close(8)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! initialise variables in mpi module:
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		call mpi_cart_initialise(nm1%kp,nm1%jp,nm1%ip)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Block until processors have synced	     						   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		call block_ring(MPI_COMM_WORLD,mp1%id,world_process,mp1%rank)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Allocate and initialise arrays									   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		call allocate_and_set(  nm1%dt,nm1%runtime,grid1%ntim, &
			grid1%x, grid1%y, grid1%z, &
			grid1%xn, grid1%yn, grid1%zn, &
			grid1%ubar, grid1%vbar, grid1%wbar,grid1%thbar, grid1%qbar, &
			grid1%dampfacn, grid1%dampfac, &
			nm1%damping_layer, &
			nm1%damping_thickness,nm1%damping_tau, &
            nm1%forcing,nm1%forcing_tau, &
            grid1%forcing_tau, grid1%u_force,grid1%v_force, &
            nm1%divergence,nm1%divergence_val,nm1%divergence_hgt,grid1%w_subs, &
			grid1%u,grid1%v,grid1%w,&
			grid1%zu,grid1%zv,grid1%zw,&
			grid1%tu,grid1%tv,grid1%tw,&
			grid1%p,grid1%th,grid1%div, &
			grid1%sth,grid1%strain,grid1%vism,grid1%vist, &
			nm1%z0,nm1%z0th, &
			grid1%q, grid1%sq, grid1%viss, &
			grid1%su,grid1%sv,grid1%sw,grid1%psrc, &
			grid1%theta,grid1%thetan, &
			grid1%rhoa,grid1%rhoan, &
			grid1%pref,grid1%prefn, grid1%tref,grid1%trefn, &
			grid1%lamsq,grid1%lamsqn, &
			nm1%cvis, &
			grid1%dx, grid1%dy, grid1%dz, &
			grid1%dxn, grid1%dyn, grid1%dzn, &
			grid1%precip, &
			grid1%ip, grid1%jp, grid1%kp,&
			grid1%ipstart, grid1%jpstart, grid1%kpstart, &
			nm1%dx, nm1%dy, nm1%dz, &
			nm1%ip, nm1%jp, nm1%kp, & 
            nm1%moisture, nm1%theta_flag, nm1%nq, grid1%nq,nm1%nprec,grid1%nprec, &
            grid1%iqv, grid1%iqc, grid1%inc,nm1%drop_num_init, nm1%drop_num, &
			nm1%n_levels, &
			q_read(1:nm1%nq,1:nm1%n_levels), &
			z_read(1:nm1%n_levels), &
			theta_read(1:nm1%n_levels), &
			nm1%psurf,nm1%tsurf, &
			nm1%adiabatic_prof, nm1%adiabatic_frac,nm1%t_cbase,nm1%t_ctop, &
			nm1%rh_above, nm1%th_jump, nm1%th_grad, &
			nm1%param_wind,nm1%param_vmax,nm1%param_z,nm1%param_sigz,nm1%param_delz, &
			nm1%radiation, & ! needed to store 
            grid1%nrad, &
            grid1%ngs,grid1%lamgs,grid1%mugs,&
			grid1%l_halo,grid1%r_halo, &
			grid1%thbase,grid1%thtop, &
			grid1%coords,mp1%dims, mp1%id, mp1%ring_comm)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Radiation                                                            !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(nm1%radiation) then
            open(8,file=nm1%rad_nmlfile,status='old', recl=80, delim='apostrophe')
            read(8,nml=rad_vars)
            close(8)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Allocate and initialise radiation arrays                         !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call allocate_and_set_radiation(nm2%start_year, nm2%start_mon, &
              nm2%start_day, nm2%start_hour,nm2%start_min, nm2%start_sec, &
              grid1%ip,grid1%jp,grid1%kp,&
              radg1%ns,radg1%nl,nm2%ns,nm2%nl, &
              radg1%ntot, &
              nm2%albedo,nm2%emissivity, radg1%albedo, radg1%emiss, &
              nm2%lat_ref,nm2%lon_ref, radg1%lat, radg1%lon, &
              nm2%lambda_read_s,nm2%lambda_read_l, &
              nm2%lambda_s_low, nm2%lambda_s_high, &
              nm2%lambda_l_low, nm2%lambda_l_high, &
              radg1%lambda,radg1%lambda_low,radg1%lambda_high, radg1%delta_lambda, &
    		  radg1%nrwbin,radg1%niwbin, &
              radg1%sflux_l, radg1%b_s_g, &
              radg1%ext_s_g, radg1%flux_u, radg1%flux_d, radg1%rad_power,  &
              grid1%l_halo, grid1%r_halo, grid1%rhoan, grid1%thetan, &
              radg1%nprocv,radg1%mvrecv, &
              radg1%tdstart,radg1%tdend,radg1%a,radg1%b,radg1%c,radg1%r,radg1%u, &
              grid1%coords,mp1%dims, mp1%id, mp1%sub_comm)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! set initial q-variable properties                                    !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(nm1%moisture) then
            ! set up microphysics
            select case(nm1%microphysics_flag)
                case(3) ! pamm

                    ! p_initialise_aerosol_3d - need to have variable nq!!!
                    call p_initialise_aerosol_3d(nm1%aero_prof_flag, grid1%nq,&
                        grid1%ncat,grid1%c_s,grid1%c_e, &
                        grid1%inc, grid1%ip,grid1%jp,grid1%kp,grid1%l_halo, &
                        grid1%x,grid1%y,grid1%z,grid1%rhoan,grid1%prefn,grid1%trefn,&
                        grid1%q)  
                    if(nm1%radiation) then
                        call calculate_gamma_params(grid1%nq,grid1%ncat,grid1%n_mode,&
                            grid1%c_s,grid1%c_e,grid1%inc,grid1%iqc, &
                            grid1%inr,grid1%iqr,grid1%ini,grid1%iqi,grid1%iai, &
                            grid1%cat_am,grid1%cat_c, grid1%cat_r, grid1%cat_i,&
                            grid1%ip,grid1%jp,grid1%kp,grid1%l_halo,grid1%r_halo,grid1%q,&
                            grid1%nrad,grid1%ngs,grid1%lamgs,grid1%mugs, &
                            grid1%rhoan,nm1%ice_flag)
                    endif                        
                case default
                    print *, 'error'
                    stop
            end select
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        
                          





        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Block until processors have synced	     						   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		call block_ring(MPI_COMM_WORLD,mp1%id,world_process,mp1%rank)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        
        if((mp1%id < mp1%dx * mp1%dy * mp1%dz).and.(nm1%radiation)) then
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Set up problem                                	   	    	   !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call set_tridag(.false.,pts1%an,pts1%bn,pts1%cn,pts1%rn,&
                pts1%un,pts1%up,pts1%a,pts1%b,pts1%c,pts1%r,pts1%xsol, &
                nm1%kp,radg1%tdend,grid1%kpstart, &
                mp1%coords,mp1%dims, mp1%id, mp1%sub_comm)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        endif

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Driver code: time-loop, advance solution, output	   				   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 		call model_driver(grid1%ntim,nm1%dt,grid1%l_halo,grid1%r_halo, &
				nm1%ip, nm1%jp, nm1%kp, &
				grid1%ip, grid1%jp, grid1%kp, &
				grid1%ipstart, grid1%jpstart, grid1%kpstart, &
				grid1%x, grid1%y, grid1%z, &
				grid1%xn, grid1%yn, grid1%zn, &
				grid1%dx, grid1%dy, grid1%dz, &
				grid1%dxn, grid1%dyn, grid1%dzn, &
                grid1%ubar, grid1%vbar, grid1%wbar,grid1%thbar, grid1%qbar, &
                grid1%dampfacn, grid1%dampfac, &
                grid1%u_force, grid1%v_force, grid1%forcing_tau, &     
                grid1%w_subs, &                        
				grid1%u,grid1%v,grid1%w,&
				grid1%zu,grid1%zv,grid1%zw,&
				grid1%tu,grid1%tv,grid1%tw,&
				grid1%th,grid1%sth,grid1%p, &
				grid1%su,grid1%sv,grid1%sw,grid1%psrc, &
				grid1%div, &
				grid1%strain, grid1%vism, grid1%vist, &
				nm1%z0,nm1%z0th, nm1%ptol, &
				grid1%c_s,grid1%c_e,grid1%cat_am,grid1%cat_c,grid1%cat_r, grid1%cat_i, &
				grid1%n_mode,grid1%inc,grid1%iqc, grid1%inr,grid1%iqr, &
				grid1%ini,grid1%iqi,grid1%iai, &
				grid1%q_name, grid1%q,grid1%sq,grid1%viss, &
				grid1%precip, &
				grid1%theta,grid1%thetan, &
				grid1%tref,grid1%trefn, &
				grid1%rhoa,grid1%rhoan, &
				grid1%pref,grid1%prefn, &
				grid1%lamsq,grid1%lamsqn, &
				grid1%thbase,grid1%thtop, &
				grid1%micro_init, &
				io1%new_file, nm1%outputfile, nm1%output_interval, &
				nm1%viscous_dissipation, &
				nm1%advection_scheme, nm1%kord, nm1%monotone, &
				nm1%moisture, nm1%microphysics_flag, &
				nm1%ice_flag, &
				nm1%hm_flag, &
				nm1%theta_flag, &
				nm1%damping_layer,  nm1%forcing, nm1%divergence, nm1%radiation, &
				nm1%j_stochastic, nm1%ice_nuc_flag, nm1%mode2_ice_flag, nm1%heyms_west, &
				grid1%nq, grid1%nprec, grid1%ncat, &
                    radg1%tdstart,radg1%tdend, &
                    radg1%a,radg1%b,radg1%c,radg1%r,radg1%u, &
                    radg1%ntot, radg1%ns, radg1%nl, &
                    radg1%flux_u, radg1%flux_d, radg1%rad_power, &
                    radg1%lambda,radg1%lambda_low,radg1%lambda_high, radg1%delta_lambda, &
                    radg1%nrwbin,radg1%niwbin, &
                    radg1%sflux_l, radg1%b_s_g, &
                    nm2%start_year, nm2%start_mon, nm2%start_day,&
                    nm2%start_hour, nm2%start_min,nm2%start_sec, &
                    radg1%lat, radg1%lon, radg1%albedo, radg1%emiss,nm2%quad_flag, &
                    nm2%asymmetry_water, &
                    grid1%nrad, grid1%ngs,grid1%lamgs,grid1%mugs, &
                    radg1%nprocv,radg1%mvrecv, &
				grid1%coords, &
				mp1%dims,mp1%id, world_process, mp1%rank, mp1%ring_comm, &
				mp1%sub_horiz_comm,mp1%sub_comm)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Terminate MPI											    		   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		call MPI_Finalize ( mp1%error )
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    end program main



