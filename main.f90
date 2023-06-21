	!> @mainpage
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@copyright 2017
	!>@brief
	!>Radiative Transfer Model (RTM): 
	!>Solve radiative transfer equations on a cartesian grid
    !>
    !>
	!> <br><br>
	!> compile using the Makefile (note requires netcdf) and then run using: <br>
	!> ./main.exe namelist.in
	!> <br><br>
	!> (namelist used for initialisation).
	!> <br><br>



	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>main programme reads in information, allocates arrays, then calls the model driver

    program main
        use variables
        use mpi
        use mpi_module
        use initialisation
        use radiation
        use drivers
        use pts
        
        implicit none
        character (len=200) :: nmlfile = ' '
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! namelist for run variables                                           !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        namelist /run_vars/ nm1
        namelist /rad_vars/ nm2
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
        ! read in namelists													   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call getarg(1,nmlfile)
        open(8,file=nmlfile,status='old', recl=80, delim='apostrophe')
        read(8,nml=run_vars)
        read(8,nml=rad_vars)
        close(8)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		grid1%l_halo=1 
		grid1%r_halo=1 

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! initialise variables in mpi module:
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		call mpi_cart_initialise(nm1%kp,nm1%jp,nm1%ip)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Block until processors have synced	     						   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		call mpi_barrier(MPI_COMM_WORLD,mp1%error)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!










        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Allocate and initialise arrays for main model						   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		call allocate_and_set(  nm1%dt,nm1%runtime,grid1%ntim, &
			grid1%x, grid1%y, grid1%z, &
			grid1%xn, grid1%yn, grid1%zn, &
			grid1%u,grid1%v,grid1%w,&
			grid1%zu,grid1%zv,grid1%zw,&
			grid1%tu,grid1%tv,grid1%tw,&
			grid1%p,grid1%th,grid1%rho, grid1%q, &
			grid1%su,grid1%sv,grid1%sw,grid1%psrc, &
			grid1%theta,grid1%thetan, &
			grid1%tref,grid1%trefn, &
			grid1%rhoa,grid1%rhoan, &
			grid1%lamsq,grid1%lamsqn, &
			nm1%cvis, &
			grid1%dx, grid1%dy, grid1%dz, &
			grid1%dxn, grid1%dyn, grid1%dzn, &
			grid1%ip, grid1%jp, grid1%kp,&
			grid1%ipstart, grid1%jpstart, grid1%kpstart, &
			nm1%dx, nm1%dy, nm1%dz, &
			nm1%ip, nm1%jp, nm1%kp, & 
			nm1%n_levels,nm1%z_read(1:nm1%n_levels), &
			nm1%theta_read(1:nm1%n_levels), nm1%q_read(1,1:nm1%n_levels), &
			nm1%psurf,nm1%tsurf, &
			grid1%l_halo,grid1%r_halo, &
			grid1%coords,mp1%dims, mp1%id, mp1%ring_comm)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Allocate and initialise arrays for radiative transfer model		   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		call allocate_and_set_radiation(nm2%start_year, nm2%start_mon, &
		  nm2%start_day, nm2%start_hour,nm2%start_min, nm2%start_sec, &
		  nm2%gases_file, &
		  radg1%probs_read, radg1%lambda_low_read, radg1%lambda_high_read, &
		  radg1%h2o_read, radg1%press_read, radg1%temp_read, radg1%bli_read, &
		  radg1%molecularWeights_read, &
		  radg1%nh2o, radg1%npress, radg1%ntemp, radg1%nweights, &
		  radg1%nmolecule, radg1%nbands, &
		  radg1%itemp, radg1%ipress, &
		  grid1%ip,grid1%jp,grid1%kp,&
		  radg1%ns,radg1%nl,nm2%ns,nm2%nl, &
		  radg1%ntot, &
		  nm2%albedo,nm2%emissivity, radg1%albedo, radg1%emiss, &
		  nm2%lat_ref,nm2%lon_ref, radg1%lat, radg1%lon, &
		  nm2%nmolecule, nm2%moleculeID(1:nm2%nmolecule), &
		  nm2%lambda_read_s,nm2%lambda_read_l, &
		  nm2%lambda_s_low, nm2%lambda_s_high, &
		  nm2%lambda_l_low, nm2%lambda_l_high, &
		  radg1%lambda,radg1%lambda_low,radg1%lambda_high, radg1%delta_lambda, &
		  radg1%nrwbin,radg1%niwbin, &
		  radg1%sflux_l, radg1%b_s_g, &
		  radg1%ext_s_g, radg1%flux_u, radg1%flux_d, radg1%rad_power,  &
		  grid1%l_halo, grid1%r_halo, grid1%rhoan, grid1%thetan, grid1%trefn, &
		  radg1%nprocv,radg1%mvrecv, &
		  radg1%tdstart,radg1%tdend,radg1%a,radg1%b,radg1%c,radg1%r,radg1%u, &
		  grid1%coords,mp1%dims, mp1%id, mp1%sub_comm)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Block until processors have synced	     						   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		call mpi_barrier(MPI_COMM_WORLD,mp1%error)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        if(mp1%id < mp1%dx * mp1%dy * mp1%dz) then
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
 		call radiation_driver(grid1%ntim,nm1%dt,grid1%l_halo,grid1%r_halo, &
				nm1%ip, nm1%jp, nm1%kp, &
				grid1%ip, grid1%jp, grid1%kp, &
				grid1%ipstart, grid1%jpstart, grid1%kpstart, &
				radg1%tdstart,radg1%tdend, &
				radg1%a,radg1%b,radg1%c,radg1%r,radg1%u, &
				radg1%ntot, radg1%ns, radg1%nl, &
				radg1%flux_u, radg1%flux_d, radg1%rad_power, &
				grid1%x, grid1%y, grid1%z, &
				grid1%dx, grid1%dy, grid1%dz, &
				grid1%dxn, grid1%dyn, grid1%dzn, &
				grid1%theta,grid1%thetan, &
				grid1%tref,grid1%trefn, &
				grid1%rhoa,grid1%rhoan, &
				grid1%q, &
				radg1%lambda,radg1%lambda_low,radg1%lambda_high, radg1%delta_lambda, &
				radg1%nrwbin,radg1%niwbin, &
				radg1%sflux_l, radg1%b_s_g, &
				nm2%start_year, nm2%start_mon, nm2%start_day,&
				nm2%start_hour, nm2%start_min,nm2%start_sec, &
				radg1%lat, radg1%lon, radg1%albedo, radg1%emiss,nm2%quad_flag, &
				nm2%gas_absorption, &
				nm2%asymmetry_water, &
                  radg1%probs_read, &
                  radg1%h2o_read, radg1%press_read, radg1%temp_read, radg1%bli_read, &
                  radg1%nh2o, radg1%npress, radg1%ntemp, radg1%nweights, &
                  radg1%nmolecule, &
                  nm2%moleculeID(1:nm2%nmolecule), nm2%moleculePPM(1:nm2%nmolecule), &
                  radg1%molecularWeights_read, &
                  radg1%itemp, radg1%ipress, &
				radg1%nprocv,radg1%mvrecv, &
				grid1%coords, &
				io1%new_file, nm1%outputfile, nm1%output_interval, &
				mp1%dims,mp1%id, world_process, mp1%rank, mp1%ring_comm, &
				mp1%sub_comm)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Terminate MPI											    		   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		call MPI_Finalize ( mp1%error )
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    end program main



