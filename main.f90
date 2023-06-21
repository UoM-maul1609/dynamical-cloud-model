	!> @mainpage
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@copyright 2020
	!>@brief
	!>Land Surface Model (LSM): 
	!>Solve land surface model on a cartesian grid
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
        use numerics_type
        use variables
        use mpi
        use mpi_module
        use initialisation
        use lsm
        use drivers
        
        implicit none
        character (len=200) :: nmlfile = ' '
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! namelist for run variables                                           !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        namelist /run_vars/ nm1
        namelist /lsm_vars/ nm3
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
        read(8,nml=lsm_vars)
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
			grid1%p,grid1%th,grid1%rho, &
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
			nm1%theta_read(1:nm1%n_levels), &
			nm1%psurf,nm1%tsurf, &
			grid1%l_halo,grid1%r_halo, &
			grid1%coords,mp1%dims, mp1%id, mp1%ring_comm)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Allocate and initialise arrays for land surface model		           !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		call allocate_and_set_lsm(grid1%ip,grid1%jp,grid1%kp,&
		  grid1%l_halo, grid1%r_halo,  &
		  nm3%nsoil_lay,nm3%soil_thickness, nm3%soil_types, &
		  lsmg1%skp, lsmg1%sz,lsmg1%szn, lsmg1%dsz,lsmg1%dszn, &
		  lsmg1%t,lsmg1%wg, lsmg1%tsurf,  &
		  lsmg1%tdend,lsmg1%a,lsmg1%b,lsmg1%c,lsmg1%r,lsmg1%u, lsmg1%pscs,  &
		  lsmg1%b1, lsmg1%wgs, lsmg1%wfc, lsmg1%phi_ps, lsmg1%kgs, &
		  nm3%tinit, nm3%wg_pc_init, nm1%tsurf, nm1%land_surface_init, &
		  grid1%coords,mp1%dims, mp1%id, mp1%sub_horiz_comm) ! no comms needed
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Block until processors have synced	     						   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		call mpi_barrier(MPI_COMM_WORLD,mp1%error)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Driver code: time-loop, advance solution, output	   				   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 		call lsm_driver(grid1%ntim,nm1%dt,grid1%l_halo,grid1%r_halo, &
				nm1%ip, nm1%jp, nm1%kp, &
				grid1%ip, grid1%jp, grid1%kp, &
				grid1%ipstart, grid1%jpstart, grid1%kpstart, &
				grid1%x, grid1%y, grid1%z, &
				grid1%dx, grid1%dy, grid1%dz, &
				grid1%dxn, grid1%dyn, grid1%dzn, &
				lsmg1%skp, lsmg1%sz,lsmg1%szn,lsmg1%dsz,lsmg1%dszn, &
                lsmg1%tdend,lsmg1%t,lsmg1%wg, lsmg1%tsurf, &
                lsmg1%a,lsmg1%b,lsmg1%c,lsmg1%r,lsmg1%u, lsmg1%pscs,  &		
                lsmg1%b1, lsmg1%wgs, lsmg1%wfc, lsmg1%phi_ps, lsmg1%kgs, &		
				grid1%coords, &
				io1%new_file, nm1%outputfile, nm1%output_interval, &
				mp1%dims,mp1%id, world_process, mp1%rank, mp1%ring_comm, &
				mp1%sub_horiz_comm) 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Terminate MPI											    		   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		call MPI_Finalize ( mp1%error )
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    end program main



