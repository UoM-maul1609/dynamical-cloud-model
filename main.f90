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
        use initialisation
        use drivers
        
        implicit none
        
        character (len=200) :: nmlfile = ' '
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! namelist for run variables                                           !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        namelist /run_vars/ nm1
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
			grid1%lamsq,grid1%lamsqn, &
			nm1%cvis, &
			grid1%dx, grid1%dy, grid1%dz, &
			grid1%dxn, grid1%dyn, grid1%dzn, &
			grid1%ip, grid1%jp, grid1%kp,&
			grid1%ipstart, grid1%jpstart, grid1%kpstart, &
			nm1%dx, nm1%dy, nm1%dz, &
			nm1%ip, nm1%jp, nm1%kp, & 
            nm1%moisture, nm1%nq, &
			nm1%n_levels,nm1%z_read(1:nm1%n_levels), &
			nm1%theta_read(1:nm1%n_levels), &
			nm1%psurf,nm1%tsurf, &
			grid1%l_halo,grid1%r_halo, &
			grid1%thbase,grid1%thtop, &
			grid1%coords,mp1%dims, mp1%id, mp1%ring_comm)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Block until processors have synced	     						   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		call block_ring(MPI_COMM_WORLD,mp1%id,world_process,mp1%rank)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




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
				grid1%u,grid1%v,grid1%w,&
				grid1%zu,grid1%zv,grid1%zw,&
				grid1%tu,grid1%tv,grid1%tw,&
				grid1%th,grid1%sth,grid1%p, &
				grid1%su,grid1%sv,grid1%sw,grid1%psrc, &
				grid1%div, &
				grid1%strain, grid1%vism, grid1%vist, &
				nm1%z0,nm1%z0th, nm1%ptol, &
				grid1%q,grid1%sq,grid1%viss, &
				grid1%theta,grid1%thetan, &
				grid1%rhoa,grid1%rhoan, &
				grid1%lamsq,grid1%lamsqn, &
				grid1%thbase,grid1%thtop, &
				grid1%coords, &
				io1%new_file, nm1%outputfile, nm1%output_interval, &
				nm1%viscous_dissipation, &
				nm1%advection_scheme, nm1%kord, nm1%monotone, &
				nm1%moisture, nm1%damping_layer, nm1%forcing, &
				nm1%nq, &
				mp1%dims,mp1%id, world_process, mp1%rank, mp1%ring_comm, &
				mp1%sub_horiz_comm,mp1%sub_comm)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Terminate MPI											    		   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		call MPI_Finalize ( mp1%error )
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    end program main



