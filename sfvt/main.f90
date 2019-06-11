	!> @mainpage
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@copyright 2018
	!>@brief
	!>Scalar Flux Vector Transport (SFVT): 
	!>3-D Cartesian advection for dealing with particles and their properties
	!> <br> <b>Transport:</b> <br>
	!>\f$	\frac{\partial \psi}{\partial t}+\nabla \cdot \psi\bar{v}=0 
    !><br><br>
    !> This is a 3-D transport model for transporting particle properties.
    !> Scalar transport calculated.
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
		grid3%l_halo=1 
		grid3%r_halo=1 




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
		call allocate_and_set(  nm1%dt,nm1%runtime,grid3%ntim, &
			grid3%x, grid3%y, grid3%z, &
			grid3%xn, grid3%yn, grid3%zn, &
			grid3%u,grid3%v,grid3%w,&
			grid3%q, &
			grid3%rhoa,grid3%rhoan, &
			grid3%lamsq,grid3%lamsqn, &
			grid3%lbc,grid3%ubc, &
			nm1%cvis, &
			grid3%dx, grid3%dy, grid3%dz, &
			grid3%dxn, grid3%dyn, grid3%dzn, &
			grid3%ip, grid3%jp, grid3%kp, grid3%nq, &
			grid3%ipstart, grid3%jpstart, grid3%kpstart, &
			nm1%dx, nm1%dy, nm1%dz, &
			nm1%ip, nm1%jp, nm1%kp, n_q,& 
			grid3%l_halo,grid3%r_halo, &
			grid3%coords,mp1%dims, mp1%id, mp1%ring_comm)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Block until processors have synced	     						   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		call block_ring(MPI_COMM_WORLD,mp1%id,world_process,mp1%rank)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Driver code: time-loop, advance solution, output	   				   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 		call model_driver(grid3%ntim,nm1%dt,grid3%nq, &
 		        grid3%l_halo,grid3%r_halo, &
				nm1%ip, nm1%jp, nm1%kp, &
				grid3%ip, grid3%jp, grid3%kp, &
				grid3%ipstart, grid3%jpstart, grid3%kpstart, &
				grid3%x, grid3%y, grid3%z, &
				grid3%dx, grid3%dy, grid3%dz, &
				grid3%dxn, grid3%dyn, grid3%dzn, &
				grid3%u,grid3%v,grid3%w,&
				grid3%q,&
				grid3%rhoa,grid3%rhoan, &
				grid3%lamsq,grid3%lamsqn, &
				grid3%lbc,grid3%ubc, &
				grid3%coords, &
				io1%new_file, nm1%outputfile, nm1%output_interval, &
				nm1%viscous_dissipation, &
				nm1%advection_scheme, nm1%kord, nm1%monotone, &
				mp1%dims,mp1%id, world_process, mp1%rank, mp1%ring_comm)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Terminate MPI											    		   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		call MPI_Finalize ( mp1%error )
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    end program main



