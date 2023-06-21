	!> @mainpage
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@copyright 2018
	!>@brief
	!>Scalar Flux Vector Transport (SFVT): 
	!>advection for dealing with particles and their properties
	!> <br> <b>Transport:</b> <br>
	!>\f$	\frac{\partial \psi}{\partial t}+\nabla \cdot \psi\bar{v}=0 
    !><br><br>
    !> This is a transport model for transporting particle properties.
    !> Scalar transport calculated.
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
        use initialisation
        use drivers_ser
        
        implicit none
        
        character (len=200) :: nmlfile = ' '
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! namelist for run variables                                           !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        namelist /run_vars/ nm1
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
        ! Allocate and initialise arrays									   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		call allocate_and_set_1d(  nm1%dt,nm1%runtime,grid1%ntim, &
			grid1%z, &
			grid1%zn, &
			grid1%w,&
			grid1%q, &
			grid1%rhoa,grid1%rhoan, &
			grid1%lamsq,grid1%lamsqn, &
			nm1%cvis, &
			grid1%dz, &
			grid1%dzn, &
			grid1%kp, grid1%nq, &
			grid1%kpstart, &
			nm1%dz, &
			nm1%kp, n_q,& 
			grid1%l_halo,grid1%r_halo)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Driver code: time-loop, advance solution, output	   				   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 		call model_driver(grid1%ntim,nm1%dt,grid1%nq, &
 		        grid1%l_halo,grid1%r_halo, &
				nm1%kp, &
				grid1%kp, &
				grid1%kpstart, &
				grid1%z, &
				grid1%dz, &
				grid1%dzn, &
				grid1%w,&
				grid1%q,&
				grid1%rhoa,grid1%rhoan, &
				grid1%lamsq,grid1%lamsqn, &
				io1%new_file, nm1%outputfile, nm1%output_interval, &
				nm1%viscous_dissipation, &
				nm1%advection_scheme, nm1%kord, nm1%monotone, nm1%neumann)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





    end program main



