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
		grid2%l_halo=1 
		grid2%r_halo=1 









        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Allocate and initialise arrays									   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		call allocate_and_set_2d(  nm1%dt,nm1%runtime,grid2%ntim, &
			grid2%x, grid2%z, &
			grid2%xn, grid2%zn, &
			grid2%u,grid2%w,&
			grid2%q, &
			grid2%rhoa,grid2%rhoan, &
			grid2%lamsq,grid2%lamsqn, &
			nm1%cvis, &
			grid2%dx, grid2%dz, &
			grid2%dxn, grid2%dzn, &
			grid2%ip, grid2%kp, grid2%nq, &
			grid2%ipstart, grid2%kpstart, &
			nm1%dx, nm1%dz, &
			nm1%ip, nm1%kp, n_q,& 
			grid2%l_halo,grid2%r_halo)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!











        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Driver code: time-loop, advance solution, output	   				   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 		call model_driver_2d(grid2%ntim,nm1%dt,grid2%nq, &
 		        grid2%l_halo,grid2%r_halo, &
				nm1%ip, nm1%kp, &
				grid2%ip, grid2%kp, &
				grid2%ipstart, grid2%kpstart, &
				grid2%x, grid2%z, &
				grid2%dx, grid2%dz, &
				grid2%dxn, grid2%dzn, &
				grid2%u,grid2%w,&
				grid2%q,&
				grid2%rhoa,grid2%rhoan, &
				grid2%lamsq,grid2%lamsqn, &
				io1%new_file, nm1%outputfile, nm1%output_interval, &
				nm1%viscous_dissipation, &
				nm1%advection_scheme, nm1%kord, nm1%monotone,nm1%neumann)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
















    end program main



