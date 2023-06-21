	!> @mainpage
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@copyright 2018
	!>@brief
	!>Sub-Grid scale Model (SGM): 
	!>Calculate sub-grid model for cartesian models
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
			grid2%u,grid2%w,grid2%th, &
			grid2%zu,grid2%zw,&
			grid2%tu,grid2%tw,&
			grid2%q, &
			grid2%vism,grid2%vist,grid2%viss,&
			grid2%tau11,grid2%tau33,&
			grid2%tau13,grid2%strain,&
			grid2%su,grid2%sw,grid2%sth,grid2%sq, & 
			grid2%rhoa,grid2%rhoan,grid2%theta, grid2%thetan, &
			grid2%lamsq,grid2%lamsqn, &
			nm1%z0,nm1%z0th, &
			grid2%lbc,grid2%ubc, &
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
				grid2%xn, grid2%zn, &
				grid2%dx, grid2%dz, &
				grid2%dxn, grid2%dzn, &
				grid2%u,grid2%w,grid2%th,&
				grid2%zu,grid2%zw, &
				grid2%tu,grid2%tw, &
				grid2%q,&
				grid2%vism,grid2%vist,grid2%viss, &
				grid2%tau11,grid2%tau33, &
				grid2%tau13,grid2%strain,&
				grid2%su,grid2%sw,grid2%sth,grid2%sq, &
				grid2%rhoa,grid2%rhoan, &
				grid2%theta,grid2%thetan, &
				grid2%lamsq,grid2%lamsqn, &
				nm1%z0,nm1%z0th, &
				grid2%lbc,grid2%ubc, &
				io1%new_file, nm1%outputfile, nm1%output_interval, &
				nm1%viscous_dissipation, &
				nm1%subgrid_scheme, nm1%kord, nm1%monotone)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
















    end program main



