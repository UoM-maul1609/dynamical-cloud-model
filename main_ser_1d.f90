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
        use nrtype
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
			grid1%u, &
			grid1%w,&
			grid1%th, &
			grid1%zu,grid1%zw,grid1%tu,grid1%tw, &
			grid1%q, &
			grid1%vism,grid1%vist,grid1%viss, &
			grid1%tau11,grid1%tau33,grid1%tau13,grid1%strain, &
			grid1%su,grid1%sw,grid1%sth,grid1%sq, &
			grid1%rhoa,grid1%rhoan, &
			grid1%theta,grid1%thetan, &
			grid1%lamsq,grid1%lamsqn, &
			nm1%z0,nm1%z0th, &
			grid1%lbc,grid1%ubc, &
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
				grid1%zn,&
				grid1%dz, &
				grid1%dzn, &
				grid1%u, &
				grid1%w,&
				grid1%th, &
				grid1%zu,grid1%zw, grid1%tu,grid1%tw, &
				grid1%q,&
				grid1%vism,grid1%vist,grid1%viss, &
				grid1%tau11,grid1%tau33, &
				grid1%tau13,grid1%strain,&
				grid1%su,grid1%sw,grid1%sth,grid1%sq, &
				grid1%rhoa,grid1%rhoan, &
				grid1%theta,grid1%thetan, &
				grid1%lamsq,grid1%lamsqn, &
				nm1%z0,nm1%z0th, &
				grid1%lbc,grid2%ubc, &
				io1%new_file, nm1%outputfile, nm1%output_interval, &
				nm1%viscous_dissipation, &
				nm1%subgrid_scheme, nm1%kord, nm1%monotone)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





    end program main



