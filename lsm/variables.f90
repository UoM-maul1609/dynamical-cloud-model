	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>variables for the land surface model
    module variables
	use numerics_type
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>variables and types for the land surface model

    implicit none
    
	integer(i4b), parameter :: n_lev=1000, n_q=9
    
		!>@brief
		!>main model prognostic variables
        type grid
            ! variables for grid
            integer(i4b) :: ip, jp, kp, ntim, l_halo, r_halo, ipstart, jpstart, kpstart
            integer(i4b), dimension(3) :: coords
            real(wp) :: f, re, g, dt
            real(wp), dimension(:,:,:), allocatable :: u,v, w, rho, th, p, &
            										su, sv, sw, psrc,zu,zv,zw,tu,tv,tw
            real(wp), dimension(:,:,:,:), allocatable :: q
            real(wp), dimension(:), allocatable ::	dx, dy, dz, dxn,dyn,dzn, &
            										x, y, z, xn,yn,zn, theta, thetan, &
            										tref,trefn, &
            										rhoa, rhoan, lamsq, lamsqn
            										 
        end type grid



    											
				
	

		!>@brief
		!>variables for namelist input
        type namelist_input
            character (len=200) :: inputfile='input'
            character (len=200) :: outputfile='output'
            logical :: add_random_height_noise, &
            			initially_geostrophic, &
            			viscous_dissipation, &
            			dissipate_h, nudge, restart, &
            			monotone
            integer(i4b) :: ip, jp, kp, subgrid_model, advection_scheme, kord, &
                            land_surface_init
            real(wp) :: vis, &
            			runtime, dt, output_interval, &
            			rotation_period_hours, &
            			nudge_timescale, &
            			cvis,  &
            			dx, dy, dz
            real(wp) :: psurf,tsurf
            integer(i4b) :: n_levels
            real(wp), dimension(n_lev) :: theta_read, z_read
            real(wp), dimension(n_q,n_lev) :: q_read
        end type namelist_input





        
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




