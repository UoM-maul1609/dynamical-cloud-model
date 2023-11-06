	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>initialisation for the dynamical cloud model
    module initialisation
    use numerics_type
    implicit none
    
    real(wp), parameter :: ra=287._wp, rv=461.5_wp,eps1=ra/rv,grav=9.81_wp,cp=1005._wp, &
                            cw=4190._wp, ttr=273.15_wp, lv=2.5e6_wp
	real(wp), dimension(:), pointer :: zr1, thr1
	real(wp) :: tsurf_glob, p_glob,theta_q_sat_glob,rv_glob,t_glob, rh_glob,th_grad_glob,&
	            z_glob, adiabatic_frac_glob, adiabatic_frac_glob2
	integer(i4b) :: nl1
	
    private
    public :: allocate_and_set, allocate_nml_qs
    contains
    
    
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>allocate arrays for q_type and q_init
	!>@param[in] nq number of q fields
	!>@param[in] n_levels number of levels for reading in sounding
	!>@param[inout] q_type: integer array
	!>@param[inout] q_init: logical array
	!>@param[inout] q_read: real array
	!>@param[inout] z_read: real array
	!>@param[inout] th_read: real array
	!>@param[inout] u_read: real array
	!>@param[inout] v_read: real array
    subroutine allocate_nml_qs(nq,n_levels,q_type,q_init,q_read, &
                z_read,th_read, u_read, v_read)
    use numerics_type
    implicit none
    integer(i4b), intent(in) :: nq, n_levels
    integer(i4b), dimension(:), allocatable, intent(inout) :: q_type
    logical, dimension(:), allocatable, intent(inout) :: q_init
    real(wp), dimension(:), allocatable, intent(inout) :: z_read,th_read, u_read, v_read
    real(wp), dimension(:,:), allocatable, intent(inout) :: q_read
    ! local variables:
    integer(i4b) :: AllocateStatus
    
    ! allocate arrays
    allocate( q_type(1:nq), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    allocate( q_init(1:nq), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    allocate( q_read(1:nq,1:n_levels), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    allocate( z_read(1:n_levels), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    allocate( th_read(1:n_levels), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    allocate( u_read(1:n_levels), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    allocate( v_read(1:n_levels), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    
    
    end subroutine allocate_nml_qs
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Allocate and set arrays                                                            !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>allocate arrays on each PE, and initialise them
	!>@param[inout] dt,runtime,ntim - time variables
	!>@param[inout] x,y,z,xn,yn,zn,u,v,w,p,th,div - grid positions and prognostics
	!>@param[inout] ubar,vbar,wbar,thbar,qbar,dampfacn,dampfac
	!>@param[in] damping_layer, damping_thickness, damping_tau
	!>@param[in] forcing, forcing_tau_nml
	!>@param[inout] forcing_tau_g, u_force, v_force
	!>@param[in] divergence,divergence_val,divergence_hgt
	!>@param[inout] w_subs
	!>@param[inout] sth,strain,vism,vist - more prognostics
	!>@param[in] z0,z0th - roughness lengths 
	!>@param[inout] q,sq,viss - q-variable
	!>@param[inout] zu,zv,zw,tu,tv,tw - previous time-step and temporary storage
	!>@param[inout] su,sv,sw,psrc - grid positions and prognostics
	!>@param[inout] theta,thetan - reference potential temperatures
	!>@param[inout] rhoa,rhoan - reference density
	!>@param[inout] pref,prefn - reference pressures
	!>@param[inout] tref,trefn - reference temperatures
	!>@param[inout] lamsq,lamsqn - mixing length
	!>@param[in] cvis - smagorinsky parameter
	!>@param[inout] dx,dy,dz - grid spacing on grid
	!>@param[inout] dxn,dyn,dzn - grid spacing on grid - staggered
	!>@param[inout] precip
	!>@param[inout] ipp,jpp,kpp,ipstart,jpstart,kpstart - number of grid / starting position
	!>@param[in] dx_nm,dy_nm,dz_nm - grid spacing from namelist
	!>@param[in] ip,jp,kp - grid points from namelist
	!>@param[in] moisture, nq,nprec - whether to have clouds and number of q-variables
	!>@param[in] theta_flag - whether latent heat is released
	!>@param[in] nqg,nprecg - number of precipitation categories
	!>@param[in] iqv,iqc,inc,drop_num_init, num_drop: moisture variables
	!>@param[in] n_levels
	!>@params[inout] q_read,z_read,theta_read, 
	!>@params[inout] u_read,v_read,psurf,tsurf - sounding variables
	!>@params[in] adiabatic_prof,adiabatic_frac,t_cbase,t_ctop, rh_above, th_jump,th_grad
	!>@param[in] bubble, param_theta, param_wind, param_vmax,param_z,param_sigz,param_delz
	!>@param[in] radiation,nrad
	!>@param[inout] ngs,lamgs,mugs
	!>@param[in] l_h,r_h - halo for arrays
	!>@param[inout] thbase, thtop
	!>@param[in] coords,dims - dimensions of cartesian topology
	!>@param[in] id - id of this PE
	!>@param[in] comm3d - communicator for cartesian topology
	subroutine allocate_and_set(dt,runtime,ntim,x, y, z, &
			xn, yn, zn, &
			ubar, vbar, wbar,thbar, qbar, &
			dampfacn, dampfac, &
			damping_layer, &
			damping_thickness,damping_tau, &
            forcing,forcing_tau_nml, &
            forcing_tau_g, u_force,v_force, &
            divergence,divergence_val,divergence_hgt,w_subs, &
			u,v,w,&
			zu,zv,zw,&
			tu,tv,tw,&
			p,th,div, &
			sth,strain,vism,vist,z0,z0th, &
			q, sq, viss, &
			su,sv,sw,psrc, &
			theta,thetan, &
			rhoa,rhoan, &
			pref,prefn,tref,trefn, &
			lamsq,lamsqn, &
			cvis, &
			dx, dy, dz, &
			dxn, dyn, dzn, &
			precip, &
			ipp, jpp, kpp,&
			ipstart, jpstart, kpstart, &
			dx_nm, dy_nm, dz_nm, &
			ip, jp, kp, &
			moisture,theta_flag,nq,nqg,nprec,nprecg,&
			iqv,iqc,inc,drop_num_init, num_drop, &
			n_levels, &
			q_read, &
			z_read, theta_read,u_read,v_read,psurf,tsurf, &
			adiabatic_prof,adiabatic_frac,t_cbase,t_ctop, rh_above, th_jump, th_grad,&
			bubble, param_theta,param_wind, param_vmax,param_z,param_sigz,param_delz, &
			radiation,nrad,ngs,lamgs,mugs, &
			l_h,r_h, &
			thbase,thtop, &
			coords,dims, id, comm3d)
				
		use numerics_type
		use mpi
		use netcdf
		use numerics, only : find_pos, poly_int, vode_integrate
		use random, only : random_normal
		use mpi_module
		
		implicit none
		real(wp), dimension(:,:,:), allocatable, intent(inout) :: &
														u,v,w,zu,zv,zw,tu,tv,tw,&
														p,th,div, &
														su,sv,sw,psrc, &
														sth,strain,vism,vist
		real(wp), dimension(:,:,:,:), allocatable, intent(inout) :: &
														q,sq,viss, precip, ngs,lamgs,mugs
		real(wp), dimension(:), allocatable, intent(inout) :: x,y,z,xn,yn,zn,dx,dy,dz, &
															dxn,dyn,dzn, theta,thetan, &
															rhoa, rhoan, lamsq, lamsqn, &
															ubar, vbar, wbar,thbar, &
															dampfac,dampfacn, &
															u_force,v_force, &
															w_subs, &
															pref,prefn,tref,trefn
															
		real(wp), dimension(:,:), allocatable, intent(inout) :: qbar

		real(wp), intent(in) :: dx_nm, dy_nm, dz_nm, cvis,z0,z0th
		real(wp), intent(in) :: dt, runtime, forcing_tau_nml,num_drop,divergence_val, &
		                        divergence_hgt
		integer(i4b), intent(inout) :: ipp, jpp, kpp, ipstart, jpstart, kpstart
		integer(i4b), intent(inout) :: ntim
		integer(i4b), intent(in) :: ip, jp, kp, l_h, r_h, nq,nprec,iqv,iqc,inc, nrad
		integer(i4b), intent(inout) :: nqg,nprecg
		logical, intent(in) :: moisture,theta_flag,damping_layer, forcing, drop_num_init, &
		    adiabatic_prof,divergence, bubble, param_theta, param_wind, radiation
		integer(i4b), intent(in) :: n_levels
		real(wp), dimension(n_levels), target, intent(inout) :: z_read,theta_read, &
		    u_read, v_read
		real(wp), dimension(nq,n_levels), target, intent(inout) :: q_read
		real(wp), intent(in) :: psurf, tsurf,damping_thickness,damping_tau, &
		    adiabatic_frac,t_cbase,t_ctop,rh_above,th_jump,th_grad, &
		    param_vmax, param_z, param_sigz, param_delz
		integer(i4b), dimension(3), intent(inout) :: coords
		integer(i4b), dimension(3), intent(in) :: dims
		integer(i4b), intent(in) :: id, comm3d
		real(wp), intent(inout) :: thbase, thtop,forcing_tau_g
		
		! locals:
		integer(i4b) :: error, AllocateStatus,i,j,k
		real(wp) :: rho_surf, htry, hmin, eps2=1.e-5_wp, t, qsat
		real(wp), dimension(1) :: psolve
		real(wp) :: var, dummy,ztop
		integer(i4b) :: iloc
		! for random number:
		real(wp) :: r
		real(wp), dimension(10,10) :: rs
		integer(i4b) :: l, nbottom, ntop, tag1
		integer(i4b), allocatable, dimension(:) :: seed
		real(wp) :: rad
		real(wp), allocatable, dimension(:) :: qv1d,qc1d
		
! if the pe is not being used in the cartesian topology, do not use here
		if(id>=dims(1)*dims(2)*dims(3)) return 
		
		
		! set pointers
		zr1 => z_read
		thr1 => theta_read
		nl1 = n_levels

		
		! scalar formulae:
		ntim=ceiling(runtime/dt)
		nqg=nq
		nprecg=nprec
		forcing_tau_g=forcing_tau_nml

		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! find the number of grid points on each PE                                      !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		call MPI_CART_COORDS(comm3d, id, 3, coords, error)
		! print *,'Coords of ',id,' are ',coords

        ! number of grid points in all but last:
		ipp = floor(real(ip,wp)/real(dims(1),wp)) 
		ipstart = ipp*(coords(1))  +1   
		if(coords(1) == (dims(1)-1)) then
			ipp=ip-(dims(1)-1)*ipp ! number of grid points in last
		endif
		! number of grid points in all but last:
		jpp = floor(real(jp,wp)/real(dims(2),wp))      
		jpstart = jpp*(coords(2))  +1 
		if(coords(2) == (dims(2)-1)) then
			jpp=jp-(dims(2)-1)*jpp ! number of grid points in last
		endif
		! number of grid points in all but last:
		kpp = floor(real(kp,wp)/real(dims(3),wp))      
		kpstart = kpp*(coords(3)) +1    
		if(coords(3) == (dims(3)-1)) then
			kpp=kp-(dims(3)-1)*kpp ! number of grid points in last
		endif
 		!print *,ip,jp,kp,ipp,jpp,kpp,ipstart, jpstart,kpstart,coords
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		
		
		
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! allocate arrays                                                                !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		allocate( u(1-r_h:kpp+r_h,1-r_h:jpp+r_h,1-l_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( v(1-r_h:kpp+r_h,1-l_h:jpp+r_h,1-r_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( w(1-l_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( zu(1-r_h:kpp+r_h,1-r_h:jpp+r_h,1-l_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( zv(1-r_h:kpp+r_h,1-l_h:jpp+r_h,1-r_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( zw(1-l_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( tu(1-r_h:kpp+r_h,1-r_h:jpp+r_h,1-l_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( tv(1-r_h:kpp+r_h,1-l_h:jpp+r_h,1-r_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( tw(1-l_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( p(1-r_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( th(1-r_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( div(1-r_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		
		allocate( su(1-r_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( sv(1-r_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( sw(1-r_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( psrc(1-r_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"

		allocate( sth(1-r_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( strain(1-r_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( vism(1-r_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( vist(1-r_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"

        if(radiation) then
            allocate( ngs(1-r_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h,nrad), &
                STAT = AllocateStatus)
            if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
            allocate( lamgs(1-r_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h,nrad), &
                STAT = AllocateStatus)
            if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
            allocate( mugs(1-r_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h,nrad), &
                STAT = AllocateStatus)
            if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
        endif
        
        ! moisture variables
        allocate( qbar(1-r_h:kpp+r_h,1:nq), &
            STAT = AllocateStatus)
        allocate( q(1-r_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h,1:nq), &
            STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
        allocate( sq(1-r_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h,1:nq), &
            STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
        allocate( viss(1-r_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h,1:nq), &
            STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
        allocate( precip(1:kpp,1-r_h:jpp+r_h,1-r_h:ipp+r_h,1:nprec), &
            STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
                        
		allocate( qv1d(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( qc1d(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
            
		
		allocate( ubar(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( vbar(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( wbar(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( thbar(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		
        allocate( dampfacn(1-l_h:kpp+r_h), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"		
        allocate( dampfac(1-l_h:kpp+r_h), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"		

        if(forcing) then
            allocate( u_force(1-l_h:kpp+r_h), STAT = AllocateStatus)
            if (AllocateStatus /= 0) STOP "*** Not enough memory ***"          
            allocate( v_force(1-l_h:kpp+r_h), STAT = AllocateStatus)
            if (AllocateStatus /= 0) STOP "*** Not enough memory ***"          
        endif
        
        if(divergence) then
            allocate( w_subs(1-l_h:kpp+r_h), STAT = AllocateStatus)
            if (AllocateStatus /= 0) STOP "*** Not enough memory ***"          
        endif
        

		allocate( x(1-l_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( y(1-l_h:jpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( z(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( theta(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( rhoa(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( pref(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( tref(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( lamsq(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		
		allocate( xn(1-l_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( yn(1-l_h:jpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( zn(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( thetan(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( rhoan(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( prefn(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( trefn(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( lamsqn(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		
		allocate( dx(1-l_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( dy(1-l_h:jpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( dz(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"		

		allocate( dxn(1-l_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( dyn(1-l_h:jpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( dzn(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
				
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! set up grid spacing arrays                                                     !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! grid spacing:
		dx(:)=dx_nm
		dy(:)=dy_nm
		dz(:)=dz_nm
		! grid spacing, staggered:
		dxn(:)=dx_nm
		dyn(:)=dy_nm
		dzn(:)=dz_nm
		! set up horizontal level array
		x=dx_nm*(/(i,i=-l_h+ipstart,ipp+r_h+ipstart-1)/) - real(ip-1,wp)/2._wp*dx_nm 
		xn=x-0.5_wp*dx_nm
		! set up horizontal level array
		y=dy_nm*(/(i,i=-l_h+jpstart,jpp+r_h+jpstart-1)/) - real(jp-1,wp)/2._wp*dy_nm
		yn=y-0.5_wp*dy_nm

		! set up vertical level array
		z=dz_nm*(/(i,i=-l_h+kpstart-1,kpp+r_h+kpstart-2)/)+1.0_wp*dz_nm
		zn=z-0.5_wp*dz_nm
		
		! set up mixing length array
		lamsq=1._wp / (1._wp/(cvis*(dx_nm+dy_nm+dz)/3._wp)**2._wp + &
				1._wp/(0.4_wp*(z + z0))**2._wp)
		lamsqn=1._wp / (1._wp/(cvis*(dx_nm+dy_nm+dzn)/3._wp)**2._wp + &
				1._wp/(0.4_wp*(zn + z0))**2._wp)	
        ! zero some arrays
		ubar=0._wp	
		vbar=0._wp	
		wbar=0._wp	
		thbar=0._wp	
		qbar=0._wp
		q=0._wp
		sq=0._wp	
		div=0._wp
		th=0._wp
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! set up damping layer                                                           !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        dampfacn=0._wp
        dampfac=0._wp
        if(damping_layer) then
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! find top of array                                                          !                                           
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!              
            call find_top(comm3d, id, kpp,l_h,r_h,zn,ztop, dims,coords)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!              
            
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! set damping factors above threshold height                                 !                                           
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!              
            do k=1-l_h,kpp+r_h
                if(zn(k) .gt. (ztop-damping_thickness)) then
                    dampfacn(k)=-1._wp/damping_tau !*(&
                        !exp((zn(k)-(ztop-damping_thickness))/damping_thickness)-1._wp)
                endif

                if(z(k) .gt. (ztop-damping_thickness)) then
                    dampfac(k)=-1._wp/damping_tau !*(&
                        !exp((z(k)-(ztop-damping_thickness))/damping_thickness)-1._wp)
                endif
            enddo
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!              
!             if(coords(3)==(dims(3)-1)) then
!                 dampfac(kpp-1:kpp+1)=0._wp
!                 dampfacn(kpp-1:kpp+1)=0._wp
!             endif

        endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! set-up density, pressure, theta, etc                                           !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        rho_surf=psurf / (ra*tsurf)
		if(.not.adiabatic_prof) then
            do i=1-l_h,kpp+r_h
    
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! solve the hydrostatic equation on staggered points					 !
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                psolve=psurf
                if( z(i) <= 0._wp ) then
                    hmin=1.e-2_wp
                    !htry=-dz(i)
                    !call vode_integrate(psolve,0._wp,z(i),eps2,htry,hmin,hydrostatic1a)
                else
                    hmin=1.e-2_wp
                    htry=dz(i)
                    call vode_integrate(psolve,0._wp,z(i),eps2,htry,hmin,hydrostatic1a)
                endif
                p(i,:,:)=psolve(1)

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! locate and interpolate to find theta on staggered points:		    	 !
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                iloc=find_pos(z_read(1:n_levels),z(i))
                iloc=min(n_levels-1,iloc)
                iloc=max(1,iloc)
                ! linear interp theta
                call poly_int(z_read(iloc:iloc+1), theta_read(iloc:iloc+1), &
                            min(z(i),z_read(n_levels)), var,dummy)
    !			th(:,:,i)=var
                theta(i)=var
                rhoa(i)=psolve(1)/(ra*theta(i)*(psolve(1)/100000.)**(ra/cp))
                pref(i)=psolve(1)
                tref(i)=(theta(i)*(psolve(1)/100000.)**(ra/cp))
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! solve the hydrostatic equation on integer points				    	 !
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                psolve=psurf
                if( zn(i) < 0._wp ) then
                    hmin=1.e-2_wp
                    htry=-dzn(i)
                    call vode_integrate(psolve,0._wp,zn(i),eps2,htry,hmin,hydrostatic1a)
                else
                    hmin=1.e-2_wp
                    htry=dzn(i)
                    call vode_integrate(psolve,0._wp,zn(i),eps2,htry,hmin,hydrostatic1a)
                endif
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! locate and interpolate to find theta on integer points			     !
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                iloc=find_pos(z_read(1:n_levels),zn(i))
                iloc=min(n_levels-1,iloc)
                iloc=max(1,iloc)
                ! linear interp theta
                call poly_int(z_read(iloc:iloc+1), theta_read(iloc:iloc+1), &
                            min(zn(i),z_read(n_levels)), var,dummy)
                thetan(i)=var
                rhoan(i)=psolve(1)/(ra*thetan(i)*(psolve(1)/100000.)**(ra/cp))
                prefn(i)=psolve(1)
                trefn(i)=(thetan(i)*(psolve(1)/100000.)**(ra/cp))
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            enddo
        
            if(moisture) then
                q=0._wp
                do i=1-l_h,kpp+r_h
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! locate and interpolate to find q-vapour:							 !
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    iloc=find_pos(z_read(1:n_levels),zn(i))
                    iloc=min(n_levels-1,iloc)
                    iloc=max(1,iloc)
                    ! linear interp q-vapour
                    call poly_int(z_read(iloc:iloc+1), q_read(1,iloc:iloc+1), &
                                min(zn(i),z_read(n_levels)), var,dummy)
                    q(i,:,:,1)=var
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
                
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! set the cloud, etc                                                 !
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    t=thetan(i)*(prefn(i)/psurf)**(ra/cp)
                    qsat=eps1*svp_liq(t)/(prefn(i)-svp_liq(t))
                    if (var >= qsat) then
                        q(i,:,:,iqc)=var-qsat
                        q(i,:,:,iqv)=qsat
                        q(i,:,:,inc)=num_drop
                    else
                        q(i,:,:,inc)=0._wp
                    endif
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!               
                enddo
            endif
            
        else
            
            ! if we are assuming adiabatic cloud, etc
            call setup_adiabatic_profile(l_h,r_h,kpp, psurf,tsurf,t_cbase,t_ctop, &
                                        rho_surf, rh_above, th_jump, th_grad,&
                                        adiabatic_frac, theta_flag, &
                                        pref,tref,theta,rhoa,qv1d,qc1d,z,dz,&
                                        prefn,trefn,thetan,rhoan,zn,dzn, &
                                        coords,dims, id, comm3d)
            if(moisture) then
                q(i,:,:,inc)=0._wp
                do i=1-l_h,kpp+r_h
                    q(i,:,:,iqc)=qc1d(i)
                    q(i,:,:,iqv)=qv1d(i)
                    if(qc1d(i).gt.0._wp) q(i,:,:,inc)=num_drop
                enddo
            endif
            
        endif
        
        		
		
		p(:,:,:)=0._wp
! 		rhoa=1._wp
! 		rhoan=1._wp
!         if(coords(3)==0) then
!             rhoa(0)=rhoa(1)
!             rhoan(0)=rhoan(1)
!             theta(0)=theta(1)
!             thetan(0)=thetan(1)
!         endif
!         if((coords(3)+1)==dims(3)) then
!             rhoa(kpp+1)=rhoa(kpp)
!             rhoan(kpp+1)=rhoan(kpp)
!             theta(kpp+1)=theta(kpp)
!             thetan(kpp+1)=thetan(kpp)
!         endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! set wind field        														 !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
		do i=1,ipp
			do j=1,jpp
				do k=1,kpp
					u(k,j,i)=-5._wp*(yn(j))/sqrt(x(i)*x(i)+yn(j)*yn(j))
					v(k,j,i)=5._wp*(xn(i))/sqrt(xn(i)*xn(i)+y(j)*y(j))
				enddo
			enddo
		enddo
		u(:,:,:)=0._wp
		v(:,:,:)=0._wp
		w(:,:,:)=0._wp
		zu(:,:,:)=0._wp
		zv(:,:,:)=0._wp
		zw(:,:,:)=0._wp
		tu(:,:,:)=0._wp
		tv(:,:,:)=0._wp
		tw(:,:,:)=0._wp
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		



		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! calculate and add noise														 !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
		call random_seed(size=l)
		allocate(seed(1:l))
		seed(:)=2
		call random_seed(put=seed)
		
! 		do i=0,ip+1
! 			do j=0,jp+1
! 				do k=0,kp+1
! 					r=random_normal() ! from the Netlib
! 					if((i >= ipstart) .and. (i <=ipstart+ipp+1) &
! 						.and. (j >= jpstart) .and. (j <= jpstart+jpp+1) &
! 						.and. (k >= kpstart) .and. (k <= kpstart+kpp+1) ) then
! 					
! 						if ( (z(k-kpstart)>4900._wp) .and. (z(k-kpstart)<5100._wp) ) &
! 							th(k-kpstart,j-jpstart,i-ipstart) = + r / 30._wp/10._wp
! 						
! 
! 					endif
! 										
! 
! 				enddo
! 			enddo
! 		enddo
		
! 		do i=0,ip+1
! 			do j=0,jp+1
! 				do k=0,kp+1
! 					call random_number(r)
! 					if((i >= ipstart) .and. (i <=ipstart+ipp+1) &
! 						.and. (j >= jpstart) .and. (j <= jpstart+jpp+1) &
! 						.and. (k >= kpstart) .and. (k <= kpstart+kpp+1) ) then
! 					
! 						if ( (z(k-kpstart)>900._wp) .and. (z(k-kpstart)<1100._wp) ) &
! 							th(k-kpstart,j-jpstart,i-ipstart) = (r-0.5_wp)/100._wp
! 						
! 
! 					endif
! 										
! 
! 				enddo
! 			enddo
! 		enddo

        if(param_theta) then
            do i=0,ip+1
                do j=0,jp+1
                    do k=0,kp+1
                        call random_number(r)
                        if((i >= ipstart) .and. (i <=ipstart+ipp+1) &
                            .and. (j >= jpstart) .and. (j <= jpstart+jpp+1) &
                            .and. (k >= kpstart) .and. (k <= kpstart+kpp+1) ) then
                    
                            if ( (z(k-kpstart)>param_z) .and. &
                                (z(k-kpstart)<=param_z+param_delz) ) &
                                th(k-kpstart,j-jpstart,i-ipstart) = &
                                    -0.001_wp+(r-0.5_wp)/20._wp
                                    !-0.001_wp+(r-0.5_wp)/100._wp

                            if ( (z(k-kpstart)>param_z-param_delz) .and. &
                                (z(k-kpstart)<=param_z) ) &
                                th(k-kpstart,j-jpstart,i-ipstart) = &
                                    0.001_wp-(r-0.5_wp)/20._wp
                                    !0.001_wp-(r-0.5_wp)/100._wp
                        

                        endif
                                        

                    enddo
                enddo
            enddo
        endif
! 		
!  		do i=1-r_h,ipp+r_h
!  			do j=1-r_h,jpp+r_h
!                 do k=0,kpp+1
!                     v(k,j,i)=1._wp*exp(-((zn(k)-5000._wp)/300._wp)**2)
!                     if((coords(3)==(dims(3)-1)).and.(k==kpp)) v(k,j,i)=0._wp
!                     if((coords(3)==0).and.(k<=1)) v(k,j,i)=0._wp
!                     zv(k,j,i)=v(k,j,i)
!                     tv(k,j,i)=v(k,j,i)
!                 enddo
!             enddo
!         enddo

        if(param_wind) then
            do i=1-r_h,ipp+r_h
                do j=1-r_h,jpp+r_h
                    do k=0,kpp+1
                        v(k,j,i)=(0.5_wp*erf((zn(k)-param_z)/param_sigz)+0.5_wp)*param_vmax
                        if((coords(3)==(dims(3)-1)).and.(k==kpp)) v(k,j,i)=0._wp
                        if((coords(3)==0).and.(k<=1)) v(k,j,i)=0._wp
                        zv(k,j,i)=v(k,j,i)
                        tv(k,j,i)=v(k,j,i)
                    enddo
                enddo
            enddo
        else
            do i=1-r_h,ipp+r_h
                do j=1-r_h,jpp+r_h
                    do k=0,kpp+1
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        ! locate and interpolate to find u,v on p-points:	             !
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        iloc=find_pos(z_read(1:n_levels),zn(k))
                        iloc=min(n_levels-1,iloc)
                        iloc=max(1,iloc)
                        ! linear interp u
                        call poly_int(z_read(iloc:iloc+1), u_read(iloc:iloc+1), &
                                    min(zn(k),z_read(n_levels)), var,dummy)
                        u(k,j,i)=var
                        ! linear interp u
                        call poly_int(z_read(iloc:iloc+1), v_read(iloc:iloc+1), &
                                    min(zn(k),z_read(n_levels)), var,dummy)
                        v(k,j,i)=var
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        if((coords(3)==(dims(3)-1)).and.(k==kpp)) then
                            u(k,j,i)=0._wp
                            v(k,j,i)=0._wp
                        endif
                        if((coords(3)==0).and.(k<=1)) then
                            u(k,j,i)=0._wp
                            v(k,j,i)=0._wp
                        endif
                        zu(k,j,i)=u(k,j,i)
                        tu(k,j,i)=u(k,j,i)
                        zv(k,j,i)=v(k,j,i)
                        tv(k,j,i)=v(k,j,i)
                    enddo
                enddo
            enddo            
        endif
        
        if (bubble) then
            th=0._wp
            do i=1-r_h,ipp+r_h
                do j=1-r_h,jpp+r_h
                    do k=1-r_h,kpp+r_h
                        !th(k,j,i)=theta(k)
                
                        rad = (zn(k)-500._wp)**2._wp
                        
                        if (ip > 1) rad=rad+xn(i)**2._wp
                        if (jp > 1) rad=rad+yn(j)**2._wp
                    
                        rad=sqrt(rad)
                        if(rad<=300._wp) then
                            th(k,j,i)=th(k,j,i)+2.5_wp
                            q(k,j,i,1)=q(k,j,i,1)*3.5_wp
                            !th(k,j,i)=th(k,j,i)+0.1_wp
                        else
                            !th(k,j,i)=0._wp
                        endif
                    enddo
                enddo
            enddo
        endif
		deallocate(seed)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		


        if(forcing) then
            do k=0,kpp+1
                u_force(k)=u(k,1,1)
                v_force(k)=v(k,1,1)
            enddo        
        endif


        if(divergence) then
            do k=0,kpp+1
                if(z(k).le.divergence_hgt) then
                    w_subs(k)=-divergence_val*z(k)
                else
                    w_subs(k)=-divergence_val*divergence_hgt
                endif
            enddo        
        endif



		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! set halos																		 !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
		call exchange_full(comm3d, id, kpp, jpp, ipp, r_h,r_h,r_h,r_h,l_h,r_h, u,&
		        0._wp,0._wp, dims,coords)
		call exchange_full(comm3d, id, kpp, jpp, ipp, r_h,r_h,l_h,r_h,r_h,r_h, v,&
		        0._wp,0._wp, dims,coords)
		call exchange_full(comm3d, id, kpp, jpp, ipp, l_h,r_h,r_h,r_h,r_h,r_h, w,&
		        0._wp,0._wp, dims,coords)
		call exchange_full(comm3d, id, kpp, jpp, ipp, r_h,r_h,r_h,r_h,r_h,r_h, th,&
		        0._wp,0._wp, dims,coords)
		call exchange_full(comm3d, id, kpp, jpp, ipp, r_h,r_h,r_h,r_h,r_h,r_h, q(:,:,:,1),&
		        0._wp,0._wp, dims,coords)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! find top and base of array                                                     !                                           
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!              
        call find_base_top(comm3d, id, kpp,l_h,r_h,thetan,thbase,thtop, dims,coords)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!              
        
        ! deallocate
        deallocate(qv1d)
        deallocate(qc1d)

		! no dangling pointers
		if (associated(zr1) ) nullify(zr1)
		if (associated(thr1) ) nullify(thr1)

	end subroutine allocate_and_set
	


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Setup the profile to be adiabatic                                                  !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>adiabatic profile:                                        
    subroutine setup_adiabatic_profile(l_h,r_h,kpp, psurf,tsurf,t_cbase,t_ctop, &
                                        rho_surf, rh_above, th_jump, th_grad,&
                                        adiabatic_frac, theta_flag, &
                                        pref,tref,theta,rhoa,qv1d,qc1d,z,dz,&
                                        prefn,trefn,thetan,rhoan,zn,dzn, &
                                        coords,dims, id, comm3d)                      
        use numerics_type
		use numerics, only : vode_integrate, zeroin
        implicit none
        logical, intent(in) :: theta_flag
        integer(i4b), intent(in) :: l_h,r_h,kpp
        real(wp), intent(in) :: psurf,tsurf,t_cbase,t_ctop,rh_above, th_jump,th_grad, &
                            adiabatic_frac
        real(wp), intent(inout) :: rho_surf
        real(wp), dimension(1-l_h:kpp+r_h), intent(in) :: z,zn,dz,dzn
        real(wp), dimension(1-l_h:kpp+r_h), intent(inout) :: pref,tref,theta,rhoa,qv1d, &
                                                    qc1d,prefn,trefn,thetan,rhoan
		integer(i4b), dimension(3), intent(inout) :: coords
		integer(i4b), dimension(3), intent(in) :: dims
		integer(i4b), intent(in) :: id, comm3d



                                                    
        real(wp) :: zcb, htry, hmin, eps2=1.e-5_wp,t,pcb,rv_sub, theta_q_sat, lmr_ctop, &
                q_tot,pct,zct,t_ctop_evap,pctop_dry
        real(wp),dimension(1) :: psolve, zsolve
        integer(i4b) :: k                        
        
                   
        ! 1. Solve hydrostatic equation up to cloud-base, 
        !    conserving dry potential temperature
        ! 2. Calculate the cloud-base pressure,
        !     by integrating hydrostatic equation
        ! 3. Solve hydrostatic equation up to cloud-top, 
        !    conserving moist potential temperature
        ! 4. Calculate the cloud-top pressure
        ! 5. Calculate the cloud-top pressure
        
        
        
        
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! 1. Solve hydrostatic equation up to cloud-base, conserving dry potential       !
        !    temperature                                                                 !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        zcb=(tsurf-t_cbase)*cp/grav   
        qv1d(:)=0._wp
        qc1d(:)=0._wp
        tsurf_glob=tsurf
        do k=1-l_h,kpp+r_h
            if(z(k).gt.0._wp) then
                hmin=1.e-2_wp
                htry=dz(k)
            else
                hmin=1.e-2_wp
                htry=-dz(k)
            endif
            psolve(1)=psurf
            call vode_integrate(psolve,0._wp,z(k),eps2,htry,hmin,hydrostatic1b)

            t=tsurf*(psolve(1)/1.e5_wp)**(ra/cp)
            if(t.lt.t_cbase) exit
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! set reference variables                                                    !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            pref(k)=psolve(1)
            tref(k)=t
            theta(k)=tsurf
            rhoa(k)=pref(k)/(ra*tref(k))
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            
        enddo
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! 1. done                                                                        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        
        
        
        
        
        
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! 2. calculate the cloud-base pressure                                           !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        psolve(1)=psurf
        call vode_integrate(psolve,0._wp,zcb,eps2,htry,hmin,hydrostatic1b)
        t=tsurf*(psolve(1)/1.e5_wp)**(ra/cp)
        pcb=psolve(1)
        rv_sub=eps1*svp_liq(t_cbase)/(pcb-svp_liq(t_cbase))
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! 2. done                                                                           !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        
        
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! 3. Solve hydrostatic equation up to cloud-top, conserving moist potential      !
        !    temperature                                                                 !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        hmin=1.e-2_wp
        adiabatic_frac_glob=adiabatic_frac
        rv_glob=rv_sub
        if(.not.theta_flag) then
            adiabatic_frac_glob2=1._wp
        else
            adiabatic_frac_glob2=0._wp
            rv_glob=0._wp
        endif
        p_glob=pcb
        theta_q_sat_glob=0._wp
        theta_q_sat=calc_theta_q(t_cbase)        
        theta_q_sat_glob=theta_q_sat
        p_glob=pcb
        t=zeroin(1.01_wp*t_cbase,t_cbase,calc_theta_q,1.e-5_wp)
        t_glob=t_cbase

        do k=1-l_h,kpp+r_h
            if(zcb.ge.z(k)) then 
                qv1d(k)=rv_sub
                cycle
            endif
            psolve(1)=pcb
            htry=dz(k)
            call vode_integrate(psolve,zcb,z(k),eps2,htry,hmin,hydrostatic2b)
            p_glob=psolve(1)
            t=zeroin(20._wp,t_cbase,calc_theta_q,1.e-5_wp)
            if(t.lt.t_ctop) exit
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! set reference variables                                                    !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            pref(k)=psolve(1)
            tref(k)=t
            theta(k)=t*(1.e5_wp/psolve(1))**(ra/cp)
            rhoa(k)=pref(k)/(ra*tref(k))
            qv1d(k)=eps1*svp_liq(t)/(psolve(1)-svp_liq(t))
            qc1d(k)=(rv_sub-qv1d(k))*adiabatic_frac
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        enddo
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! 3. done                                                                        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
      
      

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! 4. calculate the cloud-top pressure                                            !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        t_glob=t_ctop
	    pctop_dry=pcb*exp(cp/ra*log(t_ctop/t_cbase))
        pct=zeroin(pcb,pctop_dry/2._wp,calc_theta_q2,1.e-5_wp)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! 4. done                                                                        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
         
        
        
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! 5. calculate the cloud-top height                                              !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        hmin=1.e-2_wp
        htry=-10._wp
        zsolve(1)=zcb
        t_glob=t_cbase
        call vode_integrate(zsolve,pcb,pct,eps2,htry,hmin,hydrostatic2c)
        zct=zsolve(1)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! 5. done                                                                        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! 6. Solve hydrostatic equation up to domain top, with theta gradient and        !
        ! constant rh                                                                    !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        z_glob=zct
        t_glob=(t_ctop+th_jump)*(1.e5_wp/pct)**(ra/cp) ! potential temperature at cloud top
        hmin=1.e-2_wp
        th_grad_glob=th_grad
        do k=1-l_h,kpp+r_h
            if(zct.ge.z(k)) cycle
            htry=dz(k)
            psolve(1)=pct
            call vode_integrate(psolve,zct,z(k),eps2,htry,hmin,hydrostatic2e)


            t=(t_glob+th_grad_glob*(z(k)-zct))*(psolve(1)/1.e5_wp)**(ra/cp)
            
            rv_glob=rh_above*eps1*svp_liq(t)/(psolve(1)-svp_liq(t))


            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! set reference variables                                                    !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            pref(k)=psolve(1)
            tref(k)=t
            theta(k)=t*(1.e5_wp/psolve(1))**(ra/cp)
            rhoa(k)=pref(k)/(ra*tref(k))
            qv1d(k)=rh_above*eps1*svp_liq(t)/(psolve(1)-svp_liq(t))
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        enddo
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! 6. done                                                                        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                             
                                             
                                             
                                             
                                             
                       
                       
                       
                       
                                             
                                             
                                              
        ! now do the same, but for the node values
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! 1. Solve hydrostatic equation up to cloud-base, conserving dry potential       !
        !    temperature                                                                 !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        zcb=(tsurf-t_cbase)*cp/grav   
        qv1d(:)=0._wp
        qc1d(:)=0._wp
        tsurf_glob=tsurf
        do k=1-l_h,kpp+r_h
            if(zn(k).gt.0._wp) then
                hmin=1.e-2_wp
                htry=dzn(k)
            else
                hmin=1.e-2_wp
                htry=-dzn(k)
            endif
            psolve(1)=psurf
            call vode_integrate(psolve,0._wp,zn(k),eps2,htry,hmin,hydrostatic1b)

            t=tsurf*(psolve(1)/1.e5_wp)**(ra/cp)
            if(t.lt.t_cbase) exit
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! set reference variables                                                    !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            prefn(k)=psolve(1)
            trefn(k)=t
            thetan(k)=tsurf
            rhoan(k)=prefn(k)/(ra*trefn(k))
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            
        enddo
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! 1. done                                                                        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        
        
        
        
        
        
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! 2. calculate the cloud-base pressure                                           !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        psolve(1)=psurf
        call vode_integrate(psolve,0._wp,zcb,eps2,htry,hmin,hydrostatic1b)
        t=tsurf*(psolve(1)/1.e5_wp)**(ra/cp)
        pcb=psolve(1)
        rv_sub=eps1*svp_liq(t_cbase)/(pcb-svp_liq(t_cbase))
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! 2. done                                                                           !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        
        
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! 3. Solve hydrostatic equation up to cloud-top, conserving moist potential      !
        !    temperature                                                                 !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        hmin=1.e-2_wp
        rv_glob=rv_sub
        adiabatic_frac_glob=adiabatic_frac
        rv_glob=rv_sub
        if(theta_flag) then
            adiabatic_frac_glob2=1._wp
        else
            adiabatic_frac_glob2=0._wp
            rv_glob=0._wp
        endif
        p_glob=pcb
        theta_q_sat_glob=0._wp
        theta_q_sat=calc_theta_q(t_cbase)        
        theta_q_sat_glob=theta_q_sat
        p_glob=pcb
        t=zeroin(1.01_wp*t_cbase,t_cbase,calc_theta_q,1.e-5_wp)
        t_glob=t_cbase
        do k=1-l_h,kpp+r_h
            if(zcb.ge.zn(k)) then 
                qv1d(k)=rv_sub
                cycle
            endif
            psolve(1)=pcb
            htry=dzn(k)
            call vode_integrate(psolve,zcb,zn(k),eps2,htry,hmin,hydrostatic2b)
            p_glob=psolve(1)
            t=zeroin(20._wp,t_cbase,calc_theta_q,1.e-5_wp)
            if(t.lt.t_ctop) exit
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! set reference variables                                                    !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            prefn(k)=psolve(1)
            trefn(k)=t
            thetan(k)=t*(1.e5_wp/psolve(1))**(ra/cp)
            rhoan(k)=pref(k)/(ra*tref(k))
            qv1d(k)=eps1*svp_liq(t)/(psolve(1)-svp_liq(t))
            qc1d(k)=(rv_sub-qv1d(k))*adiabatic_frac
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        enddo
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! 3. done                                                                        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     
     
     
     
     
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! 4. calculate the cloud-top pressure                                            !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        t_glob=t_ctop
	    pctop_dry=pcb*exp(cp/ra*log(t_ctop/t_cbase))
        pct=zeroin(pcb,pctop_dry/2._wp,calc_theta_q2,1.e-5_wp)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! 4. done                                                                        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        
        
        
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! 5. calculate the cloud-top height                                              !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        hmin=1.e-2_wp
        htry=-10._wp
        zsolve(1)=zcb
        t_glob=t_cbase
        call vode_integrate(zsolve,pcb,pct,eps2,htry,hmin,hydrostatic2c)
        zct=zsolve(1)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! 5. done                                                                        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! 6. Solve hydrostatic equation up to domain top, with theta gradient and        !
        ! constant rh                                                                    !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        z_glob=zct
        t_glob=(t_ctop+th_jump)*(1.e5_wp/pct)**(ra/cp) ! potential temperature at cloud top
        hmin=1.e-2_wp
        th_grad_glob=th_grad
        do k=1-l_h,kpp+r_h
            if(zct.ge.zn(k)) cycle
            htry=dzn(k)
            psolve(1)=pct
            call vode_integrate(psolve,zct,zn(k),eps2,htry,hmin,hydrostatic2e)


            t=(t_glob+th_grad_glob*(zn(k)-zct))*(psolve(1)/1.e5_wp)**(ra/cp)
            
            rv_glob=rh_above*eps1*svp_liq(t)/(psolve(1)-svp_liq(t))


            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! set reference variables                                                    !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            prefn(k)=psolve(1)
            trefn(k)=t
            thetan(k)=t*(1.e5_wp/psolve(1))**(ra/cp)
            rhoan(k)=prefn(k)/(ra*trefn(k))
            qv1d(k)=rh_above*eps1*svp_liq(t)/(psolve(1)-svp_liq(t))
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        enddo
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! 6. done                                                                        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        
    end subroutine setup_adiabatic_profile
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! calculate hydrostatic equation                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculate the hydrostatic equation (version for th gradient):
	!>\f$ \frac{\partial P}{\partial z} = -\rho g \f$
	!>@param[in] z,p
	!>@param[inout] dpdz
	subroutine hydrostatic2e(z,p,dpdz)
        use numerics_type
        use numerics, only : zeroin
        implicit none
        real(wp), intent(in) :: z
        real(wp), dimension(:), intent(in) :: p
        real(wp), dimension(:), intent(out) :: dpdz
        real(wp) :: t,theta
    
        theta=t_glob+th_grad_glob*(z-z_glob)
	
		
		! calculate t from theta and p:
		t=theta*(p(1)/1.e5_wp)**(ra/cp)
		! hydrostatic equation:
		dpdz(1)=-(grav*p(1)) / (ra*t) 	
	end subroutine hydrostatic2e
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! calculate hydrostatic equation                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculate the hydrostatic equation (version for constant rh):
	!>\f$ \frac{\partial P}{\partial z} = -\rho g \f$
	!>@param[in] z,p
	!>@param[inout] dpdz
	subroutine hydrostatic2d(z,p,dpdz)
        use numerics_type
        use numerics, only : zeroin
        implicit none
        real(wp), intent(in) :: z
        real(wp), dimension(:), intent(in) :: p
        real(wp), dimension(:), intent(out) :: dpdz
        real(wp) :: t
    
        p_glob=p(1)

        t=zeroin(20._wp,t_glob*1.5_wp,rh_specified,1.e-5_wp)

        dpdz(1)=-(grav*p(1))/(ra*t)
	
	end subroutine hydrostatic2d
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! calculate hydrostatic equation                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculate the hydrostatic equation (version for constant theta_q_sat):
	!>\f$ \frac{\partial P}{\partial z} = -\rho g \f$
	!>@param[in] z,p
	!>@param[inout] dpdz
	subroutine hydrostatic2c(p,z,dzdp)
        use numerics_type
        use numerics, only : zeroin
        implicit none
        real(wp), intent(in) :: p
        real(wp), dimension(:), intent(in) :: z
        real(wp), dimension(:), intent(out) :: dzdp
        real(wp) :: t
    
        p_glob=p ! actual pressure, t_glob is cloud base temp

        t=tsurf_glob*(p_glob/1.e5_wp)**(ra/cp)
        ! find the temperature by iteration
        t=zeroin(20._wp,t_glob,calc_theta_q,1.e-5_wp)

        dzdp(1)=-(ra*t)/(grav*p)
	
	end subroutine hydrostatic2c
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! calculate hydrostatic equation                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculate the hydrostatic equation (version for constant theta_q_sat):
	!>\f$ \frac{\partial P}{\partial z} = -\rho g \f$
	!>@param[in] z,p
	!>@param[inout] dpdz
	subroutine hydrostatic2b(z,p,dpdz)
        use numerics_type
        use numerics, only : zeroin
        implicit none
        real(wp), intent(in) :: z
        real(wp), dimension(:), intent(in) :: p
        real(wp), dimension(:), intent(out) :: dpdz
        real(wp) :: t
    
        p_glob=p(1)

        t=tsurf_glob*(p_glob/1.e5_wp)**(ra/cp)
        ! find the temperature by iteration
        t=zeroin(20._wp,t_glob,calc_theta_q,1.e-5_wp)

        dpdz(1)=-(grav*p(1))/(ra*t)
	
	end subroutine hydrostatic2b
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! calculate hydrostatic equation                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculate the hydrostatic equation (version for constant theta):
	!>\f$ \frac{\partial P}{\partial z} = -\rho g \f$
	!>@param[in] z,p
	!>@param[inout] dpdz
	subroutine hydrostatic1b(z,p,dpdz)
		use numerics_type
		implicit none
		real(wp), intent(in) :: z
		real(wp), dimension(:), intent(in) :: p
		real(wp), dimension(:), intent(out) :: dpdz
		real(wp) :: t
		integer(i4b) :: iloc
		
		
		! calculate t from theta and p:
		t=tsurf_glob*(p(1)/1.e5_wp)**(ra/cp)
		! hydrostatic equation:
		dpdz(1)=-(grav*p(1)) / (ra*t) 
	
	end subroutine hydrostatic1b
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! calculate hydrostatic equation                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculate the hydrostatic equation:
	!>\f$ \frac{\partial P}{\partial z} = -\rho g \f$
	!>@param[in] z,p
	!>@param[inout] dpdz
	subroutine hydrostatic1a(z,p,dpdz)
		use numerics_type
		use variables, only : nm1
		use numerics, only : find_pos, poly_int
		implicit none
		real(wp), intent(in) :: z
		real(wp), dimension(:), intent(in) :: p
		real(wp), dimension(:), intent(out) :: dpdz
		real(wp) :: t,theta, var, dummy
		integer(i4b) :: iloc
		
		! locate and interpolate to find theta:
		iloc=find_pos(zr1(1:nl1),z)
		iloc=min(nl1-1,iloc)
		iloc=max(1,iloc)
		! linear interp theta
		call poly_int(zr1(iloc:iloc+1), thr1(iloc:iloc+1), &
					min(z,zr1(nl1)), var,dummy)
		theta=var
		
		! calculate t from theta and p:
		t=theta*(p(1)/1.e5_wp)**(ra/cp)
		! hydrostatic equation:
		dpdz(1)=-(grav*p(1)) / (ra*t) 
	
	end subroutine hydrostatic1a
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! calculate the value of temperature required to reach specified RH                  !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the value of temperature required to evaporate 
	!> liquid water and reach specified rh (note!!! p_glob, rv_glob 
	!>      and rh_glob must be set externally before calling)
	!>@param[in] t: temperature
	!>@return rh_specified: value of rh_specified
	function rh_specified(t)
        use numerics_type

        implicit none
        real(wp), intent(in) :: t
        real(wp) :: rh_specified
        real(wp) :: ws
        
    
        ws=eps1*svp_liq(t)/(p_glob-svp_liq(t))
        
        rh_specified= rv_glob/(ws)-rh_glob

	end function rh_specified
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! calculate the value of theta_q_sat                                                 !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the value of theta_q_sat (note!!! t_glob, rv_glob 
	!>      and theta_q_sat_glob must be set externally before calling)
	!>@param[in] t: temperature
	!>@return calc_theta_q2: value of theta_q_sat
	function calc_theta_q2(p)
        use numerics_type

        implicit none
        real(wp), intent(in) :: p
        real(wp) :: calc_theta_q2
        real(wp) :: ws
        
        ws=eps1*svp_liq(t_glob)/(p-svp_liq(t_glob))*adiabatic_frac_glob2
        calc_theta_q2=t_glob*(1e5_wp/p)**(ra/(cp+cw*rv_glob))* &
            exp(lv*ws/(cp+cw*rv_glob)/t_glob)-theta_q_sat_glob

	end function calc_theta_q2  
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! calculate the value of theta_q_sat                                                 !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the value of theta_q_sat (note!!! p_glob, rv_glob 
	!>      and theta_q_sat_glob must be set externally before calling)
	!>@param[in] t: temperature
	!>@return calc_theta_q: value of theta_q_sat
	function calc_theta_q(t)
        use numerics_type
        implicit none
        real(wp), intent(in) :: t
        real(wp) :: calc_theta_q
        real(wp) :: ws
        
        
        ws=eps1*svp_liq(t)/(p_glob-svp_liq(t))*adiabatic_frac_glob2
        calc_theta_q=t*(1.e5_wp/p_glob)**(ra/(cp+cw*rv_glob))* &
            exp(lv*ws/(t*(cp+cw*rv_glob)))-theta_q_sat_glob

	end function calc_theta_q    
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! saturation vapour pressure over liquid                                             !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the saturation vapour pressure over liquid water according to buck fit
	!>@param[in] t: temperature
	!>@return svp_liq: saturation vapour pressure over liquid water
	function svp_liq(t)
		use numerics_type
		implicit none
		real(wp), intent(in) :: t
		real(wp) :: svp_liq
		svp_liq = 100._wp*6.1121_wp* &
			  exp((18.678_wp - (t-ttr)/ 234.5_wp)* &
			  (t-ttr)/(257.14_wp + (t-ttr)))
	end function svp_liq


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! HELPER ROUTINE                                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine check(status)
	use netcdf
	use numerics_type
	integer(i4b), intent ( in) :: status

	if(status /= nf90_noerr) then
		print *, trim(nf90_strerror(status))
		stop "Stopped"
	end if
	end subroutine check
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	end module initialisation
	
