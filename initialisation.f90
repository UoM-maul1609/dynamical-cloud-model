	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>initialisation for the dynamical cloud model
    module initialisation
    use nrtype
    implicit none
    
    real(sp), parameter :: ra=287._sp, grav=9.81_sp,cp=1005._sp
	real(sp), dimension(:), pointer :: zr1, thr1
	integer(i4b) :: nl1
	
    private
    public :: allocate_and_set
    contains
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Allocate and set arrays                                                            !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>allocate arrays on each PE, and initialise them
	!>@param[inout] dt,runtime,ntim - time variables
	!>@param[inout] x,y,z,xn,yn,zn,u,v,w,p,th,rho - grid positions and prognostics
	!>@param[inout] ubar,vbar,wbar,thbar,qbar,dampfacn,dampfac
	!>@param[in] damping_layer, damping_thickness, damping_tau
	!>@param[inout] sth,strain,vism,vist - more prognostics
	!>@param[in] z0,z0th - roughness lengths 
	!>@param[inout] q,sq,viss - q-variable
	!>@param[inout] zu,zv,zw,tu,tv,tw - previous time-step and temporary storage
	!>@param[inout] su,sv,sw,psrc - grid positions and prognostics
	!>@param[inout] theta,thetan - reference potential temperatures
	!>@param[inout] rhoa,rhoan - reference potential temperatures
	!>@param[inout] lamsq,lamsqn - mixing length
	!>@param[in] cvis - smagorinsky parameter
	!>@param[inout] dx,dy,dz - grid spacing on grid
	!>@param[inout] dxn,dyn,dzn - grid spacing on grid - staggered
	!>@param[inout] ipp,jpp,kpp,ipstart,jpstart,kpstart - number of grid / starting position
	!>@param[in] dx_nm,dy_nm,dz_nm - grid spacing from namelist
	!>@param[in] ip,jp,kp - grid points from namelist
	!>@param[in] moisture, nq - whether to have clouds and number of q-variables
	!>@param[in] n_levels,z_read,theta_read, psurf,tsurf - sounding variables
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
			u,v,w,&
			zu,zv,zw,&
			tu,tv,tw,&
			p,th,rho, &
			sth,strain,vism,vist,z0,z0th, &
			q, sq, viss, &
			su,sv,sw,psrc, &
			theta,thetan, &
			rhoa,rhoan, &
			lamsq,lamsqn, &
			cvis, &
			dx, dy, dz, &
			dxn, dyn, dzn, &
			ipp, jpp, kpp,&
			ipstart, jpstart, kpstart, &
			dx_nm, dy_nm, dz_nm, &
			ip, jp, kp, &
			moisture,nq,&
			n_levels,z_read, theta_read,psurf,tsurf, &
			l_h,r_h, &
			thbase,thtop, &
			coords,dims, id, comm3d)
				
		use nrtype
		use mpi
		use netcdf
		use nr, only : locate, polint, rkqs, odeint
		use random, only : random_normal
		use mpi_module
		
		implicit none
		real(sp), dimension(:,:,:), allocatable, intent(inout) :: &
														u,v,w,zu,zv,zw,tu,tv,tw,&
														p,th,rho, &
														su,sv,sw,psrc, &
														sth,strain,vism,vist
		real(sp), dimension(:,:,:,:), allocatable, intent(inout) :: &
														q,sq,viss
		real(sp), dimension(:), allocatable, intent(inout) :: x,y,z,xn,yn,zn,dx,dy,dz, &
															dxn,dyn,dzn, theta,thetan, &
															rhoa, rhoan, lamsq, lamsqn, &
															ubar, vbar, wbar,thbar, &
															dampfac,dampfacn
															
		real(sp), dimension(:,:), allocatable, intent(inout) :: qbar

		real(sp), intent(in) :: dx_nm, dy_nm, dz_nm, cvis,z0,z0th
		real(sp), intent(in) :: dt, runtime
		integer(i4b), intent(inout) :: ipp, jpp, kpp, ipstart, jpstart, kpstart
		integer(i4b), intent(inout) :: ntim
		integer(i4b), intent(in) :: ip, jp, kp, l_h, r_h, nq
		logical :: moisture,damping_layer
		integer(i4b), intent(in) :: n_levels
		real(sp), dimension(n_levels), target, intent(in) :: z_read,theta_read
		real(sp), intent(in) :: psurf, tsurf,damping_thickness,damping_tau
		integer(i4b), dimension(3), intent(inout) :: coords
		integer(i4b), dimension(3), intent(in) :: dims
		integer(i4b), intent(in) :: id, comm3d
		real(sp), intent(inout) :: thbase, thtop
		
		! locals:
		integer(i4b) :: error, AllocateStatus,i,j,k
		real(sp) :: rho_surf, htry, hmin, eps2=1.e-5_sp
		real(sp), dimension(1) :: psolve
		real(sp) :: var, dummy,ztop
		integer(i4b) :: iloc
		! for random number:
		real(sp) :: r
		real(sp), dimension(10,10) :: rs
		integer(i4b) :: l, nbottom, ntop, tag1
		integer(i4b), allocatable, dimension(:) :: seed
		real(sp) :: rad
		
! if the pe is not being used in the cartesian topology, do not use here
		if(id>=dims(1)*dims(2)*dims(3)) return 
		
		
		! set pointers
		zr1 => z_read
		thr1 => theta_read
		nl1 = n_levels

		
		! scalar formulae:
		ntim=ceiling(runtime/dt)

		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! find the number of grid points on each PE                                      !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		call MPI_CART_COORDS(comm3d, id, 3, coords, error)
		! print *,'Coords of ',id,' are ',coords

        ! number of grid points in all but last:
		ipp = floor(real(ip,sp)/real(dims(1),sp)) 
		ipstart = ipp*(coords(1))  +1   
		if(coords(1) == (dims(1)-1)) then
			ipp=ip-(dims(1)-1)*ipp ! number of grid points in last
		endif
		! number of grid points in all but last:
		jpp = floor(real(jp,sp)/real(dims(2),sp))      
		jpstart = jpp*(coords(2))  +1 
		if(coords(2) == (dims(2)-1)) then
			jpp=jp-(dims(2)-1)*jpp ! number of grid points in last
		endif
		! number of grid points in all but last:
		kpp = floor(real(kp,sp)/real(dims(3),sp))      
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
		allocate( rho(1-r_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h), STAT = AllocateStatus)
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
		if (moisture) then
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
		endif	
		
		allocate( ubar(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( vbar(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( wbar(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( thbar(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		
		if(damping_layer) then
            allocate( dampfacn(1-l_h:kpp+r_h), STAT = AllocateStatus)
            if (AllocateStatus /= 0) STOP "*** Not enough memory ***"		
            allocate( dampfac(1-l_h:kpp+r_h), STAT = AllocateStatus)
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
		x=dx_nm*(/(i,i=-l_h+ipstart,ipp+r_h+ipstart-1)/) - real(ip+1,sp)/2._sp*dx_nm
		xn=x-0.5_sp*dx_nm
		! set up horizontal level array
		y=dy_nm*(/(i,i=-l_h+jpstart,jpp+r_h+jpstart-1)/) - real(jp+1,sp)/2._sp*dy_nm
		yn=y-0.5_sp*dy_nm

		! set up vertical level array
		z=dz_nm*(/(i,i=-l_h+kpstart-1,kpp+r_h+kpstart-2)/)+0.5_sp*dz_nm
		zn=z-0.5_sp*dz_nm
		
		! set up mixing length array
		lamsq=1._sp / (1._sp/(cvis*(dx_nm+dy_nm+dz)/3._sp)**2._sp + &
				1._sp/(0.4_sp*(z + z0))**2._sp)
		lamsqn=1._sp / (1._sp/(cvis*(dx_nm+dy_nm+dzn)/3._sp)**2._sp + &
				1._sp/(0.4_sp*(zn + z0))**2._sp)	
        ! zero some arrays
		ubar=0._sp	
		vbar=0._sp	
		wbar=0._sp	
		thbar=0._sp	
		qbar=0._sp	
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! set up damping layer                                                           !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(damping_layer) then
            dampfacn=0._sp
            dampfac=0._sp
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
                    dampfacn(k)=-1._sp/damping_tau*(&
                        exp((zn(k)-(ztop-damping_thickness))/damping_thickness)-1._sp)
                endif

                if(z(k) .gt. (ztop-damping_thickness)) then
                    dampfac(k)=-1._sp/damping_tau*(&
                        exp((z(k)-(ztop-damping_thickness))/damping_thickness)-1._sp)
                endif
            enddo
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!              


        endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! set-up density, pressure, theta, etc                                           !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		rho_surf=psurf / (ra*tsurf)
		do i=1-l_h,kpp+r_h
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! solve the hydrostatic equation											 !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			psolve=psurf
			if( z(i) < 0._sp ) then
				hmin=-1.e-2_sp
				htry=-dz(i)
				call odeint(psolve,0._sp,z(i),eps2,htry,hmin,hydrostatic1a,rkqs)
			else
				hmin=1.e-2_sp
				htry=dz(i)
				call odeint(psolve,0._sp,z(i),eps2,htry,hmin,hydrostatic1a,rkqs)
			endif
			p(i,:,:)=psolve(1)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! locate and interpolate to find theta:										 !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			iloc=locate(z_read(1:n_levels),z(i))
			iloc=min(n_levels-1,iloc)
			iloc=max(1,iloc)
			! linear interp theta
			call polint(z_read(iloc:iloc+1), theta_read(iloc:iloc+1), &
						min(z(i),z_read(n_levels)), var,dummy)
!			th(:,:,i)=var
			theta(i)=var
			rhoa(i)=psolve(1)/(ra*theta(i)*(psolve(1)/psurf)**(ra/cp))
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! solve the hydrostatic equation for k+1/2 levels							 !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			psolve=psurf
			if( zn(i) < 0._sp ) then
				hmin=-1.e-2_sp
				htry=-dzn(i)
				call odeint(psolve,0._sp,zn(i),eps2,htry,hmin,hydrostatic1a,rkqs)
			else
				hmin=1.e-2_sp
				htry=dzn(i)
				call odeint(psolve,0._sp,zn(i),eps2,htry,hmin,hydrostatic1a,rkqs)
			endif
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! locate and interpolate to find theta:										 !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			iloc=locate(z_read(1:n_levels),zn(i))
			iloc=min(n_levels-1,iloc)
			iloc=max(1,iloc)
			! linear interp theta
			call polint(z_read(iloc:iloc+1), theta_read(iloc:iloc+1), &
						min(zn(i),z_read(n_levels)), var,dummy)
			thetan(i)=var
			rhoan(i)=psolve(1)/(ra*thetan(i)*(psolve(1)/psurf)**(ra/cp))
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		enddo
		p(:,:,:)=0._sp
! 		rhoa=1._sp
! 		rhoan=1._sp
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! set wind field        														 !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
		do i=1,ipp
			do j=1,jpp
				do k=1,kpp
					u(k,j,i)=-5._sp*(yn(j))/sqrt(x(i)*x(i)+yn(j)*yn(j))
					v(k,j,i)=5._sp*(xn(i))/sqrt(xn(i)*xn(i)+y(j)*y(j))
				enddo
			enddo
		enddo
		u(:,:,:)=0._sp
		v(:,:,:)=0._sp
		w(:,:,:)=0._sp
		zu(:,:,:)=0._sp
		zv(:,:,:)=0._sp
		zw(:,:,:)=0._sp
		tu(:,:,:)=0._sp
		tv(:,:,:)=0._sp
		tw(:,:,:)=0._sp
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		



		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! calculate and add noise														 !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
		call random_seed(size=l)
		allocate(seed(1:l))
		seed(:)=2
		call random_seed(put=seed)
		
		do i=0,ip+1
			do j=0,jp+1
				do k=0,kp+1
					r=random_normal() ! from the Netlib
					if((i >= ipstart) .and. (i <=ipstart+ipp+1) &
						.and. (j >= jpstart) .and. (j <= jpstart+jpp+1) &
						.and. (k >= kpstart) .and. (k <= kpstart+kpp+1) ) then
					
						if ( (z(k-kpstart)>3000._sp) .and. (z(k-kpstart)<6000._sp) ) &
							th(k-kpstart,j-jpstart,i-ipstart) = + r / 30._sp
						

					endif
					
					
					

				enddo
			enddo
		enddo

 		th=0._sp
		
		
		do i=1-r_h,ipp+r_h
			do j=1-r_h,jpp+r_h
				do k=1-r_h,kpp+r_h
					!th(k,j,i)=theta(k)
				
					rad = (z(k)-3000._sp)**2._sp
						
					if (ip > 1) rad=rad+xn(i)**2._sp
					if (jp > 1) rad=rad+yn(j)**2._sp
					
					rad=sqrt(rad)
					if(rad<=1000._sp) then
						th(k,j,i)=th(k,j,i)+0.1_sp
					else
						!th(k,j,i)=0._sp
					endif
				enddo
			enddo
		enddo
		deallocate(seed)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		



		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! set halos																		 !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
		call exchange_full(comm3d, id, kpp, jpp, ipp, r_h,r_h,r_h,r_h,l_h,r_h, u,&
		        0._sp,0._sp, dims,coords)
		call exchange_full(comm3d, id, kpp, jpp, ipp, r_h,r_h,l_h,r_h,r_h,r_h, v,&
		        0._sp,0._sp, dims,coords)
		call exchange_full(comm3d, id, kpp, jpp, ipp, l_h,r_h,r_h,r_h,r_h,r_h, w,&
		        0._sp,0._sp, dims,coords)
		call exchange_full(comm3d, id, kpp, jpp, ipp, r_h,r_h,r_h,r_h,r_h,r_h, th,&
		        0._sp,0._sp, dims,coords)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! find top and base of array                                                     !                                           
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!              
        call find_base_top(comm3d, id, kpp,l_h,r_h,thetan,thbase,thtop, dims,coords)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!              
        


		! no dangling pointers
		if (associated(zr1) ) nullify(zr1)
		if (associated(thr1) ) nullify(thr1)

	end subroutine allocate_and_set
	




	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! calculate hydrostatic equation                                                            !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculate the hydrostatic equation:
	!>\f$ \frac{\partial P}{\partial z} = -\rho g \f$
	!>@param[in] z,p
	!>@param[inout] dpdz
	subroutine hydrostatic1a(z,p,dpdz)
		use nrtype
		use variables, only : nm1
		use nr, only : locate, polint
		implicit none
		real(sp), intent(in) :: z
		real(sp), dimension(:), intent(in) :: p
		real(sp), dimension(:), intent(out) :: dpdz
		real(sp) :: t,theta, var, dummy
		integer(i4b) :: iloc
		
		! locate and interpolate to find theta:
		iloc=locate(zr1(1:nl1),z)
		iloc=min(nl1-1,iloc)
		iloc=max(1,iloc)
		! linear interp theta
		call polint(zr1(iloc:iloc+1), thr1(iloc:iloc+1), &
					min(z,zr1(nl1)), var,dummy)
		theta=var
		
		! calculate t from theta and p:
		t=theta*(p(1)/1.e5_sp)**(ra/cp)
		! hydrostatic equation:
		dpdz(1)=-(grav*p(1)) / (ra*t) 
	
	end subroutine hydrostatic1a
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! HELPER ROUTINE                                                       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine check(status)
	use netcdf
	use nrtype
	integer(i4b), intent ( in) :: status

	if(status /= nf90_noerr) then
		print *, trim(nf90_strerror(status))
		stop "Stopped"
	end if
	end subroutine check
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	end module initialisation
	