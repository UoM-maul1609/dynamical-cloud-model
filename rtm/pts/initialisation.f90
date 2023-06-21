	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>initialisation for the dynamical cloud model
    module initialisation
    use numerics_type
    implicit none
    
    real(wp), parameter :: ra=287._wp, grav=9.81_wp,cp=1005._wp
	real(wp), dimension(:), pointer :: zr1, thr1
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
	!>@param[in] n_levels,z_read,theta_read, psurf,tsurf - sounding variables
	!>@param[in] l_h,r_h - halo for arrays
	!>@param[in] coords,dims - dimensions of cartesian topology
	!>@param[in] id - id of this PE
	!>@param[in] comm3d - communicator for cartesian topology
	subroutine allocate_and_set(dt,runtime,ntim,x, y, z, &
			xn, yn, zn, &
			u,v,w,&
			zu,zv,zw,&
			tu,tv,tw,&
			p,th,rho, &
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
			n_levels,z_read, theta_read,psurf,tsurf, &
			l_h,r_h, &
			coords,dims, id, comm3d)
				
		use numerics_type
		use mpi
		use numerics, only : find_pos, poly_int, vode_integrate
		use random, only : random_normal
		use mpi_module
		
		implicit none
		real(wp), dimension(:,:,:), allocatable, intent(inout) :: &
														u,v,w,zu,zv,zw,tu,tv,tw,&
														p,th,rho, &
														su,sv,sw,psrc
		real(wp), dimension(:), allocatable, intent(inout) :: x,y,z,xn,yn,zn,dx,dy,dz, &
															dxn,dyn,dzn, theta,thetan, &
															rhoa, rhoan, lamsq, lamsqn

		real(wp), intent(in) :: dx_nm, dy_nm, dz_nm, cvis
		real(wp), intent(in) :: dt, runtime
		integer(i4b), intent(inout) :: ipp, jpp, kpp, ipstart, jpstart, kpstart
		integer(i4b), intent(inout) :: ntim
		integer(i4b), intent(in) :: ip, jp, kp, l_h, r_h
		integer(i4b), intent(in) :: n_levels
		real(wp), dimension(n_levels), target, intent(in) :: z_read,theta_read
		real(wp), intent(in) :: psurf, tsurf
		integer(i4b), dimension(3), intent(inout) :: coords
		integer(i4b), dimension(3), intent(in) :: dims
		integer(i4b), intent(in) :: id, comm3d
		
		! locals:
		integer(i4b) :: error, AllocateStatus,i,j,k
		real(wp) :: rho_surf, htry, hmin, eps2=1.e-5_wp
		real(wp), dimension(1) :: psolve
		real(wp) :: var, dummy
		integer(i4b) :: iloc
		! for random number:
		real(wp) :: r
		real(wp), dimension(10,10) :: rs
		integer(i4b) :: l, nbottom, ntop, tag1
		integer(i4b), allocatable, dimension(:) :: seed
		real(wp) :: rad
		
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
		x=dx_nm*(/(i,i=-l_h+ipstart,ipp+r_h+ipstart-1)/) - real(ip+1,wp)/2._wp*dx_nm
		xn=x+0.5_wp*dx_nm
		! set up horizontal level array
		y=dy_nm*(/(i,i=-l_h+jpstart,jpp+r_h+jpstart-1)/) - real(jp+1,wp)/2._wp*dy_nm
		yn=y+0.5_wp*dy_nm

		! set up vertical level array
		z=dz_nm*(/(i,i=-l_h+kpstart-1,kpp+r_h+kpstart-2)/)+0.5_wp*dz_nm
		zn=z-0.5_wp*dz_nm
		
		! set up mixing length array
		lamsq=1._wp / (1._wp/(cvis*(dx_nm+dy_nm+dz)/3._wp)**2._wp + &
				1._wp/(0.4_wp*(z + 1.e-4_wp))**2._wp)
		lamsqn=1._wp / (1._wp/(cvis*(dx_nm+dy_nm+dzn)/3._wp)**2._wp + &
				1._wp/(0.4_wp*(zn + 1.e-4_wp))**2._wp)		
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
			if( z(i) < 0._wp ) then
				hmin=-1.e-2_wp
				htry=-dz(i)
				call vode_integrate(psolve,0._wp,z(i),eps2,htry,hmin,hydrostatic1a)
			else
				hmin=1.e-2_wp
				htry=dz(i)
				call vode_integrate(psolve,0._wp,z(i),eps2,htry,hmin,hydrostatic1a)
			endif
			p(i,:,:)=psolve(1)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! locate and interpolate to find theta:										 !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			iloc=find_pos(z_read(1:n_levels),z(i))
			iloc=min(n_levels-1,iloc)
			iloc=max(1,iloc)
			! linear interp theta
			call poly_int(z_read(iloc:iloc+1), theta_read(iloc:iloc+1), &
						min(z(i),z_read(n_levels)), var,dummy)
!			th(:,:,i)=var
			theta(i)=var
			rhoa(i)=psolve(1)/(ra*theta(i)*(psolve(1)/psurf)**(ra/cp))
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! solve the hydrostatic equation for k+1/2 levels							 !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			psolve=psurf
			if( zn(i) < 0._wp ) then
				hmin=-1.e-2_wp
				htry=-dzn(i)
				call vode_integrate(psolve,0._wp,zn(i),eps2,htry,hmin,hydrostatic1a)
			else
				hmin=1.e-2_wp
				htry=dzn(i)
				call vode_integrate(psolve,0._wp,zn(i),eps2,htry,hmin,hydrostatic1a)
			endif
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! locate and interpolate to find theta:										 !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			iloc=find_pos(z_read(1:n_levels),zn(i))
			iloc=min(n_levels-1,iloc)
			iloc=max(1,iloc)
			! linear interp theta
			call poly_int(z_read(iloc:iloc+1), theta_read(iloc:iloc+1), &
						min(zn(i),z_read(n_levels)), var,dummy)
			thetan(i)=var
			rhoan(i)=psolve(1)/(ra*thetan(i)*(psolve(1)/psurf)**(ra/cp))
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		enddo
		p(:,:,:)=0._wp
! 		rhoa=1._wp
! 		rhoan=1._wp
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! set wind field        														 !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
		do i=1,ipp
			do j=1,jpp
				do k=1,kpp
					u(k,j,i)=-5._wp*(y(j))/sqrt(xn(i)*xn(i)+y(j)*y(j))
					v(k,j,i)=5._wp*(x(i))/sqrt(x(i)*x(i)+yn(j)*yn(j))
				enddo
			enddo
		enddo
		u(:,:,:)=0._wp
		v(:,:,:)=0._wp
		w(:,:,:)=0._wp
		zu(:,:,:)=0._wp
		zv(:,:,:)=0._wp
		zw(:,:,:)=0._wp
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
					
						if ( (z(k-kpstart)>3000._wp) .and. (z(k-kpstart)<6000._wp) ) &
							th(k-kpstart,j-jpstart,i-ipstart) = + r / 30._wp
						

					endif
					
					
					

				enddo
			enddo
		enddo

 		th=0._wp
		
		
		do i=1-r_h,ipp+r_h
			do j=1-r_h,jpp+r_h
				do k=1-r_h,kpp+r_h
					!th(k,j,i)=theta(k)
				
					rad = (z(k)-3000._wp)**2._wp
						
					if (ip > 1) rad=rad+x(i)**2._wp
					if (jp > 1) rad=rad+y(j)**2._wp
					
					rad=sqrt(rad)
					if(rad<=1000._wp) then
						th(k,j,i)=th(k,j,i)+0.1_wp
					else
						!th(k,j,i)=0._wp
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
		    0._wp,0._wp, dims,coords)
		call exchange_full(comm3d, id, kpp, jpp, ipp, r_h,r_h,l_h,r_h,r_h,r_h, v,&
		    0._wp,0._wp, dims,coords)
		call exchange_full(comm3d, id, kpp, jpp, ipp, l_h,r_h,r_h,r_h,r_h,r_h, w,&
		    0._wp,0._wp, dims,coords)
		call exchange_full(comm3d, id, kpp, jpp, ipp, r_h,r_h,r_h,r_h,r_h,r_h, th,&
		    0._wp,0._wp, dims,coords)
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






	end module initialisation
	