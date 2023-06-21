	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>initialisation for the sub-grid scale model
    module initialisation
    use numerics_type
    implicit none
    
    real(wp), parameter :: ra=287._wp, grav=9.81_wp,cp=1005._wp
	real(wp), dimension(:), pointer :: zr1, thr1
	integer(i4b) :: nl1
	
    private
    public :: allocate_and_set, allocate_and_set_1d, allocate_and_set_2d
    contains
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Allocate and set arrays                                                            !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>allocate arrays on each PE, and initialise them
	!>@param[inout] dt,runtime,ntim - time variables
	!>@param[inout] x,y,z,xn,yn,zn,u,v,w,th,q - grid positions and prognostics
	!>@param[inout] zu,zv,zw,tu,tv,tw - prognostics
	!>@param[inout] vism,vist,viss,tau11,tau22,tau33,
	!>@param[inout] tau12,tau13,tau23,strain,su,sv,sw,sth,sq, 
	!>@param[inout] rhoa,rhoan, theta, thetan - reference potential temperatures
	!>@param[inout] lamsq,lamsqn - mixing length
	!>@param[in] z0, z0th: roughness length for momentum and scalars (z0th not used yet)
	!>@param[inout] lbc, ubc
	!>@param[in] cvis - smagorinsky parameter
	!>@param[inout] dx,dy,dz - grid spacing on grid
	!>@param[inout] dxn,dyn,dzn - grid spacing on grid - staggered
	!>@param[inout] ipp,jpp,kpp,ipstart,jpstart,kpstart - number of grid / starting position
	!>@param[in] dx_nm,dy_nm,dz_nm - grid spacing from namelist
	!>@param[in] ip,jp,kp - grid points from namelist
	!>@param[in] l_h,r_h - halo for arrays
	!>@param[in] coords,dims - dimensions of cartesian topology
	!>@param[in] id - id of this PE
	!>@param[in] comm3d - communicator for cartesian topology
	subroutine allocate_and_set(dt,runtime,ntim,x, y, z, &
			xn, yn, zn, &
			u,v,w,th,&
			zu,zv,zw,&
			tu,tv,tw,&
			q, &
			vism,vist,viss,tau11,tau22,tau33,tau12,tau13,tau23,strain,su,sv,sw,sth,sq, & 
			rhoa,rhoan, theta, thetan, &
			lamsq,lamsqn, &
			z0,z0th, &
			lbc,ubc, &
			cvis, &
			dx, dy, dz, &
			dxn, dyn, dzn, &
			ipp, jpp, kpp,nqg,&
			ipstart, jpstart, kpstart, &
			dx_nm, dy_nm, dz_nm, &
			ip, jp, kp, nq, &
			l_h,r_h, &
			coords,dims, id, comm3d)
				
		use numerics_type
		use mpi
		use netcdf
		use random, only : random_normal
		use mpi_module
		
		implicit none
		real(wp), dimension(:,:,:), allocatable, intent(inout) :: &
														u,v,w, th,&
														zu,zv,zw, &
														tu,tv,tw, &
														vism,vist,&
														tau11,tau22,tau33,&
														tau12,tau13,tau23,strain,&
														su,sv,sw,sth
		real(wp), dimension(:,:,:,:), allocatable, intent(inout) :: &
														q,sq,viss
		real(wp), dimension(:), allocatable, intent(inout) :: x,y,z,xn,yn,zn,dx,dy,dz, &
															dxn,dyn,dzn, &
															rhoa, rhoan, &
															theta, thetan, &
															lamsq, lamsqn, lbc, ubc

		real(wp), intent(in) :: dx_nm, dy_nm, dz_nm, cvis, z0, z0th
		real(wp), intent(in) :: dt, runtime
		integer(i4b), intent(inout) :: ipp, jpp, kpp, ipstart, jpstart, kpstart,nqg
		integer(i4b), intent(inout) :: ntim
		integer(i4b), intent(in) :: ip, jp, kp, l_h, r_h,nq
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
		
		

		
		! scalar formulae:
		nqg=nq
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
		allocate( vism(1-l_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( vist(1-l_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( tau11(1-l_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( tau22(1-l_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( tau33(1-l_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( tau12(1-l_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( tau13(1-l_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( tau23(1-l_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( strain(1-l_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( th(1-l_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( sth(1-l_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"

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

		allocate( su(1-r_h:kpp+r_h,1-r_h:jpp+r_h,1-l_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( sv(1-r_h:kpp+r_h,1-l_h:jpp+r_h,1-r_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( sw(1-l_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		

		allocate( q(1-l_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h,1:nqg), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( sq(1-l_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h,1:nqg), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( viss(1-l_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h,1:nqg), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"

		allocate( lbc(1:nqg), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( ubc(1:nqg), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"		

		allocate( x(1-l_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( y(1-l_h:jpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( z(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( rhoa(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( lamsq(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( theta(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		
		allocate( xn(1-l_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( yn(1-l_h:jpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( zn(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( rhoan(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( lamsqn(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( thetan(1-l_h:kpp+r_h), STAT = AllocateStatus)
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

		! lower and upper bcs
		lbc=0._wp
		ubc=0._wp

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! set up grid spacing arrays                                                     !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! set some arrays to zero
		su=0._wp
		sv=0._wp
		sw=0._wp
		sth=0._wp
		sq=0._wp
		vism=0._wp
		vist=0._wp
		viss=0._wp
		tau11=0._wp
		tau22=0._wp
		tau33=0._wp
		tau12=0._wp
		tau13=0._wp
		tau23=0._wp
		strain=0._wp
		q=0._wp

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
				
		! temporary density
		rhoa=1._wp
		rhoan=1._wp
		theta=300._wp-z/1000._wp
		thetan=300._wp-zn/1000._wp
		
		! set up mixing length array
		lamsq=1._wp / (1._wp/(cvis*(dx_nm+dy_nm+dz)/3._wp)**2._wp + &
				1._wp/(0.4_wp*(z + z0))**2._wp)
		lamsqn=1._wp / (1._wp/(cvis*(dx_nm+dy_nm+dzn)/3._wp)**2._wp + &
				1._wp/(0.4_wp*(zn + z0))**2._wp)		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! set wind field        														 !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
		v=0._wp	
		w=0._wp
		u=0._wp
		do i=1-l_h,ipp+r_h
			do j=1-l_h,jpp+r_h
				do k=1-l_h,kpp+r_h
! 					u(k,j,i)=-5.e-4_wp*(yn(j))  !/sqrt(xn(i)*xn(i)+y(j)*y(j))
! 					v(k,j,i)=7.5e-4_wp*(xn(i))  !/sqrt(x(i)*x(i)+yn(j)*yn(j))
					v(k,j,i)=5._wp !1.e-4_wp*(zn(k))  
					!if(zn(k)>5000._wp) v(k,j,i)=1.e-4_wp*(5000._wp) 
					u(k,j,i)=0._wp  
				enddo
			enddo
		enddo
		if(coords(3)==0) v(1,:,:)=0._wp
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! set scalar field        														 !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
        q(:,:,:,:)=0._wp
        th=0._wp
        rad=0._wp
		do i=1-r_h,ipp+r_h
			do j=1-r_h,jpp+r_h
				do k=1-r_h,kpp+r_h
				    
				    !if(kp> 1) rad = (z(k)-2000._wp)**2._wp
					rad=0._wp
					if (ip > 1) rad=rad+x(i)**2._wp
					if (jp > 1) rad=rad+2._wp*y(j)**2._wp
					
					rad=sqrt(rad)
					if(rad<=1000._wp) then
						th(k,j,i)=0.1_wp
						q(k,j,i,1)=q(k,j,i,1)+0.1_wp
                        q(k,j,i,2)=q(k,j,i,1)*2._wp
                        q(k,j,i,3)=q(k,j,i,1)*3._wp
                        q(k,j,i,4)=q(k,j,i,1)*4._wp
                        q(k,j,i,5)=q(k,j,i,1)*5._wp
                        q(k,j,i,6)=q(k,j,i,1)*6._wp
                        q(k,j,i,7)=q(k,j,i,1)*7._wp
                        q(k,j,i,8)=q(k,j,i,1)*8._wp
                        q(k,j,i,9)=q(k,j,i,1)*9._wp
					endif
				enddo
			enddo
		enddo
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		





		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! set halos																		 !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
		call exchange_full(comm3d, id, kpp, jpp, ipp, r_h,r_h,r_h,r_h,l_h,r_h, u,&
		        0._wp,0._wp,dims,coords)
		call exchange_full(comm3d, id, kpp, jpp, ipp, r_h,r_h,l_h,r_h,r_h,r_h, v,&
		        0._wp,0._wp,dims,coords)
		call exchange_full(comm3d, id, kpp, jpp, ipp, l_h,r_h,r_h,r_h,r_h,r_h, w,&
		        0._wp,0._wp,dims,coords)
		do k=1,nqg
    		call exchange_full(comm3d, id, kpp, jpp, ipp, &
    		    l_h,r_h,r_h,r_h,r_h,r_h, q(:,:,:,k),&
		        0._wp,0._wp,dims,coords)	
		enddo
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		






	end subroutine allocate_and_set
	

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Allocate and set arrays                                                            !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>allocate arrays on each PE, and initialise them
	!>@param[inout] dt,runtime,ntim - time variables
	!>@param[inout] z,zn,u,w,zu,zw,tu,tw,q - grid positions and prognostics
	!>@param[inout] vism,vist,viss,tau11,tau33,tau13, strain - subgrid
	!>@param[inout] rhoa,rhoan - reference potential temperatures
	!>@param[inout] lamsq,lamsqn - mixing length
	!>@param[in] z0, z0th - roughness lengths
	!>@param[inout] lbc, ubc - boundary conditions
	!>@param[in] cvis - smagorinsky parameter
	!>@param[inout] dz - grid spacing on grid
	!>@param[inout] dzn - grid spacing on grid - staggered
	!>@param[inout] kpp,kpstart - number of grid / starting position
	!>@param[in] dz_nm - grid spacing from namelist
	!>@param[in] kp - grid points from namelist
	!>@param[in] l_h,r_h - halo for arrays
	subroutine allocate_and_set_1d(dt,runtime,ntim,z, &
			zn, &
			u,w,th,&
			zu,zw, &
			tu,tw,&
			q, &
			vism,vist,viss,&
			tau11,tau33,&
			tau13,strain,&
			su,sw,sth,sq, & 
			rhoa,rhoan,theta, thetan, &
			lamsq,lamsqn, &
			z0,z0th, &
			lbc,ubc, &
			cvis, &
			dz, &
			dzn, &
			kpp,nqg,&
			kpstart, &
			dz_nm, &
			kp, nq, &
			l_h,r_h)


				
		use numerics_type
		use netcdf
		use random, only : random_normal
		

		implicit none
		real(wp), dimension(:), allocatable, intent(inout) :: &
														u,w, zu,zw,tu,tw,vism, vist, &
														tau11,tau33,tau13, strain, &
														su,sw,sth,th
		real(wp), dimension(:,:), allocatable, intent(inout) :: &
														q,viss,sq
		real(wp), dimension(:), allocatable, intent(inout) :: z,zn,dz, &
															dzn, &
															rhoa, rhoan, &
															theta,thetan, &
															lamsq, lamsqn,lbc,ubc

		real(wp), intent(in) :: dz_nm, cvis
		real(wp), intent(in) :: dt, runtime
		real(wp), intent(inout) :: z0th,z0
		integer(i4b), intent(inout) :: kpp, kpstart,nqg
		integer(i4b), intent(inout) :: ntim
		integer(i4b), intent(in) :: kp, l_h, r_h,nq
		
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
		
		
						
		

		
		! scalar formulae:
		nqg=nq
		ntim=ceiling(runtime/dt)

		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! find the number of grid points on each PE                                      !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		kpp=kp
		kpstart=1
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		
		
		
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! allocate arrays                                                                !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		allocate( vism(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( vist(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( tau11(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( tau33(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( tau13(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( strain(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( th(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( sth(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"

		allocate( u(1-r_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( w(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( zu(1-r_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( zw(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( tu(1-r_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( tw(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"

		allocate( su(1-r_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( sw(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		

		allocate( q(1-l_h:kpp+r_h,1:nqg), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( sq(1-l_h:kpp+r_h,1:nqg), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( viss(1-l_h:kpp+r_h,1:nqg), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"

		allocate( lbc(1:nqg), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( ubc(1:nqg), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"		

		allocate( z(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( rhoa(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( lamsq(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( theta(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		
		allocate( zn(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( rhoan(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( lamsqn(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( thetan(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		
		allocate( dz(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"		

		allocate( dzn(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




		! lower and upper bcs
		lbc=0._wp
		ubc=0._wp

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! set up grid wpacing arrays                                                     !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! set some arrays to zero
		su=0._wp
		sw=0._wp
		sth=0._wp
		sq=0._wp
		vism=0._wp
		vist=0._wp
		viss=0._wp
		tau11=0._wp
		tau33=0._wp
		tau13=0._wp
		strain=0._wp

		! grid spacing:
		dz(:)=dz_nm
		! grid spacing, staggered:
		dzn(:)=dz_nm
		! set up vertical level array
		z=dz_nm*(/(i,i=-l_h+kpstart-1,kpp+r_h+kpstart-2)/)+1.0_wp*dz_nm
		zn=z-0.5_wp*dz_nm
		
		! temporary density
		rhoa=1._wp
		rhoan=1._wp
		theta=300._wp-z/1000._wp
		thetan=300._wp-zn/1000._wp
		
		! set up mixing length array
		lamsq=1._wp / (1._wp/(cvis*(dz))**2._wp + &
				1._wp/(0.4_wp*(z + z0))**2._wp)
		lamsqn=1._wp / (1._wp/(cvis*(dzn))**2._wp + &
				1._wp/(0.4_wp*(zn + z0))**2._wp)		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! set wind field        														 !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
		w=0._wp
		u=0._wp
        do k=1-l_h,kpp+r_h
! 					u(k,j,i)=-5.e-4_wp*(yn(j))  !/sqrt(xn(i)*xn(i)+y(j)*y(j))
! 					v(k,j,i)=7.5e-4_wp*(xn(i))  !/sqrt(x(i)*x(i)+yn(j)*yn(j))
            u(k)=5._wp !1.e-4_wp*(zn(k))  
            !if(zn(k)>5000._wp) v(k,j,i)=1.e-4_wp*(5000._wp) 
        enddo
		u(1)=0._wp



		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! set scalar field        														 !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
        q(1-l_h:kpp+r_h,1:nqg)=0._wp
        rad=0._wp
        do k=1-r_h,kpp+r_h
            
!             if(kp> 1) rad = (z(k)-3000._wp)**2._wp
                
            
            rad=sqrt(rad)
            if(rad<=1000._wp) then
                q(k,1)=q(k,1)+0.1_wp
                q(k,2)=q(k,1)*2._wp
                q(k,3)=q(k,1)*3._wp
                q(k,4)=q(k,1)*4._wp
                q(k,5)=q(k,1)*5._wp
                q(k,6)=q(k,1)*6._wp
                q(k,7)=q(k,1)*7._wp
                q(k,8)=q(k,1)*8._wp
                q(k,9)=q(k,1)*9._wp
            endif
        enddo
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		



		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! set halos																		 !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
		u(0)=0._wp
		u(kpp+1)=0._wp
		w(0)=0._wp
		w(kpp+1)=0._wp
		do k=1,nqg
            q(0,k)=0._wp
            q(kpp+1,k)=0._wp
		enddo
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		







	end subroutine allocate_and_set_1d
	

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Allocate and set arrays                                                            !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>allocate arrays on each PE, and initialise them
	!>@param[inout] dt,runtime,ntim - time variables
	!>@param[inout] x,z,xn,zn,u,w,zu,zw,tu,tw,q - grid positions and prognostics
	!>@param[inout] vism,vist,viss,tau11,tau33,tau13, strain - subgrid
	!>@param[inout] rhoa,rhoan - reference potential temperatures
	!>@param[inout] lamsq,lamsqn - mixing length
	!>@param[in] z0, z0th - roughness lengths
	!>@param[inout] lbc, ubc - boundary conditions
	!>@param[in] cvis - smagorinsky parameter
	!>@param[inout] dx,dz - grid spacing on grid
	!>@param[inout] dxn,dzn - grid spacing on grid - staggered
	!>@param[inout] ipp,kpp,ipstart,jpstart,kpstart - number of grid / starting position
	!>@param[in] dx_nm,dz_nm - grid spacing from namelist
	!>@param[in] ip,kp - grid points from namelist
	!>@param[in] l_h,r_h - halo for arrays
	subroutine allocate_and_set_2d(dt,runtime,ntim,x, z, &
			xn, zn, &
			u,w,th, &
			zu,zw, &
			tu, tw, &
			q, &
			vism,vist,viss,&
			tau11,tau33,&
			tau13,strain,&
			su,sw,sth,sq, & 
			rhoa,rhoan,theta, thetan, &
			lamsq,lamsqn, &
			z0,z0th, &
			lbc,ubc, &
			cvis, &
			dx, dz, &
			dxn, dzn, &
			ipp, kpp,nqg,&
			ipstart, kpstart, &
			dx_nm, dz_nm, &
			ip, kp, nq, &
			l_h,r_h)
				
		use numerics_type
		use mpi
		use netcdf
		use random, only : random_normal
		use mpi_module
		
		implicit none
		real(wp), dimension(:,:), allocatable, intent(inout) :: &
														u,w, zu,zw,tu,tw,vism, vist, &
														tau11,tau33,tau13, strain, &
														su,sw,sth,th
		real(wp), dimension(:,:,:), allocatable, intent(inout) :: &
														q,viss,sq
		real(wp), dimension(:), allocatable, intent(inout) :: x,z,xn,zn,dx,dz, &
															dxn,dzn, &
															rhoa, rhoan, &
															theta,thetan, &
															lamsq, lamsqn,lbc,ubc

		real(wp), intent(in) :: dx_nm, dz_nm, cvis
		real(wp), intent(in) :: dt, runtime
		real(wp), intent(inout) :: z0th,z0
		integer(i4b), intent(inout) :: ipp, kpp, ipstart, kpstart,nqg
		integer(i4b), intent(inout) :: ntim
		integer(i4b), intent(in) :: ip, kp, l_h, r_h,nq
		
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
		
		
		

		
		! scalar formulae:
		nqg=nq
		ntim=ceiling(runtime/dt)

		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! find the number of grid points on each PE                                      !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		ipp=ip
		ipstart=1
		kpp=kp
		kpstart=1
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		
		
		
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! allocate arrays                                                                !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		allocate( vism(1-l_h:kpp+r_h,1-r_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( vist(1-l_h:kpp+r_h,1-r_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( tau11(1-l_h:kpp+r_h,1-r_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( tau33(1-l_h:kpp+r_h,1-r_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( tau13(1-l_h:kpp+r_h,1-r_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( strain(1-l_h:kpp+r_h,1-r_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( th(1-l_h:kpp+r_h,1-r_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( sth(1-l_h:kpp+r_h,1-r_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"

		allocate( u(1-r_h:kpp+r_h,1-l_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( w(1-l_h:kpp+r_h,1-r_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( zu(1-r_h:kpp+r_h,1-l_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( zw(1-l_h:kpp+r_h,1-r_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( tu(1-r_h:kpp+r_h,1-l_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( tw(1-l_h:kpp+r_h,1-r_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"

		allocate( su(1-r_h:kpp+r_h,1-l_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( sw(1-l_h:kpp+r_h,1-r_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		

		allocate( q(1-l_h:kpp+r_h,1-r_h:ipp+r_h,1:nqg), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( sq(1-l_h:kpp+r_h,1-r_h:ipp+r_h,1:nqg), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( viss(1-l_h:kpp+r_h,1-r_h:ipp+r_h,1:nqg), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"

		allocate( lbc(1:nqg), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( ubc(1:nqg), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"		

		allocate( x(1-l_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( z(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( rhoa(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( lamsq(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( theta(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		
		allocate( xn(1-l_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( zn(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( rhoan(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( lamsqn(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( thetan(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		
		allocate( dx(1-l_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( dz(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"		

		allocate( dxn(1-l_h:ipp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( dzn(1-l_h:kpp+r_h), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



		! lower and upper bcs
		lbc=0._wp
		ubc=0._wp

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! set up grid spacing arrays                                                     !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! set some arrays to zero
		su=0._wp
		sw=0._wp
		sth=0._wp
		sq=0._wp
		vism=0._wp
		vist=0._wp
		viss=0._wp
		tau11=0._wp
		tau33=0._wp
		tau13=0._wp
		strain=0._wp

		! grid spacing:
		dx(:)=dx_nm
		dz(:)=dz_nm
		! grid spacing, staggered:
		dxn(:)=dx_nm
		dzn(:)=dz_nm
		! set up horizontal level array
		x=dx_nm*(/(i,i=-l_h+ipstart,ipp+r_h+ipstart-1)/) - real(ip-1,wp)/2._wp*dx_nm 
		xn=x-0.5_wp*dx_nm

		! set up vertical level array
		z=dz_nm*(/(i,i=-l_h+kpstart-1,kpp+r_h+kpstart-2)/)+1.0_wp*dz_nm
		zn=z-0.5_wp*dz_nm

		! temporary density
		rhoa=1._wp
		rhoan=1._wp
		theta=300._wp-z/1000._wp
		thetan=300._wp-zn/1000._wp
		
		! set up mixing length array
		lamsq=1._wp / (1._wp/(cvis*(dx_nm+dz)/2._wp)**2._wp + &
				1._wp/(0.4_wp*(z + z0))**2._wp)
		lamsqn=1._wp / (1._wp/(cvis*(dx_nm+dzn)/2._wp)**2._wp + &
				1._wp/(0.4_wp*(zn + z0))**2._wp)		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! set wind field        														 !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
		w=0._wp
		u=0._wp
		do i=1-l_h,ipp+r_h
            do k=1-l_h,kpp+r_h
! 					u(k,j,i)=-5.e-4_wp*(yn(j))  !/sqrt(xn(i)*xn(i)+y(j)*y(j))
! 					v(k,j,i)=7.5e-4_wp*(xn(i))  !/sqrt(x(i)*x(i)+yn(j)*yn(j))
                u(k,i)=5._wp !1.e-4_wp*(zn(k))  
                !if(zn(k)>5000._wp) v(k,j,i)=1.e-4_wp*(5000._wp) 
            enddo
		enddo
		u(1,:)=0._wp
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! set scalar field        														 !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
        q(:,:,:)=0._wp
        rad=0._wp
		do i=1-r_h,ipp+r_h
            do k=1-r_h,kpp+r_h
                
                !if(kp> 1) rad = (z(k)-2000._wp)**2._wp
                rad=0._wp
                if (ip > 1) rad=rad+x(i)**2._wp
                
                rad=sqrt(rad)
                if(rad<=1000._wp) then
                    th(k,i)=0.1_wp
                    q(k,i,1)=q(k,i,1)+0.1_wp
                    q(k,i,2)=q(k,i,1)*2._wp
                    q(k,i,3)=q(k,i,1)*3._wp
                    q(k,i,4)=q(k,i,1)*4._wp
                    q(k,i,5)=q(k,i,1)*5._wp
                    q(k,i,6)=q(k,i,1)*6._wp
                    q(k,i,7)=q(k,i,1)*7._wp
                    q(k,i,8)=q(k,i,1)*8._wp
                    q(k,i,9)=q(k,i,1)*9._wp
                endif
            enddo
		enddo
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		





		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! set halos																		 !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
		u(0,:)=0._wp
		u(kpp+1,:)=0._wp
		u(:,0)=u(:,ipp)
		u(:,ipp+1)=u(:,1)
		w(0,:)=0._wp
		w(kpp+1,:)=0._wp
		w(:,0)=w(:,ipp)
		w(:,ipp+1)=w(:,1)
		do k=1,nqg
            q(0,:,k)=0._wp
            q(kpp+1,:,k)=0._wp
            q(:,0,k)=q(:,ipp,k)
            q(:,ipp+1,k)=q(:,1,k)
		    
		enddo
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		






	end subroutine allocate_and_set_2d



	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! HELPER ROUTINE                                                       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine check(status)
	use netcdf
	use numerics_type
	integer(i4b), intent ( in) :: status

	if(status /= nf90_noerr) then
		print *, trim(nf90_strerror(status))
		stop "Stopped"
	end if
	end subroutine check
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	end module initialisation
	