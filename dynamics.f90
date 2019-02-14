	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>dynamics solver routines for the dynamical cloud model
    module d_solver
    use nrtype
    implicit none
	real(sp), dimension(:,:,:), allocatable :: &
					a_e,a_w,a_n,a_s,a_u,a_d,a_p,a_c
	logical :: bicgstab_first = .true. ! true if this is the first time bicgstab called
	real(sp), parameter :: grav=9.81_sp ! gravity
	
    private
    public :: bicgstab, sources, advance_momentum
    
	contains
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>set-up matrix A so that Ax is equivalent to \f$\nabla ^2 x\f$ using finite-diff
	!>@param[in] comm3d,id, dims, coords
	!>@param[in] dt,dx,dy,dz,dxn,dyn,dzn
	!>@param[in] ip,jp,kp,halo
	!>@param[inout] a_e,a_w,a_n,a_s,a_u,a_d,a_p: terms to make up 7-point stencil
	!>@param[in] su,sv,sw
	subroutine set_mat_a(comm3d,id,dims,coords, dt,dx,dy,dz,dxn,dyn,dzn,ip,jp,kp,halo, &
						a_e,a_w,a_n,a_s,a_u,a_d,a_p,a_c,su,sv,sw)
		use nrtype
		use mpi_module
		implicit none
		integer(i4b), intent(in) :: id, comm3d
		integer(i4b), dimension(3), intent(in) :: dims,coords
		integer(i4b), intent(in) :: ip,jp,kp, halo
		real(sp), dimension(1-halo:ip+halo) :: dx,dxn
		real(sp), dimension(1-halo:jp+halo) :: dy,dyn
		real(sp), dimension(1-halo:kp+halo) :: dz,dzn
		real(sp), dimension(-halo+1:kp+halo,-halo+1:jp+halo,-halo+1:ip+halo), &
			intent(inout) :: a_e,a_w,a_n,a_s,a_u,a_d,a_p, a_c
		real(sp), dimension(-halo+1:kp+halo,-halo+1:jp+halo,-halo+1:ip+halo), &
			intent(in) :: su,sv,sw
		real(sp), intent(in) :: dt
				
		! local
		integer(i4b) :: i,j,k

		a_p=1._sp
		a_n=0._sp
		a_s=0._sp
		a_e=0._sp
		a_w=0._sp
		a_u=0._sp
		a_d=0._sp
		a_c=0._sp
!$omp simd	
		do i=1,ip
			do j=1,jp			    
				do k=1,kp
					a_e(k,j,i)=1._sp/( dx(i-1)*dxn(i) )
					a_w(k,j,i)=1._sp/( dx(i-1)*dxn(i-1) )
					a_n(k,j,i)=1._sp/( dy(j-1)*dyn(j) )
					a_s(k,j,i)=1._sp/( dy(j-1)*dyn(j-1) )
					a_u(k,j,i)=1._sp/( dz(k-1)*dzn(k) )
					a_d(k,j,i)=1._sp/( dz(k-1)*dzn(k-1) )
					a_p(k,j,i)=-1._sp/( dx(i-1)*dxn(i) ) -1._sp/( dx(i-1)*dxn(i-1) ) &
                            -1._sp/( dy(j-1)*dyn(j) ) -1._sp/( dy(j-1)*dyn(j-1) )  &
                            -1._sp/( dz(k-1)*dzn(k) ) -1._sp/( dz(k-1)*dzn(k-1) ) 
				enddo
				if(coords(3)==0) then
				    k=1
                    a_e(k,j,i)=0._sp
                    a_w(k,j,i)=0._sp
                    a_n(k,j,i)=0._sp
                    a_s(k,j,i)=0._sp
                    a_u(k,j,i)=1._sp/( dz(k-1)*dzn(k) )
                    a_d(k,j,i)=1._sp/( dz(k-1)*dzn(k-1) )
                    a_p(k,j,i)=-1._sp/( dx(i-1)*dxn(i) ) -1._sp/( dx(i-1)*dxn(i-1) ) &
                            -1._sp/( dy(j-1)*dyn(j) ) -1._sp/( dy(j-1)*dyn(j-1) )  &
                            -1._sp/( dz(k-1)*dzn(k) ) -1._sp/( dz(k-1)*dzn(k-1) ) + &
                            1._sp/( dx(i-1)*dxn(i) ) + & ! a_e
                            1._sp/( dx(i-1)*dxn(i-1) ) + & ! a_w
                            1._sp/( dy(j-1)*dyn(j) ) + & ! a_n
                            1._sp/( dy(j-1)*dyn(j-1) )  ! a_s
                    a_c(k,j,i)=1._sp/( dx(i-1)*dxn(i) )*dxn(i)*su(k,j,i)/(2._sp*dt) - &
                            1._sp/( dx(i-1)*dxn(i-1) )*dxn(i-1)*su(k,j,i-1)/(2._sp*dt) + &
                            1._sp/( dy(j-1)*dyn(j) )*dyn(j)*sv(k,j,i)/(2._sp*dt) - &
                            1._sp/( dy(j-1)*dyn(j-1) )*dyn(j-1)*sv(k,j-1,i)/(2._sp*dt) 
                endif
			enddo
		enddo
!$omp end simd

!     if(coords(3)==0) then
! 		do i=1,ip
! 			do j=1,jp
!                 a_u(1,j,i)=a_u(1,j,i)+a_d(1,j,i)
!                 a_d(1,j,i)=0._sp    
!             enddo
!         enddo
!     endif
		
!     if(coords(3)==dims(3)-1) then
! 		do i=1,ip
! 			do j=1,jp
!                 a_u(kp,j,i)=a_u(kp,j,i)+a_d(kp,j,i)
!                 a_d(kp,j,i)=0._sp    
!             enddo
!         enddo
!     endif
		
! 	call exchange_halos(comm3d, id, kp, jp, ip, &
! 		halo,halo,halo,halo,halo,halo, a_p,dims,coords)
! 	call exchange_halos(comm3d, id, kp, jp, ip, &
! 		halo,halo,halo,halo,halo,halo, a_n,dims,coords)
! 	call exchange_halos(comm3d, id, kp, jp, ip, &
! 		halo,halo,halo,halo,halo,halo, a_s,dims,coords)
! 	call exchange_halos(comm3d, id, kp, jp, ip, &
! 		halo,halo,halo,halo,halo,halo, a_e,dims,coords)
! 	call exchange_halos(comm3d, id, kp, jp, ip, &
! 		halo,halo,halo,halo,halo,halo, a_w,dims,coords)
! 	call exchange_halos(comm3d, id, kp, jp, ip, &
! 		halo,halo,halo,halo,halo,halo, a_u,dims,coords)
! 	call exchange_halos(comm3d, id, kp, jp, ip, &
! 		halo,halo,halo,halo,halo,halo, a_d,dims,coords)

	end subroutine set_mat_a	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>multiply matrix A by x
	!>@param[in] comm3d,id, dims, coords
	!>@param[in] ip,jp,kp,halo
	!>@param[in] a_e,a_w,a_n,a_s,a_u,a_d,a_p,a_c: terms to make up 7-point stencil
	!>@param[inout] ax,x: matrix ax and x
	subroutine mat_ax(comm3d,id,dims,coords, ip,jp,kp,halo, &
						a_e,a_w,a_n,a_s,a_u,a_d,a_p,a_c,ax,x)
		use nrtype
		use mpi_module, only : exchange_along_dim, exchange_full
		implicit none
		integer(i4b), intent(in) :: id, comm3d
		integer(i4b), dimension(3), intent(in) :: dims,coords
		integer(i4b), intent(in) :: ip,jp,kp, halo
		real(sp), dimension(-halo+1:kp+halo,-halo+1:jp+halo,-halo+1:ip+halo), &
			intent(in) :: a_e,a_w,a_n,a_s,a_u,a_d,a_p,a_c
		real(sp), dimension(-halo+1:kp+halo,-halo+1:jp+halo,-halo+1:ip+halo), &
			intent(inout) :: ax,x
		
		! local
		integer(i4b) :: i,j,k
		
		call exchange_along_dim(comm3d, id, kp, jp, ip, &
			halo,halo,halo,halo,halo,halo, x,0._sp,0._sp,dims,coords)

		ax=0._sp
!$omp simd		
		do i=1,ip
			do j=1,jp
				do k=1,kp
                    ax(k,j,i) = a_w(k,j,i)*x(k,j,i-1) + a_e(k,j,i)*x(k,j,i+1) + &
                                a_s(k,j,i)*x(k,j-1,i) + a_n(k,j,i)*x(k,j+1,i) + &
                                a_d(k,j,i)*x(k-1,j,i) + a_u(k,j,i)*x(k+1,j,i) + &
                                a_p(k,j,i)*x(k,j,i) +a_c(k,j,i)
				enddo
			enddo
		enddo
!$omp end simd

! 		call exchange_full(comm3d, id, kp, jp, ip, &
! 			halo,halo,halo,halo,halo,halo, ax,dims,coords)

! 		ax(0,:,:)=ax(ip,:,:)
! 		ax(ip+1,:,:)=ax(1,:,:)
! 		ax(:,0,:)=ax(:,jp,:)
! 		ax(:,jp+1,:)=ax(:,1,:)

	end subroutine mat_ax
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>precondition the matrix using Jacobi iteration
	!>@param[in] comm3d,id, dims,coords
	!>@param[in] ip,jp,kp,halo
	!>@param[in] a_e,a_w,a_n,a_s,a_u,a_d,a_p,a_c: terms to make up 7-point stencil and x
	!>@param[inout] s,r,prec,prit
	subroutine precondition(comm3d,id,dims,coords, &
				ip,jp,kp,halo, a_e,a_w,a_n,a_s,a_u,a_d,a_p,a_c, &
						s,r,prit)
		use nrtype
		use mpi_module, only : exchange_along_dim
		implicit none
		integer(i4b), intent(in) :: id, comm3d
		integer(i4b), dimension(3), intent(in) :: dims,coords
		integer(i4b), intent(in) :: ip,jp,kp, halo, prit
		real(sp), dimension(-halo+1:kp+halo,-halo+1:jp+halo,-halo+1:ip+halo), &
			intent(in) :: r,a_e,a_w,a_n,a_s,a_u,a_d,a_p,a_c
		real(sp), dimension(-halo+1:kp+halo,-halo+1:jp+halo,-halo+1:ip+halo), &
			intent(inout) :: s
		
		! local
		integer(i4b) :: it
		real(sp), dimension(-halo+1:kp+halo,-halo+1:jp+halo,-halo+1:ip+halo) :: t
		

		s = r/(a_p)
		do it=1,prit
			call mat_ax(comm3d, id, dims,coords,ip,jp,kp,halo, &
						a_e,a_w,a_n,a_s,a_u,a_d,a_p,a_c,t,s)
			s = s + (r-t)/(a_p)
			
			! MPI:
! 			call exchange_halos(comm3d, id, kp, jp, ip, &
! 				halo,halo,halo,halo,halo,halo, s,dims,coords)

		enddo


	end subroutine precondition
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>set-up and solve poisson's equation:
	!>\f$\nabla ^2 f = \rho \f$
	!>  using bi-conjugate gradient stabilised method (bicgstab)
	!>@param[in] comm3d,id, rank, dims, coords
	!>@param[in] dt
	!>@param[in] xg,yg,zg,dx,dy,dz,dxn,dyn,dzn
	!>@param[in] ip,jp,kp,l_h,r_h
	!>@param[in] su
	!>@param[in] sv
	!>@param[in] sw
	!>@param[inout] x
	!>@param[inout] b
	!>@param[in] test_solver	
	subroutine bicgstab(comm3d,id,rank, dims, coords, &
			dt,xg,yg,zg,dx,dy,dz,dxn,dyn,dzn,ip,jp,kp,l_h,r_h,su,sv,sw,x,b,test_solver)
		use nrtype
		implicit none
		integer(i4b), intent(in) :: id, comm3d, rank
		integer(i4b), dimension(3), intent(in) :: dims,coords
		real(sp), intent(in) :: dt
		integer(i4b), intent(in) :: ip, jp, kp, l_h,r_h
		real(sp), dimension(-r_h+1:kp+r_h,-r_h+1:jp+r_h,-l_h+1:ip+r_h), &
			intent(in) :: su
		real(sp), dimension(-r_h+1:kp+r_h,-l_h+1:jp+r_h,-r_h+1:ip+r_h), &
			intent(in) :: sv
		real(sp), dimension(-l_h+1:kp+r_h,-r_h+1:jp+r_h,-r_h+1:ip+r_h), &
			intent(in) :: sw
		real(sp), dimension(-r_h+1:kp+r_h,-r_h+1:jp+r_h,-r_h+1:ip+r_h), &
			intent(inout) :: x
		real(sp), dimension(-r_h+1:kp+r_h,-r_h+1:jp+r_h,-r_h+1:ip+r_h), &
			intent(inout) :: b
		real(sp), dimension(-l_h+1:ip+r_h), intent(in) :: xg,dx, dxn
		real(sp), dimension(-l_h+1:jp+r_h), intent(in) :: yg,dy, dyn
		real(sp), dimension(-l_h+1:kp+r_h), intent(in) :: zg,dz, dzn
		logical :: test_solver
	
		! locals
		real(sp), parameter :: tol=1.e-5_sp, tiny=1.e-16 !epsilon(dt)
		integer(i4b), parameter :: itmax=999	, prit=1	
		real(sp) :: sc_err, err, rho, alf, omg, bet, nrm, tt, ts
		real(sp), dimension(-r_h+1:kp+r_h,-r_h+1:jp+r_h,-r_h+1:ip+r_h) :: &
							ax, r,cr,cs,p,v1,t,s
		integer(i4b) :: it,status, i,j,k
		
		
		! allocate some arrays and set matrix A:
		if(bicgstab_first ) then
			!print *,'First call to pressure solver'
			allocate( a_n(1-r_h:kp+r_h,1-r_h:jp+r_h,1-r_h:ip+r_h), STAT = status)
			if (status /= 0) STOP "*** Not enough memory ***"
			allocate( a_s(1-r_h:kp+r_h,1-r_h:jp+r_h,1-r_h:ip+r_h), STAT = status)
			if (status /= 0) STOP "*** Not enough memory ***"
			allocate( a_e(1-r_h:kp+r_h,1-r_h:jp+r_h,1-r_h:ip+r_h), STAT = status)
			if (status /= 0) STOP "*** Not enough memory ***"
			allocate( a_w(1-r_h:kp+r_h,1-r_h:jp+r_h,1-r_h:ip+r_h), STAT = status)
			if (status /= 0) STOP "*** Not enough memory ***"
			allocate( a_u(1-r_h:kp+r_h,1-r_h:jp+r_h,1-r_h:ip+r_h), STAT = status)
			if (status /= 0) STOP "*** Not enough memory ***"
			allocate( a_d(1-r_h:kp+r_h,1-r_h:jp+r_h,1-r_h:ip+r_h), STAT = status)
			if (status /= 0) STOP "*** Not enough memory ***"
			allocate( a_p(1-r_h:kp+r_h,1-r_h:jp+r_h,1-r_h:ip+r_h), STAT = status)
			if (status /= 0) STOP "*** Not enough memory ***"
			allocate( a_c(1-r_h:kp+r_h,1-r_h:jp+r_h,1-r_h:ip+r_h), STAT = status)
			if (status /= 0) STOP "*** Not enough memory ***"
			
			call set_mat_a(comm3d,id,dims,coords,dt,dx,dy,dz,dxn,dyn,dzn,ip,jp,kp,r_h, &
						a_e,a_w,a_n,a_s,a_u,a_d,a_p,a_c,su,sv,sw)
		
			if(test_solver) then
				do i=1,ip
					do j=1,jp
						do k=1,kp
							b(k,j,i) = &
							 exp(-1._sp/5.e6_sp* &
								(xg(i)*xg(i)+yg(j)*yg(j)+(zg(k)-4000._sp)**2._sp))
						enddo
					enddo
				enddo
			endif
			
			bicgstab_first=.false.
		endif
! 		b=b +0.01
	
		sc_err = sqrt( inn_prod(comm3d,b,b,ip,jp,kp,r_h) )
		sc_err = max(sc_err, 0.01_sp)

	
		! calculate the initial residual
		call mat_ax(comm3d,id,dims,coords,ip,jp,kp,r_h, &
		    a_e,a_w,a_n,a_s,a_u,a_d,a_p,a_c,ax,x)
		r  = b - ax
		cr = r
		p  = 0.0_sp
		v1  = 0.0_sp
		
		err = sqrt( inn_prod(comm3d,r,r,ip,jp,kp,r_h) ) / sc_err
		if (err < tol ) return
		
		alf = 1.0_sp
		omg = 1.0_sp
		nrm = 1.0_sp
		
		! start iteration
		do it=1, itmax
			rho = inn_prod(comm3d,r,cr,ip,jp,kp,r_h)
			bet = ( rho/nrm )*( alf/omg )
			t   = r - bet*(omg*v1)
			call precondition(comm3d,id,dims,coords, &
						ip,jp,kp,r_h, a_e,a_w,a_n,a_s,a_u,a_d,a_p,a_c, &
						s,t,prit)
			p   = s + bet*p	
						
			call mat_ax(comm3d,id,dims,coords, &
						ip,jp,kp,r_h, a_e,a_w,a_n,a_s,a_u,a_d,a_p,a_c,v1,p) 
			nrm = inn_prod(comm3d,cr,v1,ip,jp,kp,r_h)
			
			alf = rho/nrm
			s   = r - alf*v1
			call precondition(comm3d,id,dims,coords, &
						ip,jp,kp,r_h, a_e,a_w,a_n,a_s,a_u,a_d,a_p, a_c,&
						cs,s,prit)
			call mat_ax(comm3d,id,dims,coords, &
						ip,jp,kp,r_h, a_e,a_w,a_n,a_s,a_u,a_d,a_p,a_c,t,cs)
			tt  = inn_prod(comm3d,t,t,ip,jp,kp,r_h)
			ts  = inn_prod(comm3d,t,s,ip,jp,kp,r_h)
			omg = ts/tt

			x   = x + alf*p + omg*cs
			r   = s - omg*t
			nrm = rho
			
			if( abs(omg) < tiny ) then
				print *,'convergence problem omega= ',omg,it
				stop
			endif
			
			! check residual for convergence
			err = sqrt( inn_prod(comm3d,r,r,ip,jp,kp,r_h) )/sc_err
			if( err < tol ) exit
			
		enddo
		if( err > tol ) then
			print *,'convergence failure with error,it ',err,it
			stop
		endif
		
		

	end subroutine bicgstab
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>perform an inner product 
	!>@param[in] comm3d: communicator
	!>@param[in] x,y: arrays to perform inner product over
	!>@param[in] ip,jp,kp,halo: dimensions and halos
	function inn_prod(comm3d,x,y,ip,jp,kp,halo)
		use nrtype
		use mpi
		implicit none
		integer(i4b), intent(in) :: comm3d,ip,jp,kp,halo
		real(sp), dimension(-halo+1:kp+halo,-halo+1:jp+halo,-halo+1:ip+halo), &
			intent(in) :: x,y
		real(sp) :: inn_prod,inn_prod1
		real(sp), dimension(-halo+1:kp+halo,-halo+1:jp+halo,-halo+1:ip+halo) :: b
		
		integer(i4b) :: i,j,k, error
		
		inn_prod1=0.0_sp
!$omp simd
		do i=1,ip
			do j=1,jp
				do k=1,kp
					inn_prod1 = inn_prod1 + x(k,j,i)*y(k,j,i)
				enddo
			enddo
		enddo
!$omp end simd
		call mpi_allreduce(inn_prod1,inn_prod,1,MPI_REAL8,MPI_SUM, &
			comm3d,error)
 	
	end function inn_prod
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>find momentum sources and rhs of Poisson equation:
	!>\f$\nabla ^2 f = \rho \f$
	!>@param[in] comm3d,id, rank, dims, coords
	!>@param[in] dt
	!>@param[in] xg,yg,zg,zng,dx,dy,dz,dxn,dyn,dzn
	!>@param[in] ip,jp,kp,l_h,r_h,nq
	!>@param[in] ubar,vbar,wbar,thbar,qbar, dampfacn,dampfac
	!>@param[in] zu
	!>@param[in] zv
	!>@param[in] zw
	!>@param[in] u
	!>@param[in] v
	!>@param[in] w
	!>@param[inout] su,sv,sw,rhs
	!>@param[in] q
	!>@param[inout] sq,viss
	!>@param[in] th,theta,thetan,rhoa,rhoan
	!>@param[inout] sth,strain,vism,vist
	!>@param[in] lamsq,lamsqn: smag mixing length parameter
	!>@param[in] z0,z0th: roughness lengths
	!>@param[in] viscous: logical for applying viscosity
	!>@param[in] moisture: if we have moisture
	!>@param[in] damping_layer: flag for damping layer
	subroutine sources(comm3d,id,rank, dims, coords, &
			dt,xg,yg,zg,zng,dx,dy,dz,dxn,dyn,dzn,ip,jp,kp,l_h,r_h,&
			nq, &
			ubar,vbar,wbar,thbar,qbar, dampfacn,dampfac, &
			zu,zv,zw, &
			u,v,w,su,sv,sw,&
			q,sq,viss, &
			rhs, &
			th,sth,strain,vism,vist,theta,thetan,rhoa,rhoan,lamsq,lamsqn,&
			z0,z0th,&
			viscous, &
			moisture,damping_layer)
		use nrtype
		use mpi_module, only : exchange_full, exchange_along_dim
		use subgrid_3d, only : calculate_subgrid_3d
		implicit none
		integer(i4b), intent(in) :: id, comm3d, rank
		integer(i4b), dimension(3), intent(in) :: dims, coords
		logical, intent(in) :: viscous,moisture,damping_layer
		
		real(sp), intent(in) :: dt,z0,z0th
		integer(i4b), intent(in) :: ip, jp, kp, l_h,r_h,nq
		real(sp), dimension(-r_h+1:kp+r_h,-r_h+1:jp+r_h,-l_h+1:ip+r_h), &
			intent(in) :: u,zu
		real(sp), dimension(-r_h+1:kp+r_h,-l_h+1:jp+r_h,-r_h+1:ip+r_h), &
			intent(in) :: v,zv
		real(sp), dimension(-l_h+1:kp+r_h,-r_h+1:jp+r_h,-r_h+1:ip+r_h), &
			intent(in) :: w,zw
			
		real(sp), dimension(-r_h+1:kp+r_h,-r_h+1:jp+r_h,-r_h+1:ip+r_h), &
			intent(inout) :: su
		real(sp), dimension(-r_h+1:kp+r_h,-r_h+1:jp+r_h,-r_h+1:ip+r_h), &
			intent(inout) :: sv
		real(sp), dimension(-r_h+1:kp+r_h,-r_h+1:jp+r_h,-r_h+1:ip+r_h), &
			intent(inout) :: sw
		real(sp), dimension(-r_h+1:kp+r_h,-r_h+1:jp+r_h,-r_h+1:ip+r_h), &
			intent(inout) :: rhs,sth,strain,vism,vist

		real(sp), dimension(-r_h+1:kp+r_h,-r_h+1:jp+r_h,-r_h+1:ip+r_h,1:nq), &
			intent(in) :: q
		real(sp), dimension(-r_h+1:kp+r_h,-r_h+1:jp+r_h,-r_h+1:ip+r_h,1:nq), &
			intent(inout) :: sq,viss

		real(sp), dimension(-l_h+1:ip+r_h), intent(in) :: xg,dx, dxn
		real(sp), dimension(-l_h+1:jp+r_h), intent(in) :: yg,dy, dyn
		real(sp), dimension(-l_h+1:kp+r_h), intent(in) :: zg,zng,dz, dzn, theta,thetan, &
															rhoa, rhoan,lamsq,lamsqn, &
															ubar,vbar,wbar,thbar, &
															dampfacn,dampfac
		real(sp), dimension(-l_h+1:kp+r_h,nq), intent(in) :: qbar
		real(sp), dimension(-r_h+1:kp+r_h,-r_h+1:jp+r_h,-r_h+1:ip+r_h), &
			intent(in) :: th
	
		! locals
		integer(i4b) :: status, i,j,k
		real(sp), dimension(-r_h+1:kp+r_h,-r_h+1:jp+r_h,-r_h+1:ip+r_h) :: &
			s,tau11,tau12,tau13,tau22,tau23,tau33
		
		sth=0._sp
		! using current time-step
		! sources for u:
!$omp simd	
		do i=2-l_h,ip
			do j=1,jp
				do k=1,kp
					su(k,j,i)=rhoan(k)* ( ( ( u(k,j,i-1)*(u(k,j,i)+u(k,j,i-1)) - &
								u(k,j,i+1)*(u(k,j,i)+u(k,j,i+1)) ) / &
								(2._sp*(dx(i)+dx(i-1))) + &
							  ( u(k,j-1,i)*(v(k,j-1,i)+v(k,j-1,i+1)) - &
								u(k,j+1,i)*(v(k,j,i)+v(k,j,i+1)) ) / &
								(2._sp*(dyn(j)+dyn(j-1))) + &
							  ( u(k-1,j,i)*(w(k-1,j,i)+w(k-1,j,i+1)) - &
								u(k+1,j,i)*(w(k,j,i)+w(k,j,i+1)) ) / &
								(2._sp*(dzn(k)+dzn(k-1))) ) )
				enddo
			enddo
		enddo
!$omp end simd	
		
		! sources for v:
!$omp simd	
		do i=1,ip
			do j=2-l_h,jp
				do k=1,kp
					sv(k,j,i)=rhoan(k)* ( ( ( v(k,j,i-1)*(u(k,j,i-1)+u(k,j+1,i-1)) - &
								v(k,j,i+1)*(u(k,j,i)+u(k,j+1,i)) ) / &
								(2._sp*(dxn(i)+dxn(i-1))) + &
							  ( v(k,j-1,i)*(v(k,j,i)+v(k,j-1,i)) - &
								v(k,j+1,i)*(v(k,j,i)+v(k,j+1,i)) ) / &
								(2._sp*(dy(j)+dy(j-1))) + &
							  ( v(k-1,j,i)*(w(k-1,j,i)+w(k-1,j+1,i)) - &
								v(k+1,j,i)*(w(k,j,i)+w(k,j+1,i)) ) / &
								(2._sp*(dzn(k)+dzn(k-1))) ) )
				enddo
			enddo
		enddo
!$omp end simd	

		! sources for w:
!$omp simd	
		do i=1,ip
			do j=1,jp
				do k=2-l_h,kp
					sw(k,j,i)=rhoa(k)* ( ( ( w(k,j,i-1)*(u(k,j,i-1)+u(k+1,j,i-1)) - &
								w(k,j,i+1)*(u(k,j,i)+u(k+1,j,i)) ) / &
								(2._sp*(dxn(i)+dxn(i-1))) + &
							  ( w(k,j-1,i)*(v(k,j-1,i)+v(k+1,j-1,i)) - &
								w(k,j+1,i)*(v(k,j,i)+v(k+1,j,i)) ) / &
								(2._sp*(dyn(j)+dyn(j-1))) + &
							  ( w(k-1,j,i)*(w(k,j,i)+w(k-1,j,i)) - &
								w(k+1,j,i)*(w(k,j,i)+w(k+1,j,i)) ) / &
								(2._sp*(dz(k)+dz(k-1))) ) + &
								grav*(th(k,j,i)+th(k+1,j,i) ) / &
									(2._sp*theta(k)) )
				enddo
			enddo
		enddo		
!$omp end simd	

        ! damping layer for theta, u,v,w=0 here (need to calculate horizontal means
        ! total sum / MPI_reduce) - at start make sure we know ip and jp on each proc
        if(damping_layer) then
!$omp simd	
            do i=2-l_h,ip
                do j=1,jp
                    do k=1,kp
                        su(k,j,i)=su(k,j,i)
                    enddo
                enddo
            enddo		
!$omp end simd	
!$omp simd	
            do i=1,ip
                do j=2-l_h,jp
                    do k=1,kp
                        sv(k,j,i)=sv(k,j,i)+rhoan(k)*dampfacn(k)*(v(k,j,i)-vbar(k))
                    enddo
                enddo
            enddo		
!$omp end simd	
!$omp simd	
            do i=1,ip
                do j=1,jp
                    do k=2-l_h,kp
                        sw(k,j,i)=sw(k,j,i)+rhoa(k)*dampfac(k)*(w(k,j,i)-wbar(k))
                    enddo
                enddo
            enddo		
!$omp end simd	
!$omp simd	
            do i=1,ip
                do j=1,jp
                    do k=2-l_h,kp
                        sth(k,j,i)=sth(k,j,i)+rhoan(k)*dampfacn(k)*(th(k,j,i)-thbar(k))
                    enddo
                enddo
            enddo		
!$omp end simd	
        endif

        if(viscous) then
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! sub-grid model
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call calculate_subgrid_3d(dt,zg,zng,dx,dy,dz,rhoa, theta,&
                                        dxn,dyn,dzn,rhoan, thetan,&
                                        ip,jp,kp,nq,l_h,r_h,u,v,w,zu,zv, zw, &
                                        th,q, &
                                        lamsq, z0,z0th, &
                                        vism,vist, viss,strain, &
                                        tau11,tau22,tau33,tau12,&
                                        tau13,tau23, &
                                        su,sv,sw,sth,sq, &
                                        moisture, &
                                        comm3d,id,dims,coords)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		else
            ! full exchange    
            call exchange_full(comm3d, id, kp, jp, ip, r_h,r_h,r_h,r_h,l_h,r_h,su,&
                0._sp,0._sp,dims,coords)
            call exchange_full(comm3d, id, kp, jp, ip, r_h,r_h,r_h,r_h,l_h,r_h,sv,&
                0._sp,0._sp,dims,coords)
            call exchange_full(comm3d, id, kp, jp, ip, r_h,r_h,r_h,r_h,l_h,r_h,sw,&
                0._sp,0._sp,dims,coords)		
        endif
        
        

		! su, sv, and su are now centred on i+1/2        

		! rhs of poisson (centred difference):
!$omp simd	
		do i=1,ip
			do j=1,jp
				do k=1,kp
					rhs(k,j,i)=( ( su(k,j,i)-su(k,j,i-1) )/ dx(i-1) + &
							 ( sv(k,j,i)-sv(k,j-1,i) )/ dy(j-1) + &
							 ( sw(k,j,i)-sw(k-1,j,i) )/ dz(k-1) )
				enddo
			enddo
		enddo
!$omp end simd	
		!rhs(:,:,0:1)=0._sp
		!sw(:,:,0:1)=0._sp
		
! 	call exchange_halos(comm3d, id, ip, jp, kp, &
! 			r_h,r_h,r_h,r_h,r_h,r_h, rhs,dims,coords)

	end subroutine sources	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>advance momentum using sources and pressure gradient:
	!>\f$\nabla ^2 f = \rho \f$
	!>@param[in] comm3d,id, rank
	!>@param[in] dt
	!>@param[in] dx,dy,dz,dxn,dyn,dzn
	!>@param[in] ip,jp,kp,l_h,r_h
	!>@param[inout] u,zu
	!>@param[inout] v,zv
	!>@param[inout] w,zw
	!>@param[in] p,su,sv,sw
	!>@param[in] rhoa,rhoan
	!>@param[in] dims, coords
	subroutine advance_momentum(comm3d,id,rank, &
			dt,dx,dy,dz,dxn,dyn,dzn,rhoa,rhoan,ip,jp,kp,l_h,r_h,&
			u,v,w,zu,zv,zw,su,sv,sw,p,dims,coords)
		use nrtype
		implicit none
		integer(i4b), intent(in) :: id, comm3d, rank

		real(sp), intent(in) :: dt
		integer(i4b), intent(in) :: ip, jp, kp, l_h,r_h
		real(sp), dimension(-r_h+1:kp+r_h,-r_h+1:jp+r_h,-l_h+1:ip+r_h), &
			intent(inout) :: u,zu
		real(sp), dimension(-r_h+1:kp+r_h,-l_h+1:jp+r_h,-r_h+1:ip+r_h), &
			intent(inout) :: v,zv
		real(sp), dimension(-l_h+1:kp+r_h,-r_h+1:jp+r_h,-r_h+1:ip+r_h), &
			intent(inout) :: w,zw
			
		real(sp), dimension(-r_h+1:kp+r_h,-r_h+1:jp+r_h,-r_h+1:ip+r_h), &
			intent(in) :: su
		real(sp), dimension(-r_h+1:kp+r_h,-r_h+1:jp+r_h,-r_h+1:ip+r_h), &
			intent(in) :: sv
		real(sp), dimension(-r_h+1:kp+r_h,-r_h+1:jp+r_h,-r_h+1:ip+r_h), &
			intent(in) :: sw
		real(sp), dimension(-r_h+1:kp+r_h,-r_h+1:jp+r_h,-r_h+1:ip+r_h), &
			intent(in) :: p
		real(sp), dimension(-l_h+1:ip+r_h), intent(in) :: dx, dxn
		real(sp), dimension(-l_h+1:jp+r_h), intent(in) :: dy, dyn
		real(sp), dimension(-l_h+1:kp+r_h), intent(in) :: dz, dzn, rhoa, rhoan
		integer(i4b), dimension(3), intent(in) :: dims,coords
	
		! locals
		integer(i4b) :: i,j,k
		
		! advance u-momentum:		
!$omp simd	
		do i=1,ip
			do j=1,jp
				do k=1,kp
					u(k,j,i)=(rhoan(k)*zu(k,j,i)+( su(k,j,i) - &
						( p(k,j,i+1)-p(k,j,i) ) / (dxn(i)) )*dt )/rhoan(k)
				enddo
				if(coords(3)==0) then
                    k=1
                    u(k,j,i)=(rhoan(k)*zu(k,j,i)+( su(k,j,i) - &
                        ( p(k,j,i)+dxn(i)*su(k,j,i)-p(k,j,i) ) / (dxn(i)) )*dt )/rhoan(k)
                endif
			enddo
		enddo
!$omp end simd	

		! advance v-momentum:		
!$omp simd	
		do i=1,ip
			do j=1,jp
				do k=1,kp
					v(k,j,i)=(rhoan(k)*zv(k,j,i)+( sv(k,j,i) - &
						( p(k,j+1,i)-p(k,j,i) ) / (dyn(j)) )*dt ) /rhoan(k)
				enddo
				if(coords(3)==0) then
                    k=1
                    v(k,j,i)=(rhoan(k)*zv(k,j,i)+( sv(k,j,i) - &
                        ( p(k,j,i)+dyn(j)*sv(k,j,i)-p(k,j,i) ) / (dyn(j)) )*dt ) /rhoan(k)	
                endif
			enddo
		enddo
!$omp end simd	

		! advance w-momentum:		
!$omp simd	
		do i=1,ip
			do j=1,jp
				do k=1,kp
					w(k,j,i)=(rhoa(k)*zw(k,j,i)+( sw(k,j,i) - &
						( p(k+1,j,i)-p(k,j,i) ) / (dzn(k)) )*dt) / rhoa(k)
				enddo
			enddo
		enddo
!$omp end simd	
		


	end subroutine advance_momentum	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    end module d_solver
