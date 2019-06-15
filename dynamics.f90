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
    public :: bicgstab,sources, advance_momentum
    
	contains
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>set-up matrix A so that Ax is equivalent to \f$\nabla ^2 x\f$ using finite-diff
	!>@param[in] comm3d,id, dims, coords
	!>@param[in] dt,dx,dy,dz,dxn,dyn,dzn, rhoa,rhoan
	!>@param[in] ip,jp,kp,halo
	!>@param[inout] a_e,a_w,a_n,a_s,a_u,a_d,a_p: terms to make up 7-point stencil
	!>@param[in] su,sv,sw,zu,zv,zw
	subroutine set_mat_a(comm3d,id,dims,coords, dt,dx,dy,dz,dxn,dyn,dzn, &
	                    rhoa,rhoan,ip,jp,kp,halo, &
						a_e,a_w,a_n,a_s,a_u,a_d,a_p,a_c,su,sv,sw,zu,zv,zw)
		use nrtype
		use mpi_module
		implicit none
		integer(i4b), intent(in) :: id, comm3d
		integer(i4b), dimension(3), intent(in) :: dims,coords
		integer(i4b), intent(in) :: ip,jp,kp, halo
		real(sp), dimension(1-halo:ip+halo) :: dx,dxn
		real(sp), dimension(1-halo:jp+halo) :: dy,dyn
		real(sp), dimension(1-halo:kp+halo) :: dz,dzn,rhoa,rhoan
		real(sp), dimension(-halo+1:kp+halo,-halo+1:jp+halo,-halo+1:ip+halo), &
			intent(inout) :: a_e,a_w,a_n,a_s,a_u,a_d,a_p, a_c
		real(sp), dimension(-halo+1:kp+halo,-halo+1:jp+halo,-halo+1:ip+halo), &
			intent(in) :: su,sv,sw,zu,zv,zw
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
        !> <b>Del^2 P </b>: <br>
        !>\f$	\frac{\frac{P_{i+1}-P_i}{\Delta x_{n,i}} 
        !> -\frac{P_{i}-P_{i-1}}{\Delta x_{n,i-1}}}{\Delta x_{i-1}} 
        !> +	\frac{\frac{P_{j+1}-P_j}{\Delta y_{n,j}} 
        !> -\frac{P_{j}-P_{j-1}}{\Delta y_{n,j-1}}}{\Delta y_{j-1}} 
        !> +	\frac{\frac{P_{k+1}-P_k}{\Delta z_{n,k}} 
        !> -\frac{P_{k}-P_{k-1}}{\Delta z_{n,k-1}}}{\Delta z_{k-1}} 
        !>\f$
		!> <br><br> for the upper and lower boundary assume w's are equal and opposite
		!> as w is staggered on z
		!> <br><br> From the vertical equation of motion we have:<br>
		!>\f$ w_{k+1/2}=zw_{k+1/2}+\frac{1}{\rho _{k+1/2}}\left( sw_{k+1/2} + 
		!> \frac{P_k-P_{k+1}}{\Delta z_k}\right)\times 2\Delta t
		!>\f$
		!<br> and <br>
		!>\f$ w_{k-1/2}=zw_{k-1/2}+\frac{1}{\rho _{k-1/2}}\left( sw_{k-1/2} + 
		!> \frac{P_{k-1}-P_{k}}{\Delta z_{k-1}}\right )\times 2\Delta t
		!>\f$
		!> <br> set the sum of the two equations to be equal and opposite<br>
		
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
                !> <br><br>
                !> For the bottom b.c. we have:<br>
                !> \f$ P_{-1}=P_0-\Delta z_{-1} sw_{-1/2}-
                !> \frac{\Delta z_{-1}\rho_{-1/2}}{2\Delta t}\left(zw_{1/2}+zw_{-1/2}+ 
                !> \frac{2\Delta t}{\rho _{1/2}}\left(sw_{1/2}+
                !>  \frac{P_0-P_1}{\Delta z_0}\right)\right)
                !> \f$
                !> <br>
                !> so to compute \f$P_{-1}\f$
                !> the a_c (constant term is the part that doesn't depend on P)<br>
                !> there are an additional two a_p terms (the \f$P_0\f$) terms <br>
                !> and there is an additional a_u term (the \f$P_1\f$) term
				if(coords(3)==0) then
				    k=1
					a_e(k,j,i)=0._sp !1._sp/( dx(i-1)*dxn(i) )
					a_w(k,j,i)=0._sp !1._sp/( dx(i-1)*dxn(i-1) )
					a_n(k,j,i)=0._sp !1._sp/( dy(j-1)*dyn(j) )
					a_s(k,j,i)=0._sp !1._sp/( dy(j-1)*dyn(j-1) )
					a_u(k,j,i)=1._sp/( dz(k-1)*dzn(k) )
					a_d(k,j,i)=0._sp !1._sp/( dz(k-1)*dzn(k-1) )
					a_p(k,j,i)=-1._sp/( dz(k-1)*dzn(k) ) -1._sp/( dz(k-1)*dzn(k-1) ) 

                    ! constant component of P-1
                    a_c(k,j,i)=-(dzn(k-1)*sw(k-1,j,i)+ &
                            dzn(k-1)*rhoa(k-1)/(2._sp*dt)*&
                            (zw(k,j,i)+zw(k-1,j,i)+2._sp*dt/rhoa(k)*(sw(k,j,i))) ) / &
                            ( dz(k-1)*dzn(k-1) )
                    ! in derivative format
                    a_p(k,j,i)=a_p(k,j,i)+(1._sp - &
                        dzn(k-1)*rhoa(k-1)/rhoa(k)/dzn(k))  /( dz(k-1)*dzn(k-1) )
                    ! in derivative format
                    a_u(k,j,i)= a_u(k,j,i)+ &
                        dzn(k-1)*rhoa(k-1)/rhoa(k)/dzn(k)  /( dz(k-1)*dzn(k-1) )
                                
                endif
                !> <br><br>
                !> For the top b.c. we have:<br>
                !> \f$ P_{kp+1}=P_{kp}+\Delta z_{kp} sw_{kp+1/2}-
                !> \frac{\Delta z_{kp}\rho_{kp+1/2}}{2\Delta t}\left(zw_{kp+1/2}+zw_{kp-1/2}+ 
                !> \frac{2\Delta t}{\rho _{kp-1/2}}\left(sw_{kp-1/2}+
                !>  \frac{P_{kp-1}-P_{kp}}{\Delta z_{kp-1}}\right)\right)
                !> \f$
                !> <br>
                !> so to compute \f$P_{kp+1}\f$
                !> the a_c (constant term is the part that doesn't depend on P)<br>
                !> there are an additional two a_p terms (the \f$P_{kp}\f$) terms <br>
                !> and there is an additional a_d term (the \f$P_{kp-1}\f$) term
				if(coords(3)==(dims(3)-1)) then
				    k=kp
					a_e(k,j,i)=0._sp/( dx(i-1)*dxn(i) )
					a_w(k,j,i)=0._sp/( dx(i-1)*dxn(i-1) )
					a_n(k,j,i)=0._sp/( dy(j-1)*dyn(j) )
					a_s(k,j,i)=0._sp/( dy(j-1)*dyn(j-1) )
					a_u(k,j,i)=0._sp !1._sp/( dz(k-1)*dzn(k) )
					a_d(k,j,i)=1._sp/( dz(k-1)*dzn(k-1) )
					a_p(k,j,i)=-0._sp/( dx(i-1)*dxn(i) ) -0._sp/( dx(i-1)*dxn(i-1) ) &
                            -0._sp/( dy(j-1)*dyn(j) ) -0._sp/( dy(j-1)*dyn(j-1) )  &
                            -1._sp/( dz(k-1)*dzn(k) ) -1._sp/( dz(k-1)*dzn(k-1) ) 

                    ! constant component of Pkp+1
                    a_c(k,j,i)=(dzn(k)*sw(k,j,i)+ &
                            dzn(k)*rhoa(k)/(2._sp*dt)*&
                            (zw(k,j,i)+zw(k-1,j,i)+2._sp*dt/rhoa(k-1)*(sw(k-1,j,i))) ) / &
                            ( dz(k-1)*dzn(k) )

                    ! in derivative format
                    a_p(k,j,i)=a_p(k,j,i)+(1._sp - &
                        dzn(k)*rhoa(k)/rhoa(k-1)/dzn(k-1))  /( dz(k-1)*dzn(k) )
                    ! in derivative format
                    a_d(k,j,i)= a_d(k,j,i)+ &
                        dzn(k)*rhoa(k)/rhoa(k-1)/dzn(k-1)  /( dz(k-1)*dzn(k) )                         
                endif
			enddo
		enddo
!$omp end simd



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
                                a_p(k,j,i)*x(k,j,i) + a_c(k,j,i)
				enddo
			enddo
		enddo
!$omp end simd


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
	!>@param[in] xg,yg,zg,dx,dy,dz,dxn,dyn,dzn,rhoa,rhoan
	!>@param[in] ip,jp,kp,l_h,r_h
	!>@param[in] su,zu
	!>@param[in] sv,zv
	!>@param[in] sw,zw
	!>@param[inout] x
	!>@param[inout] b
	!>@param[in] ptol
	!>@param[in] test_solver	
	subroutine bicgstab(comm3d,id,rank, dims, coords, &
			dt,xg,yg,zg,dx,dy,dz,dxn,dyn,dzn,rhoa,rhoan, &
			ip,jp,kp,l_h,r_h,su,sv,sw,zu,zv,zw,x,b,tol,&
			test_solver)
		use nrtype
		implicit none
		integer(i4b), intent(in) :: id, comm3d, rank
		integer(i4b), dimension(3), intent(in) :: dims,coords
		real(sp), intent(in) :: dt
		integer(i4b), intent(in) :: ip, jp, kp, l_h,r_h
		real(sp), dimension(-r_h+1:kp+r_h,-r_h+1:jp+r_h,-l_h+1:ip+r_h), &
			intent(in) :: su,zu
		real(sp), dimension(-r_h+1:kp+r_h,-l_h+1:jp+r_h,-r_h+1:ip+r_h), &
			intent(in) :: sv,zv
		real(sp), dimension(-l_h+1:kp+r_h,-r_h+1:jp+r_h,-r_h+1:ip+r_h), &
			intent(in) :: sw,zw
		real(sp), dimension(-r_h+1:kp+r_h,-r_h+1:jp+r_h,-r_h+1:ip+r_h), &
			intent(inout) :: x
		real(sp), dimension(-r_h+1:kp+r_h,-r_h+1:jp+r_h,-r_h+1:ip+r_h), &
			intent(inout) :: b
		real(sp), dimension(-l_h+1:ip+r_h), intent(in) :: xg,dx, dxn
		real(sp), dimension(-l_h+1:jp+r_h), intent(in) :: yg,dy, dyn
		real(sp), dimension(-l_h+1:kp+r_h), intent(in) :: zg,dz, dzn,rhoa,rhoan
		logical :: test_solver
		real(sp), intent(in) :: tol
	
		! locals
		real(sp), parameter :: tiny=1.e-16 !epsilon(dt)
		integer(i4b), parameter :: itmax=9999	, prit=1	
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
			
			call set_mat_a(comm3d,id,dims,coords,dt,dx,dy,dz,dxn,dyn,dzn, &
			            rhoa,rhoan,ip,jp,kp,r_h, &
						a_e,a_w,a_n,a_s,a_u,a_d,a_p,a_c,su,sv,sw, &
						zu,zv,zw)
		
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

			if( abs(omg) < tiny ) then
				print *,'convergence problem omega= ',omg,it
				stop
			endif

			x   = x + alf*p + omg*cs
			r   = s - omg*t

			! check residual for convergence
			err = sqrt( inn_prod(comm3d,r,r,ip,jp,kp,r_h) )/sc_err

			nrm = rho
			
			
			if( err < tol ) exit
			
		enddo
		if( err > tol ) then
			print *,'convergence failure with error,it ',err,it
			stop
		endif
		
		if(coords(3)==(dims(3)-1)) then
		    x(kp+1,:,:)=x(kp,:,:)+dzn(kp)*sw(kp,:,:)+dzn(kp)*rhoa(kp)/(2._sp*dt) * &
		        (zw(kp,:,:)+zw(kp-1,:,:)+&
		        2._sp*dt/rhoa(kp-1)*(sw(kp-1,:,:)+(x(kp-1,:,:)-x(kp,:,:))/dzn(kp-1)))
		endif
		if(coords(3)==0) then
		    x(0,:,:)=x(1,:,:)-dzn(0)*sw(0,:,:)-dzn(0)*rhoa(0)/(2._sp*dt) * &
		        (zw(1,:,:)+zw(0,:,:)+&
		        2._sp*dt/rhoa(1)*(sw(1,:,:)+(x(1,:,:)-x(2,:,:))/dzn(1)))
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
									(2._sp*theta(k) +thbar(k)+thbar(K+1)) )
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
                        su(k,j,i)=su(k,j,i)+rhoan(k)*dampfacn(k)*(zu(k,j,i)-ubar(k))
                    enddo
                enddo
            enddo		
!$omp end simd	
!$omp simd	
            do i=1,ip
                do j=2-l_h,jp
                    do k=1,kp
                        sv(k,j,i)=sv(k,j,i)+rhoan(k)*dampfacn(k)*(zv(k,j,i)-vbar(k))
                    enddo
                enddo
            enddo		
!$omp end simd	
!$omp simd	
            do i=1,ip
                do j=1,jp
                    do k=2-l_h,kp
                        sw(k,j,i)=sw(k,j,i)+rhoa(k)*dampfac(k)*(zw(k,j,i)-0._sp)
                    enddo
                enddo
            enddo		
!$omp end simd	
!$omp simd	
            do i=1,ip
                do j=1,jp
                    do k=2-l_h,kp
                        sth(k,j,i)=sth(k,j,i)+dampfacn(k)*(th(k,j,i)-0._sp)
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
                                   
            ! vertical diffusion of reference state - is this the wrong sign in the LEM?     
!             do i=1,ip
!                 do j=1,jp
!                     do k=2-l_h,kp
!                         sth(k,j,i)=sth(k,j,i)- &
!                             (vist(k,j,i)*(theta(k+1)-theta(k))*rhoa(k)/rhoan(k)/(dz(k-1)*dzn(k))  - &
!                             vist(k-1,j,i)*(theta(k)-theta(k-1))*rhoa(k-1)/rhoan(k)/(dz(k-1)*dzn(k-1)))
!                     enddo
!                 enddo
!             enddo		
!             call exchange_full(comm3d, id, kp, jp, ip, r_h,r_h,r_h,r_h,l_h,r_h,sth,&
!                 0._sp,0._sp,dims,coords)		
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
        if(coords(3)==0) sw(0,:,:)=-sw(1,:,:)/rhoa(1)*rhoa(0)
        if(coords(3)==(dims(3)-1)) then 
            sw(kp,:,:)=-sw(kp-1,:,:)/rhoa(kp-1)*rhoa(kp)
            sw(kp+1,:,:)=0._sp
        endif
!         if(coords(3)==(dims(3)-1)) su(kp:kp+1,:,:)=0._sp
!         if(coords(3)==(dims(3)-1)) sv(kp:kp+1,:,:)=0._sp
        
		! rhs of poisson (centred difference):
!$omp simd	
		do i=1,ip
			do j=1,jp
				do k=1,kp+1
					rhs(k,j,i)=( ( su(k,j,i)-su(k,j,i-1) )/ dx(i-1) + &
							 ( sv(k,j,i)-sv(k,j-1,i) )/ dy(j-1) + &
							 ( sw(k,j,i)-sw(k-1,j,i) )/ dz(k-1) )
				enddo
				if(coords(3)==0) then
                    rhs(0,j,i)=rhs(1,j,i)
                endif
			enddo
		enddo
!$omp end simd	


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
				do k=0,kp
					u(k,j,i)=zu(k,j,i)+(( su(k,j,i) - &
						( p(k,j,i+1)-p(k,j,i) ) / (dxn(i)) )*dt )/rhoan(k)
				enddo
! 				if(coords(3)==0) then
!                     k=1
!                     u(k,j,i)=zu(k,j,i)+(( su(k,j,i) - &
!                         ( p(k,j,i)+dxn(i)*su(k,j,i)-p(k,j,i) ) / (dxn(i)) )*dt )/rhoan(k)
!                 endif
			enddo
		enddo
!$omp end simd	

		! advance v-momentum:		
!$omp simd	
		do i=1,ip
			do j=1,jp
				do k=0,kp
					v(k,j,i)=zv(k,j,i)+(( sv(k,j,i) - &
						( p(k,j+1,i)-p(k,j,i) ) / (dyn(j)) )*dt ) /rhoan(k)
				enddo
! 				if(coords(3)==0) then
!                     k=1
!                     v(k,j,i)=zv(k,j,i)+(( sv(k,j,i) - &
!                         ( p(k,j,i)+dyn(j)*sv(k,j,i)-p(k,j,i) ) / (dyn(j)) )*dt ) /rhoan(k)	
!                 endif
			enddo
		enddo
!$omp end simd	

		! advance w-momentum:		
!$omp simd	
		do i=1,ip
			do j=1,jp
				do k=0,kp
					w(k,j,i)=zw(k,j,i)+(( sw(k,j,i) - &
						( p(k+1,j,i)-p(k,j,i) ) / (dzn(k)) )*dt) / rhoa(k)
				enddo
			enddo
		enddo
!$omp end simd	
		


	end subroutine advance_momentum	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    end module d_solver
