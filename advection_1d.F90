	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>advection code 
    module advection_s_1d
    use numerics_type
    
    private
    public :: mpdata_1d, mpdata_vec_1d, first_order_upstream_1d
    real(wp), parameter :: small=1e-60_wp
			
	contains
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Simple first order upstream scheme                                                 !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>perform 1 time-step of 1-d first order upstream method 
	!>\f$ \frac{\partial \psi}{\partial t} + \nabla \cdot \bar{v} \psi = 0 \f$
	!>@param[in] dt
	!>@param[in] dzn, rhoa, rhoan
	!>@param[in] kp,l_h,r_h
	!>@param[in] w
	!>@param[inout] psi
	!>@param[in] neumann: boundary condition flag for advection schemes
	subroutine first_order_upstream_1d(dt,dzn,rhoa,rhoan, kp,l_h,r_h,w,psi,neumann)
	use numerics_type
	implicit none
	real(wp), intent(in) :: dt
	integer(i4b), intent(in) :: kp, l_h, r_h
	real(wp), dimension(-l_h+1:kp+r_h), &
		intent(in) :: w
	real(wp), dimension(-r_h+1:kp+r_h), &
		intent(inout) :: psi
	real(wp), dimension(-l_h+1:kp+r_h), intent(in) :: dzn, rhoa, rhoan
	integer(i4b), intent(in) :: neumann
	
	! locals
	real(wp), dimension(kp) :: fz_r, fz_l
	integer(i4b) :: k

!$omp simd	
    do k=1,kp
        fz_r(k)=( (w(k)+abs(w(k)))*rhoa(k)*psi(k)+ &
            (w(k)-abs(w(k)))*rhoan(k+1)*psi(k+1) )*dt/ &
            (2._wp*dzn(k)*rhoa(k))

        fz_l(k)=( (w(k-1)+abs(w(k-1)))*rhoa(k-1)*psi(k-1)+ &
            (w(k-1)-abs(w(k-1)))*rhoan(k)*psi(k) )*dt/ &
            (2._wp*dzn(k)*rhoa(k))
    enddo
!$omp end simd
	
    ! neumann boundary condition top and bottom
    if(neumann==1) then
        fz_l(1)=fz_r(1)
        fz_r(kp)=fz_l(kp)    
    elseif(neumann==2) then
        fz_l(1)=-fz_r(1)
        fz_r(kp)=-fz_l(kp)    
    endif

	! could do a loop here and transport
	psi(1:kp)=psi(1:kp)-(fz_r-fz_l)
	end subroutine first_order_upstream_1d
	
	
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! multi-dimensional advection using the smolarkiewicz scheme                                                 !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>advect using 1st order upstream
	!>then re-advect using 1st order upstream with antidiffusive velocities to
	!>correct diffusiveness of the 1st order upstream and iterate
	!>nft option is also coded
	!>solves the 1-d advection equation:
	!>\f$ \frac{\partial \psi}{\partial t} + \nabla \cdot \bar{v} \psi = 0 \f$
	!>@param[in] dt
	!>@param[in] dz, dzn, rhoa, rhoan
	!>@param[in] kp,l_h,r_h
	!>@param[in] w
	!>@param[inout] psi
	!>@param[in] kord, monotone: order of MPDATA and whether it is monotone
	!>@param[in] neumann: boundary condition flag for advection schemes
	!>@param[in] boundary_cond
	!>@param[in] comm3d, id, dims, coords: mpi variables
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef MPI
	subroutine mpdata_1d(dt,dz,dzn,&
						rhoa,rhoan,kp,l_h,r_h,w,psi_in,kord,monotone,neumann, &
						boundary_cond,comm3d, id, &
						dims,coords)
#else
	subroutine mpdata_1d(dt,dz,dzn,&
						rhoa,rhoan,kp,l_h,r_h,w,psi_in,kord,monotone, neumann, &
						boundary_cond)
#endif
	use numerics_type
#ifdef MPI
	use mpi_module
	use mpi
#endif
	implicit none
	
#ifdef MPI	
	integer(i4b), intent(in) :: id, comm3d
	integer(i4b), dimension(1), intent(in) :: dims, coords
#endif
	real(wp), intent(in) :: dt
	integer(i4b), intent(in) :: kp, l_h, r_h,kord
	real(wp), dimension(-l_h+1:kp+r_h), &
		intent(in), target :: w
	real(wp), dimension(-r_h+1:kp+r_h), &
		intent(inout), target :: psi_in
	real(wp), dimension(-l_h+1:kp+r_h), intent(in) :: dz, dzn, rhoa,rhoan
	logical, intent(in) :: monotone
	integer(i4b), intent(in) :: neumann
	integer(i4b), intent(in) :: boundary_cond
	
	! locals
	real(wp) :: u_div3, u_j_bar3, &
			denom1, denom2, minlocal, minglobal, psi_local_sum, psi_sum
	integer(i4b) :: i,j,k, it, it2, error
	real(wp), dimension(:), pointer :: wt
	real(wp), dimension(:), pointer :: wt_sav
	real(wp), dimension(:), pointer :: psi, psi_old
	real(wp), dimension(-r_h+1:kp+r_h), target :: & 
		u_store1, u_store2
	real(wp), dimension(-l_h+1:kp+r_h), target :: &
		w_store1, w_store2
	real(wp), dimension(-r_h+1:kp+r_h), target :: psi_store
	real(wp), dimension(-r_h+1:kp+r_h) :: &
						psi_k_max,psi_k_min, &
						beta_k_up, beta_k_down
	

	! has to be positive definite
	minlocal=minval(psi_in(1:kp))
#ifdef MPI	
	call mpi_allreduce(minlocal,minglobal,1,MPI_REAL8,MPI_MIN, comm3d,error)
#else
    minglobal=minlocal
#endif
	psi_in=psi_in-minglobal
	
	psi_local_sum=sum(psi_in(1:kp))
#ifdef MPI
	call MPI_Allreduce(psi_local_sum, psi_sum, 1, MPI_REAL8, MPI_SUM, comm3d, error)
#else
    psi_sum=psi_local_sum
#endif
 	if(psi_sum.lt.small) then 
 	    psi_in(:)=psi_in(:)+minglobal
     	return
	endif
	w_store2=0._wp
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! associate pointers to targets                                                      !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	wt => w	
	wt_sav => w_store1
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! if kord > 1 we need a copy of the scalar field                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if (kord > 1) then
		psi_store = psi_in   ! array copy
		psi => psi_store	
	endif
	psi_old     => psi_in	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	it2=0
	do it=1,kord
		
		if(it > 1) then
			it2=it2+1
!$omp simd	
            do k=1,kp
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! calculate the anti-diffusive velocities                        !
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! for divergent flow: eq 38 smolarkiewicz 1984 
                ! last part of w wind:
                u_div3=0._wp 
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! equation 13 page 330 of smolarkiewicz (1984) 
                ! journal of computational physics
                ! second term of w wind:
                u_j_bar3 = 0._wp
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! last update of equation 13
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! w wind:
                wt_sav(k)=(abs(wt(k))*dz(k)-dt*wt(k)*wt(k) ) * &
                    (psi_old(k+1)-psi_old(k) ) / &
                    (psi_old(k+1)+psi_old(k)+small) /dzn(k) - u_j_bar3
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! last update of eq 38 smolarkiewicz 1984
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                wt_sav(k)=wt_sav(k) - 0.25_wp*dt*wt(k) * &
                 ( (wt(k+1)-wt(k-1))/(dz(k-1))+u_div3 )
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


            enddo
!$omp end simd

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! associate pointers to targets                                              !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			if(modulo(it2,2).eq.1) then  
										 
				wt => w_store1
				wt_sav => w_store2
			else if(modulo(it2,2).eq.0) then			
				wt => w_store2
				wt_sav => w_store1
			endif
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 			call exchange_full(comm3d, id, kp, jp, ip, r_h,r_h,r_h,r_h,l_h,r_h, &
! 														ut,dims,coords)
! 			call exchange_full(comm3d, id, kp, jp, ip, l_h,r_h,r_h,r_h,r_h,r_h, &
! 														wt,dims,coords)
            ! halos
!             ut(:,-l_h+1:0)=ut(:,ip-l_h+1:ip)
!             wt(:,ip+1:ip+r_h)=wt(:,1:r_h)
            wt(0)=wt(1)
            wt(kp:kp+1)=wt(kp-1)
			
		endif
		
		
		
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Non oscillatory forward in time (NFT) flux limiter -                           !
        ! Smolarkiewicz and Grabowski (1990)                                             !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if( (it > 1) .and. monotone) then
			it2=it2+1
		
		
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! equations 20a and b of Smolarkiewicz and Grabowski (1990, JCP, 86)         !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$omp simd	
            do k=1,kp			
                ! z direction - note: should the last q in the max/min be q+1?
                psi_k_max(k)=max(psi(k-1),psi(k),psi(k+1), &
                            psi_old(k-1),psi_old(k),psi_old(k+1))
            
                psi_k_min(k)=min(psi(k-1),psi(k),psi(k+1), &
                            psi_old(k-1),psi_old(k),psi_old(k+1))
            enddo
!$omp end simd
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! Equations 19a and b of Smolarkiewicz and Grabowski (1990, JCP)             !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$omp simd	
            do k=1,kp		
                denom1=(dt*((max(wt(k-1),0._wp)*psi_old(k-1)-&
                          min(wt(k),0._wp)*psi_old(k+1))/dz(k-1) &
                          +small))
                          
                denom2=(dt*((max(wt(k),0._wp)*psi_old(k)-&
                          min(wt(k-1),0._wp)*psi_old(k))/dz(k-1) &
                          +small))
                          
                beta_k_up(k)=(psi_k_max(k)-psi_old(k)) / denom1

                          
                beta_k_down(k)=(psi_old(k)-psi_k_min(k)) / denom2
            enddo
!$omp end simd
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! exchange halos for beta_i_up, down
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 			call exchange_along_dim(comm3d, id, kp, jp, ip, r_h,r_h,r_h,r_h,r_h,r_h, &
! 														beta_i_up,dims,coords)
! 			call exchange_along_dim(comm3d, id, kp, jp, ip, r_h,r_h,r_h,r_h,r_h,r_h, &
! 														beta_i_down,dims,coords)
! 			call exchange_along_dim(comm3d, id, kp, jp, ip, r_h,r_h,r_h,r_h,r_h,r_h, &
! 														beta_k_up,dims,coords)
! 			call exchange_along_dim(comm3d, id, kp, jp, ip, r_h,r_h,r_h,r_h,r_h,r_h, &
! 														beta_k_down,dims,coords)
            ! halos
!             beta_i_up(:,-l_h+1:0)       =beta_i_up(:,ip-l_h+1:ip)
!             beta_i_down(:,ip+1:ip+r_h)  =beta_i_down(:,1:r_h)
!             beta_k_up(:,-l_h+1:0)       =beta_k_up(:,ip-l_h+1:ip)
!             beta_k_down(:,ip+1:ip+r_h)  =beta_k_down(:,1:r_h)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

								
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! equation 18 of Smolarkiewicz and Grabowski (1990, JCP, 86)                 !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$omp simd	
            do k=1,kp			
                wt_sav(k)=min(1._wp,beta_k_down(k), &
                                 beta_k_up(k+1))*max(wt(k),0._wp) + &
                              min(1._wp,beta_k_up(k), &
                                 beta_k_down(k+1))*min(wt(k),0._wp)
            enddo
!$omp end simd
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! associate pointers to targets                                              !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			wt => w_store2
			wt_sav => w_store1
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 			call exchange_full(comm3d, id, kp, jp, ip, r_h,r_h,r_h,r_h,l_h,r_h, &
! 														ut,dims,coords)
! 			call exchange_full(comm3d, id, kp, jp, ip, l_h,r_h,r_h,r_h,r_h,r_h, &
! 														wt,dims,coords)
            ! halos
!             ut(:,-l_h+1:0)=ut(:,ip-l_h+1:ip)
!             wt(:,ip+1:ip+r_h)=wt(:,1:r_h)
            wt(0)=wt(1)
            wt(kp:kp+1)=wt(kp-1)

		endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! advect using first order upwind                                                !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		call first_order_upstream_1d(dt,dzn,rhoa,rhoan,&
				kp,l_h,r_h,wt,psi_old,neumann)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


		
		if((it <= kord) .and. (kord >= 1)) then
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! set halos																	 !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
! 			call exchange_full(comm3d, id, kp, jp, ip, r_h,r_h,r_h,r_h,r_h,r_h, &
! 									psi_old,dims,coords)
            ! halos
!             psi_old(:,-l_h+1:0)=psi_old(:,ip-l_h+1:ip)
            select case(boundary_cond)
                case(0) ! nothing top and bottom
                    psi_old(0)=0._wp
                    psi_old(kp+1)=0._wp
                case(1) ! same top and bottom
                    psi_old(0)=psi_old(1)
                    psi_old(kp+1)=psi_old(kp)                  
                case default
                    print *,'no such bc'
                    stop
            end select
            
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!			
		endif		
 	enddo
 

	! no dangling pointers
	if (associated(wt) ) nullify(wt)
	if (associated(wt_sav) ) nullify(wt_sav)
	if (associated(psi) ) nullify(psi)
	if (associated(psi_old) ) nullify(psi_old)

	psi_in=psi_in+minglobal
	end subroutine mpdata_1d
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	



		

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! multi-dimensional vector advection using the smolarkiewicz scheme                                                 !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>advect using 1st order upstream
	!>then re-advect using 1st order upstream with antidiffusive velocities to
	!>correct diffusiveness of the 1st order upstream and iterate
	!>nft option is also coded
	!>solves the 1-d advection equation:
	!>\f$ \frac{\partial \psi}{\partial t} + \nabla \cdot \bar{v} \psi = 0 \f$
	!>@param[in] dt
	!>@param[in] dx,dz, dxn, dyn, dzn, rhoa, rhoan
	!>@param[in] ip,kp,nq,l_h,r_h
	!>@param[in] w
	!>@param[inout] psi
	!>@param[in] kord, monotone: order of MPDATA and whether it is monotone
	!>@param[in] neumann: boundary condition flag for advection schemes
	!>@param[in] comm3d, id, dims, coords: mpi variables
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef MPI
	subroutine mpdata_vec_1d(dt,dz,dzn,&
						rhoa,rhoan,kp,nq,l_h,r_h,w,psi_in,kord,monotone, comm3d, id, &
						neumann,dims,coords)
#else
	subroutine mpdata_vec_1d(dt,dz,dzn,&
						rhoa,rhoan,kp,nq,l_h,r_h,w,psi_in,kord,monotone,neumann)
#endif
	use numerics_type
#ifdef MPI
	use mpi_module
	use mpi
#endif
	implicit none
	
#ifdef MPI
	integer(i4b), intent(in) :: id, comm3d
	integer(i4b), dimension(1), intent(in) :: dims, coords
#endif
	real(wp), intent(in) :: dt
	integer(i4b), intent(in) :: kp, nq,l_h, r_h,kord
	real(wp), dimension(-l_h+1:kp+r_h), &
		intent(in), target :: w
	real(wp), dimension(-r_h+1:kp+r_h,1:nq), &
		intent(inout), target :: psi_in
	real(wp), dimension(-l_h+1:kp+r_h), intent(in) :: dz, dzn, rhoa, rhoan
	logical, intent(in) :: monotone
	integer(i4b), intent(in) :: neumann
	
	! locals
	real(wp) :: u_div3, u_j_bar3, &
			denom1, denom2, minlocal, psi_local_sum, psi_sum
	real(wp), dimension(nq) :: minglobal
	integer(i4b) :: i,k, it, it2, n,error
	real(wp), dimension(:), pointer :: wt
	real(wp), dimension(:), pointer :: wt_sav
	real(wp), dimension(:), pointer :: psi, psi_old
	real(wp), dimension(-l_h+1:kp+r_h), target :: &
		w_store1, w_store2
	real(wp), dimension(-r_h+1:kp+r_h,1:nq), target :: psi_store
	real(wp), dimension(-r_h+1:kp+r_h) :: &
						psi_k_max,psi_k_min, &
						beta_k_up, beta_k_down
	

	! has to be positive definite
	do n=1,nq
        minlocal=minval(psi_in(1:kp,n),1)
#ifdef MPI
        call mpi_allreduce(minlocal,minglobal(n),1,MPI_REAL8,MPI_MIN, comm3d,error)
#else
        minglobal(n)=minlocal
#endif
        psi_in(:,n)=psi_in(:,n)-minglobal(n)
	enddo
	
	psi_local_sum=sum(psi_in(1:kp,1))
#ifdef MPI
	call MPI_Allreduce(psi_local_sum, psi_sum, 1, MPI_REAL8, MPI_SUM, comm3d, error)
#else
    psi_sum=psi_local_sum
#endif
 	if(psi_sum.lt.small) then 
 	    do n=1,nq
     	    psi_in(:,n)=psi_in(:,n)+minglobal(n)
     	enddo
     	return
	endif
	w_store2=0._wp
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! associate pointers to targets                                                      !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	wt => w	
	wt_sav => w_store1
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! if kord > 1 we need a copy of the scalar field                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if (kord > 1) then
		psi_store = psi_in   ! array copy
		psi(-r_h+1:kp+r_h) => psi_store(:,1)
	endif
	psi_old(-r_h+1:kp+r_h)     => psi_in(:,1)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	it2=0
	do it=1,kord
		
		if(it > 1) then
			it2=it2+1
!$omp simd	
            do k=1,kp
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! calculate the anti-diffusive velocities                        !
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! for divergent flow: eq 38 smolarkiewicz 1984 
                ! last part of w wind:
                u_div3=0._wp 
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! equation 13 page 330 of smolarkiewicz (1984) 
                ! journal of computational physics
                ! second term of w wind:
                u_j_bar3 = 0._wp
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! last update of equation 13
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! w wind:
                wt_sav(k)=(abs(wt(k))*dz(k)-dt*wt(k)*wt(k) ) * &
                    (psi_old(k+1)-psi_old(k) ) / &
                    (psi_old(k+1)+psi_old(k)+small) /dzn(k) - u_j_bar3
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! last update of eq 38 smolarkiewicz 1984
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                wt_sav(k)=wt_sav(k) - 0.25_wp*dt*wt(k) * &
                 ( (wt(k+1)-wt(k-1))/(dz(k-1))+u_div3 )
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


            enddo
!$omp end simd

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! associate pointers to targets                                              !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			if(modulo(it2,2).eq.1) then  
										 
				wt => w_store1
				wt_sav => w_store2
			else if(modulo(it2,2).eq.0) then			
				wt => w_store2
				wt_sav => w_store1
			endif
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 			call exchange_full(comm3d, id, kp, jp, ip, r_h,r_h,r_h,r_h,l_h,r_h, &
! 														ut,dims,coords)
! 			call exchange_full(comm3d, id, kp, jp, ip, l_h,r_h,r_h,r_h,r_h,r_h, &
! 														wt,dims,coords)
            ! halos
!             ut(:,-l_h+1:0)=ut(:,ip-l_h+1:ip)
!             wt(:,ip+1:ip+r_h)=wt(:,1:r_h)
            wt(0)=wt(1)
            wt(kp:kp+1)=wt(kp-1)
			
		endif
		
		
		
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Non oscillatory forward in time (NFT) flux limiter -                           !
        ! Smolarkiewicz and Grabowski (1990)                                             !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if( (it > 1) .and. monotone) then
			it2=it2+1
		
		
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! equations 20a and b of Smolarkiewicz and Grabowski (1990, JCP, 86)         !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!$omp simd	
            do k=1,kp			
                ! z direction - note: should the last q in the max/min be q+1?
                psi_k_max(k)=max(psi(k-1),psi(k),psi(k+1), &
                            psi_old(k-1),psi_old(k),psi_old(k+1))
            
                psi_k_min(k)=min(psi(k-1),psi(k),psi(k+1), &
                            psi_old(k-1),psi_old(k),psi_old(k+1))
            enddo
!$omp end simd
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! Equations 19a and b of Smolarkiewicz and Grabowski (1990, JCP)             !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$omp simd	
            do k=1,kp	
                	
                denom1=dt*((max(wt(k-1),0._wp)*psi_old(k-1)-&
                          min(wt(k),0._wp)*psi_old(k+1))/dz(k-1) &
                          +small)
                          
                denom2=dt*((max(wt(k),0._wp)*psi_old(k)-&
                          min(wt(k-1),0._wp)*psi_old(k))/dz(k-1) &
                          +small)
                          

                beta_k_up(k)=(psi_k_max(k)-psi_old(k)) / denom1
                          
                beta_k_down(k)=(psi_old(k)-psi_k_min(k)) / denom2
            enddo
!$omp end simd
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! exchange halos for beta_i_up, down
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 			call exchange_along_dim(comm3d, id, kp, jp, ip, r_h,r_h,r_h,r_h,r_h,r_h, &
! 														beta_i_up,dims,coords)
! 			call exchange_along_dim(comm3d, id, kp, jp, ip, r_h,r_h,r_h,r_h,r_h,r_h, &
! 														beta_i_down,dims,coords)
! 			call exchange_along_dim(comm3d, id, kp, jp, ip, r_h,r_h,r_h,r_h,r_h,r_h, &
! 														beta_k_up,dims,coords)
! 			call exchange_along_dim(comm3d, id, kp, jp, ip, r_h,r_h,r_h,r_h,r_h,r_h, &
! 														beta_k_down,dims,coords)
            ! halos
!             beta_i_up(:,-l_h+1:0)       =beta_i_up(:,ip-l_h+1:ip)
!             beta_i_down(:,ip+1:ip+r_h)  =beta_i_down(:,1:r_h)
!             beta_k_up(:,-l_h+1:0)       =beta_k_up(:,ip-l_h+1:ip)
!             beta_k_down(:,ip+1:ip+r_h)  =beta_k_down(:,1:r_h)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


								
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! equation 18 of Smolarkiewicz and Grabowski (1990, JCP, 86)                 !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$omp simd	
            do k=1,kp			
                wt_sav(k)=min(1._wp,beta_k_down(k), &
                                 beta_k_up(k+1))*max(wt(k),0._wp) + &
                              min(1._wp,beta_k_up(k), &
                                 beta_k_down(k+1))*min(wt(k),0._wp)
            enddo
!$omp end simd
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! associate pointers to targets                                              !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			wt => w_store2
			wt_sav => w_store1
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 			call exchange_full(comm3d, id, kp, jp, ip, r_h,r_h,r_h,r_h,l_h,r_h, &
! 														ut,dims,coords)
! 			call exchange_full(comm3d, id, kp, jp, ip, l_h,r_h,r_h,r_h,r_h,r_h, &
! 														wt,dims,coords)
            ! halos
!             ut(:,-l_h+1:0)=ut(:,ip-l_h+1:ip)
!             wt(:,ip+1:ip+r_h)=wt(:,1:r_h)
            wt(0)=wt(1)
            wt(kp:kp+1)=wt(kp-1)

		endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! advect using first order upwind                                                !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		do n=1,nq
            call first_order_upstream_1d(dt,dzn,rhoa,rhoan,&
                    kp,l_h,r_h,wt,psi_in(:,n),neumann)
            if((it <= kord) .and. (kord >= 1)) then
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! set halos													    		 !
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
    ! 			call exchange_full(comm3d, id, kp, jp, ip, r_h,r_h,r_h,r_h,r_h,r_h, &
    ! 									psi_old,dims,coords)
                ! halos
    !             psi_old(:,-l_h+1:0)=psi_old(:,ip-l_h+1:ip)
                psi_in(0,n)=0._wp
                psi_in(kp+1,n)=0._wp
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!			
    		endif		
        enddo
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


		
 	enddo
 

	! no dangling pointers
	if (associated(wt) ) nullify(wt)
	if (associated(wt_sav) ) nullify(wt_sav)
	if (associated(psi) ) nullify(psi)
	if (associated(psi_old) ) nullify(psi_old)

    do n=1,nq
    	psi_in(:,n)=psi_in(:,n)+minglobal(n)
    enddo
	end subroutine mpdata_vec_1d
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    end module advection_s_1d
