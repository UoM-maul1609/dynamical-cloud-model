	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>advection code 
    module advection_s_2d
    use nrtype
    
    private
    public :: mpdata_2d, mpdata_vec_2d, first_order_upstream_2d
    real(sp), parameter :: small=1e-60_sp
    
	contains
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Simple first order upstream scheme                                                 !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>perform 1 time-step of 2-d first order upstream method 
	!>\f$ \frac{\partial \psi}{\partial t} + \nabla \cdot \bar{v} \psi = 0 \f$
	!>@param[in] dt
	!>@param[in] dxn,dzn, rhoa, rhoan
	!>@param[in] ip,kp,l_h,r_h
	!>@param[in] u
	!>@param[in] w
	!>@param[inout] psi
	subroutine first_order_upstream_2d(dt,dxn,dzn,rhoa,rhoan,ip,kp,l_h,r_h,u,w,psi)
	use nrtype
	implicit none
	real(sp), intent(in) :: dt
	integer(i4b), intent(in) :: ip, kp, l_h, r_h
	real(sp), dimension(-r_h+1:kp+r_h,-l_h+1:ip+r_h), &
		intent(in) :: u
	real(sp), dimension(-l_h+1:kp+r_h,-r_h+1:ip+r_h), &
		intent(in) :: w
	real(sp), dimension(-r_h+1:kp+r_h,-r_h+1:ip+r_h), &
		intent(inout) :: psi
	real(sp), dimension(-l_h+1:ip+r_h), intent(in) :: dxn
	real(sp), dimension(-l_h+1:kp+r_h), intent(in) :: dzn, rhoa, rhoan
	
	! locals
	real(sp), dimension(kp,ip) :: fx_r, fx_l, fz_r, fz_l
	integer(i4b) :: i,k

!$omp simd	
	do i=1,ip
        do k=1,kp
            ! Flux going out of right cell boundary into adjacent cell-x 
            fx_r(k,i)=( (u(k,i)+abs(u(k,i)))*psi(k,i)+ &
                (u(k,i)-abs(u(k,i)))*psi(k,i+1) )*dt/ &
                (2._sp*dxn(i))
            ! Flux going through left cell boundary from adjacent cell-x
            fx_l(k,i)=( (u(k,i-1)+abs(u(k,i-1)))*psi(k,i-1)+ &
                (u(k,i-1)-abs(u(k,i-1)))*psi(k,i) )*dt/ &
                (2._sp*dxn(i))
    
    
            fz_r(k,i)=( (w(k,i)+abs(w(k,i)))*rhoa(k)*psi(k,i)+ &
                (w(k,i)-abs(w(k,i)))*rhoa(k+1)*psi(k+1,i) )*dt/ &
                (2._sp*dzn(k)*rhoa(k))
    
            fz_l(k,i)=( (w(k-1,i)+abs(w(k-1,i)))*rhoan(k-1)*psi(k-1,i)+ &
                (w(k-1,i)-abs(w(k-1,i)))*rhoan(k)*psi(k,i) )*dt/ &
                (2._sp*dzn(k)*rhoa(k))
        enddo
	enddo
!$omp end simd
	
	! could do a loop here and transport
	psi(1:kp,1:ip)=psi(1:kp,1:ip)-(fx_r-fx_l)-(fz_r-fz_l)
	end subroutine first_order_upstream_2d
	
	
	
	
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
	!>solves the 2-d advection equation:
	!>\f$ \frac{\partial \psi}{\partial t} + \nabla \cdot \bar{v} \psi = 0 \f$
	!>@param[in] dt
	!>@param[in] dx,dz, dxn, dzn, rhoa, rhoan
	!>@param[in] ip,kp,l_h,r_h
	!>@param[in] u
	!>@param[in] w
	!>@param[inout] psi
	!>@param[in] kord, monotone: order of MPDATA and whether it is monotone
	!>@param[in] comm3d, id, dims, coords: mpi variables
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef MPI
	subroutine mpdata_2d(dt,dx,dz,dxn,dzn,&
						rhoa,rhoan, ip,kp,l_h,r_h,u,w,psi_in,kord,monotone, comm3d, id, &
						dims,coords)
#else
	subroutine mpdata_2d(dt,dx,dz,dxn,dzn,&
						rhoa,rhoan, ip,kp,l_h,r_h,u,w,psi_in,kord,monotone)
#endif
	use nrtype
#ifdef MPI
	use mpi_module
	use mpi
#endif
	
	implicit none
	
#ifdef MPI	
	integer(i4b), intent(in) :: id, comm3d
	integer(i4b), dimension(2), intent(in) :: dims, coords
#endif
	real(sp), intent(in) :: dt
	integer(i4b), intent(in) :: ip, kp, l_h, r_h,kord
	real(sp), dimension(-r_h+1:kp+r_h,-l_h+1:ip+r_h), &
		intent(in), target :: u
	real(sp), dimension(-l_h+1:kp+r_h,-r_h+1:ip+r_h), &
		intent(in), target :: w
	real(sp), dimension(-r_h+1:kp+r_h,-r_h+1:ip+r_h), &
		intent(inout), target :: psi_in
	real(sp), dimension(-l_h+1:ip+r_h), intent(in) :: dx, dxn
	real(sp), dimension(-l_h+1:kp+r_h), intent(in) :: dz, dzn, rhoa, rhoan
	logical :: monotone
	
	! locals
	real(sp) :: u_div1, u_div3, u_j_bar1, u_j_bar3, &
			denom1, denom2, minlocal, minglobal, psi_local_sum, psi_sum
	integer(i4b) :: i,j,k, it, it2, error
	real(sp), dimension(:,:), pointer :: ut
	real(sp), dimension(:,:), pointer :: wt
	real(sp), dimension(:,:), pointer :: ut_sav
	real(sp), dimension(:,:), pointer :: wt_sav
	real(sp), dimension(:,:), pointer :: psi, psi_old
	real(sp), dimension(-r_h+1:kp+r_h,-l_h+1:ip+r_h), target :: & 
		u_store1, u_store2
	real(sp), dimension(-l_h+1:kp+r_h,-r_h+1:ip+r_h), target :: &
		w_store1, w_store2
	real(sp), dimension(-r_h+1:kp+r_h,-r_h+1:ip+r_h), target :: psi_store
	real(sp), dimension(-r_h+1:kp+r_h,-r_h+1:ip+r_h) :: &
						psi_i_max, psi_i_min, psi_k_max,psi_k_min, &
						beta_i_up, beta_i_down,&
						beta_k_up, beta_k_down
	

	! has to be positive definite
	minlocal=minval(psi_in(1:kp,1:ip))
#ifdef MPI
	call mpi_allreduce(minlocal,minglobal,1,MPI_REAL8,MPI_MIN, comm3d,error)
#else
    minglobal=minlocal
#endif
	psi_in=psi_in-minglobal
	
	psi_local_sum=sum(psi_in(1:kp,1:ip))
#ifdef MPI
	call MPI_Allreduce(psi_local_sum, psi_sum, 1, MPI_REAL8, MPI_SUM, comm3d, error)
#else
    psi_sum=psi_local_sum
#endif
 	if(psi_sum.lt.small) then
 	    psi_in(:,:)=psi_in(:,:)+minglobal

 	    return
	endif
	u_store2=0._sp
	w_store2=0._sp
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! associate pointers to targets                                                      !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	ut => u	
	wt => w	
	ut_sav => u_store1
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
			do i=1,ip
                do k=1,kp
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! calculate the anti-diffusive velocities                        !
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! for divergent flow: eq 38 smolarkiewicz 1984 
                    ! last part of u wind:
                    u_div1=(wt(k,i)+wt(k,i+1)-wt(k-1,i)-wt(k-1,i+1)) &
                            / dz(k-1) 
                    ! for divergent flow: eq 38 smolarkiewicz 1984 
                    ! last part of w wind:
                    u_div3=(ut(k,i)+ut(k+1,i)-ut(k,i-1)-ut(k+1,i-1)) &
                            / dx(i-1) 
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! equation 13 page 330 of smolarkiewicz (1984) 
                    ! journal of computational physics
                    ! second term of u wind:
                    u_j_bar1 = 0.5_sp*dt*ut(k,i) * ( &
                        ! equation 14:
                        0.25_sp*(wt(k,i+1)+wt(k,i)+wt(k-1,i+1)+wt(k-1,i)) * &
                        ! equation 13:
                       ( psi_old(k+1,i+1)+psi_old(k+1,i)- &
                         psi_old(k-1,i+1)-psi_old(k-1,i) ) / &
                       ( psi_old(k+1,i+1)+psi_old(k+1,i)+ &
                         psi_old(k-1,i+1)+psi_old(k-1,i)+small ) / &
                         ( 0.5_sp*(dzn(k-1)+dzn(k)) ) )
                         
                                                  
                    ! equation 13 page 330 of smolarkiewicz (1984) 
                    ! journal of computational physics
                    ! second term of w wind:
                    u_j_bar3 = 0.5_sp*dt*wt(k,i) * ( &
                        ! repeat for y dimension:
                        ! equation 14:
                        0.25_sp*(ut(k+1,i)+ut(k,i)+ut(k+1,i-1)+ut(k,i-1)) * &
                        ! equation 13:
                       ( psi_old(k+1,i+1)+psi_old(k,i+1)- &
                         psi_old(k+1,i-1)-psi_old(k,i-1) ) / &
                       ( psi_old(k+1,i+1)+psi_old(k,i+1)+ &
                         psi_old(k+1,i-1)+psi_old(k,i-1)+small ) / &
                         ( 0.5_sp*(dxn(i-1)+dxn(i)) ) )
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! last update of equation 13
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! u wind:
                    ut_sav(k,i)=(abs(ut(k,i))*dx(i)-dt*ut(k,i)*ut(k,i) ) * &
                        (psi_old(k,i+1)-psi_old(k,i) ) / &
                        (psi_old(k,i+1)+psi_old(k,i)+small) /dxn(i) - u_j_bar1
                    ! w wind:
                    wt_sav(k,i)=(abs(wt(k,i))*dz(k)-dt*wt(k,i)*wt(k,i) ) * &
                        (psi_old(k+1,i)-psi_old(k,i) ) / &
                        (psi_old(k+1,i)+psi_old(k,i)+small) /dzn(k) - u_j_bar3
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! last update of eq 38 smolarkiewicz 1984
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ut_sav(k,i)=ut_sav(k,i) - 0.25_sp*dt*ut(k,i) * &
                     ( (ut(k,i+1)-ut(k,i-1))/(dx(i-1))+u_div1 )
                    wt_sav(k,i)=wt_sav(k,i) - 0.25_sp*dt*wt(k,i) * &
                     ( (wt(k+1,i)-wt(k-1,i))/(dz(k-1))+u_div3 )
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


                enddo
			enddo
!$omp end simd

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! associate pointers to targets                                              !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			if(modulo(it2,2).eq.1) then  
										 
				ut => u_store1
				wt => w_store1
				ut_sav => u_store2
				wt_sav => w_store2
			else if(modulo(it2,2).eq.0) then			
				ut => u_store2
				wt => w_store2
				ut_sav => u_store1
				wt_sav => w_store1
			endif
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 			call exchange_full(comm3d, id, kp, jp, ip, r_h,r_h,r_h,r_h,l_h,r_h, &
! 														ut,dims,coords)
! 			call exchange_full(comm3d, id, kp, jp, ip, l_h,r_h,r_h,r_h,r_h,r_h, &
! 														wt,dims,coords)
            ! halos
            ut(:,-l_h+1:0)      =ut(:,ip-l_h+1:ip)
            ut(:,ip+1:ip+r_h)   =ut(:,1:r_h)
            ut(0,:)=0._sp
            ut(kp+1,:)=0._sp
            
            wt(:,-l_h+1:0)      =wt(:,ip-l_h+1:ip)
            wt(:,ip+1:ip+r_h)   =wt(:,1:r_h)
            wt(0,:)=0._sp
            wt(kp:kp+1,:)=0._sp
			
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
			do i=1,ip
                do k=1,kp			
                    ! x direction - note: should the last q in the max/min be q+1?
                    psi_i_max(k,i)=max(psi(k,i-1),psi(k,i),psi(k,i+1), &
                                psi_old(k,i-1),psi_old(k,i),psi_old(k,i+1))
                
                    psi_i_min(k,i)=min(psi(k,i-1),psi(k,i),psi(k,i+1), &
                                psi_old(k,i-1),psi_old(k,i),psi_old(k,i+1))
                enddo
			enddo
!$omp end simd


!$omp simd	
			do i=1,ip
                do k=1,kp			
                    ! z direction - note: should the last q in the max/min be q+1?
                    psi_k_max(k,i)=max(psi(k-1,i),psi(k,i),psi(k+1,i), &
                                psi_old(k-1,i),psi_old(k,i),psi_old(k+1,i))
                
                    psi_k_min(k,i)=min(psi(k-1,i),psi(k,i),psi(k+1,i), &
                                psi_old(k-1,i),psi_old(k,i),psi_old(k+1,i))
                enddo
			enddo
!$omp end simd
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! Equations 19a and b of Smolarkiewicz and Grabowski (1990, JCP)             !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$omp simd	
			do i=1,ip
                do k=1,kp		
                    denom1=(dt*((max(ut(k,i-1),0._sp)*psi_old(k,i-1)- &
                              min(ut(k,i),0._sp)*psi_old(k,i+1))/dx(i-1)+ &
                            (max(wt(k-1,i),0._sp)*psi_old(k-1,i)-&
                              min(wt(k,i),0._sp)*psi_old(k+1,i))/dz(k-1) &
                              +small))
                              
                    denom2=(dt*((max(ut(k,i),0._sp)*psi_old(k,i)- &
                              min(ut(k,i-1),0._sp)*psi_old(k,i))/dx(i-1) + &
                            (max(wt(k,i),0._sp)*psi_old(k,i)-&
                              min(wt(k-1,i),0._sp)*psi_old(k,i))/dz(k-1) &
                              +small))
                              
                    beta_i_up(k,i)=(psi_i_max(k,i)-psi_old(k,i)) / denom1
                        
                              
                    beta_i_down(k,i)=(psi_old(k,i)-psi_i_min(k,i)) / denom2
                                        

                    beta_k_up(k,i)=(psi_k_max(k,i)-psi_old(k,i)) / denom1

                              
                    beta_k_down(k,i)=(psi_old(k,i)-psi_k_min(k,i)) / denom2
                enddo
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
            beta_i_up(:,-l_h+1:0)       =beta_i_up(:,ip-l_h+1:ip)
            beta_i_up(:,ip+1:ip+r_h)    =beta_i_up(:,1:r_h)
            beta_i_up(0,:)=0._sp
            beta_i_up(kp+1,:)=0._sp
            
            beta_i_down(:,-l_h+1:0)     =beta_i_down(:,ip-l_h+1:ip)
            beta_i_down(:,ip+1:ip+r_h)  =beta_i_down(:,1:r_h)
            beta_i_down(0,:)=0._sp
            beta_i_down(kp+1,:)=0._sp
            
            beta_k_up(:,-l_h+1:0)       =beta_k_up(:,ip-l_h+1:ip)
            beta_k_up(:,ip+1:ip+r_h)    =beta_k_up(:,1:r_h)
            beta_k_up(0,:)=0._sp
            beta_k_up(kp+1,:)=0._sp
            
            beta_k_down(:,-l_h+1:0)     =beta_k_down(:,ip-l_h+1:ip)
            beta_k_down(:,ip+1:ip+r_h)  =beta_k_down(:,1:r_h)
            beta_k_down(0,:)=0._sp
            beta_k_down(kp+1,:)=0._sp
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

								
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! equation 18 of Smolarkiewicz and Grabowski (1990, JCP, 86)                 !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$omp simd	
			do i=1,ip
                do k=1,kp			
                    ut_sav(k,i)=min(1._sp,beta_i_down(k,i), &
                                     beta_i_up(k,i+1))*max(ut(k,i),0._sp) + &
                                  min(1._sp,beta_i_up(k,i), &
                                     beta_i_down(k,i+1))*min(ut(k,i),0._sp)
                    wt_sav(k,i)=min(1._sp,beta_k_down(k,i), &
                                     beta_k_up(k+1,i))*max(wt(k,i),0._sp) + &
                                  min(1._sp,beta_k_up(k,i), &
                                     beta_k_down(k+1,i))*min(wt(k,i),0._sp)
                enddo
			enddo 
!$omp end simd
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! associate pointers to targets                                              !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			ut => u_store2
			wt => w_store2
			ut_sav => u_store1
			wt_sav => w_store1
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 			call exchange_full(comm3d, id, kp, jp, ip, r_h,r_h,r_h,r_h,l_h,r_h, &
! 														ut,dims,coords)
! 			call exchange_full(comm3d, id, kp, jp, ip, l_h,r_h,r_h,r_h,r_h,r_h, &
! 														wt,dims,coords)
            ! halos
            ut(:,-l_h+1:0)      =ut(:,ip-l_h+1:ip)
            ut(:,ip+1:ip+r_h)   =ut(:,1:r_h)
            ut(0,:)=0._sp
            ut(kp+1,:)=0._sp
            
            wt(:,-l_h+1:0)      =wt(:,ip-l_h+1:ip)
            wt(:,ip+1:ip+r_h)   =wt(:,1:r_h)
            wt(0,:)=0._sp
            wt(kp:kp+1,:)=0._sp

		endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! advect using first order upwind                                                !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		call first_order_upstream_2d(dt,dxn,dzn,rhoa,rhoan, &
				ip,kp,l_h,r_h,ut,wt,psi_old)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


		
		if((it <= kord) .and. (kord >= 1)) then
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! set halos																	 !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
! 			call exchange_full(comm3d, id, kp, jp, ip, r_h,r_h,r_h,r_h,r_h,r_h, &
! 									psi_old,dims,coords)
            ! halos
            psi_old(:,-l_h+1:0)         =psi_old(:,ip-l_h+1:ip)
            psi_old(:,ip+1:ip+r_h)      =psi_old(:,1:r_h)
            psi_old(0,:)=0._sp
            psi_old(kp+1,:)=0._sp
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!			
		endif		
 	enddo
 

	! no dangling pointers
	if (associated(ut) ) nullify(ut)
	if (associated(wt) ) nullify(wt)
	if (associated(ut_sav) ) nullify(ut_sav)
	if (associated(wt_sav) ) nullify(wt_sav)
	if (associated(psi) ) nullify(psi)
	if (associated(psi_old) ) nullify(psi_old)

	psi_in=psi_in+minglobal
	end subroutine mpdata_2d
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
	!>solves the 2-d advection equation:
	!>\f$ \frac{\partial \psi}{\partial t} + \nabla \cdot \bar{v} \psi = 0 \f$
	!>@param[in] dt
	!>@param[in] dx,dz, dxn, dyn, dzn, rhoa, rhoan
	!>@param[in] ip,kp,nq,l_h,r_h
	!>@param[in] u
	!>@param[in] w
	!>@param[inout] psi
	!>@param[in] kord, monotone: order of MPDATA and whether it is monotone
	!>@param[in] comm3d, id, dims, coords: mpi variables
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef MPI
	subroutine mpdata_vec_2d(dt,dx,dz,dxn,dzn,&
						rhoa,rhoan, &
						ip,kp,nq,l_h,r_h,u,w,psi_in,kord,monotone, comm3d, id, &
						dims,coords)
#else
	subroutine mpdata_vec_2d(dt,dx,dz,dxn,dzn,&
						rhoa,rhoan, &
						ip,kp,nq,l_h,r_h,u,w,psi_in,kord,monotone)
#endif
	use nrtype
#ifdef MPI
	use mpi_module
	use mpi
#endif
	
	implicit none
	
#ifdef MPI	
	integer(i4b), intent(in) :: id, comm3d
	integer(i4b), dimension(2), intent(in) :: dims, coords
#endif
	real(sp), intent(in) :: dt
	integer(i4b), intent(in) :: ip, kp, nq,l_h, r_h,kord
	real(sp), dimension(-r_h+1:kp+r_h,-l_h+1:ip+r_h), &
		intent(in), target :: u
	real(sp), dimension(-l_h+1:kp+r_h,-r_h+1:ip+r_h), &
		intent(in), target :: w
	real(sp), dimension(-r_h+1:kp+r_h,-r_h+1:ip+r_h,1:nq), &
		intent(inout), target :: psi_in
	real(sp), dimension(-l_h+1:ip+r_h), intent(in) :: dx, dxn
	real(sp), dimension(-l_h+1:kp+r_h), intent(in) :: dz, dzn, rhoa,rhoan
	logical :: monotone
	
	! locals
	real(sp) :: u_div1, u_div3, u_j_bar1, u_j_bar3, &
			denom1, denom2, minlocal, psi_local_sum, psi_sum
	real(sp), dimension(nq) :: minglobal
	integer(i4b) :: i,k, it, it2, n,error
	real(sp), dimension(:,:), pointer :: ut
	real(sp), dimension(:,:), pointer :: wt
	real(sp), dimension(:,:), pointer :: ut_sav
	real(sp), dimension(:,:), pointer :: wt_sav
	real(sp), dimension(:,:), pointer :: psi, psi_old
	real(sp), dimension(-r_h+1:kp+r_h,-l_h+1:ip+r_h), target :: & 
		u_store1, u_store2
	real(sp), dimension(-r_h+1:kp+r_h,-r_h+1:ip+r_h), target :: &
		v_store1, v_store2
	real(sp), dimension(-l_h+1:kp+r_h,-r_h+1:ip+r_h), target :: &
		w_store1, w_store2
	real(sp), dimension(-r_h+1:kp+r_h,-r_h+1:ip+r_h,1:nq), target :: psi_store
	real(sp), dimension(-r_h+1:kp+r_h,-r_h+1:ip+r_h) :: &
						psi_i_max, psi_i_min, psi_k_max,psi_k_min, &
						beta_i_up, beta_i_down,&
						beta_k_up, beta_k_down
	

	! has to be positive definite
	do n=1,nq
        minlocal=minval(psi_in(1:kp,1:ip,n))
#ifdef MPI
        call mpi_allreduce(minlocal,minglobal(n),1,MPI_REAL8,MPI_MIN, comm3d,error)
#else
        minglobal(n)=minlocal
#endif
        psi_in(:,:,n)=psi_in(:,:,n)-minglobal(n)
	enddo
	
	psi_local_sum=sum(psi_in(1:kp,1:ip,1))
#ifdef MPI
	call MPI_Allreduce(psi_local_sum, psi_sum, 1, MPI_REAL8, MPI_SUM, comm3d, error)
#else
    psi_sum=psi_local_sum
#endif
 	if(psi_sum.lt.small) then
        do n=1,nq
            psi_in(:,:,n)=psi_in(:,:,n)+minglobal(n)
        enddo
 	    return
	endif
	u_store2=0._sp
	w_store2=0._sp
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! associate pointers to targets                                                      !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	ut => u	
	wt => w	
	ut_sav => u_store1
	wt_sav => w_store1
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! if kord > 1 we need a copy of the scalar field                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if (kord > 1) then
		psi_store = psi_in   ! array copy
		psi(-r_h+1:kp+r_h,-r_h+1:ip+r_h) => psi_store(:,:,1)
	endif
	psi_old(-r_h+1:kp+r_h,-r_h+1:ip+r_h)     => psi_in(:,:,1)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	it2=0
	do it=1,kord
		
		if(it > 1) then
			it2=it2+1
!$omp simd	
			do i=1,ip
                do k=1,kp
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! calculate the anti-diffusive velocities                        !
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! for divergent flow: eq 38 smolarkiewicz 1984 
                    ! last part of u wind:
                    u_div1=(wt(k,i)+wt(k,i+1)-wt(k-1,i)-wt(k-1,i+1)) &
                            / dz(k-1) 
                    ! for divergent flow: eq 38 smolarkiewicz 1984 
                    ! last part of w wind:
                    u_div3=(ut(k,i)+ut(k+1,i)-ut(k,i-1)-ut(k+1,i-1)) &
                            / dx(i-1) 
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! equation 13 page 330 of smolarkiewicz (1984) 
                    ! journal of computational physics
                    ! second term of u wind:
                    u_j_bar1 = 0.5_sp*dt*ut(k,i) * ( &
                        ! equation 14:
                        0.25_sp*(wt(k,i+1)+wt(k,i)+wt(k-1,i+1)+wt(k-1,i)) * &
                        ! equation 13:
                       ( psi_old(k+1,i+1)+psi_old(k+1,i)- &
                         psi_old(k-1,i+1)-psi_old(k-1,i) ) / &
                       ( psi_old(k+1,i+1)+psi_old(k+1,i)+ &
                         psi_old(k-1,i+1)+psi_old(k-1,i)+small ) / &
                         ( 0.5_sp*(dzn(k-1)+dzn(k)) ) )
                         
                                                  
                    ! equation 13 page 330 of smolarkiewicz (1984) 
                    ! journal of computational physics
                    ! second term of w wind:
                    u_j_bar3 = 0.5_sp*dt*wt(k,i) * ( &
                        ! repeat for y dimension:
                        ! equation 14:
                        0.25_sp*(ut(k+1,i)+ut(k,i)+ut(k+1,i-1)+ut(k,i-1)) * &
                        ! equation 13:
                       ( psi_old(k+1,i+1)+psi_old(k,i+1)- &
                         psi_old(k+1,i-1)-psi_old(k,i-1) ) / &
                       ( psi_old(k+1,i+1)+psi_old(k,i+1)+ &
                         psi_old(k+1,i-1)+psi_old(k,i-1)+small ) / &
                         ( 0.5_sp*(dxn(i-1)+dxn(i)) ) )
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! last update of equation 13
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! u wind:
                    ut_sav(k,i)=(abs(ut(k,i))*dx(i)-dt*ut(k,i)*ut(k,i) ) * &
                        (psi_old(k,i+1)-psi_old(k,i) ) / &
                        (psi_old(k,i+1)+psi_old(k,i)+small) /dxn(i) - u_j_bar1
                    ! w wind:
                    wt_sav(k,i)=(abs(wt(k,i))*dz(k)-dt*wt(k,i)*wt(k,i) ) * &
                        (psi_old(k+1,i)-psi_old(k,i) ) / &
                        (psi_old(k+1,i)+psi_old(k,i)+small) /dzn(k) - u_j_bar3
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! last update of eq 38 smolarkiewicz 1984
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ut_sav(k,i)=ut_sav(k,i) - 0.25_sp*dt*ut(k,i) * &
                     ( (ut(k,i+1)-ut(k,i-1))/(dx(i-1))+u_div1 )
                    wt_sav(k,i)=wt_sav(k,i) - 0.25_sp*dt*wt(k,i) * &
                     ( (wt(k+1,i)-wt(k-1,i))/(dz(k-1))+u_div3 )
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


                enddo
			enddo
!$omp end simd

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! associate pointers to targets                                              !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			if(modulo(it2,2).eq.1) then  
										 
				ut => u_store1
				wt => w_store1
				ut_sav => u_store2
				wt_sav => w_store2
			else if(modulo(it2,2).eq.0) then			
				ut => u_store2
				wt => w_store2
				ut_sav => u_store1
				wt_sav => w_store1
			endif
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 			call exchange_full(comm3d, id, kp, jp, ip, r_h,r_h,r_h,r_h,l_h,r_h, &
! 														ut,dims,coords)
! 			call exchange_full(comm3d, id, kp, jp, ip, l_h,r_h,r_h,r_h,r_h,r_h, &
! 														wt,dims,coords)
            ! halos
            ut(:,-l_h+1:0)      =ut(:,ip-l_h+1:ip)
            ut(:,ip+1:ip+r_h)   =ut(:,1:r_h)
            ut(0,:)=0._sp
            ut(kp+1,:)=0._sp

            wt(:,-l_h+1:0)      =wt(:,ip-l_h+1:ip)
            wt(:,ip+1:ip+r_h)   =wt(:,1:r_h)
            wt(0,:)=0._sp
            wt(kp:kp+1,:)=0._sp
			
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
			do i=1,ip
                do k=1,kp			
                    ! x direction - note: should the last q in the max/min be q+1?
                    psi_i_max(k,i)=max(psi(k,i-1),psi(k,i),psi(k,i+1), &
                                psi_old(k,i-1),psi_old(k,i),psi_old(k,i+1))
                
                    psi_i_min(k,i)=min(psi(k,i-1),psi(k,i),psi(k,i+1), &
                                psi_old(k,i-1),psi_old(k,i),psi_old(k,i+1))
                enddo
			enddo
!$omp end simd


!$omp simd	
			do i=1,ip
                do k=1,kp			
                    ! z direction - note: should the last q in the max/min be q+1?
                    psi_k_max(k,i)=max(psi(k-1,i),psi(k,i),psi(k+1,i), &
                                psi_old(k-1,i),psi_old(k,i),psi_old(k+1,i))
                
                    psi_k_min(k,i)=min(psi(k-1,i),psi(k,i),psi(k+1,i), &
                                psi_old(k-1,i),psi_old(k,i),psi_old(k+1,i))
                enddo
			enddo
!$omp end simd
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! Equations 19a and b of Smolarkiewicz and Grabowski (1990, JCP)             !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$omp simd	
			do i=1,ip
                do k=1,kp		
                    denom1=(dt*((max(ut(k,i-1),0._sp)*psi_old(k,i-1)- &
                              min(ut(k,i),0._sp)*psi_old(k,i+1))/dx(i-1)+ &
                            (max(wt(k-1,i),0._sp)*psi_old(k-1,i)-&
                              min(wt(k,i),0._sp)*psi_old(k+1,i))/dz(k-1) &
                              +small))
                              
                    denom2=(dt*((max(ut(k,i),0._sp)*psi_old(k,i)- &
                              min(ut(k,i-1),0._sp)*psi_old(k,i))/dx(i-1) + &
                            (max(wt(k,i),0._sp)*psi_old(k,i)-&
                              min(wt(k-1,i),0._sp)*psi_old(k,i))/dz(k-1) &
                              +small))
                              
                    beta_i_up(k,i)=(psi_i_max(k,i)-psi_old(k,i)) / denom1
                        
                              
                    beta_i_down(k,i)=(psi_old(k,i)-psi_i_min(k,i)) / denom2
                                        

                    beta_k_up(k,i)=(psi_k_max(k,i)-psi_old(k,i)) / denom1
                              
                    beta_k_down(k,i)=(psi_old(k,i)-psi_k_min(k,i)) / denom2
                enddo
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
            beta_i_up(:,-l_h+1:0)       =beta_i_up(:,ip-l_h+1:ip)
            beta_i_up(:,ip+1:ip+r_h)    =beta_i_up(:,1:r_h)
            beta_i_up(0,:)=0._sp
            beta_i_up(kp+1,:)=0._sp
            
            beta_i_down(:,-l_h+1:0)     =beta_i_down(:,ip-l_h+1:ip)
            beta_i_down(:,ip+1:ip+r_h)  =beta_i_down(:,1:r_h)
            beta_i_down(0,:)=0._sp
            beta_i_down(kp+1,:)=0._sp
            
            beta_k_up(:,-l_h+1:0)       =beta_k_up(:,ip-l_h+1:ip)
            beta_k_up(:,ip+1:ip+r_h)    =beta_k_up(:,1:r_h)
            beta_k_up(0,:)=0._sp
            beta_k_up(kp+1,:)=0._sp
            
            beta_k_down(:,-l_h+1:0)     =beta_k_down(:,ip-l_h+1:ip)
            beta_k_down(:,ip+1:ip+r_h)  =beta_k_down(:,1:r_h)
            beta_k_down(0,:)=0._sp
            beta_k_down(kp+1,:)=0._sp
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


								
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! equation 18 of Smolarkiewicz and Grabowski (1990, JCP, 86)                 !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$omp simd	
			do i=1,ip
                do k=1,kp			
                    ut_sav(k,i)=min(1._sp,beta_i_down(k,i), &
                                     beta_i_up(k,i+1))*max(ut(k,i),0._sp) + &
                                  min(1._sp,beta_i_up(k,i), &
                                     beta_i_down(k,i+1))*min(ut(k,i),0._sp)
                    wt_sav(k,i)=min(1._sp,beta_k_down(k,i), &
                                     beta_k_up(k+1,i))*max(wt(k,i),0._sp) + &
                                  min(1._sp,beta_k_up(k,i), &
                                     beta_k_down(k+1,i))*min(wt(k,i),0._sp)
                enddo
			enddo 
!$omp end simd
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! associate pointers to targets                                              !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			ut => u_store2
			wt => w_store2
			ut_sav => u_store1
			wt_sav => w_store1
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 			call exchange_full(comm3d, id, kp, jp, ip, r_h,r_h,r_h,r_h,l_h,r_h, &
! 														ut,dims,coords)
! 			call exchange_full(comm3d, id, kp, jp, ip, l_h,r_h,r_h,r_h,r_h,r_h, &
! 														wt,dims,coords)
            ! halos
            ut(:,-l_h+1:0)      =ut(:,ip-l_h+1:ip)
            ut(:,ip+1:ip+r_h)   =ut(:,1:r_h)
            ut(0,:)=0._sp
            ut(kp+1,:)=0._sp

            wt(:,-l_h+1:0)      =wt(:,ip-l_h+1:ip)
            wt(:,ip+1:ip+r_h)   =wt(:,1:r_h)
            wt(0,:)=0._sp
            wt(kp:kp+1,:)=0._sp


		endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! advect using first order upwind                                                !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		do n=1,nq
            call first_order_upstream_2d(dt,dxn,dzn,rhoa,rhoan, &
                    ip,kp,l_h,r_h,ut,wt,psi_in(:,:,n))
            if((it <= kord) .and. (kord >= 1)) then
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! set halos														    	 !
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
    ! 			call exchange_full(comm3d, id, kp, jp, ip, r_h,r_h,r_h,r_h,r_h,r_h, &
    ! 									psi_old,dims,coords)
                ! halos
                psi_in(:,-l_h+1:0,n)         =psi_in(:,ip-l_h+1:ip,n)
                psi_in(:,ip+1:ip+r_h,n)      =psi_in(:,1:r_h,n)
                psi_in(0,:,n)=0._sp
                psi_in(kp+1,:,n)=0._sp
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!			
            endif		
        enddo
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


		
! 		if((it <= kord) .and. (kord >= 1)) then
! 			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 			! set halos																	 !
! 			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
! ! 			call exchange_full(comm3d, id, kp, jp, ip, r_h,r_h,r_h,r_h,r_h,r_h, &
! ! 									psi_old,dims,coords)
!             ! halos
!             psi_old(:,-l_h+1:0)         =psi_old(:,ip-l_h+1:ip)
!             psi_old(:,ip+1:ip+r_h)      =psi_old(:,1:r_h)
!             psi_old(0,:)=0._sp
!             psi_old(kp,:)=0._sp
! 			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!			
! 		endif		
 	enddo
 

	! no dangling pointers
	if (associated(ut) ) nullify(ut)
	if (associated(wt) ) nullify(wt)
	if (associated(ut_sav) ) nullify(ut_sav)
	if (associated(wt_sav) ) nullify(wt_sav)
	if (associated(psi) ) nullify(psi)
	if (associated(psi_old) ) nullify(psi_old)

    do n=1,nq
    	psi_in(:,:,n)=psi_in(:,:,n)+minglobal(n)
    enddo
	end subroutine mpdata_vec_2d
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    end module advection_s_2d
