	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>diagnostic routines for the dynamical cloud model
    module diagnostics
    use nrtype
    implicit none

    private
    public :: horizontal_means, divergence_calc
    
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
	subroutine horizontal_means(comm3d,id,dims,coords, moisture, &
	    ip,jp,ipp,jpp,kpp,l_h,r_h, nq,&
	    ubar,vbar,wbar,thbar,qbar,u,v,w,th,q)
		use nrtype
		use mpi
		use mpi_module
		implicit none
		integer(i4b), intent(in) :: id, comm3d
		integer(i4b), dimension(3), intent(in) :: dims,coords
		logical, intent(in) :: moisture
		integer(i4b), intent(in) :: ip,jp,ipp,jpp,kpp, l_h,r_h, nq
		real(sp), dimension(1-l_h:kpp+r_h), intent(inout) :: ubar,vbar,wbar,thbar
		real(sp), dimension(1-l_h:kpp+r_h,1:nq), intent(inout) :: qbar
		real(sp), dimension(-l_h+1:kpp+r_h,-l_h+1:jpp+r_h,-l_h+1:ipp+r_h), &
			intent(in) :: u,v,w,th
		real(sp), dimension(-l_h+1:kpp+r_h,-l_h+1:jpp+r_h,-l_h+1:ipp+r_h,1:nq), &
			intent(in) :: q
				
		! local
		integer(i4b) :: i,j,k,n,error
		real(sp),dimension((4+nq)*(kpp+l_h+r_h+1)) :: local_sum,global_sum

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! ubar,vbar,wbar,thbar write to sum array                                        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do k=1-l_h,kpp+r_h
            local_sum(1+k)                =sum(u(k,1:jpp,1:ipp))
            local_sum(kpp+l_h+r_h+1+k)    =sum(v(k,1:jpp,1:ipp))
            local_sum((kpp+l_h+r_h+1)*2+k)=sum(w(k,1:jpp,1:ipp))
            local_sum((kpp+l_h+r_h+1)*3+k)=sum(th(k,1:jpp,1:ipp))
        enddo
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! same for moisture                                                              !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(moisture) then
            do n=1,nq
                do k=1-l_h,kpp+r_h
                     local_sum((kpp+l_h+r_h+1)*(3+n)+k)=sum(q(k,1:jpp,1:ipp,n))
                enddo        
            enddo
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! one call to allreduce                                                          !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(moisture) then
            call MPI_Allreduce(local_sum, global_sum, (4+nq)*(kpp+l_h+r_h), &
                MPI_REAL8, MPI_SUM, comm3d, error)
        else
            call MPI_Allreduce(local_sum, global_sum, 4*(kpp+l_h+r_h), &
                MPI_REAL8, MPI_SUM, comm3d, error)        
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! re-write back                                                                  !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do k=1-l_h,kpp+r_h
            ubar(k)=global_sum(1+k) / (real(ip*jp,sp))
            vbar(k)=global_sum(kpp+l_h+r_h+1+k) / (real(ip*jp,sp))
            wbar(k)=global_sum((kpp+l_h+r_h+1)*2+k) / (real(ip*jp,sp))
            thbar(k)=global_sum((kpp+l_h+r_h+1)*3+k) / (real(ip*jp,sp))
        enddo
        
        if(moisture) then
            do n=1,nq
                do k=1-l_h,kpp+r_h
                     qbar(k,n)=global_sum((kpp+l_h+r_h+1)*(3+n)+k) / (real(ip*jp,sp))
                enddo        
            enddo        
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
	end subroutine horizontal_means	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>divergence calculation
	!>@param[in] comm3d,id, dims, coords
	!>@param[in] dx,dy,dz,dxn,dyn,dzn
	!>@param[in] ipp,jpp,kpp,halo
	!>@param[inout] div
	subroutine divergence_calc(comm3d,id,dims,coords, &
	    ipp,jpp,kpp,dx,dxn,dy,dyn,dz,dzn,l_h,r_h,u,v,w,div)
		use nrtype
		use mpi
		use mpi_module
		implicit none
		integer(i4b), intent(in) :: id, comm3d
		integer(i4b), dimension(3), intent(in) :: dims,coords
		integer(i4b), intent(in) :: ipp,jpp,kpp, l_h,r_h
		real(sp), dimension(-l_h+1:ipp+r_h), intent(in) :: dx,dxn
		real(sp), dimension(-l_h+1:jpp+r_h), intent(in) :: dy,dyn
		real(sp), dimension(-l_h+1:kpp+r_h), intent(in) :: dz,dzn
		real(sp), dimension(-l_h+1:kpp+r_h,-l_h+1:jpp+r_h,-l_h+1:ipp+r_h), &
			intent(in) :: u,v,w
		real(sp), dimension(-l_h+1:kpp+r_h,-l_h+1:jpp+r_h,-l_h+1:ipp+r_h), &
			intent(inout) :: div
				
		! local
		integer(i4b) :: i,j,k,error


		do i=1,ipp
			do j=1,jpp	    
				do k=1,kpp
					div(k,j,i)=(u(k,j,i)-u(k,j,i-1))/dx(i-1) + &
					    (v(k,j,i)-v(k,j-1,i))/dy(j-1) + &
					    (w(k,j,i)-w(k-1,j,i))/dz(k-1) 
				enddo
			enddo
		enddo

        
        
	end subroutine divergence_calc	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    end module diagnostics
