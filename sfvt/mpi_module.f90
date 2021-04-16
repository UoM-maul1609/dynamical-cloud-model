	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>mpi routines for dynamical cloud model
	!> note: should be able to get significant speed up by using irecv, isend and waits
	!> need to pre-allocate correct buffer sizes for this though
	!> could probably also do away with the if statements - call isend and irecv anyway 
	!> - even if an mpi call is not required
    module mpi_module
    use nrtype
    use mpi
    implicit none
		!>@brief
		!>variables for mpi
        type mpi_faces
        	integer(i4b) :: s_west, r_west
        	integer(i4b) :: s_east, r_east
        	integer(i4b) :: s_south, r_south
        	integer(i4b) :: s_north, r_north
        	integer(i4b) :: s_bottom, r_bottom
        	integer(i4b) :: s_top, r_top
        end type mpi_faces
        
        type mpi_faces_v
        	integer(i4b) :: s_bottom, r_bottom
        	integer(i4b) :: s_top, r_top
        	integer(i4b) :: id
        end type mpi_faces_v
        
        type mpi_edges
        	integer(i4b) :: s_ws_bt, s_wn_bt, s_es_bt, s_en_bt, s_bs_we, &
        					s_bn_we, s_ts_we, s_tn_we, s_bw_sn, s_be_sn, s_tw_sn, s_te_sn
        	integer(i4b) :: r_ws_bt, r_wn_bt, r_es_bt, r_en_bt, r_bs_we, &
        					r_bn_we, r_ts_we, r_tn_we, r_bw_sn, r_be_sn, r_tw_sn, r_te_sn
        end type mpi_edges

        type mpi_corners
        	integer(i4b) :: s_wsb, s_esb, s_wnb, s_enb, s_wst, s_est, s_wnt, s_ent
        	integer(i4b) :: r_wsb, r_esb, r_wnb, r_enb, r_wst, r_est, r_wnt, r_ent
        end type mpi_corners

        type mpi_vars
        	integer(i4b) :: rank, id, error, dx, dy, dz
        	real(sp) :: wtime
        	logical, dimension(3) :: periods
        	logical :: reorder=.true.
        	integer(i4b) :: ndim=3
        	integer(i4b), dimension(3) :: dims, coords
        	logical, dimension(3) :: remain_dims=[.false.,.false.,.true.]
        	logical, dimension(3) :: remain_dims_horiz=[.true.,.true.,.false.]
        	integer(i4b) :: colour,height
        	integer(i4b) :: ring_comm
        	integer(i4b) :: sub_comm
        	integer(i4b) :: sub_horiz_comm
        	type(mpi_faces) :: face
        	type(mpi_faces_v) :: face_v
        	type(mpi_edges) :: edge
        	type(mpi_corners) :: cnr
        end type mpi_vars

		! declare an mpi type
		type(mpi_vars) :: mp1

		integer(i4b) :: world_process=0  
		integer(i4b) :: MPI_INTEGER9   

    private
    public :: mpi_define, block_ring, exchange_full, mpi_integer9, mp1, world_process, &
    			mpi_cart_initialise, exchange_along_dim, find_base_top, find_top, &
    			exchange_d_fluxes, exchange_u_fluxes, exchange_along_z, &
    			exchange_full_wo, exchange_along_dim_wo
    

	contains
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! define some types                                                                  !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>define some types to be used in the model
	!>@param[inout] MPI_INTEGER9_DEF - type to be defined as integer kind=9
	subroutine mpi_define(MPI_INTEGER9_DEF)
		implicit none
		integer(i4b), intent(inout) :: MPI_INTEGER9_DEF
		
		integer(i4b) :: error
		
		call MPI_TYPE_CREATE_F90_INTEGER (9, MPI_INTEGER9_DEF, error)
	end subroutine mpi_define
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! initialise some variables for swapping faces, edges, and corners                   !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>define some types to be used in the model
	!>@param[in] ip, jp, kp
	subroutine mpi_cart_initialise(kp,jp,ip)
		use nrtype
		implicit none
		integer(i4b), intent(in) :: ip, jp, kp
		
		integer(i4b) :: error, dummy
		integer(i4b), dimension(3) :: coords_t
		
		! define an integer:
		call mpi_define(MPI_INTEGER9)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Set-up the Cartesian topology									             	 !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! note the min is so that there is not more than 1 proc per grid point
		mp1%dx=min( floor( (real(mp1%rank,sp))**(1._sp/3._sp) ) , ip)
		mp1%dy=min( floor( &
					(real(mp1%rank,sp)/real(mp1%dx,sp))**(1._sp/2._sp) ) , jp)
		mp1%dz=min( floor( real(mp1%rank,sp) / real(mp1%dx*mp1%dy,sp) ), kp )

		if(mp1%id == world_process) then
			print *,'Cartesian topology: ',mp1%dx, mp1%dy, mp1%dz		
			if ( mp1%dx * mp1%dy * mp1%dz < mp1%rank) &
				print *, 'warning wasted processors'
		endif

		mp1%periods=[.true.,.true.,.false.]
		mp1%dims=[mp1%dx,mp1%dy,mp1%dz]
		! cart topo:
		call MPI_CART_CREATE( MPI_COMM_WORLD, mp1%ndim, mp1%dims, &
							mp1%periods, mp1%reorder,mp1%ring_comm, mp1%error )
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Find ids for sending and receiving 6 faces                                     !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if(mp1%id >= mp1%dx * mp1%dy * mp1%dz) then
			return
		endif
		call MPI_CART_SHIFT( mp1%ring_comm, 0, 1, &
							mp1%face%s_west, mp1%face%s_east, error)
		call MPI_CART_SHIFT( mp1%ring_comm, 1, 1, &
							mp1%face%s_south, mp1%face%s_north, error)
		call MPI_CART_SHIFT( mp1%ring_comm, 2, 1, &
							mp1%face%s_bottom, mp1%face%s_top, error)
							
		mp1%face%r_east  = mp1%face%s_west
		mp1%face%r_west  = mp1%face%s_east
		mp1%face%r_north = mp1%face%s_south
		mp1%face%r_south = mp1%face%s_north
		mp1%face%r_top   = mp1%face%s_bottom
		mp1%face%r_bottom= mp1%face%s_top
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Find ids for sending and receiving 12 edges                                    !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! coordinates of this PE:
		call MPI_CART_COORDS(mp1%ring_comm, mp1%id, 3, mp1%coords, error)
		mp1%colour=(mp1%coords(1)+1)*(mp1%dims(1))+(mp1%coords(2)+1)
		mp1%height=mp1%coords(3)
		
		call mpi_find_rank(mp1%coords,mp1%dims,-1,-1,0,mp1%edge%s_ws_bt,mp1%edge%r_en_bt)
		call mpi_find_rank(mp1%coords,mp1%dims,-1,+1,0, mp1%edge%s_wn_bt,mp1%edge%r_es_bt)
		call mpi_find_rank(mp1%coords,mp1%dims,+1,-1,0, mp1%edge%s_es_bt,mp1%edge%r_wn_bt)
		call mpi_find_rank(mp1%coords,mp1%dims,+1,+1,0, mp1%edge%s_en_bt,mp1%edge%r_ws_bt)
								
		call mpi_find_rank(mp1%coords,mp1%dims,0,-1,-1, mp1%edge%s_bs_we,mp1%edge%r_tn_we)
		call mpi_find_rank(mp1%coords,mp1%dims,0,+1,-1, mp1%edge%s_bn_we,mp1%edge%r_ts_we)
		call mpi_find_rank(mp1%coords,mp1%dims,0,-1,+1, mp1%edge%s_ts_we,mp1%edge%r_bn_we)
		call mpi_find_rank(mp1%coords,mp1%dims,0,+1,+1, mp1%edge%s_tn_we,mp1%edge%r_bs_we)
		
		call mpi_find_rank(mp1%coords,mp1%dims,-1,0,-1,mp1%edge%s_bw_sn,mp1%edge%r_te_sn)
		call mpi_find_rank(mp1%coords,mp1%dims,+1,0,-1,mp1%edge%s_be_sn,mp1%edge%r_tw_sn)
		call mpi_find_rank(mp1%coords,mp1%dims,-1,0,+1,mp1%edge%s_tw_sn,mp1%edge%r_be_sn)
		call mpi_find_rank(mp1%coords,mp1%dims,+1,0,+1,mp1%edge%s_te_sn,mp1%edge%r_bw_sn)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Find ids for sending and receiving 8 corners                                   !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		call mpi_find_rank(mp1%coords,mp1%dims,-1,-1,-1, mp1%cnr%s_wsb,mp1%cnr%r_ent)
		call mpi_find_rank(mp1%coords,mp1%dims,+1,-1,-1, mp1%cnr%s_esb,mp1%cnr%r_wnt)
		call mpi_find_rank(mp1%coords,mp1%dims,-1,+1,-1, mp1%cnr%s_wnb,mp1%cnr%r_est)
		call mpi_find_rank(mp1%coords,mp1%dims,+1,+1,-1, mp1%cnr%s_enb,mp1%cnr%r_wst)
		
		call mpi_find_rank(mp1%coords,mp1%dims,-1,-1,+1, mp1%cnr%s_wst,mp1%cnr%r_enb)
		call mpi_find_rank(mp1%coords,mp1%dims,+1,-1,+1, mp1%cnr%s_est,mp1%cnr%r_wnb)
		call mpi_find_rank(mp1%coords,mp1%dims,-1,+1,+1, mp1%cnr%s_wnt,mp1%cnr%r_esb)
		call mpi_find_rank(mp1%coords,mp1%dims,+1,+1,+1, mp1%cnr%s_ent,mp1%cnr%r_wsb)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! create sub communicator for vertical comms                                     !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if(mp1%id  < mp1%dx * mp1%dy * mp1%dz) then
            call MPI_Cart_sub(mp1%ring_comm,mp1%remain_dims,&
                mp1%sub_comm,mp1%error)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Find ids for sending and receiving 6 faces                                 !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call MPI_CART_SHIFT( mp1%sub_comm, 0, 0, &
                                mp1%face_v%id, dummy, error)

            call MPI_CART_SHIFT( mp1%sub_comm, 0, 1, &
                                mp1%face_v%s_bottom, mp1%face_v%s_top, error)
            mp1%face_v%id=mp1%id
            mp1%face_v%r_top   = mp1%face_v%s_bottom
            mp1%face_v%r_bottom= mp1%face_v%s_top
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! create sub communicator for horizontal comms                                   !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if(mp1%id  < mp1%dx * mp1%dy * mp1%dz) then
            call MPI_Cart_sub(mp1%ring_comm,mp1%remain_dims_horiz,&
                mp1%sub_horiz_comm,mp1%error)
        endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		

		contains
			subroutine mpi_find_rank(coords,dims,we,sn,du,send,recv)
				implicit none
				integer(i4b), dimension(3), intent(in) :: coords, dims
				integer(i4b), intent(in) :: we,sn,du
				integer(i4b), intent(inout) :: send,recv
				
				! locals:
				integer(i4b), dimension(3) :: coords_t
				integer(i4b) :: i
				
			
				coords_t(1)=mp1%coords(1)+we
				coords_t(2)=mp1%coords(2)+sn
				coords_t(3)=mp1%coords(3)+du
				
				do i=1,3
					if (mp1%periods(i)) then
						if(coords_t(i) > (mp1%dims(i)-1) ) then
							coords_t(i)=0
						else if(coords_t(i) < 0 ) then
							coords_t(i)=mp1%dims(i)-1
						endif
					else
						if(coords_t(i) > (mp1%dims(i)-1) ) then
							coords_t(i)=-1
						else if(coords_t(i) < 0 ) then
							coords_t(i)=-1
						endif
					endif
				enddo
				
				if( any(coords_t .eq. -1) ) then
					recv=-1
					send=-1
				else
					! get the rank of this coordinate
					call MPI_CART_RANK( mp1%ring_comm, coords_t, send,error)
					recv=send
				endif				 
			
			end subroutine mpi_find_rank

		

		
	end subroutine mpi_cart_initialise
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! find base and top of array using bcast                                             !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>define some types to be used in the model
	!>@param[in] comm3d, id, kpp, d_h, u_h, array, dims, coords
	!>@param[inout] base, top
    subroutine find_base_top(comm3d, id, kpp,d_h,u_h,array,base,top, dims,coords)
        implicit none
		integer(i4b), intent(in) :: comm3d, id, kpp, d_h,u_h
		real(sp), intent(in), dimension(1-d_h:u_h+kpp) :: array
		real(sp), intent(inout) :: base, top
		integer(i4b), dimension(3), intent(in) :: dims,coords
		
		! locals:
		integer(i4b), dimension(12) :: request
		integer(i4b), dimension(MPI_STATUS_SIZE, 12) :: status
		integer(i4b) :: error, tag1,num_messages,imess, tag2, root, rankl,ranku
		
		
		call MPI_CART_RANK( comm3d, [0,0,0], rankl,error)
		call MPI_CART_RANK( comm3d, [0,0,dims(3)-1], ranku,error)
		
		
        base=array(0)
        call MPI_Bcast(base,1,MPI_REAL8,rankl,comm3d,error)

        top=array(kpp+1)
        call MPI_Bcast(top,1,MPI_REAL8,ranku,comm3d,error)
    
    end subroutine find_base_top
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! find base and top of array using bcast                                             !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>define some types to be used in the model
	!>@param[in] comm3d, id, kpp, d_h, u_h, array, dims, coords
	!>@param[inout] top
    subroutine find_top(comm3d, id, kpp,d_h,u_h,array,top, dims,coords)
        implicit none
		integer(i4b), intent(in) :: comm3d, id, kpp, d_h,u_h
		real(sp), intent(in), dimension(1-d_h:u_h+kpp) :: array
		real(sp), intent(inout) :: top
		integer(i4b), dimension(3), intent(in) :: dims,coords
		
		! locals:
		integer(i4b), dimension(12) :: request
		integer(i4b), dimension(MPI_STATUS_SIZE, 12) :: status
		integer(i4b) :: error, tag1,num_messages,imess, tag2, root, rankl,ranku
		
		
		call MPI_CART_RANK( comm3d, [0,0,dims(3)-1], ranku,error)
		
		

        top=array(kpp)
        call MPI_Bcast(top,1,MPI_REAL8,ranku,comm3d,error)
    
    end subroutine find_top
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    



	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! exchange along dim for a variable using Cartesian topology                         !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>define some types to be used in the model
	!>@param[in] comm3d, id, ipp, jpp, kpp, nbands,w_h,e_h,s_h,n_h,d_h,u_h
	!>@param[inout] array: the array to exchange_halos on
	!>@param[in] lbc, ubc, dims,coords
	subroutine exchange_d_fluxes(comm3d, id, kpp, jpp, ipp, nbands,&
							d_h,u_h,s_h,n_h,w_h, e_h,  array, dims,coords)
		implicit none
		
		integer(i4b), intent(in) :: comm3d, id, ipp, jpp, kpp, nbands,&
		     w_h, e_h, s_h,n_h,d_h,u_h
		real(sp), intent(inout), &
			 dimension(1-d_h:u_h+kpp,1-s_h:n_h+jpp,1-w_h:e_h+ipp,nbands) :: &
			 array
		integer(i4b), dimension(3), intent(in) :: dims,coords
		
		! locals:
		integer(i4b), dimension(12) :: request
		integer(i4b), dimension(MPI_STATUS_SIZE, 12) :: status
		integer(i4b) :: error, tag1,num_messages,imess, tag2
		
		

			
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! message passing for adjacent cells in up / down direction                      !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		tag1=11
		! send to the bottom:
		call MPI_Issend(array(0:1,1:jpp,1:ipp,1:nbands), &
			(ipp*jpp*nbands)*(u_h+1), MPI_REAL8, mp1%face%s_bottom, &
			tag1, comm3d, request(1),error)

		! receive from the bottom of upper cell:
		call MPI_Recv(array(kpp:kpp+1,1:jpp,1:ipp,1:nbands), &
			(ipp*jpp*nbands)*(u_h+1), MPI_REAL8, mp1%face%r_bottom, &
			tag1, comm3d, status(:,1),error)
! 		call MPI_Wait(request(1), status(:,1), error)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! message passing for adjacent cells in up / down direction                      !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		tag1=12
		! send to the top:
		call MPI_Issend(array(kpp,1:jpp,1:ipp,1:nbands), &
			(ipp*jpp*nbands)*(u_h), MPI_REAL8, mp1%face%s_top, &
			tag1, comm3d, request(2),error)

		! receive from the bottom of upper cell:
		call MPI_Recv(array(0,1:jpp,1:ipp,1:nbands), &
			(ipp*jpp*nbands)*(u_h), MPI_REAL8, mp1%face%r_top, &
			tag1, comm3d, status(:,2),error)
		call MPI_Wait(request(1), status(:,1), error)
		call MPI_Wait(request(2), status(:,2), error)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	end subroutine exchange_d_fluxes
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! exchange along dim for a variable using Cartesian topology                         !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>define some types to be used in the model
	!>@param[in] comm3d, id, ipp, jpp, kpp, nbands,w_h,e_h,s_h,n_h,d_h,u_h
	!>@param[inout] array: the array to exchange_halos on
	!>@param[in] lbc, ubc, dims,coords
	subroutine exchange_u_fluxes(comm3d, id, kpp, jpp, ipp, nbands,&
							d_h,u_h,s_h,n_h,w_h, e_h,  array, dims,coords)
		implicit none
		
		integer(i4b), intent(in) :: comm3d, id, ipp, jpp, kpp, nbands,&
		     w_h, e_h, s_h,n_h,d_h,u_h
		real(sp), intent(inout), &
			 dimension(1-d_h:u_h+kpp,1-s_h:n_h+jpp,1-w_h:e_h+ipp,nbands) :: &
			 array
		integer(i4b), dimension(3), intent(in) :: dims,coords
		
		! locals:
		integer(i4b), dimension(12) :: request
		integer(i4b), dimension(MPI_STATUS_SIZE, 12) :: status
		integer(i4b) :: error, tag1,num_messages,imess, tag2
		
		

			
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! message passing for adjacent cells in up / down direction                      !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		tag1=11
		! send to the bottom:
		call MPI_Issend(array(1,1:jpp,1:ipp,1:nbands), &
			(ipp*jpp*nbands)*(u_h), MPI_REAL8, mp1%face%s_bottom, &
			tag1, comm3d, request(1),error)

		! receive from the bottom of upper cell:
		call MPI_Recv(array(kpp+1,1:jpp,1:ipp,1:nbands), &
			(ipp*jpp*nbands)*(u_h+1), MPI_REAL8, mp1%face%r_bottom, &
			tag1, comm3d, status(:,1),error)
! 		call MPI_Wait(request(1), status(:,1), error)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! message passing for adjacent cells in up / down direction                      !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		tag1=12
		! send to the top:
		call MPI_Issend(array(kpp,1:jpp,1:ipp,1:nbands), &
			(ipp*jpp*nbands)*(u_h), MPI_REAL8, mp1%face%s_top, &
			tag1, comm3d, request(2),error)

		! receive from the bottom of upper cell:
		call MPI_Recv(array(0,1:jpp,1:ipp,1:nbands), &
			(ipp*jpp*nbands)*(u_h), MPI_REAL8, mp1%face%r_top, &
			tag1, comm3d, status(:,2),error)
		call MPI_Wait(request(1), status(:,1), error)
		call MPI_Wait(request(2), status(:,2), error)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	end subroutine exchange_u_fluxes
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	
! 	subroutine exchange_fluxes(comm3d, id, kpp, jpp, ipp, nbands,&
! 							d_h,u_h,s_h,n_h,w_h, e_h,  array, dims,coords)
! 		implicit none
! 		
! 		integer(i4b), intent(in) :: comm3d, id, ipp, jpp, kpp, nbands,&
! 		     w_h, e_h, s_h,n_h,d_h,u_h
! 		real(sp), intent(inout), &
! 			 dimension(1-d_h:u_h+kpp,1-s_h:n_h+jpp,1-w_h:e_h+ipp,nbands) :: &
! 			 array
! 		integer(i4b), dimension(3), intent(in) :: dims,coords
! 		
! 		! locals:
! 		integer(i4b), dimension(2) :: request
! 		integer(i4b), dimension(MPI_STATUS_SIZE, 2) :: status
! 		integer(i4b) :: error, tag1,num_messages,imess, tag2
! 		
! 		! buffers for receiving data
! 	    integer(i4b) :: count
! 		real(sp), dimension(jpp,ipp,nbands) :: top_buf
! 		
! 
! 			
! 		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 		! message passing for adjacent cells in up / down direction                      !
! 		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 		tag1=11
! 		count=0
! 		! send to the bottom:
! 		if(mp1%face%s_bottom /= id) then
!             count=count+1
!             call MPI_Isend(array(0,1:jpp,1:ipp,1:nbands), &
!                 (ipp*jpp*nbands)*u_h, MPI_REAL8, mp1%face%s_bottom, &
!                 tag1, comm3d, request(count),error)
!         endif
!         
! 		! receive from the bottom of upper cell:
! 		if(mp1%face%r_bottom /= id) then
!             count=count+1
!             call MPI_IRecv(top_buf, &
!                 (ipp*jpp*nbands)*u_h, MPI_REAL8, mp1%face%r_bottom, &
!                 tag1, comm3d, request(count),error)
!         endif
! 		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! 
! 		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         ! wait for all messages to return before writing buffers                         !
! 		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         call MPI_Waitall(count,request(1:count), status(:,1:count), error)
! 		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        
!         
! 		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         ! write buffers                                                                  !
! 		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 		if(mp1%face%s_top /= -1) then
!             array(kpp,1:jpp,1:ipp,1:nbands)=top_buf
!         endif
! 		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! 
! 
! 
! 	end subroutine exchange_fluxes
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! exchange along z dim for a variable using Cartesian topology                       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>define some types to be used in the model
	!>@param[in] comm3d, id, ipp, jpp, kpp, w_h,e_h,s_h,n_h,d_h,u_h
	!>@param[inout] array: the array to exchange_halos on
	!>@param[in] lbc, ubc, dims,coords
	subroutine exchange_along_z(comm3d, id, kpp, jpp, ipp, &
							d_h,u_h,s_h,n_h,w_h, e_h,  array, lbc, ubc, dims,coords)
		implicit none
		
		integer(i4b), intent(in) :: comm3d, id, ipp, jpp, kpp, w_h, e_h, s_h,n_h,d_h,u_h
		real(sp), intent(inout), &
			 dimension(1-d_h:u_h+kpp,1-s_h:n_h+jpp,1-w_h:e_h+ipp) :: &
			 array
		real(sp), intent(in) :: lbc, ubc
		integer(i4b), dimension(3), intent(in) :: dims,coords
		
		! locals:
		integer(i4b), dimension(4) :: request
		integer(i4b), dimension(MPI_STATUS_SIZE, 4) :: status
		integer(i4b) :: error, tag1,num_messages,imess, tag2

		! buffers for receiving data
	    integer(i4b) :: count
		real(sp), dimension(jpp,ipp) :: top_buf, bot_buf
		
		
			
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! message passing for adjacent cells in up / down direction                      !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		tag1=10
		count=0
		! receive from the top of the lower cell:
        if(mp1%face%r_top /= mp1%face_v%id) then
		    count=count+1
            call MPI_IRecv(bot_buf, &!array(1-d_h:0,1:jpp,1:ipp), &
                (ipp*jpp)*d_h, MPI_REAL8, mp1%face_v%r_top, &
                tag1, comm3d, request(count),error)
        endif
		tag1=11
		! receive from the bottom of upper cell:
		if(mp1%face%r_bottom /= mp1%face_v%id) then
		    count=count+1
            call MPI_IRecv(top_buf, & !array(kpp+1:kpp+u_h,1:jpp,1:ipp), &
                (ipp*jpp)*u_h, MPI_REAL8, mp1%face_v%r_bottom, &
                tag1, comm3d, request(count),error)
        endif
        tag1=10
		! send to the top:
		if(mp1%face%s_top /= mp1%face_v%id) then
		    count=count+1
            call MPI_Isend(array(kpp+1-d_h:kpp,1:jpp,1:ipp), &
                (ipp*jpp)*d_h, MPI_REAL8, mp1%face_v%s_top, &
                tag1, comm3d, request(count),error)
        endif
		tag1=11
		! send to the bottom:
		if(mp1%face%s_bottom /= mp1%face_v%id) then
		    count=count+1
            call MPI_Isend(array(1:u_h,1:jpp,1:ipp), &
                (ipp*jpp)*u_h, MPI_REAL8, mp1%face_v%s_bottom, &
                tag1, comm3d, request(count),error)
        endif
        
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! case where only 1 pe in x or y directions                                      !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if ( mp1%face_v%s_top == -1) then
			! adjacent cells:
			array(kpp+1:kpp+u_h,1-s_h:jpp+n_h,1-w_h:ipp+e_h)=ubc
			! corner cells - not relevant, because below surface and above lid			
		endif
		if ( mp1%face_v%s_bottom == -1) then
			! adjacent cells:
			array(1-d_h:0,1-s_h:jpp+n_h,1-w_h:ipp+e_h)=lbc
			! corner cells - not relevant, because below surface and above lid			
		endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! wait for all messages to return before writing buffers                         !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call MPI_Waitall(count,request(1:count), status(:,1:count), error)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
        
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! write buffers                                                                  !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if(mp1%face_v%s_bottom /= -1) then
            array(0,1:jpp,1:ipp)=bot_buf
        endif
        
		if(mp1%face_v%s_top /= -1) then
            array(kpp+1,1:jpp,1:ipp)=top_buf
        endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		
	end subroutine exchange_along_z
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! exchange along dim for a variable using Cartesian topology                         !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>define some types to be used in the model
	!>@param[in] comm3d, id, ipp, jpp, kpp, w_h,e_h,s_h,n_h,d_h,u_h
	!>@param[inout] array: the array to exchange_halos on
	!>@param[in] lbc, ubc, dims,coords
	subroutine exchange_along_dim(comm3d, id, kpp, jpp, ipp, &
							d_h,u_h,s_h,n_h,w_h, e_h,  array, lbc, ubc, dims,coords)
		implicit none
		
		integer(i4b), intent(in) :: comm3d, id, ipp, jpp, kpp, w_h, e_h, s_h,n_h,d_h,u_h
		real(sp), intent(inout), &
			 dimension(1-d_h:u_h+kpp,1-s_h:n_h+jpp,1-w_h:e_h+ipp) :: &
			 array
		real(sp), intent(in) :: lbc, ubc
		integer(i4b), dimension(3), intent(in) :: dims,coords
		
		! locals:
		integer(i4b), dimension(12) :: request
		integer(i4b), dimension(MPI_STATUS_SIZE, 12) :: status
		integer(i4b) :: error, tag1,num_messages,imess, tag2
		
		! buffers for receiving data
	    integer(i4b) :: count
		real(sp), dimension(jpp,ipp) :: top_buf, bot_buf
		real(sp), dimension(kpp,jpp) :: east_buf,west_buf
		real(sp), dimension(kpp,ipp) :: south_buf,north_buf
		
		

			
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! message passing for adjacent cells in up / down direction                      !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		tag1=10
		count=0

		! receive from the top of the lower cell:
        if(mp1%face%r_top /= id) then
		    count=count+1
            call MPI_IRecv(bot_buf, &!array(1-d_h:0,1:jpp,1:ipp), &
                (ipp*jpp)*d_h, MPI_REAL8, mp1%face%r_top, &
                tag1, comm3d, request(count),error)
        endif
        tag1=11
		! receive from the bottom of upper cell:
		if(mp1%face%r_bottom /= id) then
		    count=count+1
            call MPI_IRecv(top_buf, & !array(kpp+1:kpp+u_h,1:jpp,1:ipp), &
                (ipp*jpp)*u_h, MPI_REAL8, mp1%face%r_bottom, &
                tag1, comm3d, request(count),error)
        endif
        tag1=10
		! send to the top:
		if(mp1%face%s_top /= id) then
		    count=count+1
            call MPI_Isend(array(kpp+1-d_h:kpp,1:jpp,1:ipp), &
                (ipp*jpp)*d_h, MPI_REAL8, mp1%face%s_top, &
                tag1, comm3d, request(count),error)
        endif
		tag1=11
		! send to the bottom:
		if(mp1%face%s_bottom /= id) then
		    count=count+1
            call MPI_Isend(array(1:u_h,1:jpp,1:ipp), &
                (ipp*jpp)*u_h, MPI_REAL8, mp1%face%s_bottom, &
                tag1, comm3d, request(count),error)
        endif
        
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! message passing for adjacent cells in east / west direction                    !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		tag1=12
		! receive from the east of the west cell:
		if( mp1%face%r_east /= id ) then
		    count=count+1
			call MPI_IRecv(west_buf, & !array(1:kpp,1:jpp,1-w_h:0), &
				(jpp*kpp)*w_h, MPI_REAL8, mp1%face%r_east, &
				tag1, comm3d, request(count),error)
		endif
		tag1=13
		! receive from the west of east cell:
		if( mp1%face%r_west /= id ) then
		    count=count+1
			call MPI_IRecv(east_buf, & !array(1:kpp,1:jpp,ipp+1:ipp+e_h), &
				(jpp*kpp)*e_h, MPI_REAL8, mp1%face%r_west, &
				tag1, comm3d, request(count),error)
		endif
		
		tag1=12
		! send to the east:
		if ( mp1%face%s_east /= id ) then 
		    count=count+1
			call MPI_Isend(array(1:kpp,1:jpp,ipp+1-w_h:ipp), &
				(jpp*kpp)*w_h, MPI_REAL8, mp1%face%s_east, &
				tag1, comm3d, request(count),error)
		endif

		tag1=13
		! send to the west:
		if ( mp1%face%s_west /= id ) then	
		    count=count+1
			call MPI_Isend(array(1:kpp,1:jpp,1:e_h), &
				(jpp*kpp)*e_h, MPI_REAL8, mp1%face%s_west, &
				tag1, comm3d, request(count),error)
		endif

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! message passing for adjacent cells in north / south direction                  !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		tag1=14
		! receive from north of south cell:
		if( mp1%face%r_north /= id ) then
		    count=count+1
			call MPI_IRecv(south_buf, &!array(1:kpp,1-s_h:0,1:ipp), &
				(ipp*kpp)*s_h, MPI_REAL8, mp1%face%r_north, &
				tag1, comm3d, request(count),error)
		endif
		tag1=15
		! receive from south of north cell:
		if( mp1%face%r_south /= id ) then
		    count=count+1
			call MPI_IRecv(north_buf, &!array(1:kpp,jpp+1:jpp+n_h,1:ipp), &
				(ipp*kpp)*n_h, MPI_REAL8, mp1%face%r_south, &
				tag1, comm3d, request(count),error)
		endif
		tag1=14
		! send to the north:
		if ( mp1%face%s_north /= id ) then 
		    count=count+1
			call MPI_Isend(array(1:kpp,jpp+1-s_h:jpp,1:ipp), &
				(ipp*kpp)*s_h, MPI_REAL8, mp1%face%s_north, &
				tag1, comm3d, request(count),error)
		endif
		tag1=15
		! send to the south:
		if ( mp1%face%s_south /= id ) then
		    count=count+1
			call MPI_Isend(array(1:kpp,1:n_h,1:ipp), &
				(ipp*kpp)*n_h, MPI_REAL8, mp1%face%s_south, &
				tag1, comm3d, request(count),error)
		endif

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! case where only 1 pe in x or y directions                                      !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if ( mp1%face%s_west == id) then
			! adjacent cells:
			array(1:kpp,1:jpp,1-w_h:0)=array(1:kpp,1:jpp,ipp+1-w_h:ipp)
			array(1:kpp,1:jpp,ipp+1:ipp+e_h)=array(1:kpp,1:jpp,1:e_h)
			west_buf=array(1:kpp,1:jpp,ipp)
			east_buf=array(1:kpp,1:jpp,1)
			! corner cells - not relevant, because below surface and above lid			
		endif
				
		if ( mp1%face%s_south == id) then
			! adjacent cells:
			array(1:kpp,1-s_h:0,1:ipp)=array(1:kpp,jpp+1-s_h:jpp,1:ipp)
			array(1:kpp,jpp+1:jpp+n_h,1:ipp)=array(1:kpp,1:n_h,1:ipp)
			south_buf=array(1:kpp,jpp,1:ipp)
			north_buf=array(1:kpp,1,1:ipp)
			! corner cells - not relevant, because below surface and above lid			
		endif
		if ( mp1%face%s_top == -1) then
			! adjacent cells:
			array(kpp+1:kpp+u_h,1-s_h:jpp+n_h,1-w_h:ipp+e_h)=ubc
			! corner cells - not relevant, because below surface and above lid			
		endif
		if ( mp1%face%s_bottom == -1) then
			! adjacent cells:
			array(1-d_h:0,1-s_h:jpp+n_h,1-w_h:ipp+e_h)=lbc
			! corner cells - not relevant, because below surface and above lid			
		endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! wait for all messages to return before writing buffers                         !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call MPI_Waitall(count,request(1:count), status(:,1:count), error)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! write buffers                                                                  !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if(mp1%face%s_bottom /= -1) then
            array(0,1:jpp,1:ipp)=bot_buf
        endif
        
		if(mp1%face%s_top /= -1) then
            array(kpp+1,1:jpp,1:ipp)=top_buf
        endif
        array(1:kpp,1:jpp,0)=west_buf
        array(1:kpp,1:jpp,ipp+1)=east_buf
        array(1:kpp,0,1:ipp)=south_buf
        array(1:kpp,jpp+1,1:ipp)=north_buf
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        

	end subroutine exchange_along_dim
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! exchange halos for a variable using Cartesian topology                             !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>define some types to be used in the model
	!>@param[in] comm3d, id, ipp, jpp, kpp, w_h,e_h,s_h,n_h,d_h,u_h
	!>@param[inout] array: the array to exchange_halos on
	!>@param[in] lbc, ubc, dims,coords
	subroutine exchange_full(comm3d, id, kpp, jpp, ipp, &
							d_h,u_h,s_h,n_h,w_h, e_h,  array, lbc,ubc, dims,coords)
		implicit none
		
		integer(i4b), intent(in) :: comm3d, id, ipp, jpp, kpp, w_h, e_h, s_h,n_h,d_h,u_h
		real(sp), intent(inout), &
			 dimension(1-d_h:u_h+kpp,1-s_h:n_h+jpp,1-w_h:e_h+ipp) :: &
			 array
		real(sp), intent(in) :: lbc, ubc
		integer(i4b), dimension(3), intent(in) :: dims,coords
		
		! locals:
		integer(i4b), dimension(40) :: request
		integer(i4b), dimension(MPI_STATUS_SIZE, 40) :: status
		integer(i4b) :: error, tag1,num_messages,imess, tag2, count
		
		! buffers for receiving data
		! edges
		real(sp), dimension(kpp) :: ws_bt, wn_bt, es_bt, en_bt
		real(sp), dimension(ipp) :: bs_we, bn_we, ts_we, tn_we
		real(sp), dimension(jpp) :: bw_sn, be_sn, tw_sn, te_sn
		! corners
		real(sp) :: wsb, wnb, wst, wnt, esb, enb, est, ent
				

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! exchange only only the dimensions                                              !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		call exchange_along_dim(comm3d, id, kpp, jpp, ipp, &
							d_h,u_h,s_h,n_h,w_h, e_h,  array, lbc, ubc, dims,coords)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! edges                                                                          !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		count=0
		! ws_bt
		call exchange_edges(comm3d, id,count,40,kpp,ws_bt,request,status, &
		            array,kpp,jpp,ipp,d_h,u_h,s_h,n_h,w_h, e_h, &
					1,   kpp,     1,          n_h,       1,     e_h, &
					1,   kpp,     jpp+1,      jpp+n_h,   ipp+1, ipp+e_h, &
					mp1%edge%s_ws_bt,mp1%edge%r_ws_bt)
		! wn_bt
		call exchange_edges(comm3d, id,count,40,kpp,wn_bt,request,status, &
		            array,kpp,jpp,ipp,d_h,u_h,s_h,n_h,w_h, e_h, &
					1,   kpp,     jpp-s_h+1,    jpp,       1,     e_h, &
					1,   kpp,     1-s_h,      0,        ipp+1, ipp+e_h, &
					mp1%edge%s_wn_bt,mp1%edge%r_wn_bt)
		! es_bt
		call exchange_edges(comm3d, id,count,40,kpp,es_bt,request,status, &
		            array,kpp,jpp,ipp,d_h,u_h,s_h,n_h,w_h, e_h, &
					1,   kpp,     1,          n_h,       ipp-w_h+1,     ipp, &
					1,   kpp,     jpp+1,      jpp+n_h,   1-w_h, 0, &
					mp1%edge%s_es_bt,mp1%edge%r_es_bt)
		! en_bt
		call exchange_edges(comm3d, id,count,40,kpp,en_bt,request,status, &
		            array,kpp,jpp,ipp,d_h,u_h,s_h,n_h,w_h, e_h, &
					1,   kpp,     jpp-s_h+1,    jpp,       ipp-w_h+1,     ipp, &
					1,   kpp,     1-s_h,      0,        1-w_h, 0, &
					mp1%edge%s_en_bt,mp1%edge%r_en_bt)
					
		! bs_we
		call exchange_edges(comm3d, id,count,40,ipp,bs_we,request,status, &
		            array,kpp,jpp,ipp,d_h,u_h,s_h,n_h,w_h, e_h, &
					1,   		 u_h,        1,       n_h,     1, ipp, &
					kpp+1,       kpp+u_h,    jpp+1,   jpp+n_h,   1, ipp, &
					mp1%edge%s_bs_we,mp1%edge%r_bs_we)
		! bn_we
		call exchange_edges(comm3d, id,count,40,ipp,bn_we,request,status, &
		            array,kpp,jpp,ipp,d_h,u_h,s_h,n_h,w_h, e_h, &
					1,   		 u_h,        jpp-s_h+1,       jpp,     1, ipp, &
					kpp+1,       kpp+u_h,    1-s_h,           0,       1, ipp, &
					mp1%edge%s_bn_we,mp1%edge%r_bn_we)
		! ts_we
		call exchange_edges(comm3d, id,count,40,ipp,ts_we,request,status, &
		            array,kpp,jpp,ipp,d_h,u_h,s_h,n_h,w_h, e_h, &
					kpp+1-d_h,   		 kpp,        1,       n_h,     1, ipp, &
					1-d_h,       0,    jpp+1,   jpp+n_h,   1, ipp, &
					mp1%edge%s_ts_we,mp1%edge%r_ts_we)

		! tn_we
		call exchange_edges(comm3d, id,count,40,ipp,tn_we,request,status, &
		            array,kpp,jpp,ipp,d_h,u_h,s_h,n_h,w_h, e_h, &
					kpp+1-d_h,   kpp,     jpp+1-s_h,  jpp, 1, ipp, &
					1-d_h,       0,       1-s_h,      0,   1, ipp, &
					mp1%edge%s_tn_we,mp1%edge%r_tn_we)

		! bw_sn
		call exchange_edges(comm3d, id,count,40,jpp,bw_sn,request,status, &
		            array,kpp,jpp,ipp,d_h,u_h,s_h,n_h,w_h, e_h, &
					1,   		 u_h,        1,       jpp,     1, e_h, &
					kpp+1,       kpp+u_h,    1,       jpp,   ipp+1, ipp+e_h, &
					mp1%edge%s_bw_sn,mp1%edge%r_bw_sn)
		! be_sn
		call exchange_edges(comm3d, id,count,40,jpp,be_sn,request,status, &
		            array,kpp,jpp,ipp,d_h,u_h,s_h,n_h,w_h, e_h, &
					1,   		 u_h,        1,       jpp,     ipp-w_h+1, ipp, &
					kpp+1,       kpp+u_h,    1,       jpp,   1-w_h, 0, &
					mp1%edge%s_be_sn,mp1%edge%r_be_sn)
		! tw_sn
		call exchange_edges(comm3d, id,count,40,jpp,tw_sn,request,status, &
		            array,kpp,jpp,ipp,d_h,u_h,s_h,n_h,w_h, e_h, &
					kpp+1-d_h,   		 kpp,        1,       jpp,     1, e_h, &
					1-d_h,       0,    1,       jpp,   ipp+1, ipp+e_h, &
					mp1%edge%s_tw_sn,mp1%edge%r_tw_sn)
		! te_sn
		call exchange_edges(comm3d, id,count,40,jpp,te_sn,request,status, &
		            array,kpp,jpp,ipp,d_h,u_h,s_h,n_h,w_h, e_h, &
					kpp+1-d_h,   		 kpp,      1,       jpp,     ipp-w_h+1, ipp, &
					1-d_h,       0,      1,       jpp,   1-w_h, 0, &
					mp1%edge%s_te_sn,mp1%edge%r_te_sn)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! corners                                                                        !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! wsb
		call exchange_corner(comm3d, id,count,40,wsb,request,status, &
		            array,kpp,jpp,ipp,d_h,u_h,s_h,n_h,w_h, e_h, &
					1,   u_h,     1,          n_h,       1,     e_h, &
					kpp+1,   kpp+u_h,     jpp+1,      jpp+n_h,   ipp+1, ipp+e_h, &
					mp1%cnr%s_wsb,mp1%cnr%r_wsb)
		! wnb
		call exchange_corner(comm3d, id,count,40,wnb,request,status, &
		            array,kpp,jpp,ipp,d_h,u_h,s_h,n_h,w_h, e_h, &
					1,   u_h,     jpp+1-s_h,    jpp,       1,     e_h, &
					kpp+1,   kpp+u_h,     1-s_h,      0,        ipp+1, ipp+e_h, &
					mp1%cnr%s_wnb,mp1%cnr%r_wnb)
		! wst
		call exchange_corner(comm3d, id,count,40,wst,request,status, &
		            array,kpp,jpp,ipp,d_h,u_h,s_h,n_h,w_h, e_h, &
					kpp+1-d_h,   kpp,     1,          n_h,       1,     e_h, &
					1-d_h,       0,     jpp+1,      jpp+n_h,   ipp+1, ipp+e_h, &
					mp1%cnr%s_wst,mp1%cnr%r_wst)
		! wnt
		call exchange_corner(comm3d, id,count,40,wnt,request,status, &
		            array,kpp,jpp,ipp,d_h,u_h,s_h,n_h,w_h, e_h, &
					kpp+1-d_h,   kpp,     jpp+1-s_h,    jpp,       1,     e_h, &
					1-d_h,        0,     1-s_h,      0,        ipp+1, ipp+e_h, &
					mp1%cnr%s_wnt,mp1%cnr%r_wnt)
				
					
		! esb
		call exchange_corner(comm3d, id,count,40,esb,request,status, &
		            array,kpp,jpp,ipp,d_h,u_h,s_h,n_h,w_h, e_h,&
					1,   		 u_h,        1,       n_h,     ipp+1-w_h, ipp, &
					kpp+1,       kpp+u_h,    jpp+1,   jpp+n_h,   1-w_h, 0, &
					mp1%cnr%s_esb,mp1%cnr%r_esb)
		! enb
		call exchange_corner(comm3d, id,count,40,enb,request,status, &
		            array,kpp,jpp,ipp,d_h,u_h,s_h,n_h,w_h, e_h, &
					1,   		 u_h,        jpp+1-s_h,       jpp,     ipp+1-w_h, ipp, &
					kpp+1,       kpp+u_h,    1-s_h,           0,       1-w_h, 0, &
					mp1%cnr%s_enb,mp1%cnr%r_enb)
		! est
		call exchange_corner(comm3d, id,count,40,est,request,status, &
		            array,kpp,jpp,ipp,d_h,u_h,s_h,n_h,w_h, e_h, &
					kpp+1-d_h,   		 kpp,        1,       n_h,     ipp+1-w_h, ipp, &
					1-d_h,       0,             jpp+1,   jpp+n_h,   1-w_h, 0, &
					mp1%cnr%s_est,mp1%cnr%r_est)

		! ent
		call exchange_corner(comm3d, id,count,40,ent,request,status, &
		            array,kpp,jpp,ipp,d_h,u_h,s_h,n_h,w_h, e_h, &
					kpp+1-d_h,   kpp,     jpp+1-s_h,  jpp, ipp+1-w_h, ipp, &
					1-d_h,       0,       1-s_h,      0,   1-w_h, 0, &
					mp1%cnr%s_ent,mp1%cnr%r_ent)


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! wait for all messages to return before writing buffers                         !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call MPI_Waitall(count,request(1:count), status(:,1:count), error)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! write buffers                                                                  !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! edges:
        array(1:kpp,jpp+1,ipp+1)=ws_bt
        array(1:kpp,1-s_h,ipp+1)=wn_bt
        array(1:kpp,jpp+1,1-w_h)=es_bt
        array(1:kpp,1-s_h,1-w_h)=en_bt

        array(kpp+1,jpp+1,1:ipp)=bs_we
        array(kpp+1,1-s_h,1:ipp)=bn_we
        array(1-d_h,jpp+1,1:ipp)=ts_we
        array(1-d_h,1-s_h,1:ipp)=tn_we

        array(kpp+1,1:jpp,ipp+1)=bw_sn
        array(kpp+1,1:jpp,1-w_h)=be_sn
        array(1-d_h,1:jpp,ipp+1)=tw_sn
        array(1-d_h,1:jpp,1-w_h)=te_sn

        ! corners:
        if( mp1%cnr%r_wsb /=id ) array(kpp+1,jpp+1,ipp+1)=wsb
        if( mp1%cnr%r_wnb /=id ) array(kpp+1,1-s_h,ipp+1)=wnb
        if( mp1%cnr%r_wst /=id ) array(1-d_h,jpp+1,ipp+1)=wst
        if( mp1%cnr%r_wnt /=id ) array(1-d_h,1-s_h,ipp+1)=wnt
        if( mp1%cnr%r_esb /=id ) array(kpp+1,jpp+1,1-w_h)=esb
        if( mp1%cnr%r_enb /=id ) array(kpp+1,1-s_h,1-w_h)=enb
        if( mp1%cnr%r_est /=id ) array(1-d_h,jpp+1,1-w_h)=est
        if( mp1%cnr%r_ent /=id ) array(1-d_h,1-s_h,1-w_h)=ent
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				
		contains
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! exchange corners for a variable using Cartesian topology                   !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!>@author
			!>Paul J. Connolly, The University of Manchester
			!>@brief
			!>exchange corners
			!>@param[in] comm3d, id
			!>@param[inout] array: the array to exchange_halos on
			!>@param[inout] buf: buffer for the receive
			!>@param[in] count,lenm
			!>@param[inout] request, status (mpi variables)
			!>@param[in] kp, jp, ip, d_h_, u_h, s_h, n_h, w_h, e_h, s_kl, s_ku, s_jl, 
			!            s_ju, s_il, s_iu, r_kl, r_ku, r_jl, r_ju, r_il, r_iu
			!>@param[in] send, recv
			subroutine exchange_corner(comm3d,id,count,lenm,buf,request,status, &
			                        array,kp,jp,ip,d_h,u_h,s_h,n_h,w_h, e_h, &
									s_kl,s_ku,s_jl,s_ju,s_il,s_iu, &
									r_kl,r_ku,r_jl,r_ju,r_il,r_iu, send,recv)
				implicit none
				integer(i4b), intent(in) :: comm3d, id
				integer(i4b), intent(in) :: kp,jp,ip,w_h, e_h, s_h,n_h,d_h,u_h,   &
											s_kl,s_ku,s_jl,s_ju,s_il,s_iu, &
											r_kl,r_ku,r_jl,r_ju,r_il,r_iu, send,recv, &
											lenm
				real(sp), intent(inout), &
					 dimension(1-d_h:u_h+kp,1-s_h:n_h+jp,1-w_h:e_h+ip) :: &
					 array
			    real(sp), intent(inout) :: buf
				integer(i4b), dimension(lenm), intent(inout) :: request
				integer(i4b), dimension(MPI_STATUS_SIZE,lenm), intent(inout) :: status
				integer(i4b), intent(inout) :: count
				integer(i4b) :: error,tag1=3, size1, size2 

				integer(i4b), dimension(lenm) :: request1
				integer(i4b), dimension(MPI_STATUS_SIZE,lenm) :: status1

				size1=(s_iu-s_il+1)*(s_ju-s_jl+1)*(s_ku-s_kl+1)
				size2=(r_iu-r_il+1)*(r_ju-r_jl+1)*(r_ku-r_kl+1)

                buf = array(s_kl,s_jl,s_il)
				!++++
				! receive:
				if( (recv /= id) ) then
				    count=count+1
					call MPI_IRecv(buf, &!array(r_kl:r_ku,r_jl:r_ju,r_il:r_iu), &
						size1, MPI_REAL8, recv, &
						tag1, comm3d, request(count),error)
					!call MPI_Wait(request1(1), status1(:,1), error)
				endif	
				! send:
				if ( (send /= id) ) then 
				    count=count+1
					call MPI_Isend(array(s_kl:s_ku,s_jl:s_ju,s_il:s_iu), &
						size1, MPI_REAL8, send, &
						tag1, comm3d, request(count),error)
				endif
				!----	
				
				
			end subroutine exchange_corner
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



		
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! exchange edges for a variable using Cartesian topology                     !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!>@author
			!>Paul J. Connolly, The University of Manchester
			!>@brief
			!>exchange corners
			!>@param[in] comm3d, id
			!>@param[inout] array: the array to exchange_halos on
			!>@param[inout] buf: buffer for the receive
			!>@param[in] count,lenm,lenb
			!>@param[inout] request, status (mpi variables)
			!>@param[in] kp, jp, ip, d_h_, u_h, s_h, n_h, w_h, e_h, s_kl, s_ku, s_jl, 
			!            s_ju, s_il, s_iu, r_kl, r_ku, r_jl, r_ju, r_il, r_iu
			!>@param[in] send, recv
			subroutine exchange_edges(comm3d,id,count,lenm,lenb,buf,request,status, &
			                        array,kp,jp,ip,d_h,u_h,s_h,n_h,w_h, e_h, &
									s_kl,s_ku,s_jl,s_ju,s_il,s_iu, &
									r_kl,r_ku,r_jl,r_ju,r_il,r_iu, send,recv)
				implicit none
				integer(i4b), intent(in) :: comm3d, id
				integer(i4b), intent(in) :: kp,jp,ip, w_h, e_h, s_h,n_h,d_h,u_h, &
											s_kl,s_ku,s_jl,s_ju,s_il,s_iu, &
											r_kl,r_ku,r_jl,r_ju,r_il,r_iu, send,recv, &
											lenm,lenb
				real(sp), intent(inout), &
					 dimension(1-d_h:u_h+kp,1-s_h:n_h+jp,1-w_h:e_h+ip) :: &
					 array
			    real(sp), intent(inout), dimension(lenb) :: buf
			    integer(i4b), intent(inout) :: count
				integer(i4b), dimension(lenm), intent(inout) :: request
				integer(i4b), dimension(MPI_STATUS_SIZE,lenm), intent(inout) :: status
				integer(i4b) :: error,tag1=2, size1, size2


				size1=(s_iu-s_il+1)*(s_ju-s_jl+1)*(s_ku-s_kl+1)
				size2=(r_iu-r_il+1)*(r_ju-r_jl+1)*(r_ku-r_kl+1)

                buf = &
                    reshape(array(s_kl:s_ku,s_jl:s_ju,s_il:s_iu),&
                        [(s_ku-s_kl+1)*(s_ju-s_jl+1)*(s_iu-s_il+1)])

				!++++
				! receive:
				if( (recv /= mp1%id) ) then
				    count=count+1
! 
					call MPI_IRecv(buf, &!array(r_kl:r_ku,r_jl:r_ju,r_il:r_iu), &
						size1, MPI_REAL8, recv, &
						tag1, comm3d, request(count),error)
				endif	
				! send:
				if ( (send /= mp1%id) ) then 
				    count=count+1
					call MPI_Isend(array(s_kl:s_ku,s_jl:s_ju,s_il:s_iu), &
						size1, MPI_REAL8, send, &
						tag1, comm3d, request(count),error)
				endif
				!----
		
									
									
			end subroutine exchange_edges
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
	end subroutine exchange_full
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! exchange along dim for a variable using Cartesian topology                         !
	! without top and bottom bc set                                                      !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>define some types to be used in the model
	!>@param[in] comm3d, id, ipp, jpp, kpp, w_h,e_h,s_h,n_h,d_h,u_h
	!>@param[inout] array: the array to exchange_halos on
	!>@param[in] dims,coords
	subroutine exchange_along_dim_wo(comm3d, id, kpp, jpp, ipp, &
							d_h,u_h,s_h,n_h,w_h, e_h,  array, dims,coords)
		implicit none
		
		integer(i4b), intent(in) :: comm3d, id, ipp, jpp, kpp, w_h, e_h, s_h,n_h,d_h,u_h
		real(sp), intent(inout), &
			 dimension(1-d_h:u_h+kpp,1-s_h:n_h+jpp,1-w_h:e_h+ipp) :: &
			 array
		integer(i4b), dimension(3), intent(in) :: dims,coords
		
		! locals:
		integer(i4b), dimension(12) :: request
		integer(i4b), dimension(MPI_STATUS_SIZE, 12) :: status
		integer(i4b) :: error, tag1,num_messages,imess, tag2
		
		! buffers for receiving data
	    integer(i4b) :: count
		real(sp), dimension(jpp,ipp) :: top_buf, bot_buf
		real(sp), dimension(kpp,jpp) :: east_buf,west_buf
		real(sp), dimension(kpp,ipp) :: south_buf,north_buf
		
		

			
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! message passing for adjacent cells in up / down direction                      !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		tag1=10
		count=0

		! receive from the top of the lower cell:
        if(mp1%face%r_top /= id) then
		    count=count+1
            call MPI_IRecv(bot_buf, &!array(1-d_h:0,1:jpp,1:ipp), &
                (ipp*jpp)*d_h, MPI_REAL8, mp1%face%r_top, &
                tag1, comm3d, request(count),error)
        endif
        tag1=11
		! receive from the bottom of upper cell:
		if(mp1%face%r_bottom /= id) then
		    count=count+1
            call MPI_IRecv(top_buf, & !array(kpp+1:kpp+u_h,1:jpp,1:ipp), &
                (ipp*jpp)*u_h, MPI_REAL8, mp1%face%r_bottom, &
                tag1, comm3d, request(count),error)
        endif
        tag1=10
		! send to the top:
		if(mp1%face%s_top /= id) then
		    count=count+1
            call MPI_Isend(array(kpp+1-d_h:kpp,1:jpp,1:ipp), &
                (ipp*jpp)*d_h, MPI_REAL8, mp1%face%s_top, &
                tag1, comm3d, request(count),error)
        endif
		tag1=11
		! send to the bottom:
		if(mp1%face%s_bottom /= id) then
		    count=count+1
            call MPI_Isend(array(1:u_h,1:jpp,1:ipp), &
                (ipp*jpp)*u_h, MPI_REAL8, mp1%face%s_bottom, &
                tag1, comm3d, request(count),error)
        endif
        
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! message passing for adjacent cells in east / west direction                    !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		tag1=12
		! receive from the east of the west cell:
		if( mp1%face%r_east /= id ) then
		    count=count+1
			call MPI_IRecv(west_buf, & !array(1:kpp,1:jpp,1-w_h:0), &
				(jpp*kpp)*w_h, MPI_REAL8, mp1%face%r_east, &
				tag1, comm3d, request(count),error)
		endif
		tag1=13
		! receive from the west of east cell:
		if( mp1%face%r_west /= id ) then
		    count=count+1
			call MPI_IRecv(east_buf, & !array(1:kpp,1:jpp,ipp+1:ipp+e_h), &
				(jpp*kpp)*e_h, MPI_REAL8, mp1%face%r_west, &
				tag1, comm3d, request(count),error)
		endif
		
		tag1=12
		! send to the east:
		if ( mp1%face%s_east /= id ) then 
		    count=count+1
			call MPI_Isend(array(1:kpp,1:jpp,ipp+1-w_h:ipp), &
				(jpp*kpp)*w_h, MPI_REAL8, mp1%face%s_east, &
				tag1, comm3d, request(count),error)
		endif

		tag1=13
		! send to the west:
		if ( mp1%face%s_west /= id ) then	
		    count=count+1
			call MPI_Isend(array(1:kpp,1:jpp,1:e_h), &
				(jpp*kpp)*e_h, MPI_REAL8, mp1%face%s_west, &
				tag1, comm3d, request(count),error)
		endif

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! message passing for adjacent cells in north / south direction                  !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		tag1=14
		! receive from north of south cell:
		if( mp1%face%r_north /= id ) then
		    count=count+1
			call MPI_IRecv(south_buf, &!array(1:kpp,1-s_h:0,1:ipp), &
				(ipp*kpp)*s_h, MPI_REAL8, mp1%face%r_north, &
				tag1, comm3d, request(count),error)
		endif
		tag1=15
		! receive from south of north cell:
		if( mp1%face%r_south /= id ) then
		    count=count+1
			call MPI_IRecv(north_buf, &!array(1:kpp,jpp+1:jpp+n_h,1:ipp), &
				(ipp*kpp)*n_h, MPI_REAL8, mp1%face%r_south, &
				tag1, comm3d, request(count),error)
		endif
		tag1=14
		! send to the north:
		if ( mp1%face%s_north /= id ) then 
		    count=count+1
			call MPI_Isend(array(1:kpp,jpp+1-s_h:jpp,1:ipp), &
				(ipp*kpp)*s_h, MPI_REAL8, mp1%face%s_north, &
				tag1, comm3d, request(count),error)
		endif
		tag1=15
		! send to the south:
		if ( mp1%face%s_south /= id ) then
		    count=count+1
			call MPI_Isend(array(1:kpp,1:n_h,1:ipp), &
				(ipp*kpp)*n_h, MPI_REAL8, mp1%face%s_south, &
				tag1, comm3d, request(count),error)
		endif

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! case where only 1 pe in x or y directions                                      !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if ( mp1%face%s_west == id) then
			! adjacent cells:
			array(1:kpp,1:jpp,1-w_h:0)=array(1:kpp,1:jpp,ipp+1-w_h:ipp)
			array(1:kpp,1:jpp,ipp+1:ipp+e_h)=array(1:kpp,1:jpp,1:e_h)
			west_buf=array(1:kpp,1:jpp,ipp)
			east_buf=array(1:kpp,1:jpp,1)
			! corner cells - not relevant, because below surface and above lid			
		endif
				
		if ( mp1%face%s_south == id) then
			! adjacent cells:
			array(1:kpp,1-s_h:0,1:ipp)=array(1:kpp,jpp+1-s_h:jpp,1:ipp)
			array(1:kpp,jpp+1:jpp+n_h,1:ipp)=array(1:kpp,1:n_h,1:ipp)
			south_buf=array(1:kpp,jpp,1:ipp)
			north_buf=array(1:kpp,1,1:ipp)
			! corner cells - not relevant, because below surface and above lid			
		endif
! 		if ( mp1%face%s_top == -1) then
! 			! adjacent cells:
! 			array(kpp+1:kpp+u_h,1-s_h:jpp+n_h,1-w_h:ipp+e_h)=ubc
! 			! corner cells - not relevant, because below surface and above lid			
! 		endif
! 		if ( mp1%face%s_bottom == -1) then
! 			! adjacent cells:
! 			array(1-d_h:0,1-s_h:jpp+n_h,1-w_h:ipp+e_h)=lbc
! 			! corner cells - not relevant, because below surface and above lid			
! 		endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! wait for all messages to return before writing buffers                         !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call MPI_Waitall(count,request(1:count), status(:,1:count), error)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! write buffers                                                                  !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if(mp1%face%s_bottom /= -1) then
            array(0,1:jpp,1:ipp)=bot_buf
        endif
        
		if(mp1%face%s_top /= -1) then
            array(kpp+1,1:jpp,1:ipp)=top_buf
        endif
        array(1:kpp,1:jpp,0)=west_buf
        array(1:kpp,1:jpp,ipp+1)=east_buf
        array(1:kpp,0,1:ipp)=south_buf
        array(1:kpp,jpp+1,1:ipp)=north_buf
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        

	end subroutine exchange_along_dim_wo
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


		
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! exchange halos for a variable using Cartesian topology                             !
	! without top and bottom being set                                                   !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>define some types to be used in the model
	!>@param[in] comm3d, id, ipp, jpp, kpp, w_h,e_h,s_h,n_h,d_h,u_h
	!>@param[inout] array: the array to exchange_halos on
	!>@param[in] dims,coords
	subroutine exchange_full_wo(comm3d, id, kpp, jpp, ipp, &
							d_h,u_h,s_h,n_h,w_h, e_h,  array, dims,coords)
		implicit none
		
		integer(i4b), intent(in) :: comm3d, id, ipp, jpp, kpp, w_h, e_h, s_h,n_h,d_h,u_h
		real(sp), intent(inout), &
			 dimension(1-d_h:u_h+kpp,1-s_h:n_h+jpp,1-w_h:e_h+ipp) :: &
			 array
		integer(i4b), dimension(3), intent(in) :: dims,coords
		
		! locals:
		integer(i4b), dimension(40) :: request
		integer(i4b), dimension(MPI_STATUS_SIZE, 40) :: status
		integer(i4b) :: error, tag1,num_messages,imess, tag2, count
		
		! buffers for receiving data
		! edges
		real(sp), dimension(kpp) :: ws_bt, wn_bt, es_bt, en_bt
		real(sp), dimension(ipp) :: bs_we, bn_we, ts_we, tn_we
		real(sp), dimension(jpp) :: bw_sn, be_sn, tw_sn, te_sn
		! corners
		real(sp) :: wsb, wnb, wst, wnt, esb, enb, est, ent
				

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! exchange only only the dimensions                                              !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		call exchange_along_dim_wo(comm3d, id, kpp, jpp, ipp, &
							d_h,u_h,s_h,n_h,w_h, e_h,  array, dims,coords)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! edges                                                                          !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		count=0
		! ws_bt
		call exchange_edges(comm3d, id,count,40,kpp,ws_bt,request,status, &
		            array,kpp,jpp,ipp,d_h,u_h,s_h,n_h,w_h, e_h, &
					1,   kpp,     1,          n_h,       1,     e_h, &
					1,   kpp,     jpp+1,      jpp+n_h,   ipp+1, ipp+e_h, &
					mp1%edge%s_ws_bt,mp1%edge%r_ws_bt)
		! wn_bt
		call exchange_edges(comm3d, id,count,40,kpp,wn_bt,request,status, &
		            array,kpp,jpp,ipp,d_h,u_h,s_h,n_h,w_h, e_h, &
					1,   kpp,     jpp-s_h+1,    jpp,       1,     e_h, &
					1,   kpp,     1-s_h,      0,        ipp+1, ipp+e_h, &
					mp1%edge%s_wn_bt,mp1%edge%r_wn_bt)
		! es_bt
		call exchange_edges(comm3d, id,count,40,kpp,es_bt,request,status, &
		            array,kpp,jpp,ipp,d_h,u_h,s_h,n_h,w_h, e_h, &
					1,   kpp,     1,          n_h,       ipp-w_h+1,     ipp, &
					1,   kpp,     jpp+1,      jpp+n_h,   1-w_h, 0, &
					mp1%edge%s_es_bt,mp1%edge%r_es_bt)
		! en_bt
		call exchange_edges(comm3d, id,count,40,kpp,en_bt,request,status, &
		            array,kpp,jpp,ipp,d_h,u_h,s_h,n_h,w_h, e_h, &
					1,   kpp,     jpp-s_h+1,    jpp,       ipp-w_h+1,     ipp, &
					1,   kpp,     1-s_h,      0,        1-w_h, 0, &
					mp1%edge%s_en_bt,mp1%edge%r_en_bt)
					
		! bs_we
		call exchange_edges(comm3d, id,count,40,ipp,bs_we,request,status, &
		            array,kpp,jpp,ipp,d_h,u_h,s_h,n_h,w_h, e_h, &
					1,   		 u_h,        1,       n_h,     1, ipp, &
					kpp+1,       kpp+u_h,    jpp+1,   jpp+n_h,   1, ipp, &
					mp1%edge%s_bs_we,mp1%edge%r_bs_we)
		! bn_we
		call exchange_edges(comm3d, id,count,40,ipp,bn_we,request,status, &
		            array,kpp,jpp,ipp,d_h,u_h,s_h,n_h,w_h, e_h, &
					1,   		 u_h,        jpp-s_h+1,       jpp,     1, ipp, &
					kpp+1,       kpp+u_h,    1-s_h,           0,       1, ipp, &
					mp1%edge%s_bn_we,mp1%edge%r_bn_we)
		! ts_we
		call exchange_edges(comm3d, id,count,40,ipp,ts_we,request,status, &
		            array,kpp,jpp,ipp,d_h,u_h,s_h,n_h,w_h, e_h, &
					kpp+1-d_h,   		 kpp,        1,       n_h,     1, ipp, &
					1-d_h,       0,    jpp+1,   jpp+n_h,   1, ipp, &
					mp1%edge%s_ts_we,mp1%edge%r_ts_we)

		! tn_we
		call exchange_edges(comm3d, id,count,40,ipp,tn_we,request,status, &
		            array,kpp,jpp,ipp,d_h,u_h,s_h,n_h,w_h, e_h, &
					kpp+1-d_h,   kpp,     jpp+1-s_h,  jpp, 1, ipp, &
					1-d_h,       0,       1-s_h,      0,   1, ipp, &
					mp1%edge%s_tn_we,mp1%edge%r_tn_we)

		! bw_sn
		call exchange_edges(comm3d, id,count,40,jpp,bw_sn,request,status, &
		            array,kpp,jpp,ipp,d_h,u_h,s_h,n_h,w_h, e_h, &
					1,   		 u_h,        1,       jpp,     1, e_h, &
					kpp+1,       kpp+u_h,    1,       jpp,   ipp+1, ipp+e_h, &
					mp1%edge%s_bw_sn,mp1%edge%r_bw_sn)
		! be_sn
		call exchange_edges(comm3d, id,count,40,jpp,be_sn,request,status, &
		            array,kpp,jpp,ipp,d_h,u_h,s_h,n_h,w_h, e_h, &
					1,   		 u_h,        1,       jpp,     ipp-w_h+1, ipp, &
					kpp+1,       kpp+u_h,    1,       jpp,   1-w_h, 0, &
					mp1%edge%s_be_sn,mp1%edge%r_be_sn)
		! tw_sn
		call exchange_edges(comm3d, id,count,40,jpp,tw_sn,request,status, &
		            array,kpp,jpp,ipp,d_h,u_h,s_h,n_h,w_h, e_h, &
					kpp+1-d_h,   		 kpp,        1,       jpp,     1, e_h, &
					1-d_h,       0,    1,       jpp,   ipp+1, ipp+e_h, &
					mp1%edge%s_tw_sn,mp1%edge%r_tw_sn)
		! te_sn
		call exchange_edges(comm3d, id,count,40,jpp,te_sn,request,status, &
		            array,kpp,jpp,ipp,d_h,u_h,s_h,n_h,w_h, e_h, &
					kpp+1-d_h,   		 kpp,      1,       jpp,     ipp-w_h+1, ipp, &
					1-d_h,       0,      1,       jpp,   1-w_h, 0, &
					mp1%edge%s_te_sn,mp1%edge%r_te_sn)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! corners                                                                        !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! wsb
		call exchange_corner(comm3d, id,count,40,wsb,request,status, &
		            array,kpp,jpp,ipp,d_h,u_h,s_h,n_h,w_h, e_h, &
					1,   u_h,     1,          n_h,       1,     e_h, &
					kpp+1,   kpp+u_h,     jpp+1,      jpp+n_h,   ipp+1, ipp+e_h, &
					mp1%cnr%s_wsb,mp1%cnr%r_wsb)
		! wnb
		call exchange_corner(comm3d, id,count,40,wnb,request,status, &
		            array,kpp,jpp,ipp,d_h,u_h,s_h,n_h,w_h, e_h, &
					1,   u_h,     jpp+1-s_h,    jpp,       1,     e_h, &
					kpp+1,   kpp+u_h,     1-s_h,      0,        ipp+1, ipp+e_h, &
					mp1%cnr%s_wnb,mp1%cnr%r_wnb)
		! wst
		call exchange_corner(comm3d, id,count,40,wst,request,status, &
		            array,kpp,jpp,ipp,d_h,u_h,s_h,n_h,w_h, e_h, &
					kpp+1-d_h,   kpp,     1,          n_h,       1,     e_h, &
					1-d_h,       0,     jpp+1,      jpp+n_h,   ipp+1, ipp+e_h, &
					mp1%cnr%s_wst,mp1%cnr%r_wst)
		! wnt
		call exchange_corner(comm3d, id,count,40,wnt,request,status, &
		            array,kpp,jpp,ipp,d_h,u_h,s_h,n_h,w_h, e_h, &
					kpp+1-d_h,   kpp,     jpp+1-s_h,    jpp,       1,     e_h, &
					1-d_h,        0,     1-s_h,      0,        ipp+1, ipp+e_h, &
					mp1%cnr%s_wnt,mp1%cnr%r_wnt)
				
					
		! esb
		call exchange_corner(comm3d, id,count,40,esb,request,status, &
		            array,kpp,jpp,ipp,d_h,u_h,s_h,n_h,w_h, e_h,&
					1,   		 u_h,        1,       n_h,     ipp+1-w_h, ipp, &
					kpp+1,       kpp+u_h,    jpp+1,   jpp+n_h,   1-w_h, 0, &
					mp1%cnr%s_esb,mp1%cnr%r_esb)
		! enb
		call exchange_corner(comm3d, id,count,40,enb,request,status, &
		            array,kpp,jpp,ipp,d_h,u_h,s_h,n_h,w_h, e_h, &
					1,   		 u_h,        jpp+1-s_h,       jpp,     ipp+1-w_h, ipp, &
					kpp+1,       kpp+u_h,    1-s_h,           0,       1-w_h, 0, &
					mp1%cnr%s_enb,mp1%cnr%r_enb)
		! est
		call exchange_corner(comm3d, id,count,40,est,request,status, &
		            array,kpp,jpp,ipp,d_h,u_h,s_h,n_h,w_h, e_h, &
					kpp+1-d_h,   		 kpp,        1,       n_h,     ipp+1-w_h, ipp, &
					1-d_h,       0,             jpp+1,   jpp+n_h,   1-w_h, 0, &
					mp1%cnr%s_est,mp1%cnr%r_est)

		! ent
		call exchange_corner(comm3d, id,count,40,ent,request,status, &
		            array,kpp,jpp,ipp,d_h,u_h,s_h,n_h,w_h, e_h, &
					kpp+1-d_h,   kpp,     jpp+1-s_h,  jpp, ipp+1-w_h, ipp, &
					1-d_h,       0,       1-s_h,      0,   1-w_h, 0, &
					mp1%cnr%s_ent,mp1%cnr%r_ent)


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! wait for all messages to return before writing buffers                         !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call MPI_Waitall(count,request(1:count), status(:,1:count), error)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! write buffers                                                                  !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! edges:
        array(1:kpp,jpp+1,ipp+1)=ws_bt
        array(1:kpp,1-s_h,ipp+1)=wn_bt
        array(1:kpp,jpp+1,1-w_h)=es_bt
        array(1:kpp,1-s_h,1-w_h)=en_bt

        array(kpp+1,jpp+1,1:ipp)=bs_we
        array(kpp+1,1-s_h,1:ipp)=bn_we
        array(1-d_h,jpp+1,1:ipp)=ts_we
        array(1-d_h,1-s_h,1:ipp)=tn_we

        array(kpp+1,1:jpp,ipp+1)=bw_sn
        array(kpp+1,1:jpp,1-w_h)=be_sn
        array(1-d_h,1:jpp,ipp+1)=tw_sn
        array(1-d_h,1:jpp,1-w_h)=te_sn

        ! corners:
        if( mp1%cnr%r_wsb /=id ) array(kpp+1,jpp+1,ipp+1)=wsb
        if( mp1%cnr%r_wnb /=id ) array(kpp+1,1-s_h,ipp+1)=wnb
        if( mp1%cnr%r_wst /=id ) array(1-d_h,jpp+1,ipp+1)=wst
        if( mp1%cnr%r_wnt /=id ) array(1-d_h,1-s_h,ipp+1)=wnt
        if( mp1%cnr%r_esb /=id ) array(kpp+1,jpp+1,1-w_h)=esb
        if( mp1%cnr%r_enb /=id ) array(kpp+1,1-s_h,1-w_h)=enb
        if( mp1%cnr%r_est /=id ) array(1-d_h,jpp+1,1-w_h)=est
        if( mp1%cnr%r_ent /=id ) array(1-d_h,1-s_h,1-w_h)=ent
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				
		contains
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! exchange corners for a variable using Cartesian topology                   !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!>@author
			!>Paul J. Connolly, The University of Manchester
			!>@brief
			!>exchange corners
			!>@param[in] comm3d, id
			!>@param[inout] array: the array to exchange_halos on
			!>@param[inout] buf: buffer for the receive
			!>@param[in] count,lenm
			!>@param[inout] request, status (mpi variables)
			!>@param[in] kp, jp, ip, d_h_, u_h, s_h, n_h, w_h, e_h, s_kl, s_ku, s_jl, 
			!            s_ju, s_il, s_iu, r_kl, r_ku, r_jl, r_ju, r_il, r_iu
			!>@param[in] send, recv
			subroutine exchange_corner(comm3d,id,count,lenm,buf,request,status, &
			                        array,kp,jp,ip,d_h,u_h,s_h,n_h,w_h, e_h, &
									s_kl,s_ku,s_jl,s_ju,s_il,s_iu, &
									r_kl,r_ku,r_jl,r_ju,r_il,r_iu, send,recv)
				implicit none
				integer(i4b), intent(in) :: comm3d, id
				integer(i4b), intent(in) :: kp,jp,ip,w_h, e_h, s_h,n_h,d_h,u_h,   &
											s_kl,s_ku,s_jl,s_ju,s_il,s_iu, &
											r_kl,r_ku,r_jl,r_ju,r_il,r_iu, send,recv, &
											lenm
				real(sp), intent(inout), &
					 dimension(1-d_h:u_h+kp,1-s_h:n_h+jp,1-w_h:e_h+ip) :: &
					 array
			    real(sp), intent(inout) :: buf
				integer(i4b), dimension(lenm), intent(inout) :: request
				integer(i4b), dimension(MPI_STATUS_SIZE,lenm), intent(inout) :: status
				integer(i4b), intent(inout) :: count
				integer(i4b) :: error,tag1=3, size1, size2 

				integer(i4b), dimension(lenm) :: request1
				integer(i4b), dimension(MPI_STATUS_SIZE,lenm) :: status1

				size1=(s_iu-s_il+1)*(s_ju-s_jl+1)*(s_ku-s_kl+1)
				size2=(r_iu-r_il+1)*(r_ju-r_jl+1)*(r_ku-r_kl+1)

                buf = array(s_kl,s_jl,s_il)
				!++++
				! receive:
				if( (recv /= id) ) then
				    count=count+1
					call MPI_IRecv(buf, &!array(r_kl:r_ku,r_jl:r_ju,r_il:r_iu), &
						size1, MPI_REAL8, recv, &
						tag1, comm3d, request(count),error)
					!call MPI_Wait(request1(1), status1(:,1), error)
				endif	
				! send:
				if ( (send /= id) ) then 
				    count=count+1
					call MPI_Isend(array(s_kl:s_ku,s_jl:s_ju,s_il:s_iu), &
						size1, MPI_REAL8, send, &
						tag1, comm3d, request(count),error)
				endif
				!----	
				
				
			end subroutine exchange_corner
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



		
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! exchange edges for a variable using Cartesian topology                     !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!>@author
			!>Paul J. Connolly, The University of Manchester
			!>@brief
			!>exchange corners
			!>@param[in] comm3d, id
			!>@param[inout] array: the array to exchange_halos on
			!>@param[inout] buf: buffer for the receive
			!>@param[in] count,lenm,lenb
			!>@param[inout] request, status (mpi variables)
			!>@param[in] kp, jp, ip, d_h_, u_h, s_h, n_h, w_h, e_h, s_kl, s_ku, s_jl, 
			!            s_ju, s_il, s_iu, r_kl, r_ku, r_jl, r_ju, r_il, r_iu
			!>@param[in] send, recv
			subroutine exchange_edges(comm3d,id,count,lenm,lenb,buf,request,status, &
			                        array,kp,jp,ip,d_h,u_h,s_h,n_h,w_h, e_h, &
									s_kl,s_ku,s_jl,s_ju,s_il,s_iu, &
									r_kl,r_ku,r_jl,r_ju,r_il,r_iu, send,recv)
				implicit none
				integer(i4b), intent(in) :: comm3d, id
				integer(i4b), intent(in) :: kp,jp,ip, w_h, e_h, s_h,n_h,d_h,u_h, &
											s_kl,s_ku,s_jl,s_ju,s_il,s_iu, &
											r_kl,r_ku,r_jl,r_ju,r_il,r_iu, send,recv, &
											lenm,lenb
				real(sp), intent(inout), &
					 dimension(1-d_h:u_h+kp,1-s_h:n_h+jp,1-w_h:e_h+ip) :: &
					 array
			    real(sp), intent(inout), dimension(lenb) :: buf
			    integer(i4b), intent(inout) :: count
				integer(i4b), dimension(lenm), intent(inout) :: request
				integer(i4b), dimension(MPI_STATUS_SIZE,lenm), intent(inout) :: status
				integer(i4b) :: error,tag1=2, size1, size2


				size1=(s_iu-s_il+1)*(s_ju-s_jl+1)*(s_ku-s_kl+1)
				size2=(r_iu-r_il+1)*(r_ju-r_jl+1)*(r_ku-r_kl+1)

                buf = &
                    reshape(array(s_kl:s_ku,s_jl:s_ju,s_il:s_iu),&
                        [(s_ku-s_kl+1)*(s_ju-s_jl+1)*(s_iu-s_il+1)])

				!++++
				! receive:
				if( (recv /= mp1%id) ) then
				    count=count+1
! 
					call MPI_IRecv(buf, &!array(r_kl:r_ku,r_jl:r_ju,r_il:r_iu), &
						size1, MPI_REAL8, recv, &
						tag1, comm3d, request(count),error)
				endif	
				! send:
				if ( (send /= mp1%id) ) then 
				    count=count+1
					call MPI_Isend(array(s_kl:s_ku,s_jl:s_ju,s_il:s_iu), &
						size1, MPI_REAL8, send, &
						tag1, comm3d, request(count),error)
				endif
				!----
		
									
									
			end subroutine exchange_edges
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
	end subroutine exchange_full_wo
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	
	
	
	
	
	
	
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Block via ring                                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>block by sending / receiving a message around the ring
	!>@param[in] ring_comm - comm of the cart topology
	!>@param[in] id - id of this process
	!>@param[in] world_process - id of world process
	!>@param[in] rank - rank of mpi job
	subroutine block_ring(ring_comm,id,world_process,rank)
		implicit none
		integer(i4b), intent(in) :: ring_comm, id, world_process, rank
		integer(i4b) :: error, tag1=2010
		character (len=3) :: mesg='Yo!'
		
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! essentially blocks until all processors catch up					   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call MPI_Barrier(ring_comm, error)
! 		if (id .ne. world_process ) then
! 			! processors except 0 are waiting to recv from previous pe:
! 			call MPI_Recv(mesg, len(mesg), MPI_CHARACTER, id-1, &
! 				tag1, ring_comm, MPI_STATUS_IGNORE,error)
! 		endif
! 		if ( (world_process+1) .ne. rank ) then ! so we don't send a message to ourselves!
! 			! processor 0 will send here first (as not waiting)
! 			call MPI_Send(mesg, len(mesg), MPI_CHARACTER, mod(id+1,rank), &
! 					tag1, ring_comm, error)
! 			! lastly receive message from last process
! 			if(id == world_process) then
! 				! processor 0 waiting to recv from last pe in ring
! 				call MPI_Recv(mesg, len(mesg), MPI_CHARACTER, rank-1, &
! 					tag1, ring_comm, MPI_STATUS_IGNORE,error)
! 			endif
! 		endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	end subroutine block_ring
	
	
	
	
	end module mpi_module
	
	