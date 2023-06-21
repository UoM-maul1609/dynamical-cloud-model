	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>drivers for the dynamical cloud model
    module drivers
    use numerics_type
    !use variables
    private
    public :: model_driver
    contains
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calls IO and runs one time-step of model
	!>@param[in] ntim: number of time-levels
	!>@param[in] dt: grids
	!>@param[in] ip, jp, kp: dims for full domain
	!>@param[in] ipp, jpp, kpp: dims for this block of data
	!>@param[in] l_h,r_h, ipstart,jpstart,kpstart: halo and dims
	!>@param[in] x,y,z, dx, dy, dz, dxn, dyn, dzn: grids
	!>@param[inout] ut,vt,wt, q: prognostics
	!>@param[in] rhoa, rhoan: reference variables
	!>@param[in] lamsq, lamsqn: reference variables
	!>@param[inout] lbc,ubc: boundary conditions
	!>@param[in] coords: for Cartesian topology
	!>@param[inout] new_file: flag for if this is a new file
	!>@param[in] outputfile: netcdf output
	!>@param[in] output_interval: interval for output (s)
	!>@param[in] viscous: logical for applying viscous dissipation
	!>@param[in] advection_scheme, kord, monotone: flags for advection schemes
	!>@param[in] neumann: boundary condition flag for advection schemes
	!>@param[in] dims,id, world_process, ring_comm: mpi variables
    subroutine model_driver(ntim,dt,nq,l_h,r_h, &
    			ip,jp,kp, &
    			ipp,jpp,kpp, &
				ipstart, jpstart, kpstart, &
				x,y,z, &
				dx,dy,dz, &
				dxn,dyn,dzn, &
				ut,vt,wt,&
				q, &
				rhoa,rhoan, &
				lamsq,lamsqn, &
				lbc,ubc,&
				coords, &
				new_file,outputfile, output_interval, &
				viscous, &
				advection_scheme, kord, monotone, neumann, &
				dims,id, world_process, rank, ring_comm)
		use numerics_type
		use mpi_module, only : exchange_full, exchange_along_dim
		use advection_s_3d, only : first_order_upstream_3d, mpdata_3d, adv_ref_state, &
		                    mpdata_vec_3d

		implicit none
		logical, intent(inout) :: new_file
		logical, intent(in) :: viscous, monotone
		integer(i4b), intent(in) :: ntim,nq,ip,jp,kp, ipp,jpp,kpp, &
						l_h,r_h, ipstart, jpstart, kpstart, &
						advection_scheme, kord, neumann
		integer(i4b), intent(in) :: id, world_process, ring_comm, rank
		integer(i4b), dimension(3), intent(in) :: coords, dims
		character (len=*), intent(in) :: outputfile
		real(wp), intent(in) :: output_interval, dt
		real(wp), dimension(1-l_h:ipp+r_h), intent(in) :: x,dx, dxn
		real(wp), dimension(1-l_h:jpp+r_h), intent(in) :: y,dy,dyn
		real(wp), dimension(1-l_h:kpp+r_h), intent(in) :: z,dz,dzn,&
		                                        rhoa, rhoan,lamsq, lamsqn
			
		real(wp), &
			dimension(1-r_h:kpp+r_h,1-r_h:jpp+r_h,1-l_h:ipp+r_h), target, &
			intent(inout) :: ut
		real(wp), &
			dimension(1-r_h:kpp+r_h,1-l_h:jpp+r_h,1-r_h:ipp+r_h), target, &
			intent(inout) :: vt
		real(wp), &
			dimension(1-l_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h), target, &
			intent(inout) :: wt
		real(wp), &
			dimension(1-l_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h,1:nq), target, &
			intent(inout) :: q
		real(wp), dimension(nq), intent(inout) :: lbc,ubc
		! locals:		
		integer(i4b) :: n,n2, cur=1, i,j,k, error, rank2
		real(wp) :: time, time_last_output, output_time, a
		real(wp), dimension(:,:,:), pointer :: u,zu,tu
		real(wp), dimension(:,:,:), pointer :: v,zv,tv
		real(wp), dimension(:,:,:), pointer :: w,zw,tw
! 		real(wp), &
! 			dimension(1-r_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h) :: th2

        if(id>=dims(1)*dims(2)*dims(3)) return 
		

		time_last_output=-output_interval
		output_time=output_interval
		rank2=dims(1)*dims(2)*dims(3)
		
		! associate pointers - for efficiency, when swapping arrays in leap-frog scheme
		u => ut
		v => vt
		w => wt
		

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! time-loop                                                                      !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		do n=1,ntim	
		
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! write netcdf variables                                                     !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			time=real(n-1,wp)*dt
			if (time-time_last_output >= output_interval) then
			
			
				if (id==world_process) &
					print *,'output no ',cur,' at time (s) ', &
						time,n,' steps of ',ntim

				
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! Output to NetCDF                                                       !
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				call output(new_file,outputfile,cur,nq, &
							ip,ipp,ipstart,jp,jpp,jpstart,kp,kpp,kpstart, &
							l_h,r_h, &
							time, x,y,z,rhoa, &
							u,v,w,q, &
							id, world_process, rank2, ring_comm)
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



				time_last_output=time
				cur=cur+1
			endif
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			!if(coords(3) == 0) w(1-l_h:0,:,:)=0._wp
			
						
						
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! set halos																	 !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
			do n2=1,nq
                call exchange_full(ring_comm, id, kpp, jpp, ipp, &
                                    r_h,r_h,r_h,r_h,r_h,r_h,q(:,:,:,n2),&
                                    lbc(n2),ubc(n2),dims,coords)
            enddo
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	



			

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! advect scalar fields using mid-point                                       !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			select case (advection_scheme)
				case (0)
				    do n2=1,nq
                        call first_order_upstream_3d(dt,dxn,dyn,dzn,rhoa,rhoan,&
                            ipp,jpp,kpp,l_h,r_h,u,v,w,q(:,:,:,n2),neumann,dims,coords)
                    enddo
				case (1)
				    do n2=1,nq
                        call mpdata_3d(dt,dx,dy,dz,dxn,dyn,dzn,rhoa,rhoan,&
                            ipp,jpp,kpp,l_h,r_h,u,v,w,q(:,:,:,n2),lbc(n2),ubc(n2),&
                            kord,monotone,neumann, &
                            ring_comm,id, &
                            dims,coords)
                    enddo
				case (2)
					call mpdata_vec_3d(dt,dx,dy,dz,dxn,dyn,dzn,rhoa,rhoan,&
						ipp,jpp,kpp,nq,l_h,r_h,u,v,w,q(:,:,:,:),lbc,ubc,kord, &
						monotone,neumann,ring_comm,id, dims,coords)
				case default
					print *,'not coded'
					stop
			end select
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! set halos																	 !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
			do n2=1,nq
                call exchange_full(ring_comm, id, kpp, jpp, ipp, &
                                    r_h,r_h,r_h,r_h,r_h,r_h,q(:,:,:,n2),&
                                    lbc(n2),ubc(n2),dims,coords)
            enddo
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

			

			



			

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! set halos																	 !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
			call exchange_full(ring_comm, id, kpp, jpp, ipp, r_h,r_h,r_h,r_h,l_h,r_h,u,&
								0._wp,0._wp,dims,coords)
			call exchange_full(ring_comm, id, kpp, jpp, ipp, r_h,r_h,l_h,r_h,r_h,r_h,v,&
								0._wp,0._wp,dims,coords)
			call exchange_full(ring_comm, id, kpp, jpp, ipp, l_h,r_h,r_h,r_h,r_h,r_h,w,&
								0._wp,0._wp,dims,coords)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		enddo
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        if (associated(u) ) nullify(u)
        if (associated(v) ) nullify(v)
        if (associated(w) ) nullify(w)






		
	end subroutine model_driver
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	
	
	
	


	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>outputs variables to NetCDF file using MPI
	!>@param[inout] new_file: flag if this is a new file
	!>@param[in] outputfile: outputfilename
	!>@param[in] n: time-level
	!>@param[in] nq: number of q fields
	!>@param[in] ip: number of x global grid
	!>@param[in] ipp: number of x levels on this PE
	!>@param[in] ipstart: start of i index on global grid
	!>@param[in] jp: ditto for y
	!>@param[in] jpp: ditto for y
	!>@param[in] jpstart: start of j index on global grid
	!>@param[in] kp: ditto for z
	!>@param[in] kpp: ditto for z
	!>@param[in] kpstart: start of j index on global grid
	!>@param[in] l_h,r_h: halo
	!>@param[in] time: time (s)
	!>@param[in] x,y,z, rhoa: grids
	!>@param[in] u,v,w,q: prognostic variables
	!>@param[in] id: id
	!>@param[in] world_process: world_process
	!>@param[in] rank: rank
	!>@param[in] ring_comm: ring_comm
	subroutine output(new_file,outputfile,n,nq,ip,ipp,ipstart,jp,jpp,jpstart, &
					kp,kpp,kpstart,l_h,r_h, &
					time, &
					x,y,z,rhoa, &
					u,v,w,q, &
				    id, world_process, rank, ring_comm)
	
		use netcdf
		use mpi
		use mpi_module, only : mpi_integer9
		!use variables, only : MPI_INTEGER9

		implicit none
		logical, intent(inout) :: new_file
		character (len=*), intent(in) :: outputfile
		integer(i4b), intent(in) :: n, nq,ip, ipp, ipstart, jp, jpp, jpstart, &
									kp, kpp, kpstart, l_h,r_h
		real(wp), intent(in) :: time
		real(wp), dimension(1-l_h:ipp+r_h), intent(in) :: x
		real(wp), dimension(1-l_h:jpp+r_h), intent(in) :: y
		real(wp), dimension(1-l_h:kpp+r_h), intent(in) :: z,rhoa
		real(wp), &
			dimension(1-r_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h,1:nq), &
			intent(inout) :: q
		real(wp), &
			dimension(1-r_h:kpp+r_h,1-r_h:jpp+r_h,1-l_h:ipp+r_h), &
			intent(inout) :: u
		real(wp), &
			dimension(1-r_h:kpp+r_h,1-l_h:jpp+r_h,1-r_h:ipp+r_h), &
			intent(inout) :: v
		real(wp), &
			dimension(1-l_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h), &
			intent(inout) :: w
		
		integer(i4b), intent(in) :: id ,world_process, rank, ring_comm
	
		integer(i4b) :: ncid, x_dimid, nx_dimid, ny_dimid, nz_dimid, nq_dimid,&
						error, varid,a_dimid, id_go
		integer(i4b) :: i, tag1
		logical :: var


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! perform a blocking recv to wait for message from main process, 				 !
		! before carrying on															 !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if(id .ne. world_process) then
			tag1=id
			call MPI_Recv(var,1, MPI_LOGICAL, world_process, &
				tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE,error)
		endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	
		if((id==world_process) .and. new_file) then
			! open the file
		
			call check( nf90_create(outputfile, NF90_CLOBBER, ncid) )

			! define dimensions (netcdf hands back a handle)
			call check( nf90_def_dim(ncid, "times", NF90_UNLIMITED, x_dimid) )
			call check( nf90_def_dim(ncid, "ip", ip, nx_dimid) )
			call check( nf90_def_dim(ncid, "jp", jp, ny_dimid) )
			call check( nf90_def_dim(ncid, "kp", kp, nz_dimid) )
			call check( nf90_def_dim(ncid, "nq", nq, nq_dimid) )


			! close the file, freeing up any internal netCDF resources
			! associated with the file, and flush any buffers
			call check( nf90_close(ncid) )
		
			! now define some variables, units, etc
			call check( nf90_open(outputfile, NF90_WRITE, ncid) )
			
			
			! define mode
			call check( nf90_redef(ncid) )

			! define variable: time
			call check( nf90_def_var(ncid, "time", nf90_float, &
						(/x_dimid/), varid) )
			! get id to a_dimid
			call check( nf90_inq_varid(ncid, "time", a_dimid) )
			! units
			call check( nf90_put_att(ncid, a_dimid, &
					   "units", "seconds") )
						
						
			! define variable: x
			call check( nf90_def_var(ncid, "x", nf90_float, &
						(/nx_dimid/), varid) )
			! get id to a_dimid
			call check( nf90_inq_varid(ncid, "x", a_dimid) )
			! units
			call check( nf90_put_att(ncid, a_dimid, &
					   "units", "m") )

			! define variable: y
			call check( nf90_def_var(ncid, "y", nf90_float, &
						(/ny_dimid/), varid) )
			! get id to a_dimid
			call check( nf90_inq_varid(ncid, "y", a_dimid) )
			! units
			call check( nf90_put_att(ncid, a_dimid, &
					   "units", "m") )

			! define variable: z
			call check( nf90_def_var(ncid, "z", nf90_float, &
						(/nz_dimid/), varid) )
			! get id to a_dimid
			call check( nf90_inq_varid(ncid, "z", a_dimid) )
			! units
			call check( nf90_put_att(ncid, a_dimid, &
					   "units", "m") )

			! define variable: rhoa
			call check( nf90_def_var(ncid, "rhoa", nf90_float, &
						(/nz_dimid/), varid) )
			! get id to a_dimid
			call check( nf90_inq_varid(ncid, "rhoa", a_dimid) )
			! units
			call check( nf90_put_att(ncid, a_dimid, &
					   "units", "kg/m3") )


			! define variable: q
			call check( nf90_def_var(ncid, "q", nf90_float, &
						(/nz_dimid,ny_dimid,nx_dimid,nq_dimid,x_dimid/), varid) )
			! get id to a_dimid
			call check( nf90_inq_varid(ncid, "q", a_dimid) )
			! units
			call check( nf90_put_att(ncid, a_dimid, &
					   "units", "/kg") )

			! define variable: u
			call check( nf90_def_var(ncid, "u", nf90_float, &
						(/nz_dimid,ny_dimid,nx_dimid,x_dimid/), varid) )
			! get id to a_dimid
			call check( nf90_inq_varid(ncid, "u", a_dimid) )
			! units
			call check( nf90_put_att(ncid, a_dimid, &
					   "units", "m/s") )

			! define variable: v
			call check( nf90_def_var(ncid, "v", nf90_float, &
						(/nz_dimid,ny_dimid,nx_dimid,x_dimid/), varid) )
			! get id to a_dimid
			call check( nf90_inq_varid(ncid, "v", a_dimid) )
			! units
			call check( nf90_put_att(ncid, a_dimid, &
					   "units", "m/s") )

			! define variable: w
			call check( nf90_def_var(ncid, "w", nf90_float, &
						(/nz_dimid,ny_dimid,nx_dimid,x_dimid/), varid) )
			! get id to a_dimid
			call check( nf90_inq_varid(ncid, "w", a_dimid) )
			! units
			call check( nf90_put_att(ncid, a_dimid, &
					   "units", "m/s") )




			! exit define mode
			call check( nf90_enddef(ncid) )
			
			
			call check( nf90_close(ncid) )

			new_file=.false.
		endif
	
! 	



		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! now send messages from the main process to all other processes                 !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if(id == world_process) then
			do i=1,rank-1
				tag1=i
				call MPI_Send(var, 1, MPI_LOGICAL, i, &
						tag1, MPI_COMM_WORLD, error)
			enddo
		endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	
	

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! perform a blocking recv to wait for message from main process,                 !
		! before carrying on                               								 !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if(id .ne. world_process) then
			tag1=id
			call MPI_Recv(id_go,1, MPI_INTEGER9, id-1, &
				tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE,error)
		else
			id_go=world_process ! lets us go for first run
		endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		
		
		
		

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! ****WRITE****																	 !			
		! now we can write to file - each PE writes its own segment						 !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		call check( nf90_open(outputfile, NF90_WRITE, ncid) )
		
		if(n == 1) then
			! write variable: x
			call check( nf90_inq_varid(ncid, "x", varid ) )
			call check( nf90_put_var(ncid, varid, x(1:ipp), &
						start = (/ipstart/)))	

			! write variable: y
			call check( nf90_inq_varid(ncid, "y", varid ) )
			call check( nf90_put_var(ncid, varid, y(1:jpp), &
						start = (/jpstart/)))	
						
			! write variable: z
			call check( nf90_inq_varid(ncid, "z", varid ) )
			call check( nf90_put_var(ncid, varid, z(1:kpp), &
						start = (/kpstart/)))	
			! write variable: rhoa
			call check( nf90_inq_varid(ncid, "rhoa", varid ) )
			call check( nf90_put_var(ncid, varid, rhoa(1:kpp), &
						start = (/kpstart/)))	
		endif


		if(id==world_process) then
			! write variable: time
			call check( nf90_inq_varid(ncid, "time", varid ) )
			call check( nf90_put_var(ncid, varid, time, &
						start = (/n/)))	
	    endif
	    
	    
		! write variable: q
		call check( nf90_inq_varid(ncid, "q", varid ) )
		call check( nf90_put_var(ncid, varid, q(1:kpp,1:jpp,1:ipp,1:nq), &
					start = (/kpstart,jpstart,ipstart,1,n/)))	

		! write variable: u
		call check( nf90_inq_varid(ncid, "u", varid ) )
		call check( nf90_put_var(ncid, varid, u(1:kpp,1:jpp,1:ipp), &
					start = (/kpstart,jpstart,ipstart,n/)))	

		! write variable: v
		call check( nf90_inq_varid(ncid, "v", varid ) )
		call check( nf90_put_var(ncid, varid, v(1:kpp,1:jpp,1:ipp), &
					start = (/kpstart,jpstart,ipstart,n/)))	

		! write variable: w
		call check( nf90_inq_varid(ncid, "w", varid ) )
		call check( nf90_put_var(ncid, varid, w(1:kpp,1:jpp,1:ipp), &
					start = (/kpstart,jpstart,ipstart,n/)))	


		call check( nf90_close(ncid) )
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




		

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! perform a send, to essentially allow next PE to resume and start the write     !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if((id == id_go)) then
			tag1=id+1
			if( ((id+1).lt.rank) ) then
				call MPI_Send(id+1, 1, MPI_INTEGER9, id+1, &
							tag1, MPI_COMM_WORLD, error)
			endif
			
							
		endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



		if (rank > 1 ) then
			! send to world_process to complete ring
			tag1=2010
			if( ((id+1).eq.rank) ) then
				call MPI_Send(var, 1, MPI_LOGICAL, world_process, &
							tag1, MPI_COMM_WORLD, error)			
			endif


			! receive at world_process to complete ring
			if((id==world_process) ) then
				call MPI_Recv(var,1, MPI_LOGICAL,rank-1, &
					tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE,error)
			endif
		endif	


	end subroutine output
	
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
	
	end module drivers
