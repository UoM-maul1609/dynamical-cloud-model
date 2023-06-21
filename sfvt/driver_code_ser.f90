	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>drivers for the dynamical cloud model
    module drivers_ser
    use numerics_type
    !use variables
    private
    public :: model_driver, model_driver_2d
    contains
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calls IO and runs one time-step of model
	!>@param[in] ntim: number of time-levels
	!>@param[in] dt: grids
	!>@param[in] ip, kp: dims for full domain
	!>@param[in] ipp, kpp: dims for this block of data
	!>@param[in] l_h,r_h, ipstart,jpstart,kpstart: halo and dims
	!>@param[in] x,z, dx, dz, dxn, dzn: grids
	!>@param[inout] ut,wt, q: prognostics
	!>@param[in] rhoa, rhoan: reference variables
	!>@param[in] lamsq, lamsqn: reference variables
	!>@param[inout] new_file: flag for if this is a new file
	!>@param[in] outputfile: netcdf output
	!>@param[in] output_interval: interval for output (s)
	!>@param[in] viscous: logical for applying viscous dissipation
	!>@param[in] advection_scheme, kord, monotone: flags for advection schemes
	!>@param[in] neumann: boundary condition flag for advection schemes
    subroutine model_driver_2d(ntim,dt,nq,l_h,r_h, &
    			ip,kp, &
    			ipp,kpp, &
				ipstart, kpstart, &
				x,z, &
				dx,dz, &
				dxn,dzn, &
				ut,wt,&
				q, &
				rhoa,rhoan, &
				lamsq,lamsqn, &
				new_file,outputfile, output_interval, &
				viscous, &
				advection_scheme, kord, monotone,neumann)
		use numerics_type
		use advection_s_2d, only : first_order_upstream_2d, mpdata_2d, &
		                    mpdata_vec_2d

		implicit none
		logical, intent(inout) :: new_file
		logical, intent(in) :: viscous, monotone
		integer(i4b), intent(in) :: ntim,nq,ip,kp, ipp,kpp, &
						l_h,r_h, ipstart, kpstart, &
						advection_scheme, kord, neumann
		character (len=*), intent(in) :: outputfile
		real(wp), intent(in) :: output_interval, dt
		real(wp), dimension(1-l_h:ipp+r_h), intent(in) :: x,dx, dxn
		real(wp), dimension(1-l_h:kpp+r_h), intent(in) :: z,dz,dzn,&
		                                        rhoa, rhoan,lamsq, lamsqn
			
		real(wp), &
			dimension(1-r_h:kpp+r_h,1-l_h:ipp+r_h), target, &
			intent(inout) :: ut
		real(wp), &
			dimension(1-l_h:kpp+r_h,1-r_h:ipp+r_h), target, &
			intent(inout) :: wt
		real(wp), &
			dimension(1-l_h:kpp+r_h,1-r_h:ipp+r_h,1:nq), target, &
			intent(inout) :: q
					
		! locals:		
		integer(i4b) :: n,n2, cur=1, i,j,k, error, rank2
		real(wp) :: time, time_last_output, output_time, a
		real(wp), dimension(:,:), pointer :: u,zu,tu
		real(wp), dimension(:,:), pointer :: w,zw,tw

		

		time_last_output=-output_interval
		output_time=output_interval
		
		! associate pointers - for efficiency, when swapping arrays in leap-frog scheme
		u => ut
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
			
			
                print *,'output no ',cur,' at time (s) ', &
                    time,n,' steps of ',ntim

				
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! Output to NetCDF                                                       !
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				call output_2d(new_file,outputfile,cur,nq, &
							ip,ipp,ipstart,kp,kpp,kpstart, &
							l_h,r_h, &
							time, x,z,rhoa, &
							u,w,q)
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



				time_last_output=time
				cur=cur+1
			endif
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! 			if(coords(3) == 0) w(1-l_h:0,:,:)=0._wp
			
						
						
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! set halos																	 !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
			do n2=1,nq
                q(:,-l_h+1:0,n2)=q(:,ip-l_h+1:ip,n2)
                q(:,ip+1:ip+r_h,n2)=q(:,1:r_h,n2)
            enddo
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	



			

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! advect scalar fields using mid-point                                       !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			select case (advection_scheme)
				case (0)
				    do n2=1,nq
                        call first_order_upstream_2d(dt,dxn,dzn,rhoa,rhoan, &
                            ipp,kpp,l_h,r_h,u,w,q(:,:,n2),neumann)
                    enddo
				case (1)
				    do n2=1,nq
                        call mpdata_2d(dt,dx,dz,dxn,dzn,rhoa,rhoan, &
                            ipp,kpp,l_h,r_h,u,w,q(:,:,n2),kord,monotone,neumann)
                    enddo
				case (2)
					call mpdata_vec_2d(dt,dx,dz,dxn,dzn,rhoa,rhoan, &
						ipp,kpp,nq,l_h,r_h,u,w,q(:,:,:),kord, &
						monotone,neumann)
				case default
					print *,'not coded'
					stop
			end select
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! set halos																	 !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
			do n2=1,nq
                q(:,-l_h+1:0,n2)=q(:,ip-l_h+1:ip,n2)
                q(:,ip+1:ip+r_h,n2)=q(:,1:r_h,n2)
            enddo
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

			

			



			

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! set halos																	 !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
            u(:,-l_h+1:0)    =u(:,ip-l_h+1:ip)
            u(:,ip+1:ip+r_h) =u(:,1:r_h)
            w(:,-l_h+1:0)    =w(:,ip-l_h+1:ip)
            w(:,ip+1:ip+r_h) =w(:,1:r_h)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		enddo
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        if (associated(u) ) nullify(u)
        if (associated(w) ) nullify(w)






		
	end subroutine model_driver_2d
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	
	


	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calls IO and runs one time-step of model
	!>@param[in] ntim: number of time-levels
	!>@param[in] dt: grids
	!>@param[in] kp: dims for full domain
	!>@param[in] kpp: dims for this block of data
	!>@param[in] l_h,r_h, kpstart: halo and dims
	!>@param[in] z, dz, dzn: grids
	!>@param[inout] wt, q: prognostics
	!>@param[in] rhoa, rhoan: reference variables
	!>@param[in] lamsq, lamsqn: reference variables
	!>@param[inout] new_file: flag for if this is a new file
	!>@param[in] outputfile: netcdf output
	!>@param[in] output_interval: interval for output (s)
	!>@param[in] viscous: logical for applying viscous dissipation
	!>@param[in] advection_scheme, kord, monotone: flags for advection schemes
	!>@param[in] neumann: boundary condition flag for advection schemes
    subroutine model_driver(ntim,dt,nq,l_h,r_h, &
    			kp, &
    			kpp, &
				kpstart, &
				z, &
				dz, &
				dzn, &
				wt,&
				q, &
				rhoa,rhoan, &
				lamsq,lamsqn, &
				new_file,outputfile, output_interval, &
				viscous, &
				advection_scheme, kord, monotone, neumann)
		use numerics_type
		use advection_s_1d, only : first_order_upstream_1d, mpdata_1d, &
		                    mpdata_vec_1d

		implicit none
		logical, intent(inout) :: new_file
		logical, intent(in) :: viscous, monotone
		integer(i4b), intent(in) :: ntim,nq,kp, kpp, &
						l_h,r_h, kpstart, &
						advection_scheme, kord, neumann
		character (len=*), intent(in) :: outputfile
		real(wp), intent(in) :: output_interval, dt
		real(wp), dimension(1-l_h:kpp+r_h), intent(in) :: z,dz,dzn,&
		                                        rhoa, rhoan,lamsq, lamsqn
			
		real(wp), &
			dimension(1-l_h:kpp+r_h), target, &
			intent(inout) :: wt
		real(wp), &
			dimension(1-l_h:kpp+r_h,1:nq), target, &
			intent(inout) :: q
					
		! locals:		
		integer(i4b) :: n,n2, cur=1, i,j,k, error, rank2
		real(wp) :: time, time_last_output, output_time, a
		real(wp), dimension(:), pointer :: w,zw,tw
		

		time_last_output=-output_interval
		output_time=output_interval
		
		! associate pointers - for efficiency, when swapping arrays in leap-frog scheme
		w => wt
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! time-loop                                                                      !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		do n=1,ntim	
		
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! write netcdf variables                                                     !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			time=real(n-1,wp)*dt
            w=sin(2._wp*PI/1200._wp*time)
			if (time-time_last_output >= output_interval) then
			
			
                print *,'output no ',cur,' at time (s) ', &
                    time,n,' steps of ',ntim

				
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! Output to NetCDF                                                       !
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				call output(new_file,outputfile,cur,nq, &
							kp,kpp,kpstart, &
							l_h,r_h, &
							time, z,rhoa, &
							w,q)
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



				time_last_output=time
				cur=cur+1
			endif
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			!w(1-l_h:0)=0._wp
			
						
						



			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! advect scalar fields using mid-point                                       !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			select case (advection_scheme)
				case (0)
				    do n2=1,nq
                        call first_order_upstream_1d(dt,dzn,rhoa,rhoan, &
                            kpp,l_h,r_h,w,q(:,n2), neumann)
                    enddo
				case (1)
				    do n2=1,nq
                        call mpdata_1d(dt,dz,dzn,rhoa,rhoan, &
                            kpp,l_h,r_h,w,q(:,n2),kord,monotone,neumann,0)
                    enddo
				case (2)
					call mpdata_vec_1d(dt,dz,dzn,rhoa,rhoan, &
						kpp,nq,l_h,r_h,w,q(:,1:nq),kord, &
						monotone,neumann)
				case default
					print *,'not coded'
					stop
 			end select
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			


			

					


		enddo
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


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
	!>@param[in] kp: ditto for z
	!>@param[in] kpp: ditto for z
	!>@param[in] kpstart: start of j index on global grid
	!>@param[in] l_h,r_h: halo
	!>@param[in] time: time (s)
	!>@param[in] x,z, rhoa: grids
	!>@param[in] u,w,q: prognostic variables
	subroutine output_2d(new_file,outputfile,n,nq,ip,ipp,ipstart, &
					kp,kpp,kpstart,l_h,r_h, &
					time, &
					x,z,rhoa, &
					u,w,q)
	
		use netcdf
		use mpi
		use mpi_module, only : mpi_integer9
		!use variables, only : MPI_INTEGER9

		implicit none
		logical, intent(inout) :: new_file
		character (len=*), intent(in) :: outputfile
		integer(i4b), intent(in) :: n, nq,ip, ipp, ipstart, &
									kp, kpp, kpstart, l_h,r_h
		real(wp), intent(in) :: time
		real(wp), dimension(1-l_h:ipp+r_h), intent(in) :: x
		real(wp), dimension(1-l_h:kpp+r_h), intent(in) :: z,rhoa
		real(wp), &
			dimension(1-r_h:kpp+r_h,1-r_h:ipp+r_h,1:nq), &
			intent(inout) :: q
		real(wp), &
			dimension(1-r_h:kpp+r_h,1-l_h:ipp+r_h), &
			intent(inout) :: u
		real(wp), &
			dimension(1-l_h:kpp+r_h,1-r_h:ipp+r_h), &
			intent(inout) :: w
		
	
		integer(i4b) :: ncid, x_dimid, nx_dimid, ny_dimid, nz_dimid, nq_dimid,&
						error, varid,a_dimid, id_go
		integer(i4b) :: i, tag1
		logical :: var


	
		if(new_file) then
			! open the file
		
			call check( nf90_create(outputfile, NF90_CLOBBER, ncid) )

			! define dimensions (netcdf hands back a handle)
			call check( nf90_def_dim(ncid, "times", NF90_UNLIMITED, x_dimid) )
			call check( nf90_def_dim(ncid, "ip", ip, nx_dimid) )
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
						(/nz_dimid,nx_dimid,nq_dimid,x_dimid/), varid) )
			! get id to a_dimid
			call check( nf90_inq_varid(ncid, "q", a_dimid) )
			! units
			call check( nf90_put_att(ncid, a_dimid, &
					   "units", "/kg") )

			! define variable: u
			call check( nf90_def_var(ncid, "u", nf90_float, &
						(/nz_dimid,nx_dimid,x_dimid/), varid) )
			! get id to a_dimid
			call check( nf90_inq_varid(ncid, "u", a_dimid) )
			! units
			call check( nf90_put_att(ncid, a_dimid, &
					   "units", "m/s") )

			! define variable: w
			call check( nf90_def_var(ncid, "w", nf90_float, &
						(/nz_dimid,nx_dimid,x_dimid/), varid) )
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
		! ****WRITE****																	 !			
		! now we can write to file - each PE writes its own segment						 !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		call check( nf90_open(outputfile, NF90_WRITE, ncid) )
		
		if(n == 1) then
			! write variable: x
			call check( nf90_inq_varid(ncid, "x", varid ) )
			call check( nf90_put_var(ncid, varid, x(1:ipp), &
						start = (/ipstart/)))	

			! write variable: z
			call check( nf90_inq_varid(ncid, "z", varid ) )
			call check( nf90_put_var(ncid, varid, z(1:kpp), &
						start = (/kpstart/)))	
			! write variable: rhoa
			call check( nf90_inq_varid(ncid, "rhoa", varid ) )
			call check( nf90_put_var(ncid, varid, rhoa(1:kpp), &
						start = (/kpstart/)))	
		endif


        ! write variable: time
        call check( nf90_inq_varid(ncid, "time", varid ) )
        call check( nf90_put_var(ncid, varid, time, &
                    start = (/n/)))	
	    
	    
		! write variable: q
		call check( nf90_inq_varid(ncid, "q", varid ) )
		call check( nf90_put_var(ncid, varid, q(1:kpp,1:ipp,1:nq), &
					start = (/kpstart,ipstart,1,n/)))	

		! write variable: u
		call check( nf90_inq_varid(ncid, "u", varid ) )
		call check( nf90_put_var(ncid, varid, u(1:kpp,1:ipp), &
					start = (/kpstart,ipstart,n/)))	

		! write variable: w
		call check( nf90_inq_varid(ncid, "w", varid ) )
		call check( nf90_put_var(ncid, varid, w(1:kpp,1:ipp), &
					start = (/kpstart,ipstart,n/)))	


		call check( nf90_close(ncid) )
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




		



	end subroutine output_2d
	
	


	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>outputs variables to NetCDF file 
	!>@param[inout] new_file: flag if this is a new file
	!>@param[in] outputfile: outputfilename
	!>@param[in] n: time-level
	!>@param[in] nq: number of q fields
	!>@param[in] kp: ditto for z
	!>@param[in] kpp: ditto for z
	!>@param[in] kpstart: start of j index on global grid
	!>@param[in] l_h,r_h: halo
	!>@param[in] time: time (s)
	!>@param[in] z, rhoa: grids
	!>@param[in] w,q: prognostic variables
	subroutine output(new_file,outputfile,n,nq, &
					kp,kpp,kpstart,l_h,r_h, &
					time, &
					z,rhoa, &
					w,q)
	
		use netcdf

		implicit none
		logical, intent(inout) :: new_file
		character (len=*), intent(in) :: outputfile
		integer(i4b), intent(in) :: n, nq,kp, kpp, kpstart, l_h,r_h
		real(wp), intent(in) :: time
		real(wp), dimension(1-l_h:kpp+r_h), intent(in) :: z,rhoa
		real(wp), &
			dimension(1-r_h:kpp+r_h,1:nq), &
			intent(inout) :: q
		real(wp), &
			dimension(1-l_h:kpp+r_h), &
			intent(inout) :: w
		
	
		integer(i4b) :: ncid, x_dimid, nx_dimid, ny_dimid, nz_dimid, nq_dimid,&
						error, varid,a_dimid, id_go
		integer(i4b) :: i, tag1
		logical :: var




		if(new_file) then
			! open the file
		
			call check( nf90_create(outputfile, NF90_CLOBBER, ncid) )

			! define dimensions (netcdf hands back a handle)
			call check( nf90_def_dim(ncid, "times", NF90_UNLIMITED, x_dimid) )
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
						(/nz_dimid,nq_dimid,x_dimid/), varid) )
			! get id to a_dimid
			call check( nf90_inq_varid(ncid, "q", a_dimid) )
			! units
			call check( nf90_put_att(ncid, a_dimid, &
					   "units", "/kg") )


			! define variable: w
			call check( nf90_def_var(ncid, "w", nf90_float, &
						(/nz_dimid,x_dimid/), varid) )
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
	


	
	

		
		
		

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! ****WRITE****																	 !			
		! now we can write to file - each PE writes its own segment						 !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		call check( nf90_open(outputfile, NF90_WRITE, ncid) )
		
		if(n == 1) then
			! write variable: z
			call check( nf90_inq_varid(ncid, "z", varid ) )
			call check( nf90_put_var(ncid, varid, z(1:kpp), &
						start = (/kpstart/)))	
			! write variable: rhoa
			call check( nf90_inq_varid(ncid, "rhoa", varid ) )
			call check( nf90_put_var(ncid, varid, rhoa(1:kpp), &
						start = (/kpstart/)))	
		endif


        ! write variable: time
        call check( nf90_inq_varid(ncid, "time", varid ) )
        call check( nf90_put_var(ncid, varid, time, &
                    start = (/n/)))	
	    
		! write variable: q
		call check( nf90_inq_varid(ncid, "q", varid ) )
		call check( nf90_put_var(ncid, varid, q(1:kpp,1:nq), &
					start = (/kpstart,1,n/)))	


		! write variable: w
		call check( nf90_inq_varid(ncid, "w", varid ) )
		call check( nf90_put_var(ncid, varid, w(1:kpp), &
					start = (/kpstart,n/)))	


		call check( nf90_close(ncid) )
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




		


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
			!stop "Stopped"
		end if
	end subroutine check
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	end module drivers_ser
