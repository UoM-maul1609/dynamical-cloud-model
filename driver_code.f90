	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>drivers for the dynamical cloud model
    module drivers
    use nrtype
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
	!>@param[in] x,y,z, xn,ym,zn,dx, dy, dz, dxn, dyn, dzn: grids
	!>@param[in] dampfacn,dampfac
	!>@param[inout] ubar,vbar,wbar,thbar,qbar
	!>@param[inout] ut,vt,wt: prognostics
	!>@param[inout] zut,zvt,zwt: prognostics - previous time-step
	!>@param[inout] tut,tvt,twt: prognostics - temp storage
	!>@param[inout] th,sth,p: prognostic variables
	!>@param[inout] su,sv,sw,psrc, div: more prognostic variables
	!>@param[inout] strain,vism,vist - viscosity subgrid 
	!>@param[in] z0,z0th - roughness lengths 
	!>@param[in] ptol - tolerance for pressure solver
	!>@param[inout] q,sq,viss - for clouds and subgrid
	!>@param[in] theta, thetan: reference variables
	!>@param[in] rhoa, rhoan: reference variables
	!>@param[in] lamsq, lamsqn: reference variables
	!>@param[inout] thbase, thtop: bottom and base th
	!>@param[in] coords: for Cartesian topology
	!>@param[inout] new_file: flag for if this is a new file
	!>@param[in] outputfile: netcdf output
	!>@param[in] output_interval: interval for output (s)
	!>@param[in] viscous: logical for applying viscous dissipation
	!>@param[in] advection_scheme, kord, monotone: flags for advection schemes
	!>@param[in] moisture: flag for moisture
	!>@param[in] damping_layer: flag for damping layer
	!>@param[in] nq: number of q-variables
	!>@param[in] dims,id, world_process, ring_comm, sub_horiz_comm: mpi variables
	!>@param[in] sub_vert_comm: mpi variables
    subroutine model_driver(ntim,dt,l_h,r_h, &
    			ip,jp,kp, &
    			ipp,jpp,kpp, &
				ipstart, jpstart, kpstart, &
				x,y,z, xn,yn,zn, &
				dx,dy,dz, &
				dxn,dyn,dzn, &
				ubar,vbar,wbar,thbar,qbar,dampfacn,dampfac, &
				ut,vt,wt,&
				zut,zvt,zwt,&
				tut,tvt,twt,&
				th,sth,p, &
				su,sv,sw,psrc, &
				div, &
				strain, vism, vist, z0,z0th, ptol, &
				q,sq,viss, &
				theta,thetan, &
				rhoa,rhoan, &
				lamsq,lamsqn, &
				thbase,thtop, &
				coords, &
				new_file,outputfile, output_interval, &
				viscous, &
				advection_scheme, kord, monotone, &
				moisture, damping_layer,nq, &
				dims,id, world_process, rank, ring_comm,sub_horiz_comm,sub_vert_comm)
		use nrtype
		use mpi_module, only : exchange_full, exchange_along_dim
		!use advection_3d, only : first_order_upstream_3d, mpdata_3d, adv_ref_state
        use advection_s_3d, only : first_order_upstream_3d, &
                    mpdata_3d, mpdata_vec_3d, adv_ref_state, mpdata_3d_add, &
                    mpdata_vert_3d
		use d_solver, only : bicgstab, sources, advance_momentum
		use subgrid_3d, only : advance_scalar_fields_3d
        use diagnostics

		implicit none
		logical, intent(inout) :: new_file
		logical, intent(in) :: viscous, monotone, moisture, damping_layer
		integer(i4b), intent(in) :: ntim,ip,jp,kp, ipp,jpp,kpp, &
						l_h,r_h, ipstart, jpstart, kpstart, &
						advection_scheme, kord, nq
		integer(i4b), intent(in) :: id, world_process, ring_comm, sub_horiz_comm, &
		    sub_vert_comm,rank
		integer(i4b), dimension(3), intent(in) :: coords, dims
		character (len=*), intent(in) :: outputfile
		real(sp), intent(in) :: output_interval, dt, z0,z0th, ptol
		real(sp), intent(inout) :: thbase, thtop
		real(sp), dimension(1-l_h:ipp+r_h), intent(in) :: x,xn,dx, dxn
		real(sp), dimension(1-l_h:jpp+r_h), intent(in) :: y,yn,dy,dyn
		real(sp), dimension(1-l_h:kpp+r_h), intent(in) :: z,zn,dz,dzn, theta, thetan, &
														rhoa, rhoan,lamsq, lamsqn
		real(sp), dimension(1-l_h:kpp+r_h), intent(inout) :: ubar,vbar,wbar,thbar,qbar
		real(sp), dimension(1-l_h:kpp+r_h), intent(in) :: dampfacn,dampfac
		real(sp), &
			dimension(1-r_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h), &
			intent(inout) :: th,p,su,sv,sw,psrc,div
			
		real(sp), &
			dimension(1-r_h:kpp+r_h,1-r_h:jpp+r_h,1-l_h:ipp+r_h), target, &
			intent(inout) :: ut,zut,tut
		real(sp), &
			dimension(1-r_h:kpp+r_h,1-l_h:jpp+r_h,1-r_h:ipp+r_h), target, &
			intent(inout) :: vt,zvt,tvt
		real(sp), &
			dimension(1-l_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h), target, &
			intent(inout) :: wt,zwt,twt

		real(sp), &
			dimension(1-l_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h), &
			intent(inout) ::sth,strain,vism,vist

		real(sp), &
			dimension(1-l_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h,1:nq), &
			intent(inout) :: q,sq,viss
					
		! locals:		
		integer(i4b) :: n,n2, cur=1, i,j,k, error, rank2
		real(sp) :: time, time_last_output, output_time, a,t1=0._sp,t2=0._sp
		real(sp), dimension(:,:,:), pointer :: u,zu,tu
		real(sp), dimension(:,:,:), pointer :: v,zv,tv
		real(sp), dimension(:,:,:), pointer :: w,zw,tw
! 		real(sp), &
! 			dimension(1-r_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h) :: th2
		

		time_last_output=-output_interval
		output_time=output_interval
		rank2=dims(1)*dims(2)*dims(3)
		
		! associate pointers - for efficiency, when swapping arrays in leap-frog scheme
		u => ut   ! current time-step
		v => vt
		w => wt
		zu => zut ! previous time-step
		zv => zvt
		zw => zwt
		tu => tut ! next-time-step
		tv => tvt
		tw => twt
		

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! time-loop                                                                      !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		do n=1,ntim	
		
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! write netcdf variables                                                     !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			time=real(n-1,sp)*dt
			if (time-time_last_output >= output_interval) then
			
			
				if (id==world_process) &
					print *,'output no ',cur,' at time (s) ', &
						time,n,' steps of ',ntim

								
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! Output to NetCDF                                                       !
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				call output(new_file,outputfile,cur, &
							ip,ipp,ipstart,jp,jpp,jpstart,kp,kpp,kpstart, &
							l_h,r_h, &
							time, x,y,z,rhoa, thetan, &
							u,v,w,th,p,div, &
							id, world_process, rank2, ring_comm)
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



				time_last_output=time
				cur=cur+1
			endif
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			
			
			
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! calculate horizontal averages                                              !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			call horizontal_means(sub_horiz_comm,id,dims,coords, moisture, &
	                ip,jp,ipp,jpp,kpp,l_h,r_h, nq,&
	                ubar,vbar,wbar,thbar,qbar,u,v,w,th,q)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! calculate divergence - test                                                !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			call divergence_calc(ring_comm,id,dims,coords, &
	            ipp,jpp,kpp,dx,dxn,dy,dyn,dz,dzn,rhoa,rhoan,l_h,r_h,u,v,w,div)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
				
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! calculate sources of momentum and pressure                                 !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			call sources(ring_comm, id, rank2, dims,coords, &
				dt,x,y,z,zn,dx,dy,dz,dxn,dyn,dzn,ipp,jpp,kpp,l_h,r_h,&
				nq, &
				ubar,vbar,wbar,thbar,qbar, dampfacn,dampfac, &
				zu,zv,zw, &
				u,v,w,su,sv,sw, &
				q,sq,viss, &
				psrc, &
				th,sth, strain,vism,vist, &
				theta,thetan,rhoa,rhoan,lamsq,lamsqn, &
				z0,z0th, &
				viscous, moisture,damping_layer)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! set halos																	 !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
			call exchange_along_dim(ring_comm, id, kpp, jpp, ipp, &
								r_h,r_h,r_h,r_h,r_h,r_h, psrc,0._sp, 0._sp, &
								dims,coords)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! find pressure perturbation                                                 !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			call bicgstab(ring_comm, id, rank2,dims,coords, &
			 dt,x,y,z,dx,dy,dz,dxn,dyn,dzn, &
			 rhoa,rhoan, &
			 ipp,jpp,kpp,l_h,r_h,su,sv,sw,zu,zv,zw,p,psrc,ptol, &
			.false.)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! set halos																	 !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
			call exchange_along_dim(ring_comm, id, kpp, jpp, ipp, &
								r_h,r_h,r_h,r_h,r_h,r_h, p,0._sp, 0._sp, dims,coords)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
            
!             if(coords(3)==(dims(3)-1)) p(kpp+1,:,:)=p(kpp,:,:)
            
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! advect scalar fields using mid-point                                       !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			select case (advection_scheme)
				case (0)
					call first_order_upstream_3d(dt,dxn,dyn,dzn,rhoa,rhoan, &
						ipp,jpp,kpp,l_h,r_h,u,v,w,th,.true.,dims,coords)
				case (1)
					call mpdata_3d(dt,dx,dy,dz,dxn,dyn,dzn,rhoa,rhoan, &
						ipp,jpp,kpp,l_h,r_h,u,v,w,th,t1,t2, &
						kord,monotone,.true.,ring_comm,id, &
						dims,coords)						
				case(2)
					call mpdata_3d_add(dt,dx,dy,dz,dxn,dyn,dzn,rhoa,rhoan, &
						ipp,jpp,kpp,l_h,r_h,u,v,w,th,thetan,thbase,thtop, &
						kord,monotone,.true.,ring_comm,id, &
						dims,coords)
				case default
					print *,'not coded'
					stop
			end select
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			


			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! set halos																	 !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
			call exchange_full(ring_comm, id, kpp, jpp, ipp, &
								r_h,r_h,r_h,r_h,r_h,r_h,th,t1,t2,dims,coords)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

			
            if((advection_scheme == 0) .or. (advection_scheme == 1)) then
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! advect the reference state	     									 !
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
                call adv_ref_state(dt,dx,dy,dz,dxn,dyn,dzn,rhoa,rhoan,ipp,jpp,kpp,l_h,r_h, &
                            u,v,w,th,thetan,ring_comm,id,dims,coords)
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! set halos																 !
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
                call exchange_full(ring_comm, id, kpp, jpp, ipp, &
                                    r_h,r_h,r_h,r_h,r_h,r_h,th,0._sp,0._sp,dims,coords)
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
            endif


			



			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! advance momentum 1 time-step                                               !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			call advance_momentum(ring_comm, id, rank2,&
				2._sp*dt,dx,dy,dz,dxn,dyn,dzn,rhoa,rhoan,ipp,jpp,kpp,l_h,r_h,&
				tu,tv,tw,zu,zv,zw,su,sv,sw,p,dims,coords)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! advance scalars 1 time-step                                                !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if(viscous) then
                    ! advance fields
                    call advance_scalar_fields_3d(dt,&
                        q(:,:,:,1:nq),sq(:,:,:,1:nq),&
                        th,sth,rhoa,rhoan,ipp,jpp,kpp,nq,l_h,r_h)            
            endif
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! swap pointers (for efficiency)                                             !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			if(modulo(n,3).eq.1) then
				u => tut
				v => tvt
				w => twt
				zu => ut
				zv => vt
				zw => wt
				tu => zut
				tv => zvt
				tw => zwt
			else if(modulo(n,3).eq.2) then
				u => zut
				v => zvt
				w => zwt
				zu => tut
				zv => tvt
				zw => twt
				tu => ut
				tv => vt
				tw => wt
			else if(modulo(n,3).eq.0) then
				u => ut
				v => vt
				w => wt
				zu => zut
				zv => zvt
				zw => zwt
				tu => tut
				tv => tvt
				tw => twt
			endif
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! set halos																	 !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
			call exchange_full(ring_comm, id, kpp, jpp, ipp, r_h,r_h,r_h,r_h,l_h,r_h,u,&
								0._sp,0._sp,dims,coords)
			call exchange_full(ring_comm, id, kpp, jpp, ipp, r_h,r_h,l_h,r_h,r_h,r_h,v,&
								0._sp,0._sp,dims,coords)
			call exchange_full(ring_comm, id, kpp, jpp, ipp, l_h,r_h,r_h,r_h,r_h,r_h,w,&
								0._sp,0._sp,dims,coords)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			if(coords(3)==0) then
                w(0,:,:)=-w(1,:,:)
                u(1,:,:)=0._sp
                v(1,:,:)=0._sp
            endif
            
            if(coords(3)==(dims(3)-1)) then
                w(kpp,:,:)=-w(kpp-1,:,:)
                w(kpp+1,:,:)=0._sp
                u(kpp,:,:)=0._sp
                v(kpp,:,:)=0._sp
            endif
		enddo
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








		
	end subroutine model_driver
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	
	
	
	


	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>outputs variables to NetCDF file using MPI
	!>@param[inout] new_file: flag if this is a new file
	!>@param[in] outputfile: outputfilename
	!>@param[in] n: time-level
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
	!>@param[in] x,y,z, rhoa, thetan: grids
	!>@param[in] u,v,w,th,p,div: prognostic variables
	!>@param[in] id: id
	!>@param[in] world_process: world_process
	!>@param[in] rank: rank
	!>@param[in] ring_comm: ring_comm
	subroutine output(new_file,outputfile,n,ip,ipp,ipstart,jp,jpp,jpstart, &
					kp,kpp,kpstart,l_h,r_h, &
					time, &
					x,y,z,rhoa, thetan, &
					u,v,w,th,p,div, &
				    id, world_process, rank, ring_comm)
	
		use netcdf
		use mpi
		use mpi_module, only : mpi_integer9
		!use variables, only : MPI_INTEGER9

		implicit none
		logical, intent(inout) :: new_file
		character (len=*), intent(in) :: outputfile
		integer(i4b), intent(in) :: n, ip, ipp, ipstart, jp, jpp, jpstart, &
									kp, kpp, kpstart, l_h,r_h
		real(sp), intent(in) :: time
		real(sp), dimension(1-l_h:ipp+r_h), intent(in) :: x
		real(sp), dimension(1-l_h:jpp+r_h), intent(in) :: y
		real(sp), dimension(1-l_h:kpp+r_h), intent(in) :: z,rhoa, thetan
		real(sp), &
			dimension(1-r_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h), &
			intent(inout) :: th,p,div
		real(sp), &
			dimension(1-r_h:kpp+r_h,1-r_h:jpp+r_h,1-l_h:ipp+r_h), &
			intent(inout) :: u
		real(sp), &
			dimension(1-r_h:kpp+r_h,1-l_h:jpp+r_h,1-r_h:ipp+r_h), &
			intent(inout) :: v
		real(sp), &
			dimension(1-l_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h), &
			intent(inout) :: w
		
		integer(i4b), intent(in) :: id ,world_process, rank, ring_comm
	
		integer(i4b) :: ncid, x_dimid, nx_dimid, ny_dimid, nz_dimid, &
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

			! define variable: theta
			call check( nf90_def_var(ncid, "thetan", nf90_float, &
						(/nz_dimid/), varid) )
			! get id to a_dimid
			call check( nf90_inq_varid(ncid, "thetan", a_dimid) )
			! units
			call check( nf90_put_att(ncid, a_dimid, &
					   "units", "K") )

			! define variable: th
			call check( nf90_def_var(ncid, "th", nf90_float, &
						(/nz_dimid,ny_dimid,nx_dimid,x_dimid/), varid) )
			! get id to a_dimid
			call check( nf90_inq_varid(ncid, "th", a_dimid) )
			! units
			call check( nf90_put_att(ncid, a_dimid, &
					   "units", "K") )

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

			! define variable: p
			call check( nf90_def_var(ncid, "p", nf90_float, &
						(/nz_dimid,ny_dimid,nx_dimid,x_dimid/), varid) )
			! get id to a_dimid
			call check( nf90_inq_varid(ncid, "p", a_dimid) )
			! units
			call check( nf90_put_att(ncid, a_dimid, &
					   "units", "pa") )

			! define variable: div
			call check( nf90_def_var(ncid, "div", nf90_float, &
						(/nz_dimid,ny_dimid,nx_dimid,x_dimid/), varid) )
			! get id to a_dimid
			call check( nf90_inq_varid(ncid, "div", a_dimid) )
			! units
			call check( nf90_put_att(ncid, a_dimid, &
					   "units", "m/s**2") )



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
			! write variable: theta
			call check( nf90_inq_varid(ncid, "thetan", varid ) )
			call check( nf90_put_var(ncid, varid, thetan(1:kpp), &
						start = (/kpstart/)))	
		endif


		if(id==world_process) then
			! write variable: time
			call check( nf90_inq_varid(ncid, "time", varid ) )
			call check( nf90_put_var(ncid, varid, time, &
						start = (/n/)))	
	    endif
	    
	    
		! write variable: th
		call check( nf90_inq_varid(ncid, "th", varid ) )
		call check( nf90_put_var(ncid, varid, th(1:kpp,1:jpp,1:ipp), &
					start = (/kpstart,jpstart,ipstart,n/)))	

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

		! write variable: p
		call check( nf90_inq_varid(ncid, "p", varid ) )
		call check( nf90_put_var(ncid, varid, p(1:kpp,1:jpp,1:ipp), &
					start = (/kpstart,jpstart,ipstart,n/)))	

		! write variable: div
		call check( nf90_inq_varid(ncid, "div", varid ) )
		call check( nf90_put_var(ncid, varid, div(1:kpp,1:jpp,1:ipp), &
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
		use nrtype
		integer(i4b), intent ( in) :: status

		if(status /= nf90_noerr) then
			print *, trim(nf90_strerror(status))
			stop "Stopped"
		end if
	end subroutine check
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	end module drivers
