	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>drivers for the radiative transfer model
    module drivers
    use nrtype
    !use variables
    private
    public :: radiation_driver
    contains
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calls IO and runs one time-step of model
	!>@param[in] ntim: number of time-levels
	!>@param[in] dt: grids
	!>@param[in] ip, jp, kp: dims for full domain
	!>@param[in] ipp, jpp, kpp: dims for this block of data
	!>@param[in] l_h,r_h, ipstart,jpstart,kpstart, nbands: halo and dims, number of bands
	!>@param[inout] flux_u, flux_d, rad_power
	!>@param[in] x,y,z, dx, dy, dz, dxn, dyn, dzn: grids
	!>@param[in] theta, thetan: reference variables
	!>@param[in] tref, trefn: reference variables
	!>@param[in] rhoa, rhoan: reference variables
	!>@param[in] lambda, lambda_low, lambda_high, delta_lambda, nrwbin,niwbin
	!>@param[in] sflux_l, b_s_g
	!>@param[in] start_year, start_mon,start_day,start_hour,start_min,start_sec
	!>@param[in] lat, lon, albedo, emiss,quad_flag, asymmetry_water
	!>@param[in] coords: for Cartesian topology
	!>@param[inout] new_file: flag for if this is a new file
	!>@param[in] outputfile: netcdf output
	!>@param[in] output_interval: interval for output (s)
	!>@param[in] dims,id, world_process, ring_comm: mpi variables
    subroutine radiation_driver(ntim,dt,l_h,r_h, &
    			ip,jp,kp, &
    			ipp,jpp,kpp, &
				ipstart, jpstart, kpstart, &
				tdstart,tdend, &
				a,b,c,r,usol, &
				nbands, ns,nl, flux_u, flux_d, rad_power, &
				x,y,z, &
				dx,dy,dz, &
				dxn,dyn,dzn, &
				theta,thetan, &
				tref,trefn, &
				rhoa,rhoan, &
				lambda,lambda_low,lambda_high, delta_lambda,&
				nrwbin,niwbin, &
				sflux_l, b_s_g, &
				start_year,start_mon, start_day,start_hour,start_min,start_sec, &
				lat, lon, albedo, emiss,quad_flag, &
				asymmetry_water, &
				nprocv,mvrecv, &
				coords, &
				new_file,outputfile, output_interval, &
				dims,id, world_process, rank, ring_comm,sub_comm)
		use nrtype
		use mpi_module, only : exchange_full, exchange_along_dim, exchange_d_fluxes, &
		                        exchange_u_fluxes
		use radiation, only : solve_fluxes
		use mpi

	
		implicit none
		logical, intent(inout) :: new_file
		real(sp), dimension(nbands) :: lambda, b_s_g, lambda_low,&
						lambda_high, delta_lambda, nrwbin, niwbin, sflux_l
		real(sp), intent(in), dimension(nprocv) :: mvrecv
		integer(i4b), intent(in) :: ntim,ip,jp,kp, ipp,jpp,kpp, &
						l_h,r_h, ipstart, jpstart, kpstart, nbands, ns, nl, &
						tdstart,tdend
		integer(i4b), intent(in) :: id, world_process, ring_comm, sub_comm,rank,nprocv
		integer(i4b), dimension(3), intent(in) :: coords, dims
		character (len=*), intent(in) :: outputfile
		real(sp), intent(inout), dimension(1:tdend) :: a,b,c,r,usol
		real(sp), intent(in) :: output_interval, dt
		real(sp), dimension(1-r_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h,1:nbands) :: &
						flux_d,flux_u
		real(sp), dimension(1-r_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h) :: &
						rad_power
		real(sp), dimension(1-l_h:ipp+r_h), intent(in) :: x,dx, dxn
		real(sp), dimension(1-l_h:jpp+r_h), intent(in) :: y,dy,dyn
		real(sp), dimension(1-l_h:kpp+r_h), intent(in) :: z,dz,dzn, theta, thetan, &
														rhoa, rhoan, tref, trefn
		real(sp), intent(in) :: lat, lon, albedo, emiss, asymmetry_water
		integer(i4b), intent(in) :: quad_flag, start_year, start_mon, start_day, &
									start_hour, start_min, start_sec
		! locals:		
		integer(i4b) :: n,n2, cur=1, i,j,k, error, rank2
		real(sp) :: time, time_last_output, output_time
        real(sp), dimension(1-r_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h) :: th
        real(sp), dimension(1-r_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h,0) :: &
                ngs,lamgs,mugs
		
        if(id>=dims(1)*dims(2)*dims(3)) return 
        
        th=0._sp ! dummy variable for this code


		
		time_last_output=-output_interval
		output_time=output_interval
		rank2=dims(1)*dims(2)*dims(3)
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! time-loop                                                                      !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		do n=1,ntim	
			time=real(n-1,sp)*dt


			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! main solver, solve fluxes in each band                                     !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			call solve_fluxes(start_year,start_mon,start_day, &
								start_hour,start_min,start_sec, &
								lat, lon, &
								time, nbands,ns,nl,ipp,jpp,kpp,r_h, &
								tdstart,tdend,a,b,c,r,usol, &
								lambda_low, lambda_high, lambda,nrwbin, niwbin, &
								b_s_g,sflux_l, &
								rhoan, tref, trefn, dz,dzn, albedo, emiss, &
								quad_flag, th, flux_u, flux_d, &
								asymmetry_water, 0, ngs,lamgs,mugs, &
								.false., &
								nprocv,mvrecv, &
								coords,dims, id, sub_comm)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


		    call exchange_d_fluxes(ring_comm, id, kpp, jpp, ipp, nbands,&
							r_h,r_h,r_h,r_h,r_h, r_h,  flux_d, dims,coords)
		    call exchange_u_fluxes(ring_comm, id, kpp, jpp, ipp, nbands,&
							r_h,r_h,r_h,r_h,r_h, r_h,  flux_u, dims,coords)
							
							
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! write netcdf variables                                                     !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			if (time-time_last_output >= output_interval) then
			
			
				if (id==world_process) &
					print *,'output no ',cur,' at time (s) ', &
						time,n,' steps of ',ntim

				
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! Output to NetCDF                                                       !
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				call output(new_file,outputfile,cur, &
							ip,ipp,ipstart,jp,jpp,jpstart,kp,kpp,kpstart, &
							l_h,r_h,nbands, &
							time, x,y,z,rhoa, thetan,trefn, &
							lambda_low, lambda_high,&
							lambda,sflux_l,flux_u,flux_d,rad_power, &
							id, world_process, rank2, ring_comm)
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


				time_last_output=time
				cur=cur+1
			endif
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


			! update temperatures here
			
			call mpi_barrier(ring_comm,error)

		enddo
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








		
	end subroutine radiation_driver
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
	!>@param[in] x,y,z, rhoa, thetan,trefn: grids
	!>@param[in] lambda_low_lambda_high,lambda, sflux_l
	!>@param[in] u,v,w,th,p: prognostic variables
	!>@param[in] id: id
	!>@param[in] world_process: world_process
	!>@param[in] rank: rank
	!>@param[in] ring_comm: ring_comm
	subroutine output(new_file,outputfile,n,ip,ipp,ipstart,jp,jpp,jpstart, &
					kp,kpp,kpstart,l_h,r_h,nbands, &
					time, &
					x,y,z,rhoa, thetan,trefn, &
					lambda_low, lambda_high,&
					lambda,sflux_l,flux_u,flux_d,rad_power, &
				    id, world_process, rank, ring_comm)
	
		use netcdf
		use mpi
		use mpi_module, only : mpi_integer9
		!use variables, only : MPI_INTEGER9

		implicit none
		logical, intent(inout) :: new_file
		character (len=*), intent(in) :: outputfile
		integer(i4b), intent(in) :: n, ip, ipp, ipstart, jp, jpp, jpstart, &
									kp, kpp, kpstart, l_h,r_h,nbands
		real(sp), intent(in) :: time
		real(sp), dimension(1-l_h:ipp+r_h), intent(in) :: x
		real(sp), dimension(1-l_h:jpp+r_h), intent(in) :: y
		real(sp), dimension(1-l_h:kpp+r_h), intent(in) :: z,rhoa, thetan, trefn
		real(sp), &
			dimension(1-r_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h,1:nbands), &
			intent(in) :: flux_u, flux_d
		real(sp), &
			dimension(1-r_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h), &
			intent(in) :: rad_power
		real(sp), dimension(nbands) :: lambda_low,lambda_high,lambda, sflux_l
		
		integer(i4b), intent(in) :: id ,world_process, rank, ring_comm
	
		integer(i4b) :: ncid, x_dimid, nx_dimid, ny_dimid, nz_dimid, &
						error, varid,a_dimid, nb_dimid, id_go
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
			call check( nf90_def_dim(ncid, "nbands", nbands, nb_dimid) )


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
						
			! define variable: lambda_low
			call check( nf90_def_var(ncid, "lambda_low", nf90_float, &
						(/nb_dimid/), varid) )
			! get id to a_dimid
			call check( nf90_inq_varid(ncid, "lambda_low", a_dimid) )
			! units
			call check( nf90_put_att(ncid, a_dimid, &
					   "units", "m") )
					   
			! define variable: lambda_high
			call check( nf90_def_var(ncid, "lambda_high", nf90_float, &
						(/nb_dimid/), varid) )
			! get id to a_dimid
			call check( nf90_inq_varid(ncid, "lambda_high", a_dimid) )
			! units
			call check( nf90_put_att(ncid, a_dimid, &
					   "units", "m") )
					   
			! define variable: lambda
			call check( nf90_def_var(ncid, "lambda", nf90_float, &
						(/nb_dimid/), varid) )
			! get id to a_dimid
			call check( nf90_inq_varid(ncid, "lambda", a_dimid) )
			! units
			call check( nf90_put_att(ncid, a_dimid, &
					   "units", "m") )
						
			! define variable: sflux
			call check( nf90_def_var(ncid, "sflux", nf90_float, &
						(/nb_dimid/), varid) )
			! get id to a_dimid
			call check( nf90_inq_varid(ncid, "sflux", a_dimid) )
			! units
			call check( nf90_put_att(ncid, a_dimid, &
					   "units", "W m-2") )
						
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

			! define variable: trefn
			call check( nf90_def_var(ncid, "trefn", nf90_float, &
						(/nz_dimid/), varid) )
			! get id to a_dimid
			call check( nf90_inq_varid(ncid, "trefn", a_dimid) )
			! units
			call check( nf90_put_att(ncid, a_dimid, &
					   "units", "K") )

			! define variable: flux_u
			call check( nf90_def_var(ncid, "flux_u", nf90_float, &
						(/nz_dimid,ny_dimid,nx_dimid,nb_dimid,x_dimid/), varid) )
			! get id to a_dimid
			call check( nf90_inq_varid(ncid, "flux_u", a_dimid) )
			! units
			call check( nf90_put_att(ncid, a_dimid, &
					   "units", "W m-2") )

			! define variable: flux_d
			call check( nf90_def_var(ncid, "flux_d", nf90_float, &
						(/nz_dimid,ny_dimid,nx_dimid,nb_dimid,x_dimid/), varid) )
			! get id to a_dimid
			call check( nf90_inq_varid(ncid, "flux_d", a_dimid) )
			! units
			call check( nf90_put_att(ncid, a_dimid, &
					   "units", "W m-2") )

			! define variable: heating rate
			call check( nf90_def_var(ncid, "rad_heat", nf90_float, &
						(/nz_dimid,ny_dimid,nx_dimid,x_dimid/), varid) )
			! get id to a_dimid
			call check( nf90_inq_varid(ncid, "rad_heat", a_dimid) )
			! units
			call check( nf90_put_att(ncid, a_dimid, &
					   "units", "W") )


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
			! write variable: tref
			call check( nf90_inq_varid(ncid, "trefn", varid ) )
			call check( nf90_put_var(ncid, varid, trefn(1:kpp), &
						start = (/kpstart/)))	
			if(id==world_process) then
				! write variable: lambda_low
				call check( nf90_inq_varid(ncid, "lambda_low", varid ) )
				call check( nf90_put_var(ncid, varid, lambda_low(1:nbands), &
							start = (/1/)))				
				! write variable: lambda_high
				call check( nf90_inq_varid(ncid, "lambda_high", varid ) )
				call check( nf90_put_var(ncid, varid, lambda_high(1:nbands), &
							start = (/1/)))				
				! write variable: lambda
				call check( nf90_inq_varid(ncid, "lambda", varid ) )
				call check( nf90_put_var(ncid, varid, lambda(1:nbands), &
							start = (/1/)))				
				! write variable: sflux
				call check( nf90_inq_varid(ncid, "sflux", varid ) )
				call check( nf90_put_var(ncid, varid, sflux_l(1:nbands), &
							start = (/1/)))				
			endif
		endif


		if(id==world_process) then
			! write variable: time
			call check( nf90_inq_varid(ncid, "time", varid ) )
			call check( nf90_put_var(ncid, varid, time, &
						start = (/n/)))	
	    endif
	    
	    
		! write variable: th
		call check( nf90_inq_varid(ncid, "flux_u", varid ) )
		call check( nf90_put_var(ncid, varid, flux_u(1:kpp,1:jpp,1:ipp,1:nbands), &
					start = (/kpstart,jpstart,ipstart,1,n/)))	

		! write variable: u
		call check( nf90_inq_varid(ncid, "flux_d", varid ) )
		call check( nf90_put_var(ncid, varid, flux_d(1:kpp,1:jpp,1:ipp,1:nbands), &
					start = (/kpstart,jpstart,ipstart,1,n/)))

		! write variable: v
! 		call check( nf90_inq_varid(ncid, "rad_heat", varid ) )
! 		call check( nf90_put_var(ncid, varid, rad_power(1:kpp,1:jpp,1:ipp), &
! 					start = (/kpstart,jpstart,ipstart,n/)))


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
