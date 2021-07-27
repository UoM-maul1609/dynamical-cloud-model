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
	!>@param[in] pref,prefn
	!>@param[inout] ubar,vbar,wbar,thbar,qbar
	!>@param[in] u_force, v_force, forcing_tau, w_subs
	!>@param[inout] ut,vt,wt: prognostics
	!>@param[inout] zut,zvt,zwt: prognostics - previous time-step
	!>@param[inout] tut,tvt,twt: prognostics - temp storage
	!>@param[inout] th,sth,p: prognostic variables
	!>@param[inout] su,sv,sw,psrc, div: more prognostic variables
	!>@param[inout] strain,vism,vist - viscosity subgrid 
	!>@param[in] z0,z0th - roughness lengths 
	!>@param[in] ptol - tolerance for pressure solver
	!>@param[in] c_s, c_e: start and end indices for a category
	!>@param[in] cat_am,cat_c, cat_r,cat_i: category index for cloud and rain and ice
	!>@param[in] n_mode,inc,iqc,inr,iqr,ini, iqi, iai
	!>@param[in] q_name: name of q-variables
	!>@param[inout] q,sq,viss - for clouds and subgrid
	!>@param[inout] precip - array for precipitation
	!>@param[in] theta, thetan: reference variables
	!>@param[in] tref, trefn: reference variables
	!>@param[in] rhoa, rhoan: reference variables
	!>@param[in] lamsq, lamsqn: reference variables
	!>@param[inout] thbase, thtop: bottom and base th
	!>@param[inout] micro_init: initialise microphysics
	!>@param[inout] new_file: flag for if this is a new file
	!>@param[in] outputfile: netcdf output
	!>@param[in] output_interval: interval for output (s)
	!>@param[in] viscous: logical for applying viscous dissipation
	!>@param[in] advection_scheme, kord, monotone: flags for advection schemes
	!>@param[in] moisture: flag for moisture
	!>@param[in] microphysics_flag: flag for microphysics
	!>@param[in] ice_flag: flag for ice microphysics
	!>@param[in] hm_flag: flag for hallett-mossop (not always used)
	!>@param[in] wr_flag: flag for warm rain
	!>@param[in] theta_flag: flag for adjusting theta
	!>@param[in] damping_layer: flag for damping layer
	!>@param[in] forcing: flag for large-scale forcing of horizontal winds
	!>@param[in] divergence: flag for large-scale divergence
	!>@param[in] radiation: radiation scheme
	!>@param[in] j_stochastic, ice_nuc_flag, mode2_ice_flag, coll_breakup_flag1,
	!>@param[in] heyms_west, lawson, recycle
	!>@param[in] nq,nprec,ncat: number of q-variables
	!>@param[in] tdstart, tdend: start and end array indices - parallel ts
	!>@param[inout] a, b, c, r, usol: tridag
    !>@param[in] nbands, ns, nl: number of bands
    !>@param[inout] flux_u, flux_d, rad_power
    !>@param[in] lambda, lambda_low, lambda_high, delta_lambda, nrwbin,niwbin
    !>@param[in] sflux_l, b_s_g
    !>@param[in] start_year, start_mon,start_day,start_hour,start_min,start_sec
    !>@param[in] lat, lon, albedo, emiss,quad_flag, asymmetry_water
    !>@param[in] nrad
    !>@param[inout] ngs,lamgs,mugs
	!>@param[in] coords: for Cartesian topology
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
                u_force,v_force,forcing_tau, &
                w_subs, &
				ut,vt,wt,&
				zut,zvt,zwt,&
				tut,tvt,twt,&
				th,sth,p, &
				su,sv,sw,psrc, &
				div, &
				strain, vism, vist, z0,z0th, ptol, &
				c_s, c_e, cat_am,cat_c,cat_r, cat_i,&
				n_mode,inc,iqc, inr,iqr, ini, iqi, iai, &
				q_name, q,sq,viss, &
				precip, &
				theta,thetan, &
				tref,trefn, &
				rhoa,rhoan, &
				pref,prefn, &
				lamsq,lamsqn, &
				thbase,thtop, &
				micro_init, &
				new_file,outputfile, output_interval, &
				viscous, &
				advection_scheme, kord, monotone, &
				moisture, microphysics_flag, ice_flag, hm_flag, wr_flag, theta_flag, &
				damping_layer,forcing, divergence, radiation, &
				j_stochastic, ice_nuc_flag, mode2_ice_flag, coll_breakup_flag1, &
				heyms_west, lawson, recycle, &
				nq, nprec, ncat, &
                    tdstart,tdend, &
                    a,b,c,r,usol, &
                    nbands, ns,nl, flux_u, flux_d, rad_power, &
                    lambda,lambda_low,lambda_high, delta_lambda,&
                    nrwbin,niwbin, &
                    sflux_l, b_s_g, &
                    start_year,start_mon, start_day,start_hour,start_min,start_sec, &
                    lat, lon, albedo, emiss,quad_flag, asymmetry_water, &
                    nrad,ngs,lamgs,mugs, &
                    nprocv,mvrecv, &
				coords, &
				dims,id, world_process, rank, ring_comm,sub_horiz_comm,sub_vert_comm)
		use nrtype
		use mpi_module, only : exchange_full, exchange_along_dim, exchange_along_dim_wo
        use advection_s_3d, only : first_order_upstream_3d, &
                    mpdata_3d, mpdata_vec_3d, adv_ref_state, mpdata_3d_add, &
                    mpdata_vert_3d, mpdata_vec_vert_3d
		use d_solver, only : bicgstab, sources, advance_momentum, radiative_transfer
		use subgrid_3d, only : advance_scalar_fields_3d
        use diagnostics
        use p_micro_module
        
		implicit none
				
		
		logical, intent(inout) :: new_file
		logical, intent(in) :: viscous, monotone, moisture, &
		                    damping_layer, forcing, theta_flag, ice_flag, hm_flag, &
		                    wr_flag, divergence, radiation, heyms_west, lawson, recycle
		integer(i4b), intent(in) :: ice_nuc_flag, nrad, mode2_ice_flag, coll_breakup_flag1
		logical, intent(inout) :: micro_init
		integer(i4b), intent(in) :: ntim,ip,jp,kp, ipp,jpp,kpp, &
						l_h,r_h, ipstart, jpstart, kpstart, &
						advection_scheme, kord, nq,nprec,ncat, microphysics_flag
		integer(i4b), intent(in) :: id, world_process, ring_comm, sub_horiz_comm, &
		    sub_vert_comm,rank
		integer(i4b), dimension(3), intent(in) :: coords, dims
		character (len=*), intent(in) :: outputfile
		real(sp), intent(in) :: output_interval, dt, z0,z0th, ptol, forcing_tau, &
		                        j_stochastic
        integer(i4b), dimension(ncat), intent(in) :: c_s, c_e
        integer(i4b), intent(in) :: cat_am,cat_c, cat_r, cat_i,n_mode,inc,iqc,inr,iqr, &
                                ini, iqi, iai
		
		real(sp), intent(inout) :: thbase, thtop
		real(sp), dimension(1-l_h:ipp+r_h), intent(in) :: x,xn,dx, dxn
		real(sp), dimension(1-l_h:jpp+r_h), intent(in) :: y,yn,dy,dyn
		real(sp), dimension(1-l_h:kpp+r_h), intent(in) :: z,zn,dz,dzn, theta, thetan, &
		                                                tref,trefn, &
														rhoa, rhoan,lamsq, lamsqn, &
														u_force, v_force, w_subs
		real(sp), dimension(1-l_h:kpp+r_h), intent(inout) :: ubar,vbar,wbar,thbar,qbar
		real(sp), dimension(1-l_h:kpp+r_h), intent(in) :: dampfacn,dampfac
		real(sp), dimension(1-l_h:kpp+r_h), intent(in) :: pref,prefn
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

		real(sp), &
			dimension(1:kpp,1-l_h:jpp+r_h,1-l_h:ipp+r_h,1:nprec), &
			intent(inout) :: precip
					
		real(sp), &
			dimension(1-l_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h,1:nrad), &
			intent(inout) :: ngs,lamgs,mugs

        character(len=20), intent(in), dimension(nq) :: q_name

        ! radiation variables		
		real(sp), dimension(nbands) :: lambda, b_s_g, lambda_low,&
							lambda_high, delta_lambda, nrwbin,niwbin, sflux_l
		real(sp), intent(in), dimension(nprocv) :: mvrecv
		integer(i4b), intent(in) :: nbands, ns, nl, tdstart,tdend
		integer(i4b), intent(in) :: nprocv
		real(sp), intent(inout), dimension(1:tdend) :: a,b,c,r,usol
		real(sp), dimension(1-r_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h,1:nbands) :: &
						flux_d,flux_u
		real(sp), dimension(1-r_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h) :: &
						rad_power
		real(sp), intent(in) :: lat, lon, albedo, emiss, asymmetry_water
		integer(i4b), intent(in) :: quad_flag, start_year, start_mon, start_day, &
									start_hour, start_min, start_sec
		!-
		


		! locals:		
		integer(i4b) :: n,n2, cur=1, i,j,k,nqc, error, rank2
		real(sp) :: time, time_last_output, output_time, t1=0._sp,t2=0._sp
		real(sp), dimension(nq) :: q1,q2
		real(sp), dimension(:,:,:), pointer :: u,zu,tu
		real(sp), dimension(:,:,:), pointer :: v,zv,tv
		real(sp), dimension(:,:,:), pointer :: w,zw,tw
! 		real(sp), &
! 			dimension(1-r_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h) :: th2
		

		time_last_output=-output_interval
		output_time=output_interval
		rank2=dims(1)*dims(2)*dims(3)
		if(id>=rank2) return 

		q1=0._sp
		q2=0._sp
		
		
		! associate pointers - for efficiency, when swapping arrays in leap-frog scheme
		u => ut;   v => vt;   w => wt   ! current time-step
		zu => zut; zv => zvt; zw => zwt ! previous time-step		
		tu => tut; tv => tvt; tw => twt ! next-time-step
		
		
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(coords(3)==0) then
            zw(0,:,:)=-zw(1,:,:); zu(1,:,:)=0._sp; zv(1,:,:)=0._sp
        endif
        
        if(coords(3)==(dims(3)-1)) then
            zw(kpp,:,:)=-zw(kpp-1,:,:); zw(kpp+1,:,:)=0._sp
            zu(kpp,:,:)=0._sp
            zv(kpp,:,:)=0._sp
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		



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
							time, x,y,z,rhoa, thetan, trefn, &
							u,v,w,th,p,div, &
							q_name,q,nq, precip, nprec,moisture, &
							id, world_process, rank2, ring_comm)
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



				time_last_output=time
				cur=cur+1
			endif
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! try doing radiation here                                                   !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if(radiation) then
                call radiative_transfer(ring_comm,id,rank2, dims, coords, &
                    dt,dz,dzn,ipp,jpp,kpp,l_h,r_h,&
                    th,sth,&
                    tref,trefn, rhoa,rhoan,&
                    sub_vert_comm, time,tdstart,tdend, &
                    a,b,c,r,usol, &
                    nbands, ns,nl, flux_u, flux_d, &
                    lambda,lambda_low,lambda_high, delta_lambda,nrwbin, niwbin, &
                    sflux_l, b_s_g, &
                    start_year,start_mon, start_day,start_hour,start_min,start_sec, &
                    lat, lon, albedo, emiss,quad_flag, &
                    asymmetry_water, &
                    nrad,ngs,lamgs,mugs, &
                    nprocv,mvrecv)	
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
			! calculate sources of momentum, theta, q, pressure                          !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			call sources(ring_comm, id, rank2, dims,coords, &
				dt,x,y,z,zn,dx,dy,dz,dxn,dyn,dzn,ipp,jpp,kpp,l_h,r_h,&
				nq, &
				ubar,vbar,wbar,thbar,qbar, dampfacn,dampfac, &
				u_force,v_force,forcing_tau, &
				w_subs, &
				zu,zv,zw, &
				u,v,w,su,sv,sw, &
				q,sq,viss, &
				psrc, &
				th,sth, strain,vism,vist, &
				theta,thetan,tref, trefn, rhoa,rhoan,lamsq,lamsqn, &
				z0,z0th, &
				viscous, moisture,damping_layer, forcing, divergence, radiation, &
                sub_vert_comm, time,tdstart,tdend, &
                a,b,c,r,usol, &
                nbands, ns,nl, flux_u, flux_d, rad_power, &
                lambda,lambda_low,lambda_high, delta_lambda,&
                nrwbin,niwbin, &
                sflux_l, b_s_g, &
                start_year,start_mon, start_day,start_hour,start_min,start_sec, &
                lat, lon, albedo, emiss,quad_flag, &
                asymmetry_water, &
                nrad,ngs,lamgs,mugs, &
                nprocv,mvrecv)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! set halos																	 !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
			call exchange_along_dim_wo(ring_comm, id, kpp, jpp, ipp, &
								r_h,r_h,r_h,r_h,r_h,r_h, psrc, &
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
			call exchange_along_dim_wo(ring_comm, id, kpp, jpp, ipp, &
								r_h,r_h,r_h,r_h,r_h,r_h, p, dims,coords)
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
			! advect scalar fields using mid-point                                       !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			select case (advection_scheme)
				case (0)
					call first_order_upstream_3d(dt,dxn,dyn,dzn,rhoa,rhoan, &
						ipp,jpp,kpp,l_h,r_h,u,v,w,th,0,dims,coords)
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! set halos													   		 !
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
                    call exchange_full(ring_comm, id, kpp, jpp, ipp, &
                                        r_h,r_h,r_h,r_h,r_h,r_h,th,t1,t2,dims,coords)
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
				case (1)
					call mpdata_3d(dt,dx,dy,dz,dxn,dyn,dzn,rhoa,rhoan, &
						ipp,jpp,kpp,l_h,r_h,u,v,w,th,t1,t2, &
						kord,monotone,0,ring_comm,id, &
						dims,coords)						
				case(2)
					call mpdata_3d_add(dt,dx,dy,dz,dxn,dyn,dzn,rhoa,rhoan, &
						ipp,jpp,kpp,l_h,r_h,u,v,w,th,thetan,thbase,thtop, &
						kord,monotone,0,ring_comm,id, &
						dims,coords)
				case default
					print *,'not coded'
					stop
			end select
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! advect q-fields                                                            !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if(moisture) then
                select case (advection_scheme)
                    case (0:2)
                        do nqc=1,ncat
                            call mpdata_vec_3d(dt,dx,dy,dz,dxn,dyn,dzn,rhoa,rhoan, &
                                ipp,jpp,kpp,c_e(nqc)-c_s(nqc)+1, &
                                l_h,r_h,u,v,w,q(:,:,:,c_s(nqc):c_e(nqc)), &
                                q1(c_s(nqc):c_e(nqc)),q2(c_s(nqc):c_e(nqc)), &
                                kord,monotone,0,ring_comm,id, &
                                dims,coords)						
                        enddo
                    case default
                        print *,'not coded'
                        stop
                end select            
            endif
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


			



			



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
                        q,sq,&
                        th,sth,rhoa,rhoan,ipp,jpp,kpp,nq,l_h,r_h,moisture)            
            endif
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! microphysics                                                               !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if(moisture) then
                select case (microphysics_flag)
                    case (0) ! null microphysics

                    case (1) ! not coded yet
                    			
                    case (2) ! not coded yet

                    case (3) ! pamm microphysics
                        call p_microphysics_3d(nq,ncat,n_mode,c_s,c_e, inc, iqc,&
                                        inr,iqr,ini,iqi,&
                                        iai, &
                                        cat_am,cat_c, cat_r, cat_i,&
                                        nprec, &
                                        ipp,jpp,kpp,l_h,r_h,dt,dz,&
                                        dzn,q(:,:,:,:),precip(:,:,:,:),&
                                        nrad,ngs(:,:,:,:),lamgs(:,:,:,:),mugs(:,:,:,:), &
                                        th(:,:,:),prefn, &
                                        zn(:),thetan,rhoa,rhoan,w(:,:,:), &
                                        micro_init,hm_flag,wr_flag, 1.e-14_sp, &
                                        ice_flag, theta_flag, &
                                        j_stochastic, ice_nuc_flag, mode2_ice_flag, &
                                        coll_breakup_flag1, &
                                        heyms_west, lawson, recycle, &
                                        radiation, &
                                        ring_comm,sub_vert_comm,id,dims,coords)		
                    case default
                        print *,'not coded'
                        stop
                end select
            endif
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! time-smoothing                                                             !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            u=0.99_sp*u+0.01_sp*zu
            v=0.99_sp*v+0.01_sp*zv
            w=0.99_sp*w+0.01_sp*zw
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! swap pointers (for efficiency)                                             !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			if(modulo(n,3).eq.1) then
				u => tut;  v => tvt;  w => twt
				zu => ut;  zv => vt;  zw => wt
				tu => zut; tv => zvt; tw => zwt
			else if(modulo(n,3).eq.2) then
				u => zut;  v => zvt;  w => zwt
				zu => tut; zv => tvt; zw => twt
				tu => ut;  tv => vt;  tw => wt
			else if(modulo(n,3).eq.0) then
				u => ut;   v => vt;   w => wt
				zu => zut; zv => zvt; zw => zwt
				tu => tut; tv => tvt; tw => twt
			endif
            ! u,v,w are now the values on the current time-step in the next iteration
            ! zu,zv,zw are now the values on the previous time-step in the next iteration
            ! tu,tv,tw are now the values on the future time-step in the next iteration
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
                w(0,:,:)=0._sp ! the vertical velocity at the surface is zero
                u(0,:,:)=-u(1,:,:) ! average results in zero horizontal velocity at ground
                v(0,:,:)=-v(1,:,:)
            endif
            
            if(coords(3)==(dims(3)-1)) then
                !w(kpp,:,:)=-w(kpp-1,:,:)
                w(kpp+1,:,:)=0._sp
                u(kpp+1,:,:)=u(kpp,:,:)
                v(kpp+1,:,:)=v(kpp,:,:)
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
	!>@param[in] nq: number of q-variables
	!>@param[in] l_h,r_h: halo
	!>@param[in] time: time (s)
	!>@param[in] x,y,z, rhoa, thetan, trefn: grids
	!>@param[in] u,v,w,th,p,div: prognostic variables
	!>@param[in] q_name: name of q-variables
	!>@param[in] q
	!>@param[in] precip, nprec
	!>@param[in] moisture: flag for moisture
	!>@param[in] id: id
	!>@param[in] world_process: world_process
	!>@param[in] rank: rank
	!>@param[in] ring_comm: ring_comm
	subroutine output(new_file,outputfile,n,ip,ipp,ipstart,jp,jpp,jpstart, &
					kp,kpp,kpstart,l_h,r_h, &
					time, &
					x,y,z,rhoa, thetan, trefn, &
					u,v,w,th,p,div, &
					q_name,q,nq, precip, nprec,&
					moisture, &
				    id, world_process, rank, ring_comm)
	
		use netcdf
		use mpi
		use mpi_module, only : mpi_integer9
		!use variables, only : MPI_INTEGER9

		implicit none
		logical, intent(inout) :: new_file
		character (len=*), intent(in) :: outputfile
		integer(i4b), intent(in) :: n, ip, ipp, ipstart, jp, jpp, jpstart, &
									kp, kpp, kpstart, l_h,r_h, nq,nprec
		real(sp), intent(in) :: time
		real(sp), dimension(1-l_h:ipp+r_h), intent(in) :: x
		real(sp), dimension(1-l_h:jpp+r_h), intent(in) :: y
		real(sp), dimension(1-l_h:kpp+r_h), intent(in) :: z,rhoa, thetan, trefn
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

        character(len=20), intent(in), dimension(nq) :: q_name
        real(sp), dimension(1-l_h:kpp+r_h,1-r_h:jpp+r_h,1-r_h:ipp+r_h,nq), intent(in) :: q
        real(sp), dimension(1:kpp,1-r_h:jpp+r_h,1-r_h:ipp+r_h,1:nprec), &
            intent(in) :: precip
		logical, intent(in) :: moisture
		integer(i4b), intent(in) :: id ,world_process, rank, ring_comm
	
		integer(i4b) :: ncid, x_dimid, nx_dimid, ny_dimid, nz_dimid, &
						error, varid,a_dimid, id_go, lq_dimid,nq_dimid,nprec_dimid
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
			if(moisture) then
                call check( nf90_def_dim(ncid, "nq", nq, nq_dimid) )
                call check( nf90_def_dim(ncid, "nprec", nprec, nprec_dimid) )
                call check( nf90_def_dim(ncid, "l_q_names", 20, lq_dimid) )
            endif

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

			! define variable: trefn
			call check( nf90_def_var(ncid, "trefn", nf90_float, &
						(/nz_dimid/), varid) )
			! get id to a_dimid
			call check( nf90_inq_varid(ncid, "trefn", a_dimid) )
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

            if(moisture) then
                ! define variable: q_names
                call check( nf90_def_var(ncid, "q_names", NF90_CHAR, &
                    (/lq_dimid,nq_dimid/), varid) )

                ! define variable: q
                call check( nf90_def_var(ncid, "q", nf90_float, &
                    (/nq_dimid, nz_dimid, ny_dimid, nx_dimid,x_dimid/), varid) )
                ! get id to a_dimid
                call check( nf90_inq_varid(ncid, "q", a_dimid) )
                ! units
                call check( nf90_put_att(ncid, a_dimid, &
                           "units", "kg or number per kg") )
            
                ! define variable: precip
                call check( nf90_def_var(ncid, "precip", nf90_float, &
                    (/nprec_dimid, nz_dimid, ny_dimid, nx_dimid,x_dimid/), varid) )
                ! get id to a_dimid
                call check( nf90_inq_varid(ncid, "precip", a_dimid) )
                ! units
                call check( nf90_put_att(ncid, a_dimid, &
                           "units", "mm hr-1 or number per m-2 s-1") )
            
            endif

			! exit define mode
			call check( nf90_enddef(ncid) )
						
			call check( nf90_close(ncid) )

			new_file=.false.

            call check( nf90_open(outputfile, NF90_WRITE, ncid) )
			if(moisture) then
                call check( nf90_inq_varid(ncid, "q_names", varid ) )
                call check( nf90_put_var(ncid, varid, q_name, start = (/1,1/)))
            endif			
    		call check( nf90_close(ncid) )
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
			! write variable: trefn
			call check( nf90_inq_varid(ncid, "trefn", varid ) )
			call check( nf90_put_var(ncid, varid, trefn(1:kpp), &
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
					
		if(moisture) then
            call check( nf90_inq_varid(ncid, "q", varid ) )
            call check( nf90_put_var(ncid, varid, &
                reshape( &
                reshape(q(1:kpp,1:jpp,1:ipp,1:nq),[nq,ipp,jpp,kpp],order=[4,3,2,1]), &
                    [nq,kpp,jpp,ipp],order=[1,4,3,2]), &
                        start = (/1,kpstart,jpstart,ipstart,n/)))
                        	
            call check( nf90_inq_varid(ncid, "precip", varid ) )
            call check( nf90_put_var(ncid, varid, &
                reshape( &
                reshape(precip(1:kpp,1:jpp,1:ipp,1:nprec),[nprec,ipp,jpp,kpp],&
                    order=[4,3,2,1]), &
                    [nprec,kpp,jpp,ipp],order=[1,4,3,2]), &
                        start = (/1,kpstart,jpstart,ipstart,n/)))	
		endif
		

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
