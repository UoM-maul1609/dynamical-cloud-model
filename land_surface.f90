	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>land surface model solver routines 
    module lsm
    use nrtype
    use nr, only : locate, polint, qromb, tridag
    implicit none


		!>@brief
		!>main model prognostic variables for lsm
        type lsmg
        	integer(i4b) :: year, mon, day, hour, sec
        	integer(i4b) :: doy, diy
        	integer(i4b) :: nprocv
        	real(sp) :: albedo, emiss, lat, lon
        	integer(i4b) :: ns, nl, ntot, skp
            real(sp), dimension(:), allocatable ::	b_s_g, & ! cross section for Rayleigh
            	lambda, lambda_low, lambda_high, delta_lambda, & 
            	nrwbin,niwbin ! bin averaged refractive indices
            real(sp), dimension(:), allocatable ::	ext_s_g ! extinction for Rayleigh
            real(sp), dimension(:,:,:,:), allocatable :: flux_u, flux_d
            real(sp), dimension(:,:,:), allocatable :: rad_power, t, wg
            real(sp), dimension(:,:), allocatable :: tsurf
            real(sp), dimension(:), allocatable :: sflux_l	
            real(sp), dimension(:), allocatable :: sz,szn	, dsz, dszn
            real(sp), dimension(:), allocatable :: a,b,c,r,u, pscs, b1, wgs, phi_ps, kgs
            integer(i4b) :: tdstart,tdend
            							 
        end type lsmg



		real(sp), parameter :: r_gas=8.314_sp, ma=29e-3_sp, ra=r_gas / ma, &
								t_stp=288._sp, p_stp=101300._sp, navog=6.022e23_sp, &
								sflux=1370._sp, mu1=1._sp/sqrt(3._sp)
		integer(i4b), parameter :: stl=20,charlen=15
		!>@brief
		!>variables for namelist radiation input
        type namelist_lsm
        	integer(i4b) :: start_year, start_mon,start_day, &
        					start_hour, start_min, start_sec
        	real(sp) :: lat_ref,lon_ref
        	real(sp) :: asymmetry_water
        	integer(i4b) :: quad_flag
        	character (len=200) :: albedo='Urban'
        	character (len=200) :: emissivity='Urban'
            integer(i4b) :: ns,nl,nsoil_lay
            character(len=charlen), dimension(stl) :: soil_types=['Sand','Sand','Sand',&
                'Sand','Sand','Sand','Sand','Sand','Sand','Sand','Sand','Sand',&
                'Sand','Sand','Sand','Sand','Sand','Sand','Sand','Sand']
            real(sp), dimension(20) :: lambda_read_s, soil_thickness
            real(sp), dimension(20) :: lambda_read_l
            real(sp) :: lambda_s_low, lambda_s_high, lambda_l_low, lambda_l_high
        end type namelist_lsm

		integer(i4b), parameter :: nsoil=11

		!>@brief
		!>variables for land surface
        type lsm_vars
            ! soil parameters
            real(sp), dimension(nsoil) :: b=[4.05, 4.38, 4.90, 5.30, 5.39, 7.12, &
                                            7.75, 8.52, 10.40, 10.40, 11.40]
            real(sp), dimension(nsoil) :: wgs=[0.395, 0.410, 0.435, 0.485, 0.451, &
                                            0.420, 0.477, 0.476, 0.426, 0.492, 0.482]
            real(sp), dimension(nsoil) :: wfc=[0.135, 0.150, 0.195, 0.255, 0.240, 0.255, &
                                            0.322, 0.325, 0.310, 0.370, 0.367]
            real(sp), dimension(nsoil) :: wwilt=[0.068,0.075,0.114,0.179,0.155, &
                                            0.175,0.218,0.250,0.219,0.283,0.286]
            real(sp), dimension(nsoil) :: phips=[-12.1,-9.0,-21.8,-78.6,-47.8,-29.9, &
                                            -35.6,-63.0,-15.3,-49.0,-40.5]
            real(sp), dimension(nsoil) :: kgs=[1.76e-4, 1.56e-4, 3.41e-5, 7.20e-6, & 
                        7.00e-6, 6.30e-6, 1.70e-6, 2.50e-6, 2.20e-6, 1.00e-6, 1.30e-6]
            real(sp), dimension(nsoil) :: pscs=[1.47e6, 1.41e6, 1.34e6, 1.27e6, 1.21e6, &
                            1.18e6, 1.32e6, 1.23e6, 1.18e6, 1.15e6, 1.09e6]

        	character(len=15), dimension(nsoil) :: soil_des = &
        		["Sand           ", &
        		 "Loamy sand     ", &
        		 "Sandy loam     ", &
        		 "Silt loam      ", &
        		 "Loam           ", &
        		 "Sandy clay loam", &
        		 "Silty clay loam", &
        		 "Clay loam      ", &
        		 "Sandy clay     ", &
        		 "Silty clay     ", &
        		 "Clay           "]
        
        end type lsm_vars
        
        
        type(lsm_vars) :: rsv

		! declare a namelist rad type
		type(namelist_lsm) :: nm3
		! declare a radiation grid type
		type(lsmg) :: lsmg1

		real(sp), parameter :: h_planck=6.6256e-34_sp, &
							c_light=2.9979e8_sp, &	
							k_boltz=1.38e-23_sp, temp_sun=5778._sp, &
						sigma_sb=2._sp*pi**5*k_boltz**4/(15._sp*h_planck**3*c_light**2), &
						fac_planck=2._sp*k_boltz**4/(h_planck**3*c_light**2)*pi, &
			N_a_0 = (p_stp / (t_stp*ra)) / ma * navog ! number of air molecules per m^3 at stp


	
		                    
		private
		public :: nm3, lsmg1,allocate_and_set_lsm, soil_solver

		contains
		
		
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!>@author
		!>Paul J. Connolly, The University of Manchester
		!>@brief
		!>initialise variables for radiation
		!>@param[in] 
		!>@param[inout] 
		subroutine allocate_and_set_lsm(ip,jp,kp, &
											l_h, r_h, &
											skp_nm, soil_thickness, soil_types, &
											skp, sz, szn, dsz,dszn, t, wg, tsurf, &
											tdend,a,b,c,r,u,pscs, &
											b1, wgs, phi_ps, kgs, &
											coords,dims, id, comm3d)
			use mpi
			implicit none
			integer(i4b), intent(in) :: l_h, r_h, skp_nm
			integer(i4b), intent(inout) :: ip, jp, kp, skp, tdend
			real(sp), dimension(:), allocatable, intent(inout) :: a,b,c,r,u,  sz, szn, &
			        dsz, dszn, pscs, b1, wgs, phi_ps, kgs
			real(sp), dimension(:,:,:), allocatable, intent(inout) :: t,wg
			real(sp), dimension(:,:), allocatable, intent(inout) :: tsurf
			real(sp), dimension(skp_nm), intent(in) :: soil_thickness
			character (len=charlen), dimension(stl), intent(in) :: soil_types
			
			integer(i4b), dimension(3), intent(inout) :: coords
			integer(i4b), dimension(3), intent(in) :: dims
			integer(i4b), intent(in) :: id, comm3d
			
		    integer(i4b), dimension(skp_nm) :: soil_index
			integer(i4b) :: AllocateStatus, i, j, error
			
! if the pe is not being used in the cartesian topology, do not use here
			if(id>=dims(1)*dims(2)*dims(3)) return 
		
				
            
			! find the soil layer type
			do j=1,skp_nm
                do i=1,nsoil
                    if(soil_types(j) == trim(rsv%soil_des(i))) then
                        exit
                    endif
                enddo
                soil_index(j)=i !rsv%emiss(i)
            enddo
			skp = skp_nm
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! allocate memory
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			tdend=skp
			allocate( a(1:tdend-1), STAT = AllocateStatus)
			if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
			allocate( b(1:tdend), STAT = AllocateStatus)
			if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
			allocate( c(1:tdend-1), STAT = AllocateStatus)
			if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
			allocate( r(1:tdend), STAT = AllocateStatus)
			if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
			allocate( u(1:tdend), STAT = AllocateStatus)
			if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
			
			allocate( pscs(1-l_h:skp+r_h), STAT = AllocateStatus)
			if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
			allocate( b1(1-l_h:skp+r_h), STAT = AllocateStatus)
			if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
			allocate( wgs(1-l_h:skp+r_h), STAT = AllocateStatus)
			if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
			allocate( phi_ps(1-l_h:skp+r_h), STAT = AllocateStatus)
			if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
			allocate( kgs(1-l_h:skp+r_h), STAT = AllocateStatus)
			if (AllocateStatus /= 0) STOP "*** Not enough memory ***"

			allocate( sz(1-l_h:skp+r_h), STAT = AllocateStatus)
			if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
			allocate( szn(1-l_h:skp+r_h), STAT = AllocateStatus)
			if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
			allocate( dsz(1-l_h:skp+r_h), STAT = AllocateStatus)
			if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
			allocate( dszn(1-l_h:skp+r_h), STAT = AllocateStatus)
			if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
			allocate( t(1-l_h:skp+r_h,1-l_h:jp+r_h,1-l_h:ip+r_h), STAT = AllocateStatus)
			if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
			allocate( wg(1-l_h:skp+r_h,1-l_h:jp+r_h,1-l_h:ip+r_h), STAT = AllocateStatus)
			if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
			allocate( tsurf(1-l_h:jp+r_h,1-l_h:ip+r_h), STAT = AllocateStatus)
			if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! grid soil properties                                                       !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			do i=1,skp
			    pscs(i)  =rsv%pscs(soil_index(i))
			    b1(i)    =rsv%b(soil_index(i))
			    wgs(i)   =rsv%wgs(soil_index(i))
			    phi_ps(i)=rsv%phips(soil_index(i))
			    kgs(i)   =rsv%kgs(soil_index(i))
			enddo
			pscs(0)=pscs(1)
			pscs(skp+1)=pscs(skp)
			b1(0)=b1(1)
			b1(skp+1)=b1(skp)
			wgs(0)=wgs(1)
			wgs(skp+1)=wgs(skp)
			phi_ps(0)=phi_ps(1)
			phi_ps(skp+1)=phi_ps(skp)
			kgs(0)=kgs(1)
			kgs(skp+1)=kgs(skp)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
			! position of grid edges
            sz=0.0_sp
            do i=skp,1,-1
                sz(i-1)=sz(i)-soil_thickness(skp-i+1)
            enddo
            sz(skp+1)=sz(skp)+(sz(skp)-sz(skp-1))
            !---checked

            ! position of grid nodes
            do i=1,skp+1
                ! sn1 = (s0.5+s1.5)/2
                szn(i)=0.5_sp*(sz(i)+sz(i-1))
            enddo
            szn(0)=0.5_sp*(sz(0)+(sz(0)-(szn(1)-szn(0))))
            !---checked
            
            ! thickness
            do i=0,skp
                dsz(i)=sz(i+1)-sz(i)
            enddo
            dsz(skp+1)=dsz(skp)
            ! thickness
            do i=0,skp
                dszn(i)=szn(i+1)-szn(i)
            enddo
            dszn(skp+1)=dszn(skp)
            
            t=283.15_sp
            wg=1.e-8_sp
		end subroutine allocate_and_set_lsm
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!>@author
		!>Paul J. Connolly, The University of Manchester
		!>@brief
		!>solve land surface over one time-step
		!>@param[in] 
		!>@param[inout] 
		subroutine soil_solver(ip,jp,skp,l_h, r_h, tdend,a,b,c,r,u, pscs,b1,wgs,phi_ps, &
		                        kgs, &
		                        dt,t,wg,tsurf,sz,szn,dsz,dszn, &
                                coords,dims, id, comm3d)
			use mpi
			use nr, only : tridag
			implicit none
			integer(i4b), intent(in) :: l_h, r_h
			integer(i4b), intent(in) :: ip, jp, skp, tdend
			real(sp), intent(in) :: dt
			real(sp), dimension(1-l_h:skp+r_h), intent(in) :: sz,szn,dsz,dszn, pscs, &
			                        b1, wgs, phi_ps, kgs
			real(sp), dimension(tdend), intent(inout) :: b,r,u
			real(sp), dimension(tdend-1), intent(inout) :: a,c
			real(sp), dimension(1-l_h:skp+r_h,1-l_h:jp+r_h,1-l_h:ip+r_h), &
			    intent(inout) :: t,wg
			real(sp), dimension(1-l_h:jp+r_h,1-l_h:ip+r_h), &
			    intent(inout) :: tsurf
			
			integer(i4b), dimension(3), intent(in) :: coords
			integer(i4b), dimension(3), intent(in) :: dims
			integer(i4b), intent(in) :: id, comm3d
			
			integer(i4b) :: i,j,k
			real(sp) :: alpha_i, beta_i, beta_ip1, kval, flux1, flux2, dp05, dm05, &
			    gamma_i, phip, kp05, km05
			real(sp), dimension(1-l_h:skp+r_h) :: dg, kg, ks, pgcg
			

! 			kval=100._sp
			flux1=10._sp
			flux2=2.78e-6_sp
! if the pe is not being used in the cartesian topology, do not use here
			if(id>=dims(1)*dims(2)*dims(3)) return 


            do i=1,ip
                do j=1,jp
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! Temperature solver - equation 8.91 in Jacobson                     !                    
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! first set coefficients
                    do k=1-l_h,skp+r_h
                        ! equation 8.92 on the grid
                        phip=phi_ps(k)*(wgs(k)/max(wg(k,j,i),1.e-8_sp))**b1(k)
                        ks(k)=max(418._sp*exp(-log10(abs(phip)-2.7_sp)),0.172_sp)
                        ! equation 8.94 - note that the soil has air gaps, which
                        !  are neglected in the first term
                        pgcg(k)=(1._sp-wgs(k))*pscs(k)+wg(k,j,i)*4.2e6_sp
                    enddo   
                    
                    ! a is the sub diagonal, b is diagonal, c is super diagonal
                    do k = 1,tdend-1
                        kp05=(ks(k+1)+ks(k))*0.5_sp
                        km05=(ks(k)+ks(k-1))*0.5_sp
                        alpha_i= dt*kp05/(2._sp*pgcg(k)*dsz(k-1)*dszn(k))
                        beta_i = dt*km05/(2._sp*pgcg(k)*dsz(k-1)*dszn(k-1))
                        beta_ip1 = dt*kp05/(2._sp*pscs(i+1)*dsz(k)*dszn(k))
                        ! super-diagonal
                        c(k)   = -alpha_i
                        ! diagonal
                        b(k)   = 1._sp+alpha_i+beta_i
                        ! sub-diagonal
                        a(k)   = -beta_ip1
                    enddo
            
            
                    k=tdend
                    kp05=(ks(k+1)+ks(k))*0.5_sp
                    km05=(ks(k)+ks(k-1))*0.5_sp
                    alpha_i= dt*kp05/(2._sp*pgcg(k)*dsz(k-1)*dszn(k))
                    beta_i = dt*km05/(2._sp*pgcg(k+1)*dsz(k)*dszn(k))
                    ! diagonal
                    b(k)   = 1._sp+alpha_i+beta_i
            
                    ! now set the RHS
                    do k=1,tdend
                        kp05=(ks(k+1)+ks(k))*0.5_sp
                        km05=(ks(k)+ks(k-1))*0.5_sp
                        alpha_i= dt*kp05/(2._sp*pgcg(k)*dsz(k-1)*dszn(k))
                        beta_i = dt*km05/(2._sp*pgcg(k)*dsz(k-1)*dszn(k-1))
                        beta_ip1 = dt*kp05/(2._sp*pgcg(k+1)*dsz(k)*dszn(k))
                        r(k)   = alpha_i*t(k+1,j,i)+(1._sp-alpha_i-beta_i)*t(k,j,i)+ &
                                 beta_i*t(k-1,j,i)    
                        ! uncomment for backward euler
                        ! r(k)=t(k,j,i)        
                    enddo
                    ! top BC - i.e. governed by radiation and heat fluxes
                    b(tdend)=1._sp
                    a(tdend-1)=0._sp
                    r(tdend)=t(tdend,1,1)+dt/(pgcg(tdend)*dsz(tdend-1))*(flux1)
            
            
                    ! bottom BC  - no flux
                    b(1)=1._sp
                    c(1)=0._sp
                    r(1)=t(1, j,i) ! should mean rate of change is zero
            
                    ! now solve
                    call tridag(a,b,c,r,u)
                    t(1:tdend,j,i)=u
                    t(0,j,i)=t(1,j,i)
                    ! surface temperature
                    kval=ks(tdend)
                    tsurf(j,i) = t(tdend,j,i)-flux1 * ((sz(tdend)-szn(tdend))/kval)
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    
                    
                    
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! Moisture solver - equation 8.95 in Jacobson                        !                    
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! first set coefficients
                    do k=1-l_h,skp+r_h
                        ! equation 8.97 on the grid
                        dg(k)=-b1(k)*kgs(k)*phi_ps(k)/ &
                                wgs(k)*(wg(k,j,i)/wgs(k))**(b1(k)+2._sp)
                        kg(k)=kgs(k)*(wg(k,j,i)/wgs(k))**(2._sp*b1(k)+3._sp)
                    enddo   
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

                    ! a is the sub diagonal, b is diagonal, c is super diagonal
                    do k = 1,tdend-1
                        dp05=(dg(k+1)+dg(k))*0.5_sp
                        dm05=(dg(k)+dg(k-1))*0.5_sp
                        alpha_i= dt*dp05/(2._sp*dsz(k-1)*dszn(k))
                        beta_i = dt*dm05/(2._sp*dsz(k-1)*dszn(k-1))
                        beta_ip1 = dt*dp05/(2._sp*dsz(k)*dszn(k))
                        ! super-diagonal
                        c(k)   = -alpha_i
                        ! diagonal
                        b(k)   = 1._sp+alpha_i+beta_i
                        ! sub-diagonal
                        a(k)   = -beta_ip1
                    enddo
                    k=tdend
                    dp05=(dg(k+1)+dg(k))*0.5_sp
                    dm05=(dg(k)+dg(k-1))*0.5_sp
                    alpha_i= dt*dp05/(2._sp*dsz(k-1)*dszn(k))
                    beta_i = dt*dm05/(2._sp*dsz(k)*dszn(k))
                    ! diagonal
                    b(k)   = 1._sp+alpha_i+beta_i
            
                    ! now set the RHS
                    do k=1,tdend
                        dp05=(dg(k+1)+dg(k))*0.5_sp
                        dm05=(dg(k)+dg(k-1))*0.5_sp
                        alpha_i= dt*dp05/(2._sp*dsz(k-1)*dszn(k))
                        beta_i = dt*dm05/(2._sp*dsz(k-1)*dszn(k-1))
                        beta_ip1 = dt*dp05/(2._sp*dsz(k)*dszn(k))
                        gamma_i=dt*(kg(k+1)-kg(k-1))/(dszn(k-1)+dszn(k))
                        r(k)   = alpha_i*wg(k+1,j,i)+(1._sp-alpha_i-beta_i)*wg(k,j,i)+ &
                                 beta_i*wg(k-1,j,i) + gamma_i
                        ! uncomment for backward euler
                        ! r(k)=wg(k,j,i) + gamma_i      
                    enddo
                    ! top BC - i.e. governed by radiation and heat fluxes
                    b(tdend)=1._sp
                    a(tdend-1)=0._sp
                    r(tdend)=wg(tdend,1,1)+dt/(dsz(tdend-1))*(flux2)
            
            
                    ! bottom BC  - no flux
                    b(1)=1._sp
                    c(1)=0._sp
                    r(1)=wg(1, j,i) ! should mean rate of change is zero
            
                    ! now solve
                    call tridag(a,b,c,r,u)
                    wg(1:tdend,j,i)=min(u,wgs(1:tdend))
                    wg(0,j,i)=wg(1,j,i)
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       
                enddo
            enddo
            
            
!            print *,wg(1:tdend,1,1)
            !print *,tsurf(1,1),t(tdend,1,1), wg(tdend,1,1)
            !stop
            
		end subroutine soil_solver
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





	end module lsm
	