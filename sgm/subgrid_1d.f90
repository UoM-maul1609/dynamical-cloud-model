	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>subgrid code for the simple cloud model
    module subgrid_1d
    use numerics_type
    
    private
    public :: calculate_subgrid_1d,advance_fields_1d
    
    ! variables / parameters used in subgrid model
    ! note when kvon=0.4, prn should be 0.95 - see Jacobson (page 242 - from Hogstrom)
    real(wp), parameter :: &!prn=0.7_wp, &
                        prn=0.95_wp, suba=1._wp/prn, subb=40._wp, subc=16._wp, &
                        subf=suba, subg=1.2_wp, subh=0._wp, subr=4._wp, ric=0.25_wp, &
                        grav=9.81_wp, small=1.e-10_wp, kvon=0.4_wp
    
	contains
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Calculate terms for Smagorinski model                                              !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculate viscosity, and strain-rate
	!>@param[in] dt
	!>@param[in] z,zn,dz, rhoa
	!>@param[in] dzn, rhoan
	!>@param[in] ip,kp,l_h,r_h, nq
	!>@param[in] u, zu
	!>@param[in] w, zw
	!>@param[in] th,q, lamsq, z0,z0th
	!>@param[inout] vism, vist,viss,strain, tau11, tau33, tau13, 
	!>                su, sw,sth,sq
	subroutine calculate_subgrid_1d(dt,z,zn,dz,rhoa, theta,&
	                                dzn,rhoan, thetan,&
	                                kp,nq,l_h,r_h,u,w,zu, zw, &
	                                th,q, &
	                                lamsq, z0,z0th, &
	                                vism,vist, viss,strain, &
	                                tau11,tau33,&
	                                tau13, &
	                                su,sw,sth,sq)
	use numerics_type
	
	implicit none

	real(wp), intent(in) :: dt
	integer(i4b), intent(in) :: kp, l_h, r_h, nq
	real(wp), dimension(-r_h+1:kp+r_h), &
		intent(in) :: u, zu
	real(wp), dimension(-l_h+1:kp+r_h), &
		intent(inout) :: w, zw, th
	real(wp), dimension(-r_h+1:kp+r_h,nq), &
		intent(in) :: q
	real(wp), dimension(-r_h+1:kp+r_h,nq), &
		intent(inout) :: viss,sq
	real(wp), dimension(-r_h+1:kp+r_h), &
		intent(inout) :: vism, vist, &
		                strain, tau11, tau33, tau13, &
		                su, sw, sth
	real(wp), dimension(-l_h+1:kp+r_h), intent(in) :: z, zn,&
	        dz, dzn, rhoa, rhoan, theta, thetan
    real(wp), dimension(-l_h+1:kp+r_h), intent(in) :: lamsq
    real(wp), intent(in) :: z0,z0th
	
	! locals
	integer(i4b) :: i,j,k,n, ktmp,ktmp1
	real(wp), dimension(-r_h+1:kp+r_h) :: rip, fm, fh

    ! u is on x levels, v on y levels, w on z levels. 
    ! so distance between u(i) and u(i+1) is dx(i)
    !    w is on xn levels
    !    distance between w(i) and w(i+1) is dxn(i)
    !    u is on zn levels
    !    distance between u(k) and u(k+1) is dzn(k)
    ! strain-rate on w-points
!$omp simd
    do k=1,kp
        ! 2*s33*s33
        strain(k)=((w(k)-w(k-1))/dz(k-1))**2 + &
                      ((w(k+1)-w(k))/dz(k))**2 
                      
        ! 2*s13*s13 - du/dz+dw/dx - averaging over 2 points
        strain(k)=strain(k)+ 0.5_wp * &
           (((u(k+1)-u(k))/dzn(k))**2 + &
           ((u(k+1)-u(k))/dz(k))**2)
    enddo
!$omp end simd

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Richardson number and functions                                                    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call richardson(strain,thetan,theta,dzn,th,kp,l_h,r_h,rip,fm,fh)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    

    ! viscosity on w-points
!$omp simd
    do k=1,kp
        vism(k)=lamsq(k)*sqrt(strain(k)) *fm(k)
    enddo
!$omp end simd

!$omp simd
    do k=1,kp
        ! viscosity for theta
        vist(k)=lamsq(k)*sqrt(strain(k)) *fh(k)*rhoa(k)
    enddo
!$omp end simd
!     vist=vism
    vism(kp-1:kp)=0._wp
    vist(kp-1:kp)=0._wp
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! surface layer - monin-obukhov similarity theory                                    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call monin_obukhov(z0,z0th,theta,thetan,z,zn,th,u,vism,vist,kp,l_h,r_h)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! full exchange
    !vism(0)=0._wp
    vism(kp+1)=0._wp

    !vist(0)=0._wp
    vist(kp+1)=0._wp


    ! viscosity for q-fields
    do n=1,nq
        viss(:,n)=vist(:)
    enddo

    ! calculate elements of tau on p-points (in vertical)
    do k=1,kp
        tau11(k)=0._wp
        tau33(k)=rhoan(k)*1._wp*(vism(k)+vism(k-1))* &
                    (zw(k)-zw(k-1))/dz(k-1)
                    
    enddo
    

    ! calculate elements of tau on w-points (in vertical)
    do k=1,kp
        ! on i+1/2
        tau13(k)=rhoa(k)*0.5_wp*(vism(k)+vism(k))* &
            ( (zu(k+1)-zu(k))/dzn(k) )
            
    enddo
    
    
    ! full exchange
    tau11(0)=0._wp
    tau11(kp+1)=0._wp

    tau33(0)=0._wp
    tau33(kp+1)=0._wp

    tau13(0)=0._wp
    tau13(kp+1)=0._wp

        

    ktmp=kp-1
    ktmp1=2
    ! calculate sources due to viscosity
!     tau23(1,:,:)=max(tau23(2,:,:),tau23(1,:,:))
    do k=1,ktmp
        sw(k)=sw(k)+ &
           ((tau33(k+1)-tau33(k))/dzn(k)  ) !/ rhoa(k)             
    enddo
    do k=ktmp1,ktmp
        ! note the sources are sources of u*rhoa
        su(k)=su(k)+ &
           ((tau13(k)-tau13(k-1))/dz(k-1)) !/ rhoan(k)             
            
    enddo
    
    ! full exchange    
    su(0)=0._wp
    su(kp+1)=0._wp
    
    sw(0)=0._wp
    sw(kp+1)=0._wp
    
    ! calculate th sources due to viscosity
    ! diffusion along i calculated via i+1/2 and i-1/2 terms
    sth=0._wp
    do k=1,kp
        ! note, vist has already been multiplied by rhoa
        sth(k)=sth(k)+ &
         ((vist(k-1)* (th(k-1)-th(k))/dz(k-1) - &
                 vist(k)* (th(k)-th(k+1))/dz(k) ) / dz(k-1)  ) &
                 / rhoan(k)                     
                
            
    enddo
    
    ! full exchange 
    sth(0)=0._wp
    sth(kp+1)=0._wp
    sq=0._wp
    ! q-fields
    do n=1,nq
        ! calculate q sources due to viscosity
        ! diffusion along i calculated via i+1/2 and i-1/2 terms
        do k=1,kp
            ! note, vist has already been multiplied by rhoa
            sq(k,n)=sq(k,n)+ &
             ((vist(k-1)* (q(k-1,n)-q(k,n))/dz(k-1) - &
                     vist(k)* (q(k,n)-q(k+1,n))/dz(k) ) / dz(k-1)  ) &
                     / rhoan(k)                     
                
            
        enddo
    
        ! full exchange  
        sq(0,n)=0._wp
        sq(kp+1,n)=0._wp
    enddo


    
	end subroutine calculate_subgrid_1d
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Calculate surface layer parameterisation                                           !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculate the momentum and potential temperature scales
	!> see page 243, Jacobson 
	!> "Noniterative parameterisation for momentum and potential temperature scales"
	!>@param[in] z0,z0th,theta,thetan,z,zn,th,u,ip,kp,l_h,r_h
	!>@param[inout] vism,vist
    subroutine monin_obukhov(z0,z0th,theta,thetan,z,zn,th,u,vism,vist,kp,l_h,r_h)
	use numerics_type
    
    implicit none
    integer(i4b), intent(in) :: kp, l_h, r_h
    real(wp), intent(in), dimension(-r_h+1:kp+r_h) :: &
            th, u
    real(wp), intent(inout), dimension(-r_h+1:kp+r_h) :: &
            vism,vist
    real(wp), intent(in), dimension(-r_h+1:kp+r_h) :: thetan, zn, theta,z
    real(wp), intent(in) :: z0, z0th
            
    ! locals
    integer(i4b) :: i,j
    real(wp), dimension(1) :: vmag,rib,gm,gh, &
                                                            ustar,thstar
    
    
    
    ! wind speed in first layer above surface
    vmag(1)=0.5_wp*sqrt((u(1)+u(1))**2  )+small
    
    ! bulk richardson number
    ! equation 8.39 (Jacobson pp 243)
    rib(1)=grav*(th(1)+thetan(1)-theta(0))*(z(1)-z0)**2 / &
        (theta(0)*vmag(1)**2*(z(1)-z0th)) 
    
    ! equation 8.41 (Jacobson pp 243)
    if(rib(1) .le. 0._wp) then
        gm(1)=1._wp-9.4_wp*rib(1)/ &
            (1._wp+70._wp*kvon**2*sqrt(abs(rib(1)*z(1)/z0))/(log(z(1)/z0))**2)
        gh(1)=1._wp-9.4_wp*rib(1)/ &
            (1._wp+50._wp*kvon**2*sqrt(abs(rib(1)*z(1)/z0))/(log(z(1)/z0))**2)
    else
        gm(1)=1._wp/(1._wp+4.7_wp*rib(1))**2
        gh(1)=1._wp/(1._wp+4.7_wp*rib(1))**2
    endif
    
    ! equation 8.40 (Jacobson pp 243)
    ustar(1)=kvon*vmag(1)/log(z(1)/z0)*sqrt(gm(1))
    thstar(1)=kvon**2*vmag(1)*(th(1)+thetan(1)-theta(0)) / &
        (ustar(1)*prn*log(z(1)/z0)**2)*gh(1)
        
        
    ! see equation 8.46, page 245 (Jacobson)
    vism(0)=ustar(1)**2*(z(1)-z0)/vmag(1)
    
    vist(0)=ustar(1)*thstar(1)/ &
        ((th(1)+thetan(1)-theta(0))/(z(0)-z0th) + small )
    
    
    

    end subroutine monin_obukhov
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Calculate Richardson number (used to determine viscosity)                          !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculate the richardson number and richardson number dependent functions
	!>@param[in] strain, thetan, theta, dzn, th, ip, kp, l_h, r_h
	!>@param[inout] rip, fm, fh
    subroutine richardson(strain,thetan,theta,dzn,th,kp,l_h,r_h,rip,fm,fh)
	use numerics_type
    
    implicit none
    integer(i4b), intent(in) :: kp, l_h, r_h
    real(wp), intent(in), dimension(-r_h+1:kp+r_h) :: &
            strain, th
    real(wp), intent(in), dimension(-r_h+1:kp+r_h) :: thetan, theta, dzn
    real(wp), intent(inout), dimension(-r_h+1:kp+r_h) :: &
            rip,fm, fh
            
    ! locals
    integer(i4b) :: i,j,k
    
    
    ! calculate bulk richardson number
    do k=1,kp
        rip(k)=grav*((th(k+1))/thetan(k+1) - &
                    (th(k))/thetan(k)) / dzn(k) / (strain(k)+1.e-15)
    enddo
    ! richardson number functions
    do k=1,kp
        if(rip(k) .lt. 0._wp) then
            fm(k)=(1._wp-subc*rip(k))**0.5_wp
            fh(k)=suba*(1._wp-subb*rip(k))**0.5_wp
        elseif((rip(k).ge.0._wp).and.(rip(k).lt.ric)) then
            fm(k)=(1._wp-rip(k)/ric)**subr*(1._wp-subh*rip(k))
            fh(k)=subf*(1._wp-rip(k)/ric)**subr*(1._wp-subg*rip(k))
        else
            fm(k)=0._wp
            fh(k)=0._wp
        endif

    enddo

    end subroutine richardson
    
    
    
    
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Advance fields due to diffusion term                                               !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>just add the source term to the current value of the field
	!>@param[in] su,sw,sq,sth,ipp, kpp, nq, l_h,r_h,dt, rhoa, rhoan
	!>@param[inout] q,th,u,w
    subroutine advance_fields_1d(dt,tu,tw, &
        zu,zw,su,sw,q,sq,th,sth,rhoa,rhoan,kp,nq,l_h,r_h)
	use numerics_type
    
    implicit none
    integer(i4b), intent(in) :: kp, nq, l_h, r_h
    real(wp), intent(in), dimension(-r_h+1:kp+r_h) :: &
        su,sw,sth
    real(wp), intent(in), dimension(-r_h+1:kp+r_h) :: &
        rhoa,rhoan
    real(wp), intent(inout), dimension(-r_h+1:kp+r_h) :: &
        tu,tw,zu,zw,th
    real(wp), intent(in), dimension(-r_h+1:kp+r_h,nq) :: sq
    real(wp), intent(inout), dimension(-r_h+1:kp+r_h,nq) :: q
    real(wp), intent(in) :: dt
    ! locals
    integer(i4b) :: i,j,k,n
    
    do k=1-l_h,kp+r_h
        tu(k)=zu(k)+su(k)*dt*2._wp/rhoan(k)
    enddo
    
    
    do k=1-l_h,kp+r_h
        tw(k)=zw(k)+sw(k)*dt*2._wp/rhoa(k)
    enddo
    
    do k=1-l_h,kp+r_h
        th(k)=th(k)+sth(k)*dt
    enddo
    
    do n=1,nq
        do k=1-l_h,kp+r_h
            q(k,n)=q(k,n)+sq(k,n)*dt
        enddo
    enddo    
    
    end subroutine advance_fields_1d
    
    	
	

    end module subgrid_1d
