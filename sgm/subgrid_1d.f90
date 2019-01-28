	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>subgrid code for the simple cloud model
    module subgrid_1d
    use nrtype
    
    private
    public :: calculate_subgrid_1d,advance_fields_1d
    
    ! variables / parameters used in subgrid model
    ! note when kvon=0.4, prn should be 0.95 - see Jacobson (page 242 - from Hogstrom)
    real(sp), parameter :: &!prn=0.7_sp, &
                        prn=0.95_sp, suba=1._sp/prn, subb=40._sp, subc=16._sp, &
                        subf=suba, subg=1.2_sp, subh=0._sp, subr=4._sp, ric=0.25_sp, &
                        grav=9.81_sp, small=1.e-10_sp, kvon=0.4_sp
    
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
	use nrtype
	
	implicit none

	real(sp), intent(in) :: dt
	integer(i4b), intent(in) :: kp, l_h, r_h, nq
	real(sp), dimension(-r_h+1:kp+r_h), &
		intent(in) :: u, zu
	real(sp), dimension(-l_h+1:kp+r_h), &
		intent(inout) :: w, zw, th
	real(sp), dimension(-r_h+1:kp+r_h,nq), &
		intent(in) :: q
	real(sp), dimension(-r_h+1:kp+r_h,nq), &
		intent(inout) :: viss,sq
	real(sp), dimension(-r_h+1:kp+r_h), &
		intent(inout) :: vism, vist, &
		                strain, tau11, tau33, tau13, &
		                su, sw, sth
	real(sp), dimension(-l_h+1:kp+r_h), intent(in) :: z, zn,&
	        dz, dzn, rhoa, rhoan, theta, thetan
    real(sp), dimension(-l_h+1:kp+r_h), intent(in) :: lamsq
    real(sp), intent(in) :: z0,z0th
	
	! locals
	integer(i4b) :: i,j,k,n, ktmp,ktmp1
	real(sp), dimension(-r_h+1:kp+r_h) :: rip, fm, fh

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
        strain(k)=strain(k)+ 0.5_sp * &
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
    vism(kp-1:kp)=0._sp
    vist(kp-1:kp)=0._sp
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! surface layer - monin-obukhov similarity theory                                    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call monin_obukhov(z0,z0th,thetan,zn,th,u,vism,vist,kp,l_h,r_h)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! full exchange
    vism(0)=0._sp
    vism(kp+1)=0._sp

    vist(0)=0._sp
    vist(kp+1)=0._sp


    ! viscosity for q-fields
    do n=1,nq
        viss(:,n)=vist(:)
    enddo

    ! calculate elements of tau on p-points (in vertical)
    do k=1,kp
        tau11(k)=0._sp
        tau33(k)=rhoan(k)*1._sp*(vism(k)+vism(k-1))* &
                    (zw(k)-zw(k-1))/dz(k-1)
                    
    enddo
    

    ! calculate elements of tau on w-points (in vertical)
    do k=1,kp
        ! on i+1/2
        tau13(k)=rhoa(k)*0.5_sp*(vism(k)+vism(k))* &
            ( (zu(k+1)-zu(k))/dzn(k) )
            
    enddo
    
    
    ! full exchange
    tau11(0)=0._sp
    tau11(kp+1)=0._sp

    tau33(0)=0._sp
    tau33(kp+1)=0._sp

    tau13(0)=0._sp
    tau13(kp+1)=0._sp

        

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
    su(0)=0._sp
    su(kp+1)=0._sp
    
    sw(0)=0._sp
    sw(kp+1)=0._sp
    
    ! calculate th sources due to viscosity
    ! diffusion along i calculated via i+1/2 and i-1/2 terms
    sth=0._sp
    do k=1,kp
        ! note, vist has already been multiplied by rhoa
        sth(k)=sth(k)+ &
         ((vist(k-1)* (th(k-1)-th(k))/dz(k-1) - &
                 vist(k)* (th(k)-th(k+1))/dz(k) ) / dz(k-1)  ) &
                 / rhoan(k)                     
                
            
    enddo
    
    ! full exchange 
    sth(0)=0._sp
    sth(kp+1)=0._sp
    sq=0._sp
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
        sq(0,n)=0._sp
        sq(kp+1,n)=0._sp
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
	!>@param[in] z0,z0th,thetan,zn,th,u,ip,kp,l_h,r_h
	!>@param[inout] vism,vist
    subroutine monin_obukhov(z0,z0th,thetan,zn,th,u,vism,vist,kp,l_h,r_h)
	use nrtype
    
    implicit none
    real(sp), intent(in), dimension(-r_h+1:kp+r_h) :: &
            th, u
    real(sp), intent(inout), dimension(-r_h+1:kp+r_h) :: &
            vism,vist
    real(sp), intent(in), dimension(-r_h+1:kp+r_h) :: thetan, zn
    integer(i4b), intent(in) :: kp, l_h, r_h
    real(sp), intent(in) :: z0, z0th
            
    ! locals
    integer(i4b) :: i,j
    real(sp), dimension(1) :: vmag,rib,gm,gh, &
                                                            ustar,thstar
    
    
    
    ! wind speed in first layer above surface
    vmag(1)=0.5_sp*sqrt((u(2)+u(2))**2  )+small
    
    ! bulk richardson number
    ! equation 8.39 (Jacobson pp 243)
    rib(1)=grav*(th(2)+thetan(2)-th(1)-thetan(1))*(zn(2)-z0)**2 / &
        (thetan(1)*vmag(1)**2*(zn(2)-z0th)) 
    
    ! equation 8.41 (Jacobson pp 243)
    if(rib(1) .le. 0._sp) then
        gm(1)=1._sp-9.4_sp*rib(1)/ &
            (1._sp+70._sp*kvon**2*sqrt(abs(rib(1)*zn(2)/z0))/(log(zn(2)/z0))**2)
        gh(1)=1._sp-9.4_sp*rib(1)/ &
            (1._sp+50._sp*kvon**2*sqrt(abs(rib(1)*zn(2)/z0))/(log(zn(2)/z0))**2)
    else
        gm(1)=1._sp/(1._sp+4.7_sp*rib(1))**2
        gh(1)=1._sp/(1._sp+4.7_sp*rib(1))**2
    endif
    
    ! equation 8.40 (Jacobson pp 243)
    ustar(1)=kvon*vmag(1)/log(zn(2)/z0)*sqrt(gm(1))
    thstar(1)=kvon**2*vmag(1)*(th(2)+thetan(2)-th(1)-thetan(1)) / &
        (ustar(1)*prn*log(zn(2)/z0)**2)*gh(1)
        
        
    ! see equation 8.46, page 245 (Jacobson)
    vism(1)=ustar(1)**2*(zn(2)-z0)/vmag(1)
    
    vist(1)=ustar(1)*thstar(1)/ &
        ((th(2)+thetan(2)-th(1)-thetan(1))/(zn(2)-z0th) + small )
    
    
    

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
	use nrtype
    
    implicit none
    real(sp), intent(in), dimension(-r_h+1:kp+r_h) :: &
            strain, th
    real(sp), intent(in), dimension(-r_h+1:kp+r_h) :: thetan, theta, dzn
    real(sp), intent(inout), dimension(-r_h+1:kp+r_h) :: &
            rip,fm, fh
    integer(i4b), intent(in) :: kp, l_h, r_h
            
    ! locals
    integer(i4b) :: i,j,k
    
    
    ! calculate bulk richardson number
    do k=1,kp
        rip(k)=grav*((th(k+1))/thetan(k+1) - &
                    (th(k))/thetan(k)) / dzn(k) / (strain(k)+1.e-15)
    enddo
    ! richardson number functions
    do k=1,kp
        if(rip(k) .lt. 0._sp) then
            fm(k)=(1._sp-subc*rip(k))**0.5_sp
            fh(k)=suba*(1._sp-subb*rip(k))**0.5_sp
        elseif((rip(k).ge.0._sp).and.(rip(k).lt.ric)) then
            fm(k)=(1._sp-rip(k)/ric)**subr*(1._sp-subh*rip(k))
            fh(k)=subf*(1._sp-rip(k)/ric)**subr*(1._sp-subg*rip(k))
        else
            fm(k)=0._sp
            fh(k)=0._sp
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
	use nrtype
    
    implicit none
    real(sp), intent(in), dimension(-r_h+1:kp+r_h) :: &
        su,sw,sth
    real(sp), intent(in), dimension(-r_h+1:kp+r_h) :: &
        rhoa,rhoan
    real(sp), intent(inout), dimension(-r_h+1:kp+r_h) :: &
        tu,tw,zu,zw,th
    real(sp), intent(in), dimension(-r_h+1:kp+r_h,nq) :: sq
    real(sp), intent(inout), dimension(-r_h+1:kp+r_h,nq) :: q
    integer(i4b), intent(in) :: kp, nq, l_h, r_h
    real(sp), intent(in) :: dt
    ! locals
    integer(i4b) :: i,j,k,n
    
    do k=1-l_h,kp+r_h
        tu(k)=zu(k)+su(k)*dt*2._sp/rhoan(k)
    enddo
    
    
    do k=1-l_h,kp+r_h
        tw(k)=zw(k)+sw(k)*dt*2._sp/rhoa(k)
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
