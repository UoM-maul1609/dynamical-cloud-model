	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>subgrid code for the thermal cloud model
    module subgrid_2d
    use nrtype
    
    private
    public :: calculate_subgrid_2d,advance_fields_2d
    
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
	!>@param[in] z,zn,dx,dz, rhoa
	!>@param[in] dxn,dzn, rhoan
	!>@param[in] ip,kp,l_h,r_h, nq
	!>@param[in] u, zu
	!>@param[in] w, zw
	!>@param[in] th,q, lamsq, z0,z0th
	!>@param[inout] vism, vist,viss,strain, tau11, tau33, tau13, 
	!>                su, sw,sth,sq
	subroutine calculate_subgrid_2d(dt,z,zn,dx,dz,rhoa, theta,&
	                                dxn,dzn,rhoan, thetan,&
	                                ip,kp,nq,l_h,r_h,u,w,zu, zw, &
	                                th,q, &
	                                lamsq, z0,z0th, &
	                                vism,vist, viss,strain, &
	                                tau11,tau33,&
	                                tau13, &
	                                su,sw,sth,sq)
	use nrtype
	
	implicit none

	real(sp), intent(in) :: dt
	integer(i4b), intent(in) :: ip, kp, l_h, r_h, nq
	real(sp), dimension(-r_h+1:kp+r_h,-l_h+1:ip+r_h), &
		intent(in) :: u, zu
	real(sp), dimension(-l_h+1:kp+r_h,-r_h+1:ip+r_h), &
		intent(inout) :: w, zw, th
	real(sp), dimension(-r_h+1:kp+r_h,-r_h+1:ip+r_h,nq), &
		intent(in) :: q
	real(sp), dimension(-r_h+1:kp+r_h,-r_h+1:ip+r_h,nq), &
		intent(inout) :: viss,sq
	real(sp), dimension(-r_h+1:kp+r_h,-r_h+1:ip+r_h), &
		intent(inout) :: vism, vist, &
		                strain, tau11, tau33, tau13, &
		                su, sw, sth
	real(sp), dimension(-l_h+1:ip+r_h), intent(in) :: dx, dxn
	real(sp), dimension(-l_h+1:kp+r_h), intent(in) :: z, zn,&
	        dz, dzn, rhoa, rhoan, theta, thetan
    real(sp), dimension(-l_h+1:kp+r_h), intent(in) :: lamsq
    real(sp), intent(in) :: z0,z0th
	
	! locals
	integer(i4b) :: i,j,k,n, ktmp,ktmp1
	real(sp), dimension(-r_h+1:kp+r_h,-r_h+1:ip+r_h) :: rip, fm, fh

    ! u is on x levels, v on y levels, w on z levels. 
    ! so distance between u(i) and u(i+1) is dx(i)
    !    w is on xn levels
    !    distance between w(i) and w(i+1) is dxn(i)
    !    u is on zn levels
    !    distance between u(k) and u(k+1) is dzn(k)
    ! strain-rate on w-points
!$omp simd
    do i=1,ip
        do k=1,kp
            ! 2*s11*s11
            strain(k,i)=((u(k+1,i)-u(k+1,i-1))/dx(i-1))**2 + &
                          ((u(k,i)-u(k,i-1))/dx(i-1))**2
            ! 2*s33*s33
            strain(k,i)=strain(k,i)+ &
                          ((w(k,i)-w(k-1,i))/dz(k-1))**2 + &
                          ((w(k+1,i)-w(k,i))/dz(k))**2 
                          
            ! 2*s13*s13 - du/dz+dw/dx - averaging over 2 points
            strain(k,i)=strain(k,i)+ 0.5_sp * &
               (((u(k+1,i)-u(k,i))/dzn(k)+(w(k,i+1)-w(k,i))/dxn(i))**2 + &
               ((u(k+1,i-1)-u(k,i-1))/dz(k)+(w(k,i)-w(k,i-1))/dxn(i-1))**2)
        enddo
    enddo
!$omp end simd

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Richardson number and functions                                                    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call richardson(strain,thetan,theta,dzn,th,ip,kp,l_h,r_h,rip,fm,fh)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    

    ! viscosity on w-points
!$omp simd
    do i=1,ip
        do k=1,kp
            vism(k,i)=lamsq(k)*sqrt(strain(k,i)) *fm(k,i)
        enddo
    enddo
!$omp end simd

!$omp simd
    do i=1,ip
        do k=1,kp
            ! viscosity for theta
            vist(k,i)=lamsq(k)*sqrt(strain(k,i)) *fh(k,i)*rhoa(k)
        enddo
    enddo
!$omp end simd
!     vist=vism
    vism(kp-1:kp,:)=0._sp
    vist(kp-1:kp,:)=0._sp
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! surface layer - monin-obukhov similarity theory                                    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call monin_obukhov(z0,z0th,thetan,zn,th,u,vism,vist,ip,kp,l_h,r_h)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! full exchange
    vism(0,:)=0._sp
    vism(kp+1,:)=0._sp
    vism(:,0)=vism(:,ip)
    vism(:,ip+1)=vism(:,1)

    vist(0,:)=0._sp
    vist(kp+1,:)=0._sp
    vist(:,0)=vist(:,ip)
    vist(:,ip+1)=vist(:,1)


    ! viscosity for q-fields
    do n=1,nq
        viss(:,:,n)=vist(:,:)
    enddo

    ! calculate elements of tau on p-points (in vertical)
    do i=1,ip
        do k=1,kp
            tau11(k,i)=rhoan(k)*1._sp*(vism(k,i)+vism(k-1,i))* &
                        (zu(k,i)-zu(k,i-1))/dx(i-1)
            tau33(k,i)=rhoan(k)*1._sp*(vism(k,i)+vism(k-1,i))* &
                        (zw(k,i)-zw(k-1,i))/dz(k-1)
                        
        enddo
    enddo
    

    ! calculate elements of tau on w-points (in vertical)
    do i=1,ip
        do k=1,kp
            ! on i+1/2
            tau13(k,i)=rhoa(k)*0.5_sp*(vism(k,i)+vism(k,i+1))* &
                ( (zu(k+1,i)-zu(k,i))/dzn(k) + (zw(k,i+1)-zw(k,i))/dxn(i))
                
        enddo
    enddo
    
    
    ! full exchange
    tau11(0,:)=0._sp
    tau11(kp+1,:)=0._sp
    tau11(:,0)=tau11(:,ip)
    tau11(:,ip+1)=tau11(:,1)

    tau33(0,:)=0._sp
    tau33(kp+1,:)=0._sp
    tau33(:,0)=tau33(:,ip)
    tau33(:,ip+1)=tau33(:,1)

    tau13(0,:)=0._sp
    tau13(kp+1,:)=0._sp
    tau13(:,0)=tau13(:,ip)
    tau13(:,ip+1)=tau13(:,1)

        

    ktmp=kp-1
    ktmp1=2
    ! calculate sources due to viscosity
!     tau23(1,:,:)=max(tau23(2,:,:),tau23(1,:,:))
    do i=1,ip
        do k=1,ktmp
            sw(k,i)=sw(k,i)+ &
               ((tau33(k+1,i)-tau33(k,i))/dzn(k) + &
                (tau13(k,i)-tau13(k,i-1))/dx(i-1) ) !/ rhoa(k)             
        enddo
        do k=ktmp1,ktmp
            ! note the sources are sources of u*rhoa
            su(k,i)=su(k,i)+ &
               ((tau11(k,i+1)-tau11(k,i))/dxn(i) + &
                (tau13(k,i)-tau13(k-1,i))/dz(k-1)) !/ rhoan(k)             
                
        enddo
    enddo
    
    ! full exchange    
    su(0,:)=0._sp
    su(kp+1,:)=0._sp
    su(:,0)=su(:,ip)
    su(:,ip+1)=su(:,1)
    
    sw(0,:)=0._sp
    sw(kp+1,:)=0._sp
    sw(:,0)=sw(:,ip)
    sw(:,ip+1)=sw(:,1)
    
    ! calculate th sources due to viscosity
    ! diffusion along i calculated via i+1/2 and i-1/2 terms
    sth=0._sp
    do i=1,ip
        do k=1,kp
            ! note, vist has already been multiplied by rhoa
            sth(k,i)=sth(k,i)+ &
             (0.25_sp*((vist(k,i)+vist(k,i-1)+vist(k-1,i)+vist(k-1,i-1))*&
                      (th(k,i-1)-th(k,i))/dxn(i-1) - &
                      (vist(k,i+1)+vist(k,i)+vist(k-1,i+1)+vist(k-1,i))*&
                      (th(k,i)-th(k,i+1))/dxn(i) ) / dx(i-1)  &
             + &
                     (vist(k-1,i)* (th(k-1,i)-th(k,i))/dz(k-1) - &
                     vist(k,i)* (th(k,i)-th(k+1,i))/dz(k) ) / dz(k-1)  ) &
                     / rhoan(k)                     
                    
                
        enddo
    enddo
    
    ! full exchange 
    sth(0,:)=0._sp
    sth(kp+1,:)=0._sp
    sth(:,0)=sth(:,ip)
    sth(:,ip+1)=sth(:,1)
    sq=0._sp
    ! q-fields
    do n=1,nq
        ! calculate q sources due to viscosity
        ! diffusion along i calculated via i+1/2 and i-1/2 terms
        do i=1,ip
            do k=1,kp
                ! note, vist has already been multiplied by rhoa
                sq(k,i,n)=sq(k,i,n)+ &
                 (0.25_sp*((vist(k,i)+vist(k,i-1)+vist(k-1,i)+vist(k-1,i-1))*&
                          (q(k,i-1,n)-q(k,i,n))/dxn(i-1) - &
                          (vist(k,i+1)+vist(k,i)+vist(k-1,i+1)+vist(k-1,i))*&
                          (q(k,i,n)-q(k,i+1,n))/dxn(i) ) / dx(i-1)  &
                 + &
                         (vist(k-1,i)* (q(k-1,i,n)-q(k,i,n))/dz(k-1) - &
                         vist(k,i)* (q(k,i,n)-q(k+1,i,n))/dz(k) ) / dz(k-1)  ) &
                         / rhoan(k)                     
                    
                
            enddo
        enddo
    
        ! full exchange  
        sq(0,:,n)=0._sp
        sq(kp+1,:,n)=0._sp
        sq(:,0,n)=sq(:,ip,n)
        sq(:,ip+1,n)=sq(:,1,n)
    enddo


    
	end subroutine calculate_subgrid_2d
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
    subroutine monin_obukhov(z0,z0th,thetan,zn,th,u,vism,vist,ip,kp,l_h,r_h)
	use nrtype
    
    integer(i4b), intent(in) :: ip, kp, l_h, r_h
    implicit none
    real(sp), intent(in), dimension(-r_h+1:kp+r_h,-l_h+1:ip+r_h) :: &
            th, u
    real(sp), intent(inout), dimension(-r_h+1:kp+r_h,-l_h+1:ip+r_h) :: &
            vism,vist
    real(sp), intent(in), dimension(-r_h+1:kp+r_h) :: thetan, zn
    real(sp), intent(in) :: z0, z0th
            
    ! locals
    integer(i4b) :: i,j
    real(sp), dimension(-l_h+1:ip+r_h) :: vmag,rib,gm,gh, &
                                                            ustar,thstar
    
    
    
    do i=1,ip
        ! wind speed in first layer above surface
        vmag(i)=0.5_sp*sqrt((u(2,i)+u(2,i-1))**2  )+small
        
        ! bulk richardson number
        ! equation 8.39 (Jacobson pp 243)
        rib(i)=grav*(th(2,i)+thetan(2)-th(1,i)-thetan(1))*(zn(2)-z0)**2 / &
            (thetan(1)*vmag(i)**2*(zn(2)-z0th)) 
        
        ! equation 8.41 (Jacobson pp 243)
        if(rib(i) .le. 0._sp) then
            gm(i)=1._sp-9.4_sp*rib(i)/ &
                (1._sp+70._sp*kvon**2*sqrt(abs(rib(i)*zn(2)/z0))/(log(zn(2)/z0))**2)
            gh(i)=1._sp-9.4_sp*rib(i)/ &
                (1._sp+50._sp*kvon**2*sqrt(abs(rib(i)*zn(2)/z0))/(log(zn(2)/z0))**2)
        else
            gm(i)=1._sp/(1._sp+4.7_sp*rib(i))**2
            gh(i)=1._sp/(1._sp+4.7_sp*rib(i))**2
        endif
        
        ! equation 8.40 (Jacobson pp 243)
        ustar(i)=kvon*vmag(i)/log(zn(2)/z0)*sqrt(gm(i))
        thstar(i)=kvon**2*vmag(i)*(th(2,i)+thetan(2)-th(1,i)-thetan(1)) / &
            (ustar(i)*prn*log(zn(2)/z0)**2)*gh(i)
            
            
        ! see equation 8.46, page 245 (Jacobson)
        vism(1,i)=ustar(i)**2*(zn(2)-z0)/vmag(i)
        
        vist(1,i)=ustar(i)*thstar(i)/ &
            ((th(2,i)+thetan(2)-th(1,i)-thetan(1))/(zn(2)-z0th) + small )
    enddo
    
    
    

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
    subroutine richardson(strain,thetan,theta,dzn,th,ip,kp,l_h,r_h,rip,fm,fh)
	use nrtype
    
    implicit none
    integer(i4b), intent(in) :: ip, kp, l_h, r_h
    real(sp), intent(in), dimension(-r_h+1:kp+r_h,-l_h+1:ip+r_h) :: &
            strain, th
    real(sp), intent(in), dimension(-r_h+1:kp+r_h) :: thetan, theta, dzn
    real(sp), intent(inout), dimension(-r_h+1:kp+r_h,-l_h+1:ip+r_h) :: &
            rip,fm, fh
            
    ! locals
    integer(i4b) :: i,j,k
    
    
    ! calculate bulk richardson number
    do i=1,ip
        do k=1,kp
            rip(k,i)=grav*((th(k+1,i))/thetan(k+1) - &
                        (th(k,i))/thetan(k)) / dzn(k) / (strain(k,i)+1.e-15)
        enddo
    enddo
    ! richardson number functions
    do i=1,ip
        do k=1,kp
            if(rip(k,i) .lt. 0._sp) then
                fm(k,i)=(1._sp-subc*rip(k,i))**0.5_sp
                fh(k,i)=suba*(1._sp-subb*rip(k,i))**0.5_sp
            elseif((rip(k,i).ge.0._sp).and.(rip(k,i).lt.ric)) then
                fm(k,i)=(1._sp-rip(k,i)/ric)**subr*(1._sp-subh*rip(k,i))
                fh(k,i)=subf*(1._sp-rip(k,i)/ric)**subr*(1._sp-subg*rip(k,i))
            else
                fm(k,i)=0._sp
                fh(k,i)=0._sp
            endif

        enddo
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
    subroutine advance_fields_2d(dt,tu,tw, &
        zu,zw,su,sw,q,sq,th,sth,rhoa,rhoan,ip,kp,nq,l_h,r_h)
	use nrtype
    
    implicit none
    integer(i4b), intent(in) :: ip, kp, nq, l_h, r_h
    real(sp), intent(in), dimension(-r_h+1:kp+r_h,-l_h+1:ip+r_h) :: &
        su,sw,sth
    real(sp), intent(in), dimension(-r_h+1:kp+r_h) :: &
        rhoa,rhoan
    real(sp), intent(inout), dimension(-r_h+1:kp+r_h,-l_h+1:ip+r_h) :: &
        tu,tw,zu,zw,th
    real(sp), intent(in), dimension(-r_h+1:kp+r_h,-l_h+1:ip+r_h,nq) :: sq
    real(sp), intent(inout), dimension(-r_h+1:kp+r_h,-l_h+1:ip+r_h,nq) :: q
    real(sp), intent(in) :: dt
    ! locals
    integer(i4b) :: i,j,k,n
    
    do i=1-l_h,ip+r_h
        do k=1-l_h,kp+r_h
            tu(k,i)=zu(k,i)+su(k,i)*dt*2._sp/rhoan(k)
        enddo
    enddo
    
    
    do i=1-l_h,ip+r_h
        do k=1-l_h,kp+r_h
            tw(k,i)=zw(k,i)+sw(k,i)*dt*2._sp/rhoa(k)
        enddo
    enddo
    
    do i=1-l_h,ip+r_h
        do k=1-l_h,kp+r_h
            th(k,i)=th(k,i)+sth(k,i)*dt
        enddo
    enddo
    
    do n=1,nq
        do i=1-l_h,ip+r_h
            do k=1-l_h,kp+r_h
                q(k,i,n)=q(k,i,n)+sq(k,i,n)*dt
            enddo
        enddo
    enddo    
    
    end subroutine advance_fields_2d
    
    	
	

    end module subgrid_2d
