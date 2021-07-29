	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>subgrid code for the dynamical cloud model
    module subgrid_3d
    use nrtype
    
    private
    public :: calculate_subgrid_3d,advance_fields_3d,advance_scalar_fields_3d
    
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
	!>@param[in] z,zn,dx,dy,dz, rhoa
	!>@param[in] dxn,dyn,dzn, rhoan
	!>@param[in] ip,jp,kp,l_h,r_h, nq
	!>@param[in] u, zu
	!>@param[in] v, zv
	!>@param[in] w, zw
	!>@param[in] th,q, lamsq, z0,z0th
	!>@param[inout] vism, vist,viss,strain, tau11, tau22, tau33, tau12, tau13, tau23, 
	!>                su, sv, sw,sth,sq
	!>@param[in] moisture - flag for moisture
    !>@param[in] comm3d, id, dims, coords: mpi variables
	subroutine calculate_subgrid_3d(dt,z,zn,dx,dy,dz,rhoa, theta,&
	                                dxn,dyn,dzn,rhoan, thetan,&
	                                ip,jp,kp,nq,l_h,r_h,u,v,w,zu,zv, zw, &
	                                th,q, &
	                                lamsq, z0,z0th, &
	                                vism,vist, viss,strain, &
	                                tau11,tau22,tau33,tau12,&
	                                tau13,tau23, &
	                                su,sv,sw,sth,sq, moisture, &
	                                comm3d,id,dims,coords)
	use nrtype
	use mpi_module
	use mpi
	
	implicit none

    integer(i4b), intent(in) :: id, comm3d
    integer(i4b), dimension(3), intent(in) :: dims, coords
	real(sp), intent(in) :: dt
	integer(i4b), intent(in) :: ip, jp, kp, l_h, r_h, nq
	real(sp), dimension(-r_h+1:kp+r_h,-r_h+1:jp+r_h,-l_h+1:ip+r_h), &
		intent(in) :: u, zu
	real(sp), dimension(-r_h+1:kp+r_h,-l_h+1:jp+r_h,-r_h+1:ip+r_h), &
		intent(in) :: v, zv
	real(sp), dimension(-l_h+1:kp+r_h,-r_h+1:jp+r_h,-r_h+1:ip+r_h), &
		intent(in) :: w, zw, th
	real(sp), dimension(-r_h+1:kp+r_h,-r_h+1:jp+r_h,-r_h+1:ip+r_h,nq), &
		intent(in) :: q
	real(sp), dimension(-r_h+1:kp+r_h,-r_h+1:jp+r_h,-r_h+1:ip+r_h,nq), &
		intent(inout) :: viss,sq
	real(sp), dimension(-r_h+1:kp+r_h,-r_h+1:jp+r_h,-r_h+1:ip+r_h), &
		intent(inout) :: vism, vist, &
		                strain, tau11, tau22, tau33, tau12, tau13, tau23, &
		                su, sv, sw, sth
	real(sp), dimension(-l_h+1:ip+r_h), intent(in) :: dx, dxn
	real(sp), dimension(-l_h+1:jp+r_h), intent(in) :: dy, dyn
	real(sp), dimension(-l_h+1:kp+r_h), intent(in) :: z, zn,&
	        dz, dzn, rhoa, rhoan, theta, thetan
    real(sp), dimension(-l_h+1:kp+r_h), intent(in) :: lamsq
    real(sp), intent(in) :: z0,z0th
    logical, intent(in) :: moisture
	
	! locals
	integer(i4b) :: i,j,k,n, ktmp,ktmp1
	real(sp), dimension(-r_h+1:kp+r_h,-r_h+1:jp+r_h,-r_h+1:ip+r_h) :: rip, fm, fh

    ! u is on x levels, v on y levels, w on z levels. 
    ! so distance between u(i) and u(i+1) is dx(i)
    !    w is on xn levels
    !    distance between w(i) and w(i+1) is dxn(i)
    !    u is on zn levels
    !    distance between u(k) and u(k+1) is dzn(k)
    ! strain-rate on w-points
!$omp simd
    do i=1,ip
        do j=1,jp
            do k=1,kp
                ! 2*s11*s11
                strain(k,j,i)=((u(k+1,j,i)-u(k+1,j,i-1))/dx(i-1))**2 + &
                              ((u(k,j,i)-u(k,j,i-1))/dx(i-1))**2
                ! 2*s22*s22
                strain(k,j,i)=strain(k,j,i)+ &
                              ((v(k+1,j,i)-v(k+1,j-1,i))/dy(j-1))**2 + &
                              ((v(k,j,i)-v(k,j-1,i))/dy(j-1))**2 
                ! 2*s33*s33
                strain(k,j,i)=strain(k,j,i)+ &
                              ((w(k,j,i)-w(k-1,j,i))/dz(k-1))**2 + &
                              ((w(k+1,j,i)-w(k,j,i))/dz(k))**2 
                              
                ! 2*s13*s13 - du/dz+dw/dx - averaging over 2 points
                strain(k,j,i)=strain(k,j,i)+ 0.5_sp * &
                   (((u(k+1,j,i)-u(k,j,i))/dzn(k)+(w(k,j,i+1)-w(k,j,i))/dxn(i))**2 + &
                   ((u(k+1,j,i-1)-u(k,j,i-1))/dz(k)+(w(k,j,i)-w(k,j,i-1))/dxn(i-1))**2)
                ! 2*s23*s23 - dw/dy+dv/dz - averaging over 2 points
                strain(k,j,i)=strain(k,j,i)+ 0.5_sp * &
                   (((w(k,j,i)-w(k,j-1,i))/dyn(j-1)+(v(k+1,j-1,i)-v(k,j-1,i))/dzn(k))**2 + &
                   ((w(k,j+1,i)-w(k,j,i))/dyn(j)+(v(k+1,j,i)-v(k,j,i))/dzn(k))**2)
                   
                ! 2*s12*s12 - du/dy+dv/dx - averaging over 8 points
                strain(k,j,i)=strain(k,j,i)+ 0.125_sp * &
                   (((u(k,j,i-1)-u(k,j-1,i-1))/dyn(j-1)+(v(k,j-1,i)-v(k,j-1,i-1))/dxn(i-1))**2 + &
                   ((u(k,j+1,i-1)-u(k,j,i-1))/dyn(j)+(v(k,j,i)-v(k,j,i-1))/dxn(i-1))**2  &
                   + &
                   ((u(k,j,i)-u(k,j-1,i))/dyn(j-1)+(v(k,j-1,i+1)-v(k,j-1,i))/dxn(i))**2 + &
                   ((u(k,j+1,i)-u(k,j,i))/dyn(j)+(v(k,j,i+1)-v(k,j,i))/dxn(i))**2  &
                   + &
                   ((u(k+1,j,i-1)-u(k+1,j-1,i-1))/dyn(j-1)+(v(k+1,j-1,i)-v(k+1,j-1,i-1))/dxn(i-1))**2 + &
                   ((u(k,j+1,i-1)-u(k,j,i-1))/dyn(j)+(v(k,j,i)-v(k,j,i-1))/dxn(i-1))**2  &
                   + &
                   ((u(k+1,j,i)-u(k+1,j-1,i))/dyn(j-1)+(v(k+1,j-1,i+1)-v(k+1,j-1,i))/dxn(i))**2 + &
                   ((u(k+1,j+1,i)-u(k+1,j,i))/dyn(j)+(v(k+1,j,i+1)-v(k+1,j,i))/dxn(i))**2 )
            enddo
        enddo
    enddo
!$omp end simd

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Richardson number and functions                                                    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call richardson(strain,thetan,theta,dzn,th,ip,jp,kp,l_h,r_h,rip,fm,fh)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    

    ! viscosity on w-points
!$omp simd
    do i=1,ip
        do j=1,jp
            do k=1,kp
                vism(k,j,i)=lamsq(k)*sqrt(strain(k,j,i)) *fm(k,j,i)
            enddo
        enddo
    enddo
!$omp end simd

!$omp simd
    do i=1,ip
        do j=1,jp
            do k=1,kp
                ! viscosity for theta
                vist(k,j,i)=lamsq(k)*sqrt(strain(k,j,i)) *fh(k,j,i)*rhoa(k)
            enddo
        enddo
    enddo
!$omp end simd
!     vist=vism
    if(coords(3)==(dims(3)-1)) then
        vism(kp-1:kp,:,:)=0._sp
        vist(kp-1:kp,:,:)=0._sp
    endif
    if(coords(3)==0) then
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! surface layer - monin-obukhov similarity theory                                !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call monin_obukhov(z0,z0th,theta,thetan,z,zn,th,u,v,vism,vist,ip,jp,kp,l_h,r_h)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    endif

    ! full exchange
    call exchange_full_wo(comm3d, id, kp, jp, ip, r_h,r_h,r_h,r_h,l_h,r_h,vism,&
        dims,coords)
    call exchange_full_wo(comm3d, id, kp, jp, ip, r_h,r_h,r_h,r_h,l_h,r_h,vist,&
        dims,coords)


    if(moisture) then
        ! viscosity for q-fields
        do n=1,nq
            viss(:,:,:,n)=vist(:,:,:)
        enddo
    endif
    
    
    ! calculate elements of tau on p-points (in vertical)
    do i=1,ip
        do j=1,jp
            do k=1,kp
                tau11(k,j,i)=rhoan(k)*1._sp*(vism(k,j,i)+vism(k-1,j,i))* &
                            (zu(k,j,i)-zu(k,j,i-1))/dx(i-1)
                tau22(k,j,i)=rhoan(k)*1._sp*(vism(k,j,i)+vism(k-1,j,i))* &
                            (zv(k,j,i)-zv(k,j-1,i))/dy(j-1)
                tau33(k,j,i)=rhoan(k)*1._sp*(vism(k,j,i)+vism(k-1,j,i))* &
                            (zw(k,j,i)-zw(k-1,j,i))/dz(k-1)
                            
                !i+1/2, j+1/2
                tau12(k,j,i)=rhoan(k)* &
                    0.125_sp* (vism(k,j,i)+vism(k,j+1,i)+vism(k-1,j,i)+vism(k-1,j+1,i) + &
                    vism(k,j,i+1)+vism(k,j+1,i+1)+vism(k-1,j,i+1)+vism(k-1,j+1,i+1) ) * &
                    ((zu(k,j+1,i)-zu(k,j,i))/dyn(j)+(zv(k,j,i+1)-zv(k,j,i))/dxn(i))
            enddo
        enddo
    enddo
    

    ! calculate elements of tau on w-points (in vertical)
    do i=1,ip
        do j=1,jp
            do k=1,kp
                ! on i+1/2
                tau13(k,j,i)=rhoa(k)*0.5_sp*(vism(k,j,i)+vism(k,j,i+1))* &
                    ( (zu(k+1,j,i)-zu(k,j,i))/dzn(k) + (zw(k,j,i+1)-zw(k,j,i))/dxn(i))
                    
                ! on j+1/2
                tau23(k,j,i)=rhoa(k)*0.5_sp*(vism(k,j,i)+vism(k,j+1,i))* &
                    ( (zw(k,j+1,i)-zw(k,j,i))/dyn(j) + (zv(k+1,j,i)-zv(k,j,i))/dzn(k))
            enddo
        enddo
    enddo
    
    
    ! full exchange
    call exchange_along_dim(comm3d, id, kp, jp, ip, r_h,r_h,r_h,r_h,l_h,r_h,tau11,&
        0._sp,0._sp,dims,coords)
    call exchange_along_dim(comm3d, id, kp, jp, ip, r_h,r_h,r_h,r_h,l_h,r_h,tau22,&
        0._sp,0._sp,dims,coords)
    call exchange_along_dim(comm3d, id, kp, jp, ip, r_h,r_h,r_h,r_h,l_h,r_h,tau33,&
        0._sp,0._sp,dims,coords)
    call exchange_along_dim(comm3d, id, kp, jp, ip, r_h,r_h,r_h,r_h,l_h,r_h,tau12,&
        0._sp,0._sp,dims,coords)
    call exchange_along_dim(comm3d, id, kp, jp, ip, r_h,r_h,r_h,r_h,l_h,r_h,tau13,&
        0._sp,0._sp,dims,coords)
    call exchange_along_dim(comm3d, id, kp, jp, ip, r_h,r_h,r_h,r_h,l_h,r_h,tau23,&
        0._sp,0._sp,dims,coords)
        

    if(coords(3) == (dims(3)-1)) then
        ktmp=kp-1
    else 
        ktmp=kp
    endif
    if(coords(3)==0) then 
        ktmp1=2
    else
        ktmp1=1
    endif
    ! calculate sources due to viscosity
!     tau23(1,:,:)=max(tau23(2,:,:),tau23(1,:,:))
    do i=1,ip
        do j=1,jp
            do k=1,ktmp
                sw(k,j,i)=sw(k,j,i)+ &
                   ((tau33(k+1,j,i)-tau33(k,j,i))/dzn(k) + &
                    (tau13(k,j,i)-tau13(k,j,i-1))/dx(i-1) + &
                    (tau23(k,j,i)-tau23(k,j-1,i))/dy(j-1)) !/ rhoa(k)             
            enddo
            do k=ktmp1,ktmp
                ! note the sources are sources of u*rhoa
                su(k,j,i)=su(k,j,i)+ &
                   ((tau11(k,j,i+1)-tau11(k,j,i))/dxn(i) + &
                    (tau12(k,j,i)-tau12(k,j-1,i))/dy(j-1) + &
                    (tau13(k,j,i)-tau13(k-1,j,i))/dz(k-1)) !/ rhoan(k)             
                    
                sv(k,j,i)=sv(k,j,i)+ &
                   ((tau22(k,j+1,i)-tau22(k,j,i))/dyn(j) + &
                    (tau12(k,j,i)-tau12(k,j,i-1))/dx(i-1) + &
                    (tau23(k,j,i)-tau23(k-1,j,i))/dz(k-1)) !/ rhoan(k)             

            enddo
        enddo
    enddo
    
    ! full exchange    
    call exchange_full(comm3d, id, kp, jp, ip, r_h,r_h,r_h,r_h,l_h,r_h,su,&
        0._sp,0._sp,dims,coords)
    call exchange_full(comm3d, id, kp, jp, ip, r_h,r_h,r_h,r_h,l_h,r_h,sv,&
        0._sp,0._sp,dims,coords)
    call exchange_full(comm3d, id, kp, jp, ip, r_h,r_h,r_h,r_h,l_h,r_h,sw,&
        0._sp,0._sp,dims,coords)
    
    ! calculate th sources due to viscosity
    ! diffusion along i calculated via i+1/2 and i-1/2 terms
!     sth=0._sp
    do i=1,ip
        do j=1,jp
            do k=1,kp
                ! note, vist has already been multiplied by rhoa
                sth(k,j,i)=sth(k,j,i)+ &
                 (0.25_sp*((vist(k,j,i)+vist(k,j,i-1)+vist(k-1,j,i)+vist(k-1,j,i-1))*&
                          (th(k,j,i-1)-th(k,j,i))/dxn(i-1) - &
                          (vist(k,j,i+1)+vist(k,j,i)+vist(k-1,j,i+1)+vist(k-1,j,i))*&
                          (th(k,j,i)-th(k,j,i+1))/dxn(i) ) / dx(i-1)  &
                 + &
                 0.25_sp*((vist(k,j,i)+vist(k,j-1,i)+vist(k-1,j,i)+vist(k-1,j-1,i))*&
                          (th(k,j-1,i)-th(k,j,i))/dyn(j-1) - &
                          (vist(k,j+1,i)+vist(k,j,i)+vist(k-1,j+1,i)+vist(k-1,j,i))*&
                          (th(k,j,i)-th(k,j+1,i))/dyn(j) ) / dy(j-1)  &
                 + &
                         (vist(k-1,j,i)* (th(k-1,j,i)-th(k,j,i))/dz(k-1) - &
                         vist(k,j,i)* (th(k,j,i)-th(k+1,j,i))/dz(k) ) / dz(k-1)  ) &
                         / rhoan(k)                     
                        
                    
            enddo
        enddo
    enddo
    
    ! full exchange    
    call exchange_full(comm3d, id, kp, jp, ip, r_h,r_h,r_h,r_h,l_h,r_h,sth,&
        0._sp,0._sp,dims,coords)
        
        
    if(moisture) then
        sq=0._sp
        ! q-fields
        do n=1,nq
            ! calculate q sources due to viscosity
            ! diffusion along i calculated via i+1/2 and i-1/2 terms
            do i=1,ip
                do j=1,jp
                    do k=1,kp
                        ! note, vist has already been multiplied by rhoa
                        sq(k,j,i,n)=sq(k,j,i,n)+ &
                         (0.25_sp*&
                         ((vist(k,j,i)+vist(k,j,i-1)+vist(k-1,j,i)+vist(k-1,j,i-1))*&
                                  (q(k,j,i-1,n)-q(k,j,i,n))/dxn(i-1) - &
                        (vist(k,j,i+1)+vist(k,j,i)+vist(k-1,j,i+1)+vist(k-1,j,i))*&
                                  (q(k,j,i,n)-q(k,j,i+1,n))/dxn(i) ) / dx(i-1)  &
                         + &
                         0.25_sp*&
                         ((vist(k,j,i)+vist(k,j-1,i)+vist(k-1,j,i)+vist(k-1,j-1,i))*&
                                  (q(k,j-1,i,n)-q(k,j,i,n))/dyn(j-1) - &
                        (vist(k,j+1,i)+vist(k,j,i)+vist(k-1,j+1,i)+vist(k-1,j,i))*&
                                  (q(k,j,i,n)-q(k,j+1,i,n))/dyn(j) ) / dy(j-1)  &
                         + &
                                 (vist(k-1,j,i)* (q(k-1,j,i,n)-q(k,j,i,n))/dz(k-1) - &
                            vist(k,j,i)* (q(k,j,i,n)-q(k+1,j,i,n))/dz(k) ) / dz(k-1)  ) &
                                 / rhoan(k)                     
                        
                    
                    enddo
                enddo
            enddo
    
            ! full exchange    
            call exchange_full(comm3d, id, kp, jp, ip, r_h,r_h,r_h,r_h,l_h,r_h,&
                sq(:,:,:,n),0._sp,0._sp,dims,coords)
        enddo
    endif

    
	end subroutine calculate_subgrid_3d
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
	!>@param[in] z0,z0th,theta,thetan,z,zn,th,u,v,ip,jp,kp,l_h,r_h
	!>@param[inout] vism,vist
    subroutine monin_obukhov(z0,z0th,theta,thetan,z,zn,th,u,v,vism,vist,ip,jp,kp,l_h,r_h)
	use nrtype
	use mpi_module
	use mpi
    
    implicit none
    integer(i4b), intent(in) :: ip, jp, kp, l_h, r_h
    real(sp), intent(in), dimension(-r_h+1:kp+r_h,-r_h+1:jp+r_h,-l_h+1:ip+r_h) :: &
            th, u, v
    real(sp), intent(inout), dimension(-r_h+1:kp+r_h,-r_h+1:jp+r_h,-l_h+1:ip+r_h) :: &
            vism,vist
    real(sp), intent(in), dimension(-r_h+1:kp+r_h) :: thetan, zn,z,theta
    real(sp), intent(in) :: z0, z0th
            
    ! locals
    integer(i4b) :: i,j
    real(sp), dimension(-r_h+1:jp+r_h,-l_h+1:ip+r_h) :: vmag,rib,gm,gh, &
                                                            ustar,thstar
    
    
    
    do i=1,ip
        do j=1,jp
            ! wind speed in first layer above surface 
            vmag(j,i)=0.5_sp*sqrt((u(1,j,i)+u(1,j,i-1))**2 + &
                                  (v(1,j,i)+v(1,j-1,i))**2  )+small
            
            ! bulk richardson number
            ! equation 8.39 (Jacobson pp 243) - theta(0) is the surface there is no pertubation here
            rib(j,i)=grav*(th(1,j,i)+thetan(1)-theta(0))*(z(1)-z0)**2 / &
                (theta(0)*vmag(j,i)**2*(z(1)-z0th)) 
            
            ! equation 8.41 (Jacobson pp 243) - need to use z(1)
            if(rib(j,i) .le. 0._sp) then
                gm(j,i)=1._sp-9.4_sp*rib(j,i)/ &
                    (1._sp+70._sp*kvon**2*sqrt(abs(rib(j,i)*z(1)/z0))/(log(z(1)/z0))**2)
                gh(j,i)=1._sp-9.4_sp*rib(j,i)/ &
                    (1._sp+50._sp*kvon**2*sqrt(abs(rib(j,i)*z(1)/z0))/(log(z(1)/z0))**2)
            else
                gm(j,i)=1._sp/(1._sp+4.7_sp*rib(j,i))**2
                gh(j,i)=1._sp/(1._sp+4.7_sp*rib(j,i))**2
            endif
            
            ! equation 8.40 (Jacobson pp 243) - need to use z(1)
            ustar(j,i)=kvon*vmag(j,i)/log(z(1)/z0)*sqrt(gm(j,i))
            thstar(j,i)=kvon**2*vmag(j,i)*(th(1,j,i)+thetan(1)-theta(0)) / &
                (ustar(j,i)*prn*log(z(1)/z0)**2)*gh(j,i)
                
                
            ! see equation 8.46, page 245 (Jacobson) - need to use z(1) / look
            vism(0,j,i)=ustar(j,i)**2*(z(1)-z0)/vmag(j,i)
            
            vist(0,j,i)=ustar(j,i)*thstar(j,i)/ &
                ((th(1,j,i)+thetan(1)-theta(0))/(z(1)-z0th) + small )
        enddo
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
	!>@param[in] strain, thetan, theta, dzn, th, ip, jp, kp, l_h, r_h
	!>@param[inout] rip, fm, fh
    subroutine richardson(strain,thetan,theta,dzn,th,ip,jp,kp,l_h,r_h,rip,fm,fh)
	use nrtype
	use mpi_module
	use mpi
    
    implicit none
    integer(i4b), intent(in) :: ip, jp, kp, l_h, r_h
    real(sp), intent(in), dimension(-r_h+1:kp+r_h,-r_h+1:jp+r_h,-l_h+1:ip+r_h) :: &
            strain, th
    real(sp), intent(in), dimension(-r_h+1:kp+r_h) :: thetan, theta, dzn
    real(sp), intent(inout), dimension(-r_h+1:kp+r_h,-r_h+1:jp+r_h,-l_h+1:ip+r_h) :: &
            rip,fm, fh
            
    ! locals
    integer(i4b) :: i,j,k
    
    
    ! calculate bulk richardson number
    do i=1,ip
        do j=1,jp
            do k=1,kp
                rip(k,j,i)=grav*((th(k+1,j,i))/thetan(k+1) - &
                            (th(k,j,i))/thetan(k)) / dzn(k) / (strain(k,j,i)+1.e-15)
            enddo
        enddo
    enddo
    ! richardson number functions
    do i=1,ip
        do j=1,jp
            do k=1,kp
                if(rip(k,j,i) .lt. 0._sp) then
                    fm(k,j,i)=(1._sp-subc*rip(k,j,i))**0.5_sp
                    fh(k,j,i)=suba*(1._sp-subb*rip(k,j,i))**0.5_sp
                elseif((rip(k,j,i).ge.0._sp).and.(rip(k,j,i).lt.ric)) then
                    fm(k,j,i)=(1._sp-rip(k,j,i)/ric)**subr*(1._sp-subh*rip(k,j,i))
                    fh(k,j,i)=subf*(1._sp-rip(k,j,i)/ric)**subr*(1._sp-subg*rip(k,j,i))
                else
                    fm(k,j,i)=0._sp
                    fh(k,j,i)=0._sp
                endif

            enddo
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
	!>@param[in] su,sv,sw,sq,sth,ipp, jpp, kpp, nq, l_h,r_h,dt, rhoa, rhoan
	!>@param[inout] q,th,u,v,w
    subroutine advance_fields_3d(dt,tu,tv,tw, &
        zu,zv,zw,su,sv,sw,q,sq,th,sth,rhoa,rhoan,ip,jp,kp,nq,l_h,r_h)
	use nrtype
	use mpi_module
	use mpi
    
    implicit none
    integer(i4b), intent(in) :: ip, jp, kp, nq, l_h, r_h
    real(sp), intent(in), dimension(-r_h+1:kp+r_h,-r_h+1:jp+r_h,-l_h+1:ip+r_h) :: &
        su,sv,sw,sth
    real(sp), intent(in), dimension(-r_h+1:kp+r_h) :: &
        rhoa,rhoan
    real(sp), intent(inout), dimension(-r_h+1:kp+r_h,-r_h+1:jp+r_h,-l_h+1:ip+r_h) :: &
        tu,tv,tw,zu,zv,zw,th
    real(sp), intent(in), dimension(-r_h+1:kp+r_h,-r_h+1:jp+r_h,-l_h+1:ip+r_h,nq) :: sq
    real(sp), intent(inout), dimension(-r_h+1:kp+r_h,-r_h+1:jp+r_h,-l_h+1:ip+r_h,nq) :: q
    real(sp), intent(in) :: dt
    ! locals
    integer(i4b) :: i,j,k,n
    
    do i=1-l_h,ip+r_h
        do j=1-l_h,jp+r_h
            do k=1-l_h,kp+r_h
                tu(k,j,i)=zu(k,j,i)+su(k,j,i)*dt*2._sp/rhoan(k)
            enddo
        enddo
    enddo
    
    do i=1-l_h,ip+r_h
        do j=1-l_h,jp+r_h
            do k=1-l_h,kp+r_h
                tv(k,j,i)=zv(k,j,i)+sv(k,j,i)*dt*2._sp/rhoan(k)
            enddo
        enddo
    enddo
    
    do i=1-l_h,ip+r_h
        do j=1-l_h,jp+r_h
            do k=1-l_h,kp+r_h
                tw(k,j,i)=zw(k,j,i)+sw(k,j,i)*dt*2._sp/rhoa(k)
            enddo
        enddo
    enddo
    
    do i=1-l_h,ip+r_h
        do j=1-l_h,jp+r_h
            do k=1-l_h,kp+r_h
                th(k,j,i)=th(k,j,i)+sth(k,j,i)*dt
            enddo
        enddo
    enddo
    
    do n=1,nq
        do i=1-l_h,ip+r_h
            do j=1-l_h,jp+r_h
                do k=1-l_h,kp+r_h
                    q(k,j,i,n)=q(k,j,i,n)+sq(k,j,i,n)*dt
                enddo
            enddo
        enddo
    enddo    
    
    end subroutine advance_fields_3d
    
    	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Advance fields due to diffusion term                                               !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>just add the source term to the current value of the field
	!>@param[in] sq,sth,ipp, jpp, kpp, nq, l_h,r_h,dt, rhoa, rhoan
	!>@param[inout] q,th
	!>@param[in] moisture
    subroutine advance_scalar_fields_3d(dt, &
        q,sq,th,sth,rhoa,rhoan,ip,jp,kp,nq,l_h,r_h,moisture)
	use nrtype
	use mpi_module
	use mpi
    
    implicit none
    integer(i4b), intent(in) :: ip, jp, kp, nq, l_h, r_h
    real(sp), intent(in), dimension(-r_h+1:kp+r_h,-r_h+1:jp+r_h,-l_h+1:ip+r_h) :: &
        sth
    real(sp), intent(in), dimension(-r_h+1:kp+r_h) :: &
        rhoa,rhoan
    real(sp), intent(inout), dimension(-r_h+1:kp+r_h,-r_h+1:jp+r_h,-l_h+1:ip+r_h) :: &
        th
    real(sp), intent(in), dimension(-r_h+1:kp+r_h,-r_h+1:jp+r_h,-l_h+1:ip+r_h,nq) :: sq
    real(sp), intent(inout), dimension(-r_h+1:kp+r_h,-r_h+1:jp+r_h,-l_h+1:ip+r_h,nq) :: q
    real(sp), intent(in) :: dt
    logical, intent(in) :: moisture
    ! locals
    integer(i4b) :: i,j,k,n
    
    do i=1-l_h,ip+r_h
        do j=1-l_h,jp+r_h
            do k=1-l_h,kp+r_h
                th(k,j,i)=th(k,j,i)+sth(k,j,i)*dt
            enddo
        enddo
    enddo
    if(moisture) then
        do n=1,nq
!             if((n<2).or.(n>14)) then
                do i=1-l_h,ip+r_h
                    do j=1-l_h,jp+r_h
                        do k=1-l_h,kp+r_h
                            q(k,j,i,n)=q(k,j,i,n)+sq(k,j,i,n)*dt
                        enddo
                    enddo
                enddo
!             endif
        enddo    
    endif    
    end subroutine advance_scalar_fields_3d
    
	

    end module subgrid_3d
