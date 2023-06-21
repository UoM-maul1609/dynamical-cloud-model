!>@author
!>Paul Connolly, The University of Manchester
!>@brief
!>advection routines
    module advection_1d
    use numerics_type
    private
    public :: bott_scheme_1d
    contains
!>@author
!>Paul J. Connolly, The University of Manchester
!>@brief
!>advects a scalar field in 1-d (see A. Bott, MWR, 1992)
!>@param[in] kp:  number of grid points
!>@param[in] ord: order of interpolation scheme
!>@param[in] o_halo: halos required for advection scheme
!>@param[in] dt:  timestep
!>@param[in] dx:  grid spacing
!>@param[in] x:   grid
!>@param[in] u:   velocity
!>@param[inout] psi: field being advected
!>@param[in] monotone: flag for monotonic advection
!>solves the 1-d advection equation:
!>\f$ \frac{\partial \psi}{\partial t} + \frac{\partial u \psi}{\partial x} = 0 \f$
    subroutine bott_scheme_1d(kp,ord,o_halo,dt,dx,x,u,psi,monotone)

    use numerics_type
    implicit none
    ! arguments:
    integer(i4b), intent(in) :: kp, ord, o_halo
    real(wp), intent(in) :: dt,dx
    real(wp), dimension(-o_halo+1:kp+o_halo), intent(in) :: x,u
    logical, intent(in) :: monotone
    real(wp), dimension(-o_halo+1:kp+o_halo), intent(inout) :: psi

    ! locals:
	real(wp), dimension(ord+1) :: a_coeff01, a_coeff02
	real(wp), dimension(-o_halo+1:kp+o_halo) :: f_plus, f_plus_mon,&
									f_minus,f_minus_mon, fj, psi_old, u1
	real(wp) :: cp_j, cp_jp1, cp_jm1, cp_jm2, cm_j, cm_jm1, cm_jp1, &
	dp_jm05=0._wp, dp_jp05=0._wp, dp_jp15=0._wp, dm_jp05=0._wp, dm_jm05=0._wp, &
	i_j, i_jp1, f_pm1, f_m
	real(wp) :: dummy1, dummy2,small=1e-60_wp
	integer(i4b) :: j, k

	if(sum(psi).lt.small) then
		return
	endif
	! zero arrays
	f_plus=0._wp
	f_plus_mon=0._wp
	f_minus=0._wp
	f_minus_mon=0._wp
	fj=0._wp

	psi_old(:)=psi(:) ! boundary conditions should be considered
	u1=u          ! boundary conditions should be considered


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate flux out of the right boundary                                 !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do j=0,kp+1
        ! f_plus jp+1 is referenced in the flux limiter 2nd do loop
        ! f_plus 0 is referenced because it it used to update field
		cp_j  =(u1(j)+abs(u1(j)))*dt/dx*0.5_wp
		cp_jp1=(u1(j+1)+abs(u1(j+1)))*dt/dx*0.5_wp
		cp_jm1=(u1(j-1)+abs(u1(j-1)))*dt/dx*0.5_wp
		cp_jm2=(u1(j-2)+abs(u1(j-2)))*dt/dx*0.5_wp
		cm_j  =-(u1(j)-abs(u1(j)))*dt/dx*0.5_wp
		cm_jm1=-(u1(j-1)-abs(u1(j-1)))*dt/dx*0.5_wp
		cm_jp1=-(u1(j+1)-abs(u1(j+1)))*dt/dx*0.5_wp

		if(monotone) then
			! deformation: eq. 17 bott, mwr (1992)
			dp_jm05=dx/dt*psi_old(j-1)*(cp_jm1-cp_jm2)
			dp_jp05=dx/dt*psi_old(j)*(cp_j-cp_jm1)
			dp_jp15=dx/dt*psi_old(j+1)*(cp_jp1-cp_j)
			dm_jp05=dx/dt*psi_old(j+1)*(cm_j-cm_jp1)
			dm_jm05=dx/dt*psi_old(j)*(cm_jm1-cm_j)
		endif

		! get coefficients for interpolation
		call coeff_bott_scheme_1d(psi_old,a_coeff01,j,ord,o_halo)
		call coeff_bott_scheme_1d(psi_old,a_coeff02,j-1,ord,o_halo)

		! step 1: calculate the monotone fluxes
		f_plus(j)=0._wp
		f_pm1=0._wp
		if(u1(j).gt.0._wp) then
			do k=0,ord
				f_plus(j)=f_plus(j)+dx/dt*a_coeff01(k+1)/ &
					((1._wp+real(k,wp))*(2._wp**(1._wp+real(k,wp))))* &
					(1._wp-(1._wp-2._wp*cp_j)**(1._wp+real(k,wp)) )
				f_pm1=f_pm1+dx/dt*a_coeff02(k+1)/ &
					((1._wp+real(k,wp))*(2._wp**(1._wp+real(k,wp))))* &
					(1._wp-(1._wp-2._wp*cp_jm1)**(1._wp+real(k,wp)) )
			enddo
			if(j.eq.0) then
				
				f_pm1=max(0._wp,f_pm1)
				f_pm1=min(dx/dt*psi_old(j-1),f_pm1)
				f_plus(j-1)=f_pm1
			endif
		endif

		if(monotone) then
			! step 2: apply the monotone flux limiter - eq 20, Bott 1992
			if((u1(j).gt.0._wp).and.(u1(j-1).gt.0._wp)) then
				dummy1=dx/dt * &
					(psi_old(j)-max(psi_old(j-1),psi_old(j)) )+f_plus(j-1)
				dummy2=dx/dt * &
					(psi_old(j)-min(psi_old(j-1),psi_old(j)) )+f_plus(j-1)
				  
				f_plus_mon(j)=max(dummy1,f_plus(j))
				f_plus_mon(j)=min(dummy2,f_plus_mon(j))
				f_plus_mon(j)=max(dummy1,f_plus_mon(j))
				!if(j.eq.0)  f_plus_mon(j)=f_plus(j)
			endif
			! step 3: add the deformation terms
			f_plus(j)=f_plus_mon(j)+dp_jp05
			! step 4: apply positive definite flux limiter (as in bott, 1989)
			f_plus(j)=max(0._wp,f_plus(j))
			f_plus(j)=min(dx/dt*psi_old(j),f_plus(j))
		else
			! step 4: apply positive definite flux limiter (as in bott, 1989)
			f_plus(j)=max(0._wp,f_plus(j))
			f_plus(j)=min(dx/dt*psi_old(j),f_plus(j))		
		endif
		
			
	enddo
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate flux into the boundary                                         !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	do j=kp+1,0,-1 ! flux into the boundary
		cp_j  =(u1(j)+abs(u1(j)))*dt/dx*0.5_wp
		cp_jp1=(u1(j+1)+abs(u1(j+1)))*dt/dx*0.5_wp
		cp_jm1=(u1(j-1)+abs(u1(j-1)))*dt/dx*0.5_wp
		cp_jm2=(u1(j-2)+abs(u1(j-2)))*dt/dx*0.5_wp
		cm_j  =-(u1(j)-abs(u1(j)))*dt/dx*0.5_wp
		cm_jm1=-(u1(j-1)-abs(u1(j-1)))*dt/dx*0.5_wp
		cm_jp1=-(u1(j+1)-abs(u1(j+1)))*dt/dx*0.5_wp

		if(monotone) then
			! deformation: eq. 17 bott, mwr (1992)
			dp_jp05=dx/dt*psi_old(j)*(cp_j-cp_jm1)
			dp_jp15=dx/dt*psi_old(j+1)*(cp_jp1-cp_j)
			dm_jp05=dx/dt*psi_old(j+1)*(cm_j-cm_jp1)
			dm_jm05=dx/dt*psi_old(j)*(cm_jm1-cm_j)
		endif

		! get coefficients for interpolation
		call coeff_bott_scheme_1d(psi_old,a_coeff01,j,ord,o_halo)
		call coeff_bott_scheme_1d(psi_old,a_coeff02,j+1,ord,o_halo)
		f_minus(j-1)=0._wp
		f_m=0._wp
		! step 1: calculate the monotone fluxes
		if(u1(j-1).lt.0._wp) then
			do k=0,ord
				f_minus(j-1)=f_minus(j-1)+dx/dt*a_coeff01(k+1)/ &
					((1._wp+real(k,wp))*(2._wp**(1._wp+real(k,wp))))*& 
					((-1._wp)**real(k,wp))*&
					(1._wp-(1._wp-2._wp*cm_jm1)**(1._wp+real(k,wp)))
				f_m=f_m+dx/dt*a_coeff02(k+1)/ &
					((1._wp+real(k,wp))*(2._wp**(1._wp+real(k,wp))))*& 
					((-1._wp)**real(k,wp))*&
					(1._wp-(1._wp-2._wp*cm_j)**(1._wp+real(k,wp)))
			enddo
			if(j.eq.(kp+1)) then
				f_m=max(0._wp,f_m)
				f_m=min(dx/dt*psi_old(j),f_m)
				f_minus(j)=f_m
			endif
		endif

		if(monotone) then
			! step 2: apply the monotone flux limiter - eq 20, Bott 1992
			if((u1(j).lt.0._wp).and.(u1(j-1).lt.0._wp)) then
				dummy1=dx/dt * &
					(psi_old(j)-max(psi_old(j+1),psi_old(j)) )+f_minus(j)
				dummy2=dx/dt * &
					(psi_old(j)-min(psi_old(j+1),psi_old(j)) )+f_minus(j)
				f_minus_mon(j-1)=min(dummy2,f_minus(j-1))
				f_minus_mon(j-1)=max(dummy1,f_minus_mon(j-1))
				f_minus_mon(j-1)=min(dummy2,f_minus_mon(j-1))
			endif
			! step 3: add the deformation terms
			f_minus(j-1)=f_minus_mon(j-1)+dm_jm05
			! step 4: apply positive definite flux limiter (as in bott, 1989)
			f_minus(j-1)=max(0._wp,f_minus(j-1))
			f_minus(j-1)=min(dx/dt*psi_old(j),f_minus(j-1))		
		else
			! step 4: apply positive definite flux limiter (as in bott, 1989)
			f_minus(j-1)=max(0._wp,f_minus(j-1))
			f_minus(j-1)=min(dx/dt*psi_old(j),f_minus(j-1))		
		endif
		
	enddo
	
	do j=kp+1,0,-1 
		! step 4: continued... apply the second condition of positive definite
		! flux limiter for divergent flows
		call coeff_bott_scheme_1d(psi_old,a_coeff01,j,ord,o_halo)
		call coeff_bott_scheme_1d(psi_old,a_coeff02,j+1,ord,o_halo)
		i_jp1=0._wp
		i_j=0._wp
		do k=0,ord
			i_j=i_j+ &
                dx/dt*a_coeff01(k+1) / &
                ( (1._wp+real(k,wp))* &
                (2._wp**(1._wp+real(k,wp))))*(((-1._wp)**real(k,wp))+1._wp)
			i_jp1=i_jp1+ &
                dx/dt*a_coeff02(k+1) / &
                ( (1._wp+real(k,wp))* &
                (2._wp**(1._wp+real(k,wp))))*(((-1._wp)**real(k,wp))+1._wp)
		enddo
		i_j=max(i_j,f_plus(j)+f_minus(j-1)+small)
		i_jp1=max(i_jp1,f_plus(j+1)+f_minus(j)+small)

		fj(j)=dx/dt*(f_plus(j)/i_j*psi_old(j)-f_minus(j)/i_jp1*psi_old(j+1))
! 		fj(j)=(f_plus(j)-f_minus(j))
	enddo

	! update the field!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	psi(1:kp) = psi(1:kp) - dt/dx*(fj(1:kp)-fj(0:kp-1)) ! eq. 2 Bott 1992

    where(psi.lt.0._wp)
		psi=1.e-60_wp
	end where
	! done!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	contains
	subroutine coeff_bott_scheme_1d(q,a_coeff,j,ord,o_halo)
		use numerics_type
		implicit none
		real(wp), intent(in), dimension(-o_halo+1:kp+1+o_halo) :: q
		real(wp), intent(inout), dimension(1:ord+1) :: a_coeff
		integer(i4b), intent(in) :: j,ord,o_halo

		! table 1: bott, mwr 1989.
		a_coeff(1)=q(j) ! same for all cases
		select case (ord)
            case (0) ! jp+1 and 0

            case (1) ! jp+2 and 0
                a_coeff(2) = q(j+1)-q(j)
            case (2) ! jp+2 and -1
                a_coeff(2) = 0.5_wp*(q(j+1)-q(j-1))
                a_coeff(3) = 0.5_wp*(q(j+1)-2._wp*q(j)+q(j-1))
                ! after bott's reply to smolarkiewicz's comment, 
                ! the zeroth term changed
                a_coeff(1) = -1._wp/24._wp*(q(j+1)-26._wp*q(j)+q(j-1))
            case (3) ! jp+3 and -1
                a_coeff(2) = (-q(j+2)+6._wp*q(j+1)-3._wp*q(j)-2._wp*q(j-1))/6._wp
                a_coeff(3) = (q(j+1)-2._wp*q(j)+q(j-1))/2._wp
                a_coeff(4) = (q(j+2)-3._wp*q(j+1)+3._wp*q(j)-q(j-1))/6._wp
            case (4) ! jp+3 and -2
                a_coeff(2) = (-q(j+2)+8._wp*q(j+1)-8._wp*q(j-1)+q(j-2))/12._wp
                a_coeff(3) = (-q(j+2)+16._wp*q(j+1)-30._wp*q(j) &
                                +16.*q(j-1)-q(j-2))/24._wp
                a_coeff(4) = (q(j+2)-2._wp*q(j+1)+2._wp*q(j-1)-q(j-2))/12._wp
                a_coeff(5) = (q(j+2)-4._wp*q(j+1)+6._wp*q(j) &
                                -4._wp*q(j-1)+q(j-2))/24._wp
                ! after bott's reply, the 0th, 1st, 2nd term changed
                a_coeff(1) = 1._wp/1920._wp*(9.*q(j+2)-116.*q(j+1) &
                            +2134.*q(j)-116.*q(j-1)+9.*q(j-2))
                a_coeff(2) = 1._wp/48._wp*(-5._wp*q(j+2)+34._wp*q(j+1) &
                            -34._wp*q(j-1)+5._wp*q(j-2))
                a_coeff(3) = 1._wp/48._wp*(-3._wp*q(j+2)+36._wp*q(j+1) &
                            -66._wp*q(j)+36._wp*q(j-1)-3._wp*q(j-2))
                !!!!!!! see costa and sampaio mwr 1997: higher order bott coefficients
            case (5) ! jp+4 and -3
                a_coeff(1) = 1._wp/1920._wp*(9._wp*q(j+2)-116._wp*q(j+1)+2134._wp*q(j)-116._wp*q(j-1)+9._wp*q(j-2))
                a_coeff(2) = 1._wp/11520._wp*(259._wp*q(j+3)-2236._wp*q(j+2) &
                            +9455._wp*q(j+1)-9455.*q(j-1)+2236.*q(j-2)-259.*q(j-3))

                a_coeff(3) = 1._wp/16._wp*(-q(j+2)+12._wp*q(j+1) &
                            -22._wp*q(j)+12._wp*q(j-1)-q(j-2))
                a_coeff(4) = 1._wp/288._wp*(-7*q(j+3)-52._wp*q(j+2) &
                            -83._wp*q(j+1)+83._wp*q(j-1)-52._wp*q(j-2)+7._wp*q(j-3))

                a_coeff(5) = 1._wp/24._wp*(q(j+2)-4._wp*q(j+1) &
                            +6._wp*q(j)-4._wp*q(j-1)+q(j-2))
                a_coeff(6) = 1._wp/240._wp*(q(j+3)-4._wp*q(j+2) &
                            +5._wp*q(j+1)-5._wp*q(j-1)+4._wp*q(j-2)-q(j-3))

            case (6) ! jp+4 and -3
                a_coeff(1) = 1._wp/107520._wp*(-75._wp*q(j+3)+&
                            954._wp*q(j+2)-7621._wp*q(j+1)+121004._wp*q(j)-&
                            7621._wp*q(j-1)+954._wp*q(j-2)-75._wp*q(j-3))

                a_coeff(2) = 1._wp/11520._wp*(259._wp*q(j+3)-2236._wp*q(j+2)+ &
                            9455._wp*q(j+1)-9455._wp*q(j-1)+ &
                            2236._wp*q(j-2)-259._wp*q(j-3))

                a_coeff(3) = 1._wp/3840._wp*(37._wp*q(j+3)-462._wp*q(j+2)+ &
                            3435._wp*q(j+1)-6020._wp*q(j)+&
                            3435._wp*q(j-1)-462._wp*q(j-2)+37._wp*q(j-3))

                a_coeff(4) = 1._wp/288._wp*(-7._wp*q(j+3)+52._wp*q(j+2)- &
                            83._wp*q(j+1)+83._wp*q(j-1)-52._wp*q(j-2)+7._wp*q(j-3))

                a_coeff(5) = 1._wp/576._wp*(-5._wp*q(j+3)+54._wp*q(j+2)- &
                            171._wp*q(j+1)+244._wp*q(j)-&
                            171._wp*q(j-1)+54._wp*q(j-2)-5._wp*q(j-3))

                a_coeff(6) = 1._wp/240._wp*(q(j+3)-4._wp*q(j+2)+ &
                            5._wp*q(j+1)-5._wp*q(j-1)+4._wp*q(j-2)-q(j-3))

                a_coeff(7) = 1._wp/720._wp*(q(j+3)-6._wp*q(j+2)+&
                            15._wp*q(j+1)-20._wp*q(j)+&
                            15._wp*q(j-1)-6._wp*q(j-2)+q(j-3))
            case (7) ! jp+5 and -4
                a_coeff(1) = 1._wp/107520._wp*(-75._wp*q(j+3)+ &
                            954._wp*q(j+2)-7621._wp*q(j+1)+ &
                            121004._wp*q(j)-7621._wp*q(j-1)+ &
                            954._wp*q(j-2)-75._wp*q(j-3))

                a_coeff(2) = 1._wp/645120._wp*(-3229._wp*q(j+4)+&
                            33878._wp*q(j+3)-170422._wp*q(j+2)+ &
                            574686._wp*q(j+1)-574686._wp*q(j-1)+170433._wp*q(j-2)-&
                            33878._wp*q(j-3)+3229._wp*q(j-4))

                a_coeff(3) = 1._wp/3840._wp*(37*q(j+3)-462._wp*q(j+2)+ &
                            3435._wp*q(j+1)-6020._wp*q(j)+ &
                            3435._wp*q(j-1)-462._wp*q(j-2)+37._wp*q(j-3))

                a_coeff(4) = 1._wp/23040._wp*(141._wp*q(j+4)-&
                            1406*q(j+3)+6134._wp*q(j+2)-&
                            8614._wp*q(j+1)+8614._wp*q(j-1)-6134._wp*q(j-2)&
                            +1406._wp*q(j-3)-141._wp*q(j-4))

                a_coeff(5) = 1._wp/576._wp*(-5._wp*q(j+3)+&
                            54._wp*q(j+2)-171._wp*q(j+1)+ &
                            244._wp*q(j)-171._wp*q(j-1)+54._wp*q(j-2)-5._wp*q(j-3))

                a_coeff(6) = 1._wp/1920._wp*(-3._wp*q(j+4)+26._wp*q(j+3)-&
                                74._wp*q(j+2)+82._wp*q(j+1)- &
                            82._wp*q(j-1)+74._wp*q(j-2)-26._wp*q(j-3)+3._wp*q(j-4))

                a_coeff(7) = 1._wp/720._wp*(q(j+3)-6._wp*q(j+2)+&
                            15._wp*q(j+1)-20._wp*q(j)+ &
                            15._wp*q(j-1)-6._wp*q(j-2)+q(j-3))
                a_coeff(8) = 1._wp/10080._wp*(q(j+4)-6._wp*q(j+3)+&
                            14._wp*q(j+2)-14._wp*q(j+1)+&
                            14._wp*q(j-1)-14._wp*q(j-2)+6._wp*q(j-3)-q(j-4))
            case (8) ! jp+5 and -4
                a_coeff(1) = 1._wp/10321920._wp*(1225._wp*q(j+4)-&
                            17000._wp*q(j+3)+125884._wp*q(j+2)- &
                            800216._wp*q(j+1)+11702134._wp*q(j)-&
                            800216._wp*q(j-1)+125884._wp*q(j-2)-&
                            17000._wp*q(j-3)+1225._wp*q(j-4))

                a_coeff(2) = 1._wp/645120._wp*(-3229._wp*q(j+4)+ &
                            33878._wp*q(j+3)-170422._wp*q(j+2)+ &
                            574686._wp*q(j+1)-574686._wp*q(j-1)+ &
                            170433._wp*q(j-2)-33878._wp*q(j-3)+3229._wp*q(j-4))

                a_coeff(3) = 1._wp/1935360._wp*(-3229._wp*q(j+4)+ &
                            44480._wp*q(j+3)-323260._wp*q(j+2)+&
                            1912064._wp*q(j+1)-3260110._wp*q(j)+ &
                            1912064._wp*q(j-1)-323260._wp*q(j-2)+&
                            44480._wp*q(j-3)-3229._wp*q(j-4))

                a_coeff(4) = 1._wp/23040._wp*(141._wp*q(j+4)- &
                            1406*q(j+3)+6134._wp*q(j+2)-&
                            8614._wp*q(j+1)+8614._wp*q(j-1)-6134._wp*q(j-2)+&
                            1406._wp*q(j-3)-141._wp*q(j-4))

                a_coeff(5) = 1._wp/27648._wp*(47._wp*q(j+4)-&
                            616._wp*q(j+3)+3908._wp*q(j+2)-&
                            10840._wp*q(j+1)+15002._wp*q(j)-10840._wp*q(j-1)+&
                            3908._wp*q(j-2)-616._wp*q(j-3)+47._wp*q(j-4))

                a_coeff(6) = 1._wp/1920._wp*(-3._wp*q(j+4)+26._wp*q(j+3)-&
                            74._wp*q(j+2)+82._wp*q(j+1)-&
                            82._wp*q(j-1)+74._wp*q(j-2)-26._wp*q(j-3)+3._wp*q(j-4))

                a_coeff(7) = 1._wp/17280._wp*(-7._wp*q(j+4)+80._wp*q(j+3)-&
                            340._wp*q(j+2)+752._wp*q(j+1)-&
                            970._wp*q(j)+752._wp*q(j-1)-340._wp*q(j-2)+&
                            80._wp*q(j-3)-7._wp*q(j-4))

                a_coeff(8) = 1._wp/10080._wp*(q(j+4)-6._wp*q(j+3)+&
                            14._wp*q(j+2)-14._wp*q(j+1)+&
                            14._wp*q(j-1)-14._wp*q(j-2)+6._wp*q(j-3)-q(j-4))
                a_coeff(9) = 1._wp/40320._wp*(q(j+4)-8._wp*q(j+3)+&
                            28._wp*q(j+2)-56._wp*q(j+1)+&
                            70._wp*q(j)-56._wp*q(j-1)+28._wp*q(j-2)-&
                            8._wp*q(j-3)+q(j-4))
            case default
                print *,'not defined for positive definite apf'
                stop
		end select


	end subroutine coeff_bott_scheme_1d
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    end subroutine bott_scheme_1d

    end module advection_1d
