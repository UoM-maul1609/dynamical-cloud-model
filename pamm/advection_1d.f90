!>@author
!>Paul Connolly, The University of Manchester
!>@brief
!>advection routines
    module advection_1d
    use nrtype
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

    use nrtype
    implicit none
    ! arguments:
    integer(i4b), intent(in) :: kp, ord, o_halo
    real(sp), intent(in) :: dt,dx
    real(sp), dimension(-o_halo+1:kp+o_halo), intent(in) :: x,u
    logical, intent(in) :: monotone
    real(sp), dimension(-o_halo+1:kp+o_halo), intent(inout) :: psi

    ! locals:
	real(sp), dimension(ord+1) :: a_coeff01, a_coeff02
	real(sp), dimension(-o_halo+1:kp+o_halo) :: f_plus, f_plus_mon,&
									f_minus,f_minus_mon, fj, psi_old, u1
	real(sp) :: cp_j, cp_jp1, cp_jm1, cp_jm2, cm_j, cm_jm1, cm_jp1, &
	dp_jm05=0._sp, dp_jp05=0._sp, dp_jp15=0._sp, dm_jp05=0._sp, dm_jm05=0._sp, &
	i_j, i_jp1, f_pm1, f_m
	real(sp) :: dummy1, dummy2,small=1e-60_sp
	integer(i4b) :: j, k

	if(sum(psi).lt.small) then
		return
	endif
	! zero arrays
	f_plus=0._sp
	f_plus_mon=0._sp
	f_minus=0._sp
	f_minus_mon=0._sp
	fj=0._sp

	psi_old(:)=psi(:) ! boundary conditions should be considered
	u1=u          ! boundary conditions should be considered


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate flux out of the right boundary                                 !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do j=0,kp+1
        ! f_plus jp+1 is referenced in the flux limiter 2nd do loop
        ! f_plus 0 is referenced because it it used to update field
		cp_j  =(u1(j)+abs(u1(j)))*dt/dx*0.5_sp
		cp_jp1=(u1(j+1)+abs(u1(j+1)))*dt/dx*0.5_sp
		cp_jm1=(u1(j-1)+abs(u1(j-1)))*dt/dx*0.5_sp
		cp_jm2=(u1(j-2)+abs(u1(j-2)))*dt/dx*0.5_sp
		cm_j  =-(u1(j)-abs(u1(j)))*dt/dx*0.5_sp
		cm_jm1=-(u1(j-1)-abs(u1(j-1)))*dt/dx*0.5_sp
		cm_jp1=-(u1(j+1)-abs(u1(j+1)))*dt/dx*0.5_sp

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
		f_plus(j)=0._sp
		f_pm1=0._sp
		if(u1(j).gt.0._sp) then
			do k=0,ord
				f_plus(j)=f_plus(j)+dx/dt*a_coeff01(k+1)/ &
					((1._sp+real(k,sp))*(2._sp**(1._sp+real(k,sp))))* &
					(1._sp-(1._sp-2._sp*cp_j)**(1._sp+real(k,sp)) )
				f_pm1=f_pm1+dx/dt*a_coeff02(k+1)/ &
					((1._sp+real(k,sp))*(2._sp**(1._sp+real(k,sp))))* &
					(1._sp-(1._sp-2._sp*cp_jm1)**(1._sp+real(k,sp)) )
			enddo
			if(j.eq.0) then
				
				f_pm1=max(0._sp,f_pm1)
				f_pm1=min(dx/dt*psi_old(j-1),f_pm1)
				f_plus(j-1)=f_pm1
			endif
		endif

		if(monotone) then
			! step 2: apply the monotone flux limiter - eq 20, Bott 1992
			if((u1(j).gt.0._sp).and.(u1(j-1).gt.0._sp)) then
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
			f_plus(j)=max(0._sp,f_plus(j))
			f_plus(j)=min(dx/dt*psi_old(j),f_plus(j))
		else
			! step 4: apply positive definite flux limiter (as in bott, 1989)
			f_plus(j)=max(0._sp,f_plus(j))
			f_plus(j)=min(dx/dt*psi_old(j),f_plus(j))		
		endif
		
			
	enddo
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate flux into the boundary                                         !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	do j=kp+1,0,-1 ! flux into the boundary
		cp_j  =(u1(j)+abs(u1(j)))*dt/dx*0.5_sp
		cp_jp1=(u1(j+1)+abs(u1(j+1)))*dt/dx*0.5_sp
		cp_jm1=(u1(j-1)+abs(u1(j-1)))*dt/dx*0.5_sp
		cp_jm2=(u1(j-2)+abs(u1(j-2)))*dt/dx*0.5_sp
		cm_j  =-(u1(j)-abs(u1(j)))*dt/dx*0.5_sp
		cm_jm1=-(u1(j-1)-abs(u1(j-1)))*dt/dx*0.5_sp
		cm_jp1=-(u1(j+1)-abs(u1(j+1)))*dt/dx*0.5_sp

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
		f_minus(j-1)=0._sp
		f_m=0._sp
		! step 1: calculate the monotone fluxes
		if(u1(j-1).lt.0._sp) then
			do k=0,ord
				f_minus(j-1)=f_minus(j-1)+dx/dt*a_coeff01(k+1)/ &
					((1._sp+real(k,sp))*(2._sp**(1._sp+real(k,sp))))*& 
					((-1._sp)**real(k,sp))*&
					(1._sp-(1._sp-2._sp*cm_jm1)**(1._sp+real(k,sp)))
				f_m=f_m+dx/dt*a_coeff02(k+1)/ &
					((1._sp+real(k,sp))*(2._sp**(1._sp+real(k,sp))))*& 
					((-1._sp)**real(k,sp))*&
					(1._sp-(1._sp-2._sp*cm_j)**(1._sp+real(k,sp)))
			enddo
			if(j.eq.(kp+1)) then
				f_m=max(0._sp,f_m)
				f_m=min(dx/dt*psi_old(j),f_m)
				f_minus(j)=f_m
			endif
		endif

		if(monotone) then
			! step 2: apply the monotone flux limiter - eq 20, Bott 1992
			if((u1(j).lt.0._sp).and.(u1(j-1).lt.0._sp)) then
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
			f_minus(j-1)=max(0._sp,f_minus(j-1))
			f_minus(j-1)=min(dx/dt*psi_old(j),f_minus(j-1))		
		else
			! step 4: apply positive definite flux limiter (as in bott, 1989)
			f_minus(j-1)=max(0._sp,f_minus(j-1))
			f_minus(j-1)=min(dx/dt*psi_old(j),f_minus(j-1))		
		endif
		
	enddo
	
	do j=kp+1,0,-1 
		! step 4: continued... apply the second condition of positive definite
		! flux limiter for divergent flows
		call coeff_bott_scheme_1d(psi_old,a_coeff01,j,ord,o_halo)
		call coeff_bott_scheme_1d(psi_old,a_coeff02,j+1,ord,o_halo)
		i_jp1=0._sp
		i_j=0._sp
		do k=0,ord
			i_j=i_j+ &
                dx/dt*a_coeff01(k+1) / &
                ( (1._sp+real(k,sp))* &
                (2._sp**(1._sp+real(k,sp))))*(((-1._sp)**real(k,sp))+1._sp)
			i_jp1=i_jp1+ &
                dx/dt*a_coeff02(k+1) / &
                ( (1._sp+real(k,sp))* &
                (2._sp**(1._sp+real(k,sp))))*(((-1._sp)**real(k,sp))+1._sp)
		enddo
		i_j=max(i_j,f_plus(j)+f_minus(j-1)+small)
		i_jp1=max(i_jp1,f_plus(j+1)+f_minus(j)+small)

		fj(j)=dx/dt*(f_plus(j)/i_j*psi_old(j)-f_minus(j)/i_jp1*psi_old(j+1))
! 		fj(j)=(f_plus(j)-f_minus(j))
	enddo

	! update the field!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	psi(1:kp) = psi(1:kp) - dt/dx*(fj(1:kp)-fj(0:kp-1)) ! eq. 2 Bott 1992

    where(psi.lt.0._sp)
		psi=1.e-60_sp
	end where
	! done!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	contains
	subroutine coeff_bott_scheme_1d(q,a_coeff,j,ord,o_halo)
		use nrtype
		implicit none
		real(sp), intent(in), dimension(-o_halo+1:kp+1+o_halo) :: q
		real(sp), intent(inout), dimension(1:ord+1) :: a_coeff
		integer(i4b), intent(in) :: j,ord,o_halo

		! table 1: bott, mwr 1989.
		a_coeff(1)=q(j) ! same for all cases
		select case (ord)
            case (0) ! jp+1 and 0

            case (1) ! jp+2 and 0
                a_coeff(2) = q(j+1)-q(j)
            case (2) ! jp+2 and -1
                a_coeff(2) = 0.5_sp*(q(j+1)-q(j-1))
                a_coeff(3) = 0.5_sp*(q(j+1)-2._sp*q(j)+q(j-1))
                ! after bott's reply to smolarkiewicz's comment, 
                ! the zeroth term changed
                a_coeff(1) = -1._sp/24._sp*(q(j+1)-26._sp*q(j)+q(j-1))
            case (3) ! jp+3 and -1
                a_coeff(2) = (-q(j+2)+6._sp*q(j+1)-3._sp*q(j)-2._sp*q(j-1))/6._sp
                a_coeff(3) = (q(j+1)-2._sp*q(j)+q(j-1))/2._sp
                a_coeff(4) = (q(j+2)-3._sp*q(j+1)+3._sp*q(j)-q(j-1))/6._sp
            case (4) ! jp+3 and -2
                a_coeff(2) = (-q(j+2)+8._sp*q(j+1)-8._sp*q(j-1)+q(j-2))/12._sp
                a_coeff(3) = (-q(j+2)+16._sp*q(j+1)-30._sp*q(j) &
                                +16.*q(j-1)-q(j-2))/24._sp
                a_coeff(4) = (q(j+2)-2._sp*q(j+1)+2._sp*q(j-1)-q(j-2))/12._sp
                a_coeff(5) = (q(j+2)-4._sp*q(j+1)+6._sp*q(j) &
                                -4._sp*q(j-1)+q(j-2))/24._sp
                ! after bott's reply, the 0th, 1st, 2nd term changed
                a_coeff(1) = 1._sp/1920._sp*(9.*q(j+2)-116.*q(j+1) &
                            +2134.*q(j)-116.*q(j-1)+9.*q(j-2))
                a_coeff(2) = 1._sp/48._sp*(-5._sp*q(j+2)+34._sp*q(j+1) &
                            -34._sp*q(j-1)+5._sp*q(j-2))
                a_coeff(3) = 1._sp/48._sp*(-3._sp*q(j+2)+36._sp*q(j+1) &
                            -66._sp*q(j)+36._sp*q(j-1)-3._sp*q(j-2))
                !!!!!!! see costa and sampaio mwr 1997: higher order bott coefficients
            case (5) ! jp+4 and -3
                a_coeff(1) = 1._sp/1920._sp*(9._sp*q(j+2)-116._sp*q(j+1)+2134._sp*q(j)-116._sp*q(j-1)+9._sp*q(j-2))
                a_coeff(2) = 1._sp/11520._sp*(259._sp*q(j+3)-2236._sp*q(j+2) &
                            +9455._sp*q(j+1)-9455.*q(j-1)+2236.*q(j-2)-259.*q(j-3))

                a_coeff(3) = 1._sp/16._sp*(-q(j+2)+12._sp*q(j+1) &
                            -22._sp*q(j)+12._sp*q(j-1)-q(j-2))
                a_coeff(4) = 1._sp/288._sp*(-7*q(j+3)-52._sp*q(j+2) &
                            -83._sp*q(j+1)+83._sp*q(j-1)-52._sp*q(j-2)+7._sp*q(j-3))

                a_coeff(5) = 1._sp/24._sp*(q(j+2)-4._sp*q(j+1) &
                            +6._sp*q(j)-4._sp*q(j-1)+q(j-2))
                a_coeff(6) = 1._sp/240._sp*(q(j+3)-4._sp*q(j+2) &
                            +5._sp*q(j+1)-5._sp*q(j-1)+4._sp*q(j-2)-q(j-3))

            case (6) ! jp+4 and -3
                a_coeff(1) = 1._sp/107520._sp*(-75._sp*q(j+3)+&
                            954._sp*q(j+2)-7621._sp*q(j+1)+121004._sp*q(j)-&
                            7621._sp*q(j-1)+954._sp*q(j-2)-75._sp*q(j-3))

                a_coeff(2) = 1._sp/11520._sp*(259._sp*q(j+3)-2236._sp*q(j+2)+ &
                            9455._sp*q(j+1)-9455._sp*q(j-1)+ &
                            2236._sp*q(j-2)-259._sp*q(j-3))

                a_coeff(3) = 1._sp/3840._sp*(37._sp*q(j+3)-462._sp*q(j+2)+ &
                            3435._sp*q(j+1)-6020._sp*q(j)+&
                            3435._sp*q(j-1)-462._sp*q(j-2)+37._sp*q(j-3))

                a_coeff(4) = 1._sp/288._sp*(-7._sp*q(j+3)+52._sp*q(j+2)- &
                            83._sp*q(j+1)+83._sp*q(j-1)-52._sp*q(j-2)+7._sp*q(j-3))

                a_coeff(5) = 1._sp/576._sp*(-5._sp*q(j+3)+54._sp*q(j+2)- &
                            171._sp*q(j+1)+244._sp*q(j)-&
                            171._sp*q(j-1)+54._sp*q(j-2)-5._sp*q(j-3))

                a_coeff(6) = 1._sp/240._sp*(q(j+3)-4._sp*q(j+2)+ &
                            5._sp*q(j+1)-5._sp*q(j-1)+4._sp*q(j-2)-q(j-3))

                a_coeff(7) = 1._sp/720._sp*(q(j+3)-6._sp*q(j+2)+&
                            15._sp*q(j+1)-20._sp*q(j)+&
                            15._sp*q(j-1)-6._sp*q(j-2)+q(j-3))
            case (7) ! jp+5 and -4
                a_coeff(1) = 1._sp/107520._sp*(-75._sp*q(j+3)+ &
                            954._sp*q(j+2)-7621._sp*q(j+1)+ &
                            121004._sp*q(j)-7621._sp*q(j-1)+ &
                            954._sp*q(j-2)-75._sp*q(j-3))

                a_coeff(2) = 1._sp/645120._sp*(-3229._sp*q(j+4)+&
                            33878._sp*q(j+3)-170422._sp*q(j+2)+ &
                            574686._sp*q(j+1)-574686._sp*q(j-1)+170433._sp*q(j-2)-&
                            33878._sp*q(j-3)+3229._sp*q(j-4))

                a_coeff(3) = 1._sp/3840._sp*(37*q(j+3)-462._sp*q(j+2)+ &
                            3435._sp*q(j+1)-6020._sp*q(j)+ &
                            3435._sp*q(j-1)-462._sp*q(j-2)+37._sp*q(j-3))

                a_coeff(4) = 1._sp/23040._sp*(141._sp*q(j+4)-&
                            1406*q(j+3)+6134._sp*q(j+2)-&
                            8614._sp*q(j+1)+8614._sp*q(j-1)-6134._sp*q(j-2)&
                            +1406._sp*q(j-3)-141._sp*q(j-4))

                a_coeff(5) = 1._sp/576._sp*(-5._sp*q(j+3)+&
                            54._sp*q(j+2)-171._sp*q(j+1)+ &
                            244._sp*q(j)-171._sp*q(j-1)+54._sp*q(j-2)-5._sp*q(j-3))

                a_coeff(6) = 1._sp/1920._sp*(-3._sp*q(j+4)+26._sp*q(j+3)-&
                                74._sp*q(j+2)+82._sp*q(j+1)- &
                            82._sp*q(j-1)+74._sp*q(j-2)-26._sp*q(j-3)+3._sp*q(j-4))

                a_coeff(7) = 1._sp/720._sp*(q(j+3)-6._sp*q(j+2)+&
                            15._sp*q(j+1)-20._sp*q(j)+ &
                            15._sp*q(j-1)-6._sp*q(j-2)+q(j-3))
                a_coeff(8) = 1._sp/10080._sp*(q(j+4)-6._sp*q(j+3)+&
                            14._sp*q(j+2)-14._sp*q(j+1)+&
                            14._sp*q(j-1)-14._sp*q(j-2)+6._sp*q(j-3)-q(j-4))
            case (8) ! jp+5 and -4
                a_coeff(1) = 1._sp/10321920._sp*(1225._sp*q(j+4)-&
                            17000._sp*q(j+3)+125884._sp*q(j+2)- &
                            800216._sp*q(j+1)+11702134._sp*q(j)-&
                            800216._sp*q(j-1)+125884._sp*q(j-2)-&
                            17000._sp*q(j-3)+1225._sp*q(j-4))

                a_coeff(2) = 1._sp/645120._sp*(-3229._sp*q(j+4)+ &
                            33878._sp*q(j+3)-170422._sp*q(j+2)+ &
                            574686._sp*q(j+1)-574686._sp*q(j-1)+ &
                            170433._sp*q(j-2)-33878._sp*q(j-3)+3229._sp*q(j-4))

                a_coeff(3) = 1._sp/1935360._sp*(-3229._sp*q(j+4)+ &
                            44480._sp*q(j+3)-323260._sp*q(j+2)+&
                            1912064._sp*q(j+1)-3260110._sp*q(j)+ &
                            1912064._sp*q(j-1)-323260._sp*q(j-2)+&
                            44480._sp*q(j-3)-3229._sp*q(j-4))

                a_coeff(4) = 1._sp/23040._sp*(141._sp*q(j+4)- &
                            1406*q(j+3)+6134._sp*q(j+2)-&
                            8614._sp*q(j+1)+8614._sp*q(j-1)-6134._sp*q(j-2)+&
                            1406._sp*q(j-3)-141._sp*q(j-4))

                a_coeff(5) = 1._sp/27648._sp*(47._sp*q(j+4)-&
                            616._sp*q(j+3)+3908._sp*q(j+2)-&
                            10840._sp*q(j+1)+15002._sp*q(j)-10840._sp*q(j-1)+&
                            3908._sp*q(j-2)-616._sp*q(j-3)+47._sp*q(j-4))

                a_coeff(6) = 1._sp/1920._sp*(-3._sp*q(j+4)+26._sp*q(j+3)-&
                            74._sp*q(j+2)+82._sp*q(j+1)-&
                            82._sp*q(j-1)+74._sp*q(j-2)-26._sp*q(j-3)+3._sp*q(j-4))

                a_coeff(7) = 1._sp/17280._sp*(-7._sp*q(j+4)+80._sp*q(j+3)-&
                            340._sp*q(j+2)+752._sp*q(j+1)-&
                            970._sp*q(j)+752._sp*q(j-1)-340._sp*q(j-2)+&
                            80._sp*q(j-3)-7._sp*q(j-4))

                a_coeff(8) = 1._sp/10080._sp*(q(j+4)-6._sp*q(j+3)+&
                            14._sp*q(j+2)-14._sp*q(j+1)+&
                            14._sp*q(j-1)-14._sp*q(j-2)+6._sp*q(j-3)-q(j-4))
                a_coeff(9) = 1._sp/40320._sp*(q(j+4)-8._sp*q(j+3)+&
                            28._sp*q(j+2)-56._sp*q(j+1)+&
                            70._sp*q(j)-56._sp*q(j-1)+28._sp*q(j-2)-&
                            8._sp*q(j-3)+q(j-4))
            case default
                print *,'not defined for positive definite apf'
                stop
		end select


	end subroutine coeff_bott_scheme_1d
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    end subroutine bott_scheme_1d

    end module advection_1d
