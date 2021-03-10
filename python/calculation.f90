! subroutine sieve(nus,kabs, len1,len2,A,MOLW,S,GAMMA_Q,NUQ,DELQ,PATM, &
!                 nulow,deltanu)
! 
! ! =====================================================
! ! Uses the sieve of Eratosthenes to compute a logical
! ! array of size n_max, where .true. in element i
! ! indicates that i is a prime.
! ! =====================================================
!     implicit none
!     integer, intent(in)   :: len1,len2
!     real(8), intent(out), dimension(len2)  :: nus
! !f2py intent(out) :: nus
!     real(8), intent(out), dimension(len2)  :: kabs
! !f2py intent(out) :: kabs
!     real(8), intent(in), dimension(len1)  :: MOLW,S,GAMMA_Q,NUQ,DELQ
!     real(8), intent(in) :: PATM, A, nulow,deltanu
!     integer :: i,j
!     real(8) :: pi=4.*atan(1.),var,var2,var3,var4,gami2
!     integer :: nthreads, tid, omp_get_num_threads, omp_get_thread_num
!     real(8), dimension(len2) :: kabspartial
!     
!     var=A*pi
!     var2=1.*PATM
! 
!     !print *,len1,len2
!     !kabspartial=0.
!     kabs=0.
!     do j = 1, len2
!         nus(j)=nulow+deltanu*real(j)
!     end do
! 
! 
!     !$OMP PARALLEL private(i,j,nthreads,tid,nus,var3,var,gami2,var4,kabspartial) 
!     tid = omp_get_thread_num()
!     !write(*,*) 'Hello from thread = ', tid
!     if (tid.eq.0) then
!         nthreads = omp_get_num_threads()
!         write(*,*) 'Number of threads = ', nthreads
!     end if
!     !$OMP DO
!     do i = 1, len1 ! over all bands
!         var3=var/MOLW(i)*S(i)*GAMMA_Q(i)
!         gami2=GAMMA_Q(i)**2
!         var4=NUQ(i)-DELQ(i)*var2
!         
!         ! var4 is the center
!         ! 2*GAMMA_Q(i) is FWHM
!         ! find 
!         
!         do j = 1, 1 ! over all nu
!             kabs(j)=kabs(j)+ &
!                 (var3/(gami2+(NUS(j)-(var4))**2))
!         end do
!     end do
!     !$OMP END DO
!     
!     
!     !$OMP CRITICAL
!     !kabs = kabs+kabspartial
!     !$OMP END CRITICAL
!     
!     !$OMP END PARALLEL
!     
!     
!     
!     return
! end subroutine

subroutine absorption_calc(nus,kabs, len1,len2,A,&
                T,Tref,c2,nq,gamma_air,gamma_self,mr,LI,Q1,Q2,LSE,&
                MOLW,NUQ,DELQ,PATM, &
                nulow,deltanu)

! =====================================================
! Uses the sieve of Eratosthenes to compute a logical
! array of size n_max, where .true. in element i
! indicates that i is a prime.
! =====================================================
    implicit none
    integer, intent(in)   :: len1,len2
    real(8), intent(out), dimension(len2)  :: nus
!f2py intent(out) :: nus
    real(8), intent(out), dimension(len2)  :: kabs
!f2py intent(out) :: kabs
    real(8), intent(in), dimension(len1)  :: MOLW,NUQ,DELQ, nq, gamma_air, gamma_self, &
                                            mr, LI, Q1, Q2, LSE
    real(8), intent(in) :: PATM, A, nulow,deltanu, T, Tref,c2
    integer :: i,j
    real(8) :: pi=4.*atan(1.),var,var2,var3,var4,gami2
    real(8), dimension(len1) :: gamma_q,S
    
    do i = 1, len1
        gamma_q(i)=(Tref/T)**nq(i)* &
            (gamma_air(i)*(Patm-Patm*mr(i))+gamma_self(i)*Patm*mr(i) )
        S(i)=LI(i)*Q1(i)/Q2(i)*exp(-c2*LSE(i)/T)*(1-exp(-c2*nuq(i)/T)) / &
            (exp(-c2*LSE(i)/Tref)*(1-exp(-c2*nuq(i)/Tref)))
    enddo
    
    
    
    var=A/pi
    var2=PATM

    kabs=0.
    do j = 1, len2
        nus(j)=nulow+deltanu*real(j)
    end do



    do i = 1, len1 ! over all bands
        var3=var/MOLW(i)*S(i)*GAMMA_Q(i)
        gami2=GAMMA_Q(i)**2
        var4=NUQ(i)-DELQ(i)*var2
        
        ! var4 is the center
        ! 2*GAMMA_Q(i)^2 is FWHM
        ! find nus where within 3xFWHM to try and speed up?
!         where ((nus.lt.(var4+abs(1.e6*GAMMA_Q(i)))).and. &
!             (nus.gt.(var4-abs(1.e6*GAMMA_Q(i)))))
!             
!             kabs(:)=kabs(:)+ &
!                 (var3/(gami2+(NUS(:)-(var4))**2))
!         end where
        
        do j = 1, len2 ! over all nu
            kabs(j)=kabs(j)+ &
                (var3/(gami2+(NUS(j)-(var4))**2))
        end do
    end do
        
    return
end subroutine