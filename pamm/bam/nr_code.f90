	! module for nrtype
	MODULE nrtype1
		INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
		INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
		INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)
		INTEGER, PARAMETER :: SP = KIND(1.0D0) ! change to the same as double?
		INTEGER, PARAMETER :: DP = KIND(1.0D0)
		INTEGER, PARAMETER :: SPC = KIND((1.0,1.0))
		INTEGER, PARAMETER :: DPC = KIND((1.0D0,1.0D0))
		INTEGER, PARAMETER :: LGT = KIND(.true.)
		REAL(SP), PARAMETER :: PI=3.141592653589793238462643383279502884197_sp
		REAL(SP), PARAMETER :: PIO2=1.57079632679489661923132169163975144209858_sp
		REAL(SP), PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_sp
		REAL(SP), PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967_sp
		REAL(SP), PARAMETER :: EULER=0.5772156649015328606065120900824024310422_sp
		REAL(DP), PARAMETER :: PI_D=3.141592653589793238462643383279502884197_dp
		REAL(DP), PARAMETER :: PIO2_D=1.57079632679489661923132169163975144209858_dp
		REAL(DP), PARAMETER :: TWOPI_D=6.283185307179586476925286766559005768394_dp
		TYPE sprs2_sp
			INTEGER(I4B) :: n,len
			REAL(SP), DIMENSION(:), POINTER :: val
			INTEGER(I4B), DIMENSION(:), POINTER :: irow
			INTEGER(I4B), DIMENSION(:), POINTER :: jcol
		END TYPE sprs2_sp
		TYPE sprs2_dp
			INTEGER(I4B) :: n,len
			REAL(DP), DIMENSION(:), POINTER :: val
			INTEGER(I4B), DIMENSION(:), POINTER :: irow
			INTEGER(I4B), DIMENSION(:), POINTER :: jcol
		END TYPE sprs2_dp
	END MODULE nrtype1


	! module for nrutil1
	MODULE nrutil1
		USE nrtype1
		IMPLICIT NONE
		INTEGER(I4B), PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8
		INTEGER(I4B), PARAMETER :: NPAR_GEOP=4,NPAR2_GEOP=2
		INTEGER(I4B), PARAMETER :: NPAR_CUMSUM=16
		INTEGER(I4B), PARAMETER :: NPAR_CUMPROD=8
		INTEGER(I4B), PARAMETER :: NPAR_POLY=8
		INTEGER(I4B), PARAMETER :: NPAR_POLYTERM=8
		INTERFACE array_copy
	!		MODULE PROCEDURE array_copy_r, array_copy_d, array_copy_i
			MODULE PROCEDURE array_copy_r, array_copy_i
		END INTERFACE
		INTERFACE swap
			MODULE PROCEDURE swap_i,swap_r,swap_rv,swap_c, &
				swap_cv,swap_cm,swap_z,swap_zv,swap_zm, &
				masked_swap_rs,masked_swap_rv,masked_swap_rm
		END INTERFACE
		INTERFACE reallocate
			MODULE PROCEDURE reallocate_rv,reallocate_rm,&
				reallocate_iv,reallocate_im,reallocate_hv
		END INTERFACE
		INTERFACE imaxloc
			MODULE PROCEDURE imaxloc_r,imaxloc_i
		END INTERFACE
		INTERFACE assert
			MODULE PROCEDURE assert1,assert2,assert3,assert4,assert_v
		END INTERFACE
		INTERFACE assert_eq
			MODULE PROCEDURE assert_eq2,assert_eq3,assert_eq4,assert_eqn
		END INTERFACE
		INTERFACE arth
	!		MODULE PROCEDURE arth_r, arth_d, arth_i
			MODULE PROCEDURE arth_r, arth_i
		END INTERFACE
		INTERFACE geop
	!		MODULE PROCEDURE geop_r, geop_d, geop_i, geop_c, geop_dv
			MODULE PROCEDURE geop_r, geop_i, geop_c, geop_dv
		END INTERFACE
		INTERFACE cumsum
			MODULE PROCEDURE cumsum_r,cumsum_i
		END INTERFACE
		INTERFACE poly
	!		MODULE PROCEDURE poly_rr,poly_rrv,poly_dd,poly_ddv,&
	!			poly_rc,poly_cc,poly_msk_rrv,poly_msk_ddv
			MODULE PROCEDURE poly_rr,poly_rrv,&
				poly_rc,poly_cc,poly_msk_rrv
		END INTERFACE
		INTERFACE poly_term
			MODULE PROCEDURE poly_term_rr,poly_term_cc
		END INTERFACE
		INTERFACE outerprod
	!		MODULE PROCEDURE outerprod_r,outerprod_d
			MODULE PROCEDURE outerprod_r
		END INTERFACE
		INTERFACE outerdiff
	!		MODULE PROCEDURE outerdiff_r,outerdiff_d,outerdiff_i
			MODULE PROCEDURE outerdiff_r,outerdiff_i
		END INTERFACE
		INTERFACE scatter_add
	!		MODULE PROCEDURE scatter_add_r,scatter_add_d
			MODULE PROCEDURE scatter_add_r
		END INTERFACE
		INTERFACE scatter_max
	!		MODULE PROCEDURE scatter_max_r,scatter_max_d
			MODULE PROCEDURE scatter_max_r
		END INTERFACE
		INTERFACE diagadd
			MODULE PROCEDURE diagadd_rv,diagadd_r
		END INTERFACE
		INTERFACE diagmult
			MODULE PROCEDURE diagmult_rv,diagmult_r
		END INTERFACE
		INTERFACE get_diag
	!		MODULE PROCEDURE get_diag_rv, get_diag_dv
			MODULE PROCEDURE get_diag_rv
		END INTERFACE
		INTERFACE put_diag
			MODULE PROCEDURE put_diag_rv, put_diag_r
		END INTERFACE
	CONTAINS
	!BL
		SUBROUTINE array_copy_r(src,dest,n_copied,n_not_copied)
		REAL(SP), DIMENSION(:), INTENT(IN) :: src
		REAL(SP), DIMENSION(:), INTENT(OUT) :: dest
		INTEGER(I4B), INTENT(OUT) :: n_copied, n_not_copied
		n_copied=min(size(src),size(dest))
		n_not_copied=size(src)-n_copied
		dest(1:n_copied)=src(1:n_copied)
		END SUBROUTINE array_copy_r
	!BL
		SUBROUTINE array_copy_d(src,dest,n_copied,n_not_copied)
		REAL(DP), DIMENSION(:), INTENT(IN) :: src
		REAL(DP), DIMENSION(:), INTENT(OUT) :: dest
		INTEGER(I4B), INTENT(OUT) :: n_copied, n_not_copied
		n_copied=min(size(src),size(dest))
		n_not_copied=size(src)-n_copied
		dest(1:n_copied)=src(1:n_copied)
		END SUBROUTINE array_copy_d
	!BL
		SUBROUTINE array_copy_i(src,dest,n_copied,n_not_copied)
		INTEGER(I4B), DIMENSION(:), INTENT(IN) :: src
		INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: dest
		INTEGER(I4B), INTENT(OUT) :: n_copied, n_not_copied
		n_copied=min(size(src),size(dest))
		n_not_copied=size(src)-n_copied
		dest(1:n_copied)=src(1:n_copied)
		END SUBROUTINE array_copy_i
	!BL
	!BL
		SUBROUTINE swap_i(a,b)
		INTEGER(I4B), INTENT(INOUT) :: a,b
		INTEGER(I4B) :: dum
		dum=a
		a=b
		b=dum
		END SUBROUTINE swap_i
	!BL
		SUBROUTINE swap_r(a,b)
		REAL(SP), INTENT(INOUT) :: a,b
		REAL(SP) :: dum
		dum=a
		a=b
		b=dum
		END SUBROUTINE swap_r
	!BL
		SUBROUTINE swap_rv(a,b)
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
		REAL(SP), DIMENSION(SIZE(a)) :: dum
		dum=a
		a=b
		b=dum
		END SUBROUTINE swap_rv
	!BL
		SUBROUTINE swap_c(a,b)
		COMPLEX(SPC), INTENT(INOUT) :: a,b
		COMPLEX(SPC) :: dum
		dum=a
		a=b
		b=dum
		END SUBROUTINE swap_c
	!BL
		SUBROUTINE swap_cv(a,b)
		COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: a,b
		COMPLEX(SPC), DIMENSION(SIZE(a)) :: dum
		dum=a
		a=b
		b=dum
		END SUBROUTINE swap_cv
	!BL
		SUBROUTINE swap_cm(a,b)
		COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
		COMPLEX(SPC), DIMENSION(size(a,1),size(a,2)) :: dum
		dum=a
		a=b
		b=dum
		END SUBROUTINE swap_cm
	!BL
		SUBROUTINE swap_z(a,b)
		COMPLEX(DPC), INTENT(INOUT) :: a,b
		COMPLEX(DPC) :: dum
		dum=a
		a=b
		b=dum
		END SUBROUTINE swap_z
	!BL
		SUBROUTINE swap_zv(a,b)
		COMPLEX(DPC), DIMENSION(:), INTENT(INOUT) :: a,b
		COMPLEX(DPC), DIMENSION(SIZE(a)) :: dum
		dum=a
		a=b
		b=dum
		END SUBROUTINE swap_zv
	!BL
		SUBROUTINE swap_zm(a,b)
		COMPLEX(DPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
		COMPLEX(DPC), DIMENSION(size(a,1),size(a,2)) :: dum
		dum=a
		a=b
		b=dum
		END SUBROUTINE swap_zm
	!BL
		SUBROUTINE masked_swap_rs(a,b,mask)
		REAL(SP), INTENT(INOUT) :: a,b
		LOGICAL(LGT), INTENT(IN) :: mask
		REAL(SP) :: swp
		if (mask) then
			swp=a
			a=b
			b=swp
		end if
		END SUBROUTINE masked_swap_rs
	!BL
		SUBROUTINE masked_swap_rv(a,b,mask)
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
		LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
		REAL(SP), DIMENSION(size(a)) :: swp
		where (mask)
			swp=a
			a=b
			b=swp
		end where
		END SUBROUTINE masked_swap_rv
	!BL
		SUBROUTINE masked_swap_rm(a,b,mask)
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a,b
		LOGICAL(LGT), DIMENSION(:,:), INTENT(IN) :: mask
		REAL(SP), DIMENSION(size(a,1),size(a,2)) :: swp
		where (mask)
			swp=a
			a=b
			b=swp
		end where
		END SUBROUTINE masked_swap_rm
	!BL
	!BL
		FUNCTION reallocate_rv(p,n)
		REAL(SP), DIMENSION(:), POINTER :: p, reallocate_rv
		INTEGER(I4B), INTENT(IN) :: n
		INTEGER(I4B) :: nold,ierr
		allocate(reallocate_rv(n),stat=ierr)
		if (ierr /= 0) call &
			nrerror('reallocate_rv: problem in attempt to allocate memory')
		if (.not. associated(p)) RETURN
		nold=size(p)
		reallocate_rv(1:min(nold,n))=p(1:min(nold,n))
		deallocate(p)
		END FUNCTION reallocate_rv
	!BL
		FUNCTION reallocate_iv(p,n)
		INTEGER(I4B), DIMENSION(:), POINTER :: p, reallocate_iv
		INTEGER(I4B), INTENT(IN) :: n
		INTEGER(I4B) :: nold,ierr
		allocate(reallocate_iv(n),stat=ierr)
		if (ierr /= 0) call &
			nrerror('reallocate_iv: problem in attempt to allocate memory')
		if (.not. associated(p)) RETURN
		nold=size(p)
		reallocate_iv(1:min(nold,n))=p(1:min(nold,n))
		deallocate(p)
		END FUNCTION reallocate_iv
	!BL
		FUNCTION reallocate_hv(p,n)
		CHARACTER(1), DIMENSION(:), POINTER :: p, reallocate_hv
		INTEGER(I4B), INTENT(IN) :: n
		INTEGER(I4B) :: nold,ierr
		allocate(reallocate_hv(n),stat=ierr)
		if (ierr /= 0) call &
			nrerror('reallocate_hv: problem in attempt to allocate memory')
		if (.not. associated(p)) RETURN
		nold=size(p)
		reallocate_hv(1:min(nold,n))=p(1:min(nold,n))
		deallocate(p)
		END FUNCTION reallocate_hv
	!BL
		FUNCTION reallocate_rm(p,n,m)
		REAL(SP), DIMENSION(:,:), POINTER :: p, reallocate_rm
		INTEGER(I4B), INTENT(IN) :: n,m
		INTEGER(I4B) :: nold,mold,ierr
		allocate(reallocate_rm(n,m),stat=ierr)
		if (ierr /= 0) call &
			nrerror('reallocate_rm: problem in attempt to allocate memory')
		if (.not. associated(p)) RETURN
		nold=size(p,1)
		mold=size(p,2)
		reallocate_rm(1:min(nold,n),1:min(mold,m))=&
			p(1:min(nold,n),1:min(mold,m))
		deallocate(p)
		END FUNCTION reallocate_rm
	!BL
		FUNCTION reallocate_im(p,n,m)
		INTEGER(I4B), DIMENSION(:,:), POINTER :: p, reallocate_im
		INTEGER(I4B), INTENT(IN) :: n,m
		INTEGER(I4B) :: nold,mold,ierr
		allocate(reallocate_im(n,m),stat=ierr)
		if (ierr /= 0) call &
			nrerror('reallocate_im: problem in attempt to allocate memory')
		if (.not. associated(p)) RETURN
		nold=size(p,1)
		mold=size(p,2)
		reallocate_im(1:min(nold,n),1:min(mold,m))=&
			p(1:min(nold,n),1:min(mold,m))
		deallocate(p)
		END FUNCTION reallocate_im
	!BL
		FUNCTION ifirstloc(mask)
		LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
		INTEGER(I4B) :: ifirstloc
		INTEGER(I4B), DIMENSION(1) :: loc
		loc=maxloc(merge(1,0,mask))
		ifirstloc=loc(1)
		if (.not. mask(ifirstloc)) ifirstloc=size(mask)+1
		END FUNCTION ifirstloc
	!BL
		FUNCTION imaxloc_r(arr)
		REAL(SP), DIMENSION(:), INTENT(IN) :: arr
		INTEGER(I4B) :: imaxloc_r
		INTEGER(I4B), DIMENSION(1) :: imax
		imax=maxloc(arr(:))
		imaxloc_r=imax(1)
		END FUNCTION imaxloc_r
	!BL
		FUNCTION imaxloc_i(iarr)
		INTEGER(I4B), DIMENSION(:), INTENT(IN) :: iarr
		INTEGER(I4B), DIMENSION(1) :: imax
		INTEGER(I4B) :: imaxloc_i
		imax=maxloc(iarr(:))
		imaxloc_i=imax(1)
		END FUNCTION imaxloc_i
	!BL
		FUNCTION iminloc(arr)
		REAL(SP), DIMENSION(:), INTENT(IN) :: arr
		INTEGER(I4B), DIMENSION(1) :: imin
		INTEGER(I4B) :: iminloc
		imin=minloc(arr(:))
		iminloc=imin(1)
		END FUNCTION iminloc
	!BL
		SUBROUTINE assert1(n1,string)
		CHARACTER(LEN=*), INTENT(IN) :: string
		LOGICAL, INTENT(IN) :: n1
		if (.not. n1) then
			write (*,*) 'nrerror: an assertion failed with this tag:', &
				string
			STOP 'program terminated by assert1'
		end if
		END SUBROUTINE assert1
	!BL
		SUBROUTINE assert2(n1,n2,string)
		CHARACTER(LEN=*), INTENT(IN) :: string
		LOGICAL, INTENT(IN) :: n1,n2
		if (.not. (n1 .and. n2)) then
			write (*,*) 'nrerror: an assertion failed with this tag:', &
				string
			STOP 'program terminated by assert2'
		end if
		END SUBROUTINE assert2
	!BL
		SUBROUTINE assert3(n1,n2,n3,string)
		CHARACTER(LEN=*), INTENT(IN) :: string
		LOGICAL, INTENT(IN) :: n1,n2,n3
		if (.not. (n1 .and. n2 .and. n3)) then
			write (*,*) 'nrerror: an assertion failed with this tag:', &
				string
			STOP 'program terminated by assert3'
		end if
		END SUBROUTINE assert3
	!BL
		SUBROUTINE assert4(n1,n2,n3,n4,string)
		CHARACTER(LEN=*), INTENT(IN) :: string
		LOGICAL, INTENT(IN) :: n1,n2,n3,n4
		if (.not. (n1 .and. n2 .and. n3 .and. n4)) then
			write (*,*) 'nrerror: an assertion failed with this tag:', &
				string
			STOP 'program terminated by assert4'
		end if
		END SUBROUTINE assert4
	!BL
		SUBROUTINE assert_v(n,string)
		CHARACTER(LEN=*), INTENT(IN) :: string
		LOGICAL, DIMENSION(:), INTENT(IN) :: n
		if (.not. all(n)) then
			write (*,*) 'nrerror: an assertion failed with this tag:', &
				string
			STOP 'program terminated by assert_v'
		end if
		END SUBROUTINE assert_v
	!BL
		FUNCTION assert_eq2(n1,n2,string)
		CHARACTER(LEN=*), INTENT(IN) :: string
		INTEGER, INTENT(IN) :: n1,n2
		INTEGER :: assert_eq2
		if (n1 == n2) then
			assert_eq2=n1
		else
			write (*,*) 'nrerror: an assert_eq failed with this tag:', &
				string
			STOP 'program terminated by assert_eq2'
		end if
		END FUNCTION assert_eq2
	!BL
		FUNCTION assert_eq3(n1,n2,n3,string)
		CHARACTER(LEN=*), INTENT(IN) :: string
		INTEGER, INTENT(IN) :: n1,n2,n3
		INTEGER :: assert_eq3
		if (n1 == n2 .and. n2 == n3) then
			assert_eq3=n1
		else
			write (*,*) 'nrerror: an assert_eq failed with this tag:', &
				string
			STOP 'program terminated by assert_eq3'
		end if
		END FUNCTION assert_eq3
	!BL
		FUNCTION assert_eq4(n1,n2,n3,n4,string)
		CHARACTER(LEN=*), INTENT(IN) :: string
		INTEGER, INTENT(IN) :: n1,n2,n3,n4
		INTEGER :: assert_eq4
		if (n1 == n2 .and. n2 == n3 .and. n3 == n4) then
			assert_eq4=n1
		else
			write (*,*) 'nrerror: an assert_eq failed with this tag:', &
				string
			STOP 'program terminated by assert_eq4'
		end if
		END FUNCTION assert_eq4
	!BL
		FUNCTION assert_eqn(nn,string)
		CHARACTER(LEN=*), INTENT(IN) :: string
		INTEGER, DIMENSION(:), INTENT(IN) :: nn
		INTEGER :: assert_eqn
		if (all(nn(2:) == nn(1))) then
			assert_eqn=nn(1)
		else
			write (*,*) 'nrerror: an assert_eq failed with this tag:', &
				string
			STOP 'program terminated by assert_eqn'
		end if
		END FUNCTION assert_eqn
	!BL
		SUBROUTINE nrerror(string)
		CHARACTER(LEN=*), INTENT(IN) :: string
		write (*,*) 'nrerror: ',string
		STOP 'program terminated by nrerror'
		END SUBROUTINE nrerror
	!BL
		FUNCTION arth_r(first,increment,n)
		REAL(SP), INTENT(IN) :: first,increment
		INTEGER(I4B), INTENT(IN) :: n
		REAL(SP), DIMENSION(n) :: arth_r
		INTEGER(I4B) :: k,k2
		REAL(SP) :: temp
		if (n > 0) arth_r(1)=first
		if (n <= NPAR_ARTH) then
			do k=2,n
				arth_r(k)=arth_r(k-1)+increment
			end do
		else
			do k=2,NPAR2_ARTH
				arth_r(k)=arth_r(k-1)+increment
			end do
			temp=increment*NPAR2_ARTH
			k=NPAR2_ARTH
			do
				if (k >= n) exit
				k2=k+k
				arth_r(k+1:min(k2,n))=temp+arth_r(1:min(k,n-k))
				temp=temp+temp
				k=k2
			end do
		end if
		END FUNCTION arth_r
	!BL
		FUNCTION arth_d(first,increment,n)
		REAL(DP), INTENT(IN) :: first,increment
		INTEGER(I4B), INTENT(IN) :: n
		REAL(DP), DIMENSION(n) :: arth_d
		INTEGER(I4B) :: k,k2
		REAL(DP) :: temp
		if (n > 0) arth_d(1)=first
		if (n <= NPAR_ARTH) then
			do k=2,n
				arth_d(k)=arth_d(k-1)+increment
			end do
		else
			do k=2,NPAR2_ARTH
				arth_d(k)=arth_d(k-1)+increment
			end do
			temp=increment*NPAR2_ARTH
			k=NPAR2_ARTH
			do
				if (k >= n) exit
				k2=k+k
				arth_d(k+1:min(k2,n))=temp+arth_d(1:min(k,n-k))
				temp=temp+temp
				k=k2
			end do
		end if
		END FUNCTION arth_d
	!BL
		FUNCTION arth_i(first,increment,n)
		INTEGER(I4B), INTENT(IN) :: first,increment,n
		INTEGER(I4B), DIMENSION(n) :: arth_i
		INTEGER(I4B) :: k,k2,temp
		if (n > 0) arth_i(1)=first
		if (n <= NPAR_ARTH) then
			do k=2,n
				arth_i(k)=arth_i(k-1)+increment
			end do
		else
			do k=2,NPAR2_ARTH
				arth_i(k)=arth_i(k-1)+increment
			end do
			temp=increment*NPAR2_ARTH
			k=NPAR2_ARTH
			do
				if (k >= n) exit
				k2=k+k
				arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
				temp=temp+temp
				k=k2
			end do
		end if
		END FUNCTION arth_i
	!BL
	!BL
		FUNCTION geop_r(first,factor,n)
		REAL(SP), INTENT(IN) :: first,factor
		INTEGER(I4B), INTENT(IN) :: n
		REAL(SP), DIMENSION(n) :: geop_r
		INTEGER(I4B) :: k,k2
		REAL(SP) :: temp
		if (n > 0) geop_r(1)=first
		if (n <= NPAR_GEOP) then
			do k=2,n
				geop_r(k)=geop_r(k-1)*factor
			end do
		else
			do k=2,NPAR2_GEOP
				geop_r(k)=geop_r(k-1)*factor
			end do
			temp=factor**NPAR2_GEOP
			k=NPAR2_GEOP
			do
				if (k >= n) exit
				k2=k+k
				geop_r(k+1:min(k2,n))=temp*geop_r(1:min(k,n-k))
				temp=temp*temp
				k=k2
			end do
		end if
		END FUNCTION geop_r
	!BL
		FUNCTION geop_d(first,factor,n)
		REAL(DP), INTENT(IN) :: first,factor
		INTEGER(I4B), INTENT(IN) :: n
		REAL(DP), DIMENSION(n) :: geop_d
		INTEGER(I4B) :: k,k2
		REAL(DP) :: temp
		if (n > 0) geop_d(1)=first
		if (n <= NPAR_GEOP) then
			do k=2,n
				geop_d(k)=geop_d(k-1)*factor
			end do
		else
			do k=2,NPAR2_GEOP
				geop_d(k)=geop_d(k-1)*factor
			end do
			temp=factor**NPAR2_GEOP
			k=NPAR2_GEOP
			do
				if (k >= n) exit
				k2=k+k
				geop_d(k+1:min(k2,n))=temp*geop_d(1:min(k,n-k))
				temp=temp*temp
				k=k2
			end do
		end if
		END FUNCTION geop_d
	!BL
		FUNCTION geop_i(first,factor,n)
		INTEGER(I4B), INTENT(IN) :: first,factor,n
		INTEGER(I4B), DIMENSION(n) :: geop_i
		INTEGER(I4B) :: k,k2,temp
		if (n > 0) geop_i(1)=first
		if (n <= NPAR_GEOP) then
			do k=2,n
				geop_i(k)=geop_i(k-1)*factor
			end do
		else
			do k=2,NPAR2_GEOP
				geop_i(k)=geop_i(k-1)*factor
			end do
			temp=factor**NPAR2_GEOP
			k=NPAR2_GEOP
			do
				if (k >= n) exit
				k2=k+k
				geop_i(k+1:min(k2,n))=temp*geop_i(1:min(k,n-k))
				temp=temp*temp
				k=k2
			end do
		end if
		END FUNCTION geop_i
	!BL
		FUNCTION geop_c(first,factor,n)
		COMPLEX(SP), INTENT(IN) :: first,factor
		INTEGER(I4B), INTENT(IN) :: n
		COMPLEX(SP), DIMENSION(n) :: geop_c
		INTEGER(I4B) :: k,k2
		COMPLEX(SP) :: temp
		if (n > 0) geop_c(1)=first
		if (n <= NPAR_GEOP) then
			do k=2,n
				geop_c(k)=geop_c(k-1)*factor
			end do
		else
			do k=2,NPAR2_GEOP
				geop_c(k)=geop_c(k-1)*factor
			end do
			temp=factor**NPAR2_GEOP
			k=NPAR2_GEOP
			do
				if (k >= n) exit
				k2=k+k
				geop_c(k+1:min(k2,n))=temp*geop_c(1:min(k,n-k))
				temp=temp*temp
				k=k2
			end do
		end if
		END FUNCTION geop_c
	!BL
		FUNCTION geop_dv(first,factor,n)
		REAL(DP), DIMENSION(:), INTENT(IN) :: first,factor
		INTEGER(I4B), INTENT(IN) :: n
		REAL(DP), DIMENSION(size(first),n) :: geop_dv
		INTEGER(I4B) :: k,k2
		REAL(DP), DIMENSION(size(first)) :: temp
		if (n > 0) geop_dv(:,1)=first(:)
		if (n <= NPAR_GEOP) then
			do k=2,n
				geop_dv(:,k)=geop_dv(:,k-1)*factor(:)
			end do
		else
			do k=2,NPAR2_GEOP
				geop_dv(:,k)=geop_dv(:,k-1)*factor(:)
			end do
			temp=factor**NPAR2_GEOP
			k=NPAR2_GEOP
			do
				if (k >= n) exit
				k2=k+k
				geop_dv(:,k+1:min(k2,n))=geop_dv(:,1:min(k,n-k))*&
					spread(temp,2,size(geop_dv(:,1:min(k,n-k)),2))
				temp=temp*temp
				k=k2
			end do
		end if
		END FUNCTION geop_dv
	!BL
	!BL
		RECURSIVE FUNCTION cumsum_r(arr,seed) RESULT(ans)
		REAL(SP), DIMENSION(:), INTENT(IN) :: arr
		REAL(SP), OPTIONAL, INTENT(IN) :: seed
		REAL(SP), DIMENSION(size(arr)) :: ans
		INTEGER(I4B) :: n,j
		REAL(SP) :: sd
		n=size(arr)
		if (n == 0_i4b) RETURN
		sd=0.0_sp
		if (present(seed)) sd=seed
		ans(1)=arr(1)+sd
		if (n < NPAR_CUMSUM) then
			do j=2,n
				ans(j)=ans(j-1)+arr(j)
			end do
		else
			ans(2:n:2)=cumsum_r(arr(2:n:2)+arr(1:n-1:2),sd)
			ans(3:n:2)=ans(2:n-1:2)+arr(3:n:2)
		end if
		END FUNCTION cumsum_r
	!BL
		RECURSIVE FUNCTION cumsum_i(arr,seed) RESULT(ans)
		INTEGER(I4B), DIMENSION(:), INTENT(IN) :: arr
		INTEGER(I4B), OPTIONAL, INTENT(IN) :: seed
		INTEGER(I4B), DIMENSION(size(arr)) :: ans
		INTEGER(I4B) :: n,j,sd
		n=size(arr)
		if (n == 0_i4b) RETURN
		sd=0_i4b
		if (present(seed)) sd=seed
		ans(1)=arr(1)+sd
		if (n < NPAR_CUMSUM) then
			do j=2,n
				ans(j)=ans(j-1)+arr(j)
			end do
		else
			ans(2:n:2)=cumsum_i(arr(2:n:2)+arr(1:n-1:2),sd)
			ans(3:n:2)=ans(2:n-1:2)+arr(3:n:2)
		end if
		END FUNCTION cumsum_i
	!BL
	!BL
		RECURSIVE FUNCTION cumprod(arr,seed) RESULT(ans)
		REAL(SP), DIMENSION(:), INTENT(IN) :: arr
		REAL(SP), OPTIONAL, INTENT(IN) :: seed
		REAL(SP), DIMENSION(size(arr)) :: ans
		INTEGER(I4B) :: n,j
		REAL(SP) :: sd
		n=size(arr)
		if (n == 0_i4b) RETURN
		sd=1.0_sp
		if (present(seed)) sd=seed
		ans(1)=arr(1)*sd
		if (n < NPAR_CUMPROD) then
			do j=2,n
				ans(j)=ans(j-1)*arr(j)
			end do
		else
			ans(2:n:2)=cumprod(arr(2:n:2)*arr(1:n-1:2),sd)
			ans(3:n:2)=ans(2:n-1:2)*arr(3:n:2)
		end if
		END FUNCTION cumprod
	!BL
	!BL
		FUNCTION poly_rr(x,coeffs)
		REAL(SP), INTENT(IN) :: x
		REAL(SP), DIMENSION(:), INTENT(IN) :: coeffs
		REAL(SP) :: poly_rr
		REAL(SP) :: pow
		REAL(SP), DIMENSION(:), ALLOCATABLE :: vec
		INTEGER(I4B) :: i,n,nn
		n=size(coeffs)
		if (n <= 0) then
			poly_rr=0.0_sp
		else if (n < NPAR_POLY) then
			poly_rr=coeffs(n)
			do i=n-1,1,-1
				poly_rr=x*poly_rr+coeffs(i)
			end do
		else
			allocate(vec(n+1))
			pow=x
			vec(1:n)=coeffs
			do
				vec(n+1)=0.0_sp
				nn=ishft(n+1,-1)
				vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
				if (nn == 1) exit
				pow=pow*pow
				n=nn
			end do
			poly_rr=vec(1)
			deallocate(vec)
		end if
		END FUNCTION poly_rr
	!BL
		FUNCTION poly_dd(x,coeffs)
		REAL(DP), INTENT(IN) :: x
		REAL(DP), DIMENSION(:), INTENT(IN) :: coeffs
		REAL(DP) :: poly_dd
		REAL(DP) :: pow
		REAL(DP), DIMENSION(:), ALLOCATABLE :: vec
		INTEGER(I4B) :: i,n,nn
		n=size(coeffs)
		if (n <= 0) then
			poly_dd=0.0_dp
		else if (n < NPAR_POLY) then
			poly_dd=coeffs(n)
			do i=n-1,1,-1
				poly_dd=x*poly_dd+coeffs(i)
			end do
		else
			allocate(vec(n+1))
			pow=x
			vec(1:n)=coeffs
			do
				vec(n+1)=0.0_dp
				nn=ishft(n+1,-1)
				vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
				if (nn == 1) exit
				pow=pow*pow
				n=nn
			end do
			poly_dd=vec(1)
			deallocate(vec)
		end if
		END FUNCTION poly_dd
	!BL
		FUNCTION poly_rc(x,coeffs)
		COMPLEX(SPC), INTENT(IN) :: x
		REAL(SP), DIMENSION(:), INTENT(IN) :: coeffs
		COMPLEX(SPC) :: poly_rc
		COMPLEX(SPC) :: pow
		COMPLEX(SPC), DIMENSION(:), ALLOCATABLE :: vec
		INTEGER(I4B) :: i,n,nn
		n=size(coeffs)
		if (n <= 0) then
			poly_rc=0.0_sp
		else if (n < NPAR_POLY) then
			poly_rc=coeffs(n)
			do i=n-1,1,-1
				poly_rc=x*poly_rc+coeffs(i)
			end do
		else
			allocate(vec(n+1))
			pow=x
			vec(1:n)=coeffs
			do
				vec(n+1)=0.0_sp
				nn=ishft(n+1,-1)
				vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
				if (nn == 1) exit
				pow=pow*pow
				n=nn
			end do
			poly_rc=vec(1)
			deallocate(vec)
		end if
		END FUNCTION poly_rc
	!BL
		FUNCTION poly_cc(x,coeffs)
		COMPLEX(SPC), INTENT(IN) :: x
		COMPLEX(SPC), DIMENSION(:), INTENT(IN) :: coeffs
		COMPLEX(SPC) :: poly_cc
		COMPLEX(SPC) :: pow
		COMPLEX(SPC), DIMENSION(:), ALLOCATABLE :: vec
		INTEGER(I4B) :: i,n,nn
		n=size(coeffs)
		if (n <= 0) then
			poly_cc=0.0_sp
		else if (n < NPAR_POLY) then
			poly_cc=coeffs(n)
			do i=n-1,1,-1
				poly_cc=x*poly_cc+coeffs(i)
			end do
		else
			allocate(vec(n+1))
			pow=x
			vec(1:n)=coeffs
			do
				vec(n+1)=0.0_sp
				nn=ishft(n+1,-1)
				vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
				if (nn == 1) exit
				pow=pow*pow
				n=nn
			end do
			poly_cc=vec(1)
			deallocate(vec)
		end if
		END FUNCTION poly_cc
	!BL
		FUNCTION poly_rrv(x,coeffs)
		REAL(SP), DIMENSION(:), INTENT(IN) :: coeffs,x
		REAL(SP), DIMENSION(size(x)) :: poly_rrv
		INTEGER(I4B) :: i,n,m
		m=size(coeffs)
		n=size(x)
		if (m <= 0) then
			poly_rrv=0.0_sp
		else if (m < n .or. m < NPAR_POLY) then
			poly_rrv=coeffs(m)
			do i=m-1,1,-1
				poly_rrv=x*poly_rrv+coeffs(i)
			end do
		else
			do i=1,n
				poly_rrv(i)=poly_rr(x(i),coeffs)
			end do
		end if
		END FUNCTION poly_rrv
	!BL
		FUNCTION poly_ddv(x,coeffs)
		REAL(DP), DIMENSION(:), INTENT(IN) :: coeffs,x
		REAL(DP), DIMENSION(size(x)) :: poly_ddv
		INTEGER(I4B) :: i,n,m
		m=size(coeffs)
		n=size(x)
		if (m <= 0) then
			poly_ddv=0.0_dp
		else if (m < n .or. m < NPAR_POLY) then
			poly_ddv=coeffs(m)
			do i=m-1,1,-1
				poly_ddv=x*poly_ddv+coeffs(i)
			end do
		else
			do i=1,n
				poly_ddv(i)=poly_dd(x(i),coeffs)
			end do
		end if
		END FUNCTION poly_ddv
	!BL
		FUNCTION poly_msk_rrv(x,coeffs,mask)
		REAL(SP), DIMENSION(:), INTENT(IN) :: coeffs,x
		LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
		REAL(SP), DIMENSION(size(x)) :: poly_msk_rrv
		poly_msk_rrv=unpack(poly_rrv(pack(x,mask),coeffs),mask,0.0_sp)
		END FUNCTION poly_msk_rrv
	!BL
		FUNCTION poly_msk_ddv(x,coeffs,mask)
		REAL(DP), DIMENSION(:), INTENT(IN) :: coeffs,x
		LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
		REAL(DP), DIMENSION(size(x)) :: poly_msk_ddv
		poly_msk_ddv=unpack(poly_ddv(pack(x,mask),coeffs),mask,0.0_dp)
		END FUNCTION poly_msk_ddv
	!BL
	!BL
		RECURSIVE FUNCTION poly_term_rr(a,b) RESULT(u)
		REAL(SP), DIMENSION(:), INTENT(IN) :: a
		REAL(SP), INTENT(IN) :: b
		REAL(SP), DIMENSION(size(a)) :: u
		INTEGER(I4B) :: n,j
		n=size(a)
		if (n <= 0) RETURN
		u(1)=a(1)
		if (n < NPAR_POLYTERM) then
			do j=2,n
				u(j)=a(j)+b*u(j-1)
			end do
		else
			u(2:n:2)=poly_term_rr(a(2:n:2)+a(1:n-1:2)*b,b*b)
			u(3:n:2)=a(3:n:2)+b*u(2:n-1:2)
		end if
		END FUNCTION poly_term_rr
	!BL
		RECURSIVE FUNCTION poly_term_cc(a,b) RESULT(u)
		COMPLEX(SPC), DIMENSION(:), INTENT(IN) :: a
		COMPLEX(SPC), INTENT(IN) :: b
		COMPLEX(SPC), DIMENSION(size(a)) :: u
		INTEGER(I4B) :: n,j
		n=size(a)
		if (n <= 0) RETURN
		u(1)=a(1)
		if (n < NPAR_POLYTERM) then
			do j=2,n
				u(j)=a(j)+b*u(j-1)
			end do
		else
			u(2:n:2)=poly_term_cc(a(2:n:2)+a(1:n-1:2)*b,b*b)
			u(3:n:2)=a(3:n:2)+b*u(2:n-1:2)
		end if
		END FUNCTION poly_term_cc
	!BL
	!BL
		FUNCTION zroots_unity(n,nn)
		INTEGER(I4B), INTENT(IN) :: n,nn
		COMPLEX(SPC), DIMENSION(nn) :: zroots_unity
		INTEGER(I4B) :: k
		REAL(SP) :: theta
		zroots_unity(1)=1.0
		theta=TWOPI/n
		k=1
		do
			if (k >= nn) exit
			zroots_unity(k+1)=cmplx(cos(k*theta),sin(k*theta),SPC)
			zroots_unity(k+2:min(2*k,nn))=zroots_unity(k+1)*&
				zroots_unity(2:min(k,nn-k))
			k=2*k
		end do
		END FUNCTION zroots_unity
	!BL
		FUNCTION outerprod_r(a,b)
		REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
		REAL(SP), DIMENSION(size(a),size(b)) :: outerprod_r
		outerprod_r = spread(a,dim=2,ncopies=size(b)) * &
			spread(b,dim=1,ncopies=size(a))
		END FUNCTION outerprod_r
	!BL
		FUNCTION outerprod_d(a,b)
		REAL(DP), DIMENSION(:), INTENT(IN) :: a,b
		REAL(DP), DIMENSION(size(a),size(b)) :: outerprod_d
		outerprod_d = spread(a,dim=2,ncopies=size(b)) * &
			spread(b,dim=1,ncopies=size(a))
		END FUNCTION outerprod_d
	!BL
		FUNCTION outerdiv(a,b)
		REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
		REAL(SP), DIMENSION(size(a),size(b)) :: outerdiv
		outerdiv = spread(a,dim=2,ncopies=size(b)) / &
			spread(b,dim=1,ncopies=size(a))
		END FUNCTION outerdiv
	!BL
		FUNCTION outersum(a,b)
		REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
		REAL(SP), DIMENSION(size(a),size(b)) :: outersum
		outersum = spread(a,dim=2,ncopies=size(b)) + &
			spread(b,dim=1,ncopies=size(a))
		END FUNCTION outersum
	!BL
		FUNCTION outerdiff_r(a,b)
		REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
		REAL(SP), DIMENSION(size(a),size(b)) :: outerdiff_r
		outerdiff_r = spread(a,dim=2,ncopies=size(b)) - &
			spread(b,dim=1,ncopies=size(a))
		END FUNCTION outerdiff_r
	!BL
		FUNCTION outerdiff_d(a,b)
		REAL(DP), DIMENSION(:), INTENT(IN) :: a,b
		REAL(DP), DIMENSION(size(a),size(b)) :: outerdiff_d
		outerdiff_d = spread(a,dim=2,ncopies=size(b)) - &
			spread(b,dim=1,ncopies=size(a))
		END FUNCTION outerdiff_d
	!BL
		FUNCTION outerdiff_i(a,b)
		INTEGER(I4B), DIMENSION(:), INTENT(IN) :: a,b
		INTEGER(I4B), DIMENSION(size(a),size(b)) :: outerdiff_i
		outerdiff_i = spread(a,dim=2,ncopies=size(b)) - &
			spread(b,dim=1,ncopies=size(a))
		END FUNCTION outerdiff_i
	!BL
		FUNCTION outerand(a,b)
		LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: a,b
		LOGICAL(LGT), DIMENSION(size(a),size(b)) :: outerand
		outerand = spread(a,dim=2,ncopies=size(b)) .and. &
			spread(b,dim=1,ncopies=size(a))
		END FUNCTION outerand
	!BL
		SUBROUTINE scatter_add_r(dest,source,dest_index)
		REAL(SP), DIMENSION(:), INTENT(OUT) :: dest
		REAL(SP), DIMENSION(:), INTENT(IN) :: source
		INTEGER(I4B), DIMENSION(:), INTENT(IN) :: dest_index
		INTEGER(I4B) :: m,n,j,i
		n=assert_eq2(size(source),size(dest_index),'scatter_add_r')
		m=size(dest)
		do j=1,n
			i=dest_index(j)
			if (i > 0 .and. i <= m) dest(i)=dest(i)+source(j)
		end do
		END SUBROUTINE scatter_add_r
		SUBROUTINE scatter_add_d(dest,source,dest_index)
		REAL(DP), DIMENSION(:), INTENT(OUT) :: dest
		REAL(DP), DIMENSION(:), INTENT(IN) :: source
		INTEGER(I4B), DIMENSION(:), INTENT(IN) :: dest_index
		INTEGER(I4B) :: m,n,j,i
		n=assert_eq2(size(source),size(dest_index),'scatter_add_d')
		m=size(dest)
		do j=1,n
			i=dest_index(j)
			if (i > 0 .and. i <= m) dest(i)=dest(i)+source(j)
		end do
		END SUBROUTINE scatter_add_d
		SUBROUTINE scatter_max_r(dest,source,dest_index)
		REAL(SP), DIMENSION(:), INTENT(OUT) :: dest
		REAL(SP), DIMENSION(:), INTENT(IN) :: source
		INTEGER(I4B), DIMENSION(:), INTENT(IN) :: dest_index
		INTEGER(I4B) :: m,n,j,i
		n=assert_eq2(size(source),size(dest_index),'scatter_max_r')
		m=size(dest)
		do j=1,n
			i=dest_index(j)
			if (i > 0 .and. i <= m) dest(i)=max(dest(i),source(j))
		end do
		END SUBROUTINE scatter_max_r
		SUBROUTINE scatter_max_d(dest,source,dest_index)
		REAL(DP), DIMENSION(:), INTENT(OUT) :: dest
		REAL(DP), DIMENSION(:), INTENT(IN) :: source
		INTEGER(I4B), DIMENSION(:), INTENT(IN) :: dest_index
		INTEGER(I4B) :: m,n,j,i
		n=assert_eq2(size(source),size(dest_index),'scatter_max_d')
		m=size(dest)
		do j=1,n
			i=dest_index(j)
			if (i > 0 .and. i <= m) dest(i)=max(dest(i),source(j))
		end do
		END SUBROUTINE scatter_max_d
	!BL
		SUBROUTINE diagadd_rv(mat,diag)
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
		REAL(SP), DIMENSION(:), INTENT(IN) :: diag
		INTEGER(I4B) :: j,n
		n = assert_eq2(size(diag),min(size(mat,1),size(mat,2)),'diagadd_rv')
		do j=1,n
			mat(j,j)=mat(j,j)+diag(j)
		end do
		END SUBROUTINE diagadd_rv
	!BL
		SUBROUTINE diagadd_r(mat,diag)
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
		REAL(SP), INTENT(IN) :: diag
		INTEGER(I4B) :: j,n
		n = min(size(mat,1),size(mat,2))
		do j=1,n
			mat(j,j)=mat(j,j)+diag
		end do
		END SUBROUTINE diagadd_r
	!BL
		SUBROUTINE diagmult_rv(mat,diag)
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
		REAL(SP), DIMENSION(:), INTENT(IN) :: diag
		INTEGER(I4B) :: j,n
		n = assert_eq2(size(diag),min(size(mat,1),size(mat,2)),'diagmult_rv')
		do j=1,n
			mat(j,j)=mat(j,j)*diag(j)
		end do
		END SUBROUTINE diagmult_rv
	!BL
		SUBROUTINE diagmult_r(mat,diag)
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
		REAL(SP), INTENT(IN) :: diag
		INTEGER(I4B) :: j,n
		n = min(size(mat,1),size(mat,2))
		do j=1,n
			mat(j,j)=mat(j,j)*diag
		end do
		END SUBROUTINE diagmult_r
	!BL
		FUNCTION get_diag_rv(mat)
		REAL(SP), DIMENSION(:,:), INTENT(IN) :: mat
		REAL(SP), DIMENSION(size(mat,1)) :: get_diag_rv
		INTEGER(I4B) :: j
		j=assert_eq2(size(mat,1),size(mat,2),'get_diag_rv')
		do j=1,size(mat,1)
			get_diag_rv(j)=mat(j,j)
		end do
		END FUNCTION get_diag_rv
	!BL
		FUNCTION get_diag_dv(mat)
		REAL(DP), DIMENSION(:,:), INTENT(IN) :: mat
		REAL(DP), DIMENSION(size(mat,1)) :: get_diag_dv
		INTEGER(I4B) :: j
		j=assert_eq2(size(mat,1),size(mat,2),'get_diag_dv')
		do j=1,size(mat,1)
			get_diag_dv(j)=mat(j,j)
		end do
		END FUNCTION get_diag_dv
	!BL
		SUBROUTINE put_diag_rv(diagv,mat)
		REAL(SP), DIMENSION(:), INTENT(IN) :: diagv
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
		INTEGER(I4B) :: j,n
		n=assert_eq2(size(diagv),min(size(mat,1),size(mat,2)),'put_diag_rv')
		do j=1,n
			mat(j,j)=diagv(j)
		end do
		END SUBROUTINE put_diag_rv
	!BL
		SUBROUTINE put_diag_r(scal,mat)
		REAL(SP), INTENT(IN) :: scal
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
		INTEGER(I4B) :: j,n
		n = min(size(mat,1),size(mat,2))
		do j=1,n
			mat(j,j)=scal
		end do
		END SUBROUTINE put_diag_r
	!BL
		SUBROUTINE unit_matrix(mat)
		REAL(SP), DIMENSION(:,:), INTENT(OUT) :: mat
		INTEGER(I4B) :: i,n
		n=min(size(mat,1),size(mat,2))
		mat(:,:)=0.0_sp
		do i=1,n
			mat(i,i)=1.0_sp
		end do
		END SUBROUTINE unit_matrix
	!BL
		FUNCTION upper_triangle(j,k,extra)
		INTEGER(I4B), INTENT(IN) :: j,k
		INTEGER(I4B), OPTIONAL, INTENT(IN) :: extra
		LOGICAL(LGT), DIMENSION(j,k) :: upper_triangle
		INTEGER(I4B) :: n
		n=0
		if (present(extra)) n=extra
		upper_triangle=(outerdiff(arth_i(1,1,j),arth_i(1,1,k)) < n)
		END FUNCTION upper_triangle
	!BL
		FUNCTION lower_triangle(j,k,extra)
		INTEGER(I4B), INTENT(IN) :: j,k
		INTEGER(I4B), OPTIONAL, INTENT(IN) :: extra
		LOGICAL(LGT), DIMENSION(j,k) :: lower_triangle
		INTEGER(I4B) :: n
		n=0
		if (present(extra)) n=extra
		lower_triangle=(outerdiff(arth_i(1,1,j),arth_i(1,1,k)) > -n)
		END FUNCTION lower_triangle
	!BL
		FUNCTION vabs(v)
		REAL(SP), DIMENSION(:), INTENT(IN) :: v
		REAL(SP) :: vabs
		vabs=sqrt(dot_product(v,v))
		END FUNCTION vabs
	!BL
	END MODULE nrutil1


	! module for nr
	MODULE nr1
	!       PRIVATE
	!        PUBLIC

			CONTAINS

	! all remaining code for numerical recipes
		SUBROUTINE trapzd(func,a,b,s,n)
		USE nrtype1; USE nrutil1, ONLY : arth
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: a,b
		REAL(SP), INTENT(INOUT) :: s
		INTEGER(I4B), INTENT(IN) :: n
		INTERFACE
			FUNCTION func(x)
			USE nrtype1
			REAL(SP), DIMENSION(:), INTENT(IN) :: x
			REAL(SP), DIMENSION(size(x)) :: func
			END FUNCTION func
		END INTERFACE
		REAL(SP) :: del,fsum
		INTEGER(I4B) :: it
		if (n == 1) then
			s=0.5_sp*(b-a)*sum(func( (/ a,b /) ))
		else
			it=2**(n-2)
			del=(b-a)/it
			fsum=sum(func(arth(a+0.5_sp*del,del,it)))
			s=0.5_sp*(s+del*fsum)
		end if
		END SUBROUTINE trapzd
		FUNCTION zbrent(func,x1,x2,tol)
		USE nrtype1; USE nrutil1, ONLY : nrerror
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x1,x2,tol
		REAL(SP) :: zbrent
		INTERFACE
			FUNCTION func(x)
			USE nrtype1
			IMPLICIT NONE
			REAL(SP), INTENT(IN) :: x
			REAL(SP) :: func
			END FUNCTION func
		END INTERFACE
		INTEGER(I4B), PARAMETER :: ITMAX=100
		REAL(SP), PARAMETER :: EPS=epsilon(x1)
		INTEGER(I4B) :: iter
		REAL(SP) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
		a=x1
		b=x2
		fa=func(a)
		fb=func(b)
		if ((fa > 0.0 .and. fb > 0.0) .or. (fa < 0.0 .and. fb < 0.0)) &
			call nrerror('root must be bracketed for zbrent')
		c=b
		fc=fb
		do iter=1,ITMAX
			if ((fb > 0.0 .and. fc > 0.0) .or. (fb < 0.0 .and. fc < 0.0)) then
				c=a
				fc=fa
				d=b-a
				e=d
			end if
			if (abs(fc) < abs(fb)) then
				a=b
				b=c
				c=a
				fa=fb
				fb=fc
				fc=fa
			end if
			tol1=2.0_sp*EPS*abs(b)+0.5_sp*tol
			xm=0.5_sp*(c-b)
			if (abs(xm) <= tol1 .or. fb == 0.0) then
				zbrent=b
				RETURN
			end if
			if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
				s=fb/fa
				if (a == c) then
					p=2.0_sp*xm*s
					q=1.0_sp-s
				else
					q=fa/fc
					r=fb/fc
					p=s*(2.0_sp*xm*q*(q-r)-(b-a)*(r-1.0_sp))
					q=(q-1.0_sp)*(r-1.0_sp)*(s-1.0_sp)
				end if
				if (p > 0.0) q=-q
				p=abs(p)
				if (2.0_sp*p  <  min(3.0_sp*xm*q-abs(tol1*q),abs(e*q))) then
					e=d
					d=p/q
				else
					d=xm
					e=d
				end if
			else
				d=xm
				e=d
			end if
			a=b
			fa=fb
			b=b+merge(d,sign(tol1,xm), abs(d) > tol1 )
			fb=func(b)
		end do
		call nrerror('zbrent: exceeded maximum iterations')
		zbrent=b
		END FUNCTION zbrent
		FUNCTION qsimp(func,a,b)
		USE nrtype1; USE nrutil1, ONLY : nrerror
	!	USE nr1, ONLY : trapzd
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: a,b
		REAL(SP) :: qsimp
		INTERFACE
			FUNCTION func(x)
			USE nrtype1
			REAL(SP), DIMENSION(:), INTENT(IN) :: x
			REAL(SP), DIMENSION(size(x)) :: func
			END FUNCTION func
		END INTERFACE
		INTEGER(I4B), PARAMETER :: JMAX=100 ! 20
		REAL(SP), PARAMETER :: EPS=1.0e-6_sp ! 6_sp
		INTEGER(I4B) :: j
		REAL(SP) :: os,ost,st
		ost=0.0
		os= 0.0
		do j=1,JMAX
			call trapzd(func,a,b,st,j)
			qsimp=(4.0_sp*st-ost)/3.0_sp
			if (j > 5) then
				if (abs(qsimp-os) < EPS*abs(os) .or. &
					(qsimp == 0.0 .and. os == 0.0)) RETURN
			end if
			os=qsimp
			ost=st
		end do
		call nrerror('qsimp: too many steps')
		END FUNCTION qsimp
		SUBROUTINE polint(xa,ya,x,y,dy)
		USE nrtype1; USE nrutil1, ONLY : assert_eq,iminloc,nrerror
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: xa,ya
		REAL(SP), INTENT(IN) :: x
		REAL(SP), INTENT(OUT) :: y,dy
		INTEGER(I4B) :: m,n,ns
		REAL(SP), DIMENSION(size(xa)) :: c,d,den,ho
		n=assert_eq(size(xa),size(ya),'polint')
		c=ya
		d=ya
		ho=xa-x
		ns=iminloc(abs(x-xa))
		y=ya(ns)
		ns=ns-1
		do m=1,n-1
			den(1:n-m)=ho(1:n-m)-ho(1+m:n)
			if (any(den(1:n-m) == 0.0)) &
				call nrerror('polint: calculation failure')
			den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
			d(1:n-m)=ho(1+m:n)*den(1:n-m)
			c(1:n-m)=ho(1:n-m)*den(1:n-m)
			if (2*ns < n-m) then
				dy=c(ns+1)
			else
				dy=d(ns)
				ns=ns-1
			end if
			y=y+dy
		end do
		END SUBROUTINE polint
		FUNCTION qromb(func,a,b)
		USE nrtype1; USE nrutil1, ONLY : nrerror
	!	USE nr1, ONLY : polint,trapzd
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: a,b
		REAL(SP) :: qromb
		INTERFACE
			FUNCTION func(x)
			USE nrtype1
			REAL(SP), DIMENSION(:), INTENT(IN) :: x
			REAL(SP), DIMENSION(size(x)) :: func
			END FUNCTION func
		END INTERFACE
		INTEGER(I4B), PARAMETER :: JMAX=100,JMAXP=JMAX+1,K=5,KM=K-1
		REAL(SP), PARAMETER :: EPS=1.0e-4_sp ! 10_sp
		REAL(SP), DIMENSION(JMAXP) :: h,s
		REAL(SP) :: dqromb
		INTEGER(I4B) :: j
		h(1)=1.0
		do j=1,JMAX
			call trapzd(func,a,b,s(j),j)
			if (j >= K) then
				call polint(h(j-KM:j),s(j-KM:j),0.0_sp,qromb,dqromb)
				if (abs(dqromb) <= EPS*abs(qromb)) RETURN
			end if
			s(j+1)=s(j)
			h(j+1)=0.25_sp*h(j)
		end do
		call nrerror('qromb: too many steps')
		END FUNCTION qromb
		FUNCTION brent(ax,bx,cx,func,tol,xmin)
		USE nrtype1; USE nrutil1, ONLY : nrerror
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: ax,bx,cx,tol
		REAL(SP), INTENT(OUT) :: xmin
		REAL(SP) :: brent
		INTERFACE
			FUNCTION func(x)
			USE nrtype1
			IMPLICIT NONE
			REAL(SP), INTENT(IN) :: x
			REAL(SP) :: func
			END FUNCTION func
		END INTERFACE
		INTEGER(I4B), PARAMETER :: ITMAX=100
		REAL(SP), PARAMETER :: CGOLD=0.3819660_sp,ZEPS=1.0e-3_sp*epsilon(ax)
		INTEGER(I4B) :: iter
		REAL(SP) :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
		a=min(ax,cx)
		b=max(ax,cx)
		v=bx
		w=v
		x=v
		e=0.0
		fx=func(x)
		fv=fx
		fw=fx
		do iter=1,ITMAX
			xm=0.5_sp*(a+b)
			tol1=tol*abs(x)+ZEPS
			tol2=2.0_sp*tol1
			if (abs(x-xm) <= (tol2-0.5_sp*(b-a))) then
				xmin=x
				brent=fx
				RETURN
			end if
			if (abs(e) > tol1) then
				r=(x-w)*(fx-fv)
				q=(x-v)*(fx-fw)
				p=(x-v)*q-(x-w)*r
				q=2.0_sp*(q-r)
				if (q > 0.0) p=-p
				q=abs(q)
				etemp=e
				e=d
				if (abs(p) >= abs(0.5_sp*q*etemp) .or. &
					p <= q*(a-x) .or. p >= q*(b-x)) then
					e=merge(a-x,b-x, x >= xm )
					d=CGOLD*e
				else
					d=p/q
					u=x+d
					if (u-a < tol2 .or. b-u < tol2) d=sign(tol1,xm-x)
				end if
			else
				e=merge(a-x,b-x, x >= xm )
				d=CGOLD*e
			end if
			u=merge(x+d,x+sign(tol1,d), abs(d) >= tol1 )
			fu=func(u)
			if (fu <= fx) then
				if (u >= x) then
					a=x
				else
					b=x
				end if
				call shft(v,w,x,u)
				call shft(fv,fw,fx,fu)
			else
				if (u < x) then
					a=u
				else
					b=u
				end if
				if (fu <= fw .or. w == x) then
					v=w
					fv=fw
					w=u
					fw=fu
				else if (fu <= fv .or. v == x .or. v == w) then
					v=u
					fv=fu
				end if
			end if
		end do
		call nrerror('brent: exceed maximum iterations')
		CONTAINS
	!BL
		SUBROUTINE shft(a,b,c,d)
		REAL(SP), INTENT(OUT) :: a
		REAL(SP), INTENT(INOUT) :: b,c
		REAL(SP), INTENT(IN) :: d
		a=b
		b=c
		c=d
		END SUBROUTINE shft
		END FUNCTION brent
		SUBROUTINE midpnt(func,a,b,s,n)
		USE nrtype1; USE nrutil1, ONLY : arth
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: a,b
		REAL(SP), INTENT(INOUT) :: s
		INTEGER(I4B), INTENT(IN) :: n
		INTERFACE
			FUNCTION func(x)
			USE nrtype1
			REAL(SP), DIMENSION(:), INTENT(IN) :: x
			REAL(SP), DIMENSION(size(x)) :: func
			END FUNCTION func
		END INTERFACE
		REAL(SP) :: del
		INTEGER(I4B) :: it
		REAL(SP), DIMENSION(2*3**(n-2)) :: x
		if (n == 1) then
			s=(b-a)*sum(func( (/0.5_sp*(a+b)/) ))
		else
			it=3**(n-2)
			del=(b-a)/(3.0_sp*it)
			x(1:2*it-1:2)=arth(a+0.5_sp*del,3.0_sp*del,it)
			x(2:2*it:2)=x(1:2*it-1:2)+2.0_sp*del
			s=s/3.0_sp+del*sum(func(x))
		end if
		END SUBROUTINE midpnt
	END MODULE nr1

