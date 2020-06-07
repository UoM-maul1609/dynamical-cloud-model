	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>parallel tridag solver
	!> see https://web.alcf.anl.gov/~zippy/publications/partrid/partrid.html
	!> but note there are two errors in part 1 of the algorithm (corrected here)
    module pts
    use nrtype
    use mpi
    implicit none
    !>variables for mpi
    type pts_vars
        integer(i4b), dimension(:), allocatable :: kps, ids
		integer(i4b), allocatable, dimension(:) :: request
		integer(i4b), dimension(MPI_STATUS_SIZE, 12) :: status
		integer(i4b) :: error, tag1,num_messages,imess, tag2,id2,nprocv,sz,thispe
        real(sp), allocatable, dimension(:) :: mvrecv
        real(sp), allocatable, dimension(:) :: outdata,reduca,reducb,reducc,reducr, &
                        coeffs,a,b,c,r,an,bn,cn,rn,un,up,xsol

    end type pts_vars
    
    ! declare a pts type
    type(pts_vars) :: pts1
    

    private
    public :: set_tridag,parallel_tridag,re_write_full,pts1
    
    contains
    
		
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>@author
    !>Paul J. Connolly, The University of Manchester
    !>@brief
    !>allocate and initialise tridiagonal solver 
    !>@param[in] init_flag: if true set to a standard problem, if false do not set yet
    !>@param[inout] an,bn,cn,rn,un,up,a,b,c,r,x: variables for matrices
    !>@param[in] kp,kpp,kpstart: dimensions for problem
    !>@param[in] comm3d,id, dims, coords: mpi stuff
	subroutine set_tridag(init_flag,an,bn,cn,rn,un,up,a,b,c,r,x, &
	    kp,kpp,kpstart,coords,dims,id,comm3d)
        use nr, only : tridag
        use mpi
        use mpi_module
        implicit none
        logical :: init_flag
        integer(i4b), intent(in) :: kp,kpp,kpstart
        integer(i4b), dimension(3), intent(inout) :: coords
        integer(i4b), dimension(3), intent(in) :: dims
        integer(i4b), intent(in) :: id, comm3d

        real(sp), dimension(:),allocatable, intent(inout) :: bn,rn,un,up
        real(sp), dimension(:),allocatable, intent(inout) :: an,cn
        
        real(sp), dimension(:),allocatable, intent(inout) :: a,b,c,r,x
	
	    !local variables
	    integer(i4b) :: AllocateStatus,error,i
	

        call MPI_CART_RANK( comm3d, coords(3), pts1%id2,error)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! allocate arrays                                                                !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		allocate( a(1:kpp), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( c(1:kpp), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( b(1:kpp), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( r(1:kpp), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( x(1:kpp), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"

        ! lower diag:
        a(:) = 0.5_sp
        ! upper diag:
        c(:) = 1._sp
        ! diagonal:
        b(:)=2._sp
        ! rhs:
        r(:) = 1._sp

		if(init_flag) then
            allocate( an(1:kp-1), STAT = AllocateStatus)
            if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
            allocate( cn(1:kp-1), STAT = AllocateStatus)
            if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
            allocate( bn(1:kp), STAT = AllocateStatus)
            if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
            allocate( rn(1:kp), STAT = AllocateStatus)
            if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
            allocate( un(1:kp), STAT = AllocateStatus)
            if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
            allocate( up(1:kp), STAT = AllocateStatus)
            if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! set the problem                                                            !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! lower diag:
            an(:) = 0.5_sp
            ! upper diag:
            cn(:) = 1._sp
            ! diagonal:
            bn(:)=2._sp
            ! rhs:
            rn(:) = 1._sp
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! split the problem onto PEs                                                 !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            b=bn(kpstart:kpstart+kpp-1) ! set submatrix diagonal
            r=rn(kpstart:kpstart+kpp-1) ! set submatrix rhs 
            ! set submatrix lower diagonal
            if(coords(3)==0) then
                a(1)=0._sp
                a(2:kpp)=an(kpstart:kpstart+kpp-2) 
            else
                a=an(kpstart-1:kpstart+kpp-2)
            endif
        
            ! set submatrix upper diagonal
            if(coords(3)==(dims(3)-1)) then
                c(kpp)=0._sp
                c(1:kpp-1)=cn(kpstart:kpstart+kpp-2) 
            else
                c=cn(kpstart:kpstart+kpp-1)
            endif
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        endif  
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! allocate variables needed for solver
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        pts1%nprocv=dims(3)
        pts1%sz=2*pts1%nprocv-2
        pts1%thispe=coords(3)+1


        ! allocate outdata array for message passing
        allocate( pts1%outdata(1:pts1%nprocv*8), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"

        allocate( pts1%mvrecv(1:8*pts1%nprocv), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"

        
        allocate(pts1%reduca(1:pts1%sz), STAT=AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
        allocate(pts1%reducb(1:pts1%sz), STAT=AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
        allocate(pts1%reducc(1:pts1%sz), STAT=AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
        allocate(pts1%reducr(1:pts1%sz), STAT=AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
        allocate(pts1%coeffs(1:pts1%sz), STAT=AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"

        call MPI_COMM_SIZE(comm3d,i,error)

        if(pts1%nprocv /= i) stop

        allocate(pts1%kps(1:i), STAT=AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
        allocate(pts1%ids(1:i), STAT=AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
        allocate(pts1%request(1:i), STAT=AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! serial problem solver                                                      !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if((pts1%id2==0) .and. init_flag) then
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! solve tridiagonal matrix:								    		     !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call tridag(an,bn,cn,rn,un)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		
	end subroutine set_tridag
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>@author
    !>Paul J. Connolly, The University of Manchester
    !>@brief
    !>write elements of solution to 1-d array on pe 1
    !>@param[inout] up,un,x: solutions
    !>@param[in] kp,kpp,kpstart: dimensions for problem
    !>@param[in] comm3d,id, dims, coords: mpi stuff
	subroutine re_write_full(up,un,x, &
	    kp,kpp,kpstart,coords,dims,id,comm3d)
        use nr, only : tridag
        use mpi
        use mpi_module
        implicit none
        real(sp), intent(inout),dimension(kp) :: up
        real(sp), intent(in),dimension(kp) :: un
        real(sp), intent(in),dimension(kpp) :: x
        integer(i4b), intent(in) :: kp,kpp,kpstart
        integer(i4b), dimension(3), intent(inout) :: coords
        integer(i4b), dimension(3), intent(in) :: dims
        integer(i4b), intent(in) :: id, comm3d

	    !local variables
	    integer(i4b) :: tag1, error,i
	
	    if(pts1%nprocv .gt. 1) then
            tag1=10
            call MPI_allgather(kpp,1,MPI_INTEGER,pts1%kps,1,MPI_INTEGER,&
                comm3d,pts1%status(:,1),error)

            ! id2 is the id in this communicator
            ! nprocv is number of processors in this topology
            if(pts1%id2/=0) then
                call MPI_issend(x,kpp, MPI_REAL8, 0, tag1, comm3d, pts1%request(1),error)
            endif

            if(pts1%id2 == 0) then
                up(1:kpp)=x
                do i=2,pts1%nprocv
                    call MPI_recv(up(sum(pts1%kps(1:i-1))+1:pts1%kps(i)),&
                        pts1%kps(i), MPI_REAL8, i-1, tag1, comm3d, &
                        pts1%status(:,1),error)
                enddo
            endif
            if(pts1%id2/=0) &
            call MPI_Wait(pts1%request(1), pts1%status(:,1), error)
        endif          
          
                
        if(pts1%id2==0) then 
            print *,up -un
            print *,' '
        endif
	end subroutine re_write_full
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>@author
    !>Paul J. Connolly, The University of Manchester
    !>@brief
    !>tridiagonal solver 
    !>@param[inout] a,b,c,r,x: variables for matrices
    !>@param[in] kp,kpp,kpstart: dimensions for problem
    !>@param[in] comm3d,id, dims, coords: mpi stuff
    subroutine parallel_tridag(a,b,c,r,x, &
        kpp,coords,dims, id, comm3d)
        use nr, only : tridag
        use mpi
        use mpi_module
        implicit none
        integer(i4b), intent(in) :: kpp
        integer(i4b), dimension(3), intent(in) :: coords
        integer(i4b), dimension(3), intent(in) :: dims
        integer(i4b), intent(in) :: id, comm3d
        
        real(sp), dimension(kpp), intent(inout) :: a,b,c,r,x
        
        ! local variables:        
        real(sp), dimension(kpp) :: xuh,xlh,xr
        real(sp), dimension(8) :: msend
        real(sp) :: denom, uhcoeff,lhcoeff
        integer(i4b) :: i,j,k, ibase,ifirst,nsig,error
        integer(i4b), dimension(3) :: coords_t

            

        
       
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! part 1: solve sub problems on each PE                                      !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        xuh(1)=c(1)/b(1)
        xlh(1)=r(1)/b(1)
        ! forward elimination
        do i=2,kpp
            denom=b(i)-a(i)*xuh(i-1)
            if(denom .eq. 0._sp) stop
            xuh(i)=c(i)/denom
            xlh(i)=(r(i)-a(i)*xlh(i-1))/denom
        enddo
        

        ! back substitution
        xr(kpp)=xlh(kpp)
        xlh(kpp)=-xuh(kpp)
        xuh(kpp)=a(kpp)/b(kpp)
        do i=kpp-1,1,-1
            xr(i)=xlh(i)-xuh(i)*xr(i+1)
            xlh(i)=-xuh(i)*xlh(i+1)
            denom=b(i)-c(i)*xuh(i+1)
            denom=sign(1._sp,denom+1.e-60_sp)*max(abs(denom),1.e-15)
            if(denom .eq. 0._sp) stop
            xuh(i)=a(i)/denom
        enddo
        
        ! forward substitution
        xuh(1)=-xuh(1)
        do i=2,kpp
            xuh(i)=-xuh(i)*xuh(i-1)
        enddo
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        if(pts1%nprocv==1) then
            x=xr
            return
        endif
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! part 2: find data for reduced tridiagonal form                             !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        pts1%outdata=0._sp
        ! set outdata array for message passing
        pts1%outdata(1+(pts1%thispe-1)*8)=-1._sp
        pts1%outdata(2+(pts1%thispe-1)*8)=xuh(1)
        pts1%outdata(3+(pts1%thispe-1)*8)=xlh(1)
        pts1%outdata(4+(pts1%thispe-1)*8)=-xr(1)
        pts1%outdata(5+(pts1%thispe-1)*8)=xuh(kpp)
        pts1%outdata(6+(pts1%thispe-1)*8)=xlh(kpp)
        pts1%outdata(7+(pts1%thispe-1)*8)=-1._sp
        pts1%outdata(8+(pts1%thispe-1)*8)=-xr(kpp) 
        ! http://mpitutorial.com/tutorials/mpi-scatter-gather-and-allgather/
        msend=pts1%outdata(1+(pts1%thispe-1)*8:8+(pts1%thispe-1)*8)  ! message to send
        
        call MPI_allgather(msend,8,MPI_REAL8,pts1%mvrecv,8,MPI_REAL8,comm3d,error)
        pts1%outdata=pts1%mvrecv
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! part 3: put outdata into reduced tridiagonal form and solve                !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        nsig=8*pts1%nprocv
        ifirst=5
        do i=1,2*pts1%nprocv-2
            ibase=modulo(ifirst+4*(i-1),nsig)
            pts1%reduca(i)=pts1%outdata(ibase)
            pts1%reducb(i)=pts1%outdata(ibase+1)
            pts1%reducc(i)=pts1%outdata(ibase+2)
            pts1%reducr(i)=pts1%outdata(ibase+3)
        enddo
        ! solve reduced system
        call tridag(pts1%reduca(2:pts1%sz),pts1%reducb,pts1%reducc(1:pts1%sz-1),&
            pts1%reducr,pts1%coeffs)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !print *,xr
        !pause

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! part 4: pick out the appropriate elements of coeffs                        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(coords(3)==0) then
            uhcoeff=0._sp
        else
            uhcoeff=pts1%coeffs(2*pts1%thispe-2)
        endif
        
        if(coords(3)==(dims(3)-1)) then
            lhcoeff=0._sp
        else
            lhcoeff=pts1%coeffs(2*pts1%thispe-1)
        endif
        do i=1,kpp
            x(i)=xr(i)+uhcoeff*xuh(i)+lhcoeff*xlh(i)
        enddo
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        

                    
    end subroutine parallel_tridag
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	end module pts
	