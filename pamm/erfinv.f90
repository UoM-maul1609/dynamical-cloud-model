	!> @mainpage
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@copyright 2019
	!>@brief
	!>A module for calculation of the inverse erf function
    module erfinv_m
    use nrtype
    implicit none
    
    real(sp), parameter :: LN2 = 6.931471805599453094172321214581e-1_sp
    real(sp), parameter :: A0 = 1.1975323115670912564578e0_sp
    real(sp), parameter :: A1 = 4.7072688112383978012285e1_sp
    real(sp), parameter :: A2 = 6.9706266534389598238465e2_sp
    real(sp), parameter :: A3 = 4.8548868893843886794648e3_sp
    real(sp), parameter :: A4 = 1.6235862515167575384252e4_sp
    real(sp), parameter :: A5 = 2.3782041382114385731252e4_sp
    real(sp), parameter :: A6 = 1.1819493347062294404278e4_sp
    real(sp), parameter :: A7 = 8.8709406962545514830200e2_sp
    real(sp), parameter :: B0 = 1.0000000000000000000e0_sp
    real(sp), parameter :: B1 = 4.2313330701600911252e1_sp
    real(sp), parameter :: B2 = 6.8718700749205790830e2_sp
    real(sp), parameter :: B3 = 5.3941960214247511077e3_sp
    real(sp), parameter :: B4 = 2.1213794301586595867e4_sp
    real(sp), parameter :: B5 = 3.9307895800092710610e4_sp
    real(sp), parameter :: B6 = 2.8729085735721942674e4_sp
    real(sp), parameter :: B7 = 5.2264952788528545610e3_sp
    real(sp), parameter :: C0 = 1.42343711074968357734e0_sp
    real(sp), parameter :: C1 = 4.63033784615654529590e0_sp
    real(sp), parameter :: C2 = 5.76949722146069140550e0_sp
    real(sp), parameter :: C3 = 3.64784832476320460504e0_sp
    real(sp), parameter :: C4 = 1.27045825245236838258e0_sp
    real(sp), parameter :: C5 = 2.41780725177450611770e-1_sp
    real(sp), parameter :: C6 = 2.27238449892691845833e-2_sp
    real(sp), parameter :: C7 = 7.74545014278341407640e-4_sp
    real(sp), parameter :: D0 = 1.4142135623730950488016887e0_sp
    real(sp), parameter :: D1 = 2.9036514445419946173133295e0_sp
    real(sp), parameter :: D2 = 2.3707661626024532365971225e0_sp
    real(sp), parameter :: D3 = 9.7547832001787427186894837e-1_sp
    real(sp), parameter :: D4 = 2.0945065210512749128288442e-1_sp
    real(sp), parameter :: D5 = 2.1494160384252876777097297e-2_sp
    real(sp), parameter :: D6 = 7.7441459065157709165577218e-4_sp
    real(sp), parameter :: D7 = 1.4859850019840355905497876e-9_sp
    real(sp), parameter :: E0 = 6.65790464350110377720e0_sp
    real(sp), parameter :: E1 = 5.46378491116411436990e0_sp
    real(sp), parameter :: E2 = 1.78482653991729133580e0_sp
    real(sp), parameter :: E3 = 2.96560571828504891230e-1_sp
    real(sp), parameter :: E4 = 2.65321895265761230930e-2_sp
    real(sp), parameter :: E5 = 1.24266094738807843860e-3_sp
    real(sp), parameter :: E6 = 2.71155556874348757815e-5_sp
    real(sp), parameter :: E7 = 2.01033439929228813265e-7_sp
    real(sp), parameter :: F0 = 1.414213562373095048801689e0_sp
    real(sp), parameter :: F1 = 8.482908416595164588112026e-1_sp
    real(sp), parameter :: F2 = 1.936480946950659106176712e-1_sp
    real(sp), parameter :: F3 = 2.103693768272068968719679e-2_sp
    real(sp), parameter :: F4 = 1.112800997078859844711555e-3_sp
    real(sp), parameter :: F5 = 2.611088405080593625138020e-5_sp
    real(sp), parameter :: F6 = 2.010321207683943062279931e-7_sp
    real(sp), parameter :: F7 = 2.891024605872965461538222e-15_sp
    
    private
    public :: erfinv
    contains
    
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculate the inverse of erf
	!>@param[in] x: input integer
	!>@param[inout] ans: the output
	subroutine erfinv(x,ans) 
	    use, intrinsic :: IEEE_ARITHMETIC
        implicit none
        real(sp), intent(in) :: x
        real(sp), intent(inout) :: ans
        complex(sp) :: num, den
        real(sp) :: abs_x, r
        
        if((x>1._sp) .or. (x<-1._sp)) then
            ans=IEEE_VALUE(0._sp,IEEE_QUIET_NAN)
            return
        elseif((x==1._sp) ) then
            ans=IEEE_VALUE(0._sp,IEEE_POSITIVE_INF)
            return
        elseif((x==-1._sp) ) then
            ans=IEEE_VALUE(0._sp,IEEE_POSITIVE_INF)
        endif
        
        
        abs_x=abs(x)
        if(abs_x <= 0.85_sp) then
            r =  0.180625_sp - 0.25_sp * x * x
            num = (((((((A7 * r + A6) * r + A5) * r + A4) * r + A3) &
                * r + A2) * r + A1) * r + A0)
            den = (((((((B7 * r + B6) * r + B5) * r + B4) * r + B3) &
                * r + B2) * r + B1) * r + B0)
            ans = x * num / den
            return
        endif
        
        r = sqrt(LN2 - log(1.0_sp - abs_x))
        if (r <= 5.0_sp) then
            r = r - 1.6;
            num = (((((((C7 * r + C6) * r + C5) * r + C4) * r + C3) * &
                r + C2) * r + C1) * r + C0)
            den = (((((((D7 * r + D6) * r + D5) * r + D4) * r + D3) * &
                r + D2) * r + D1) * r + D0)
        else 
            r = r - 5.0_sp
            num = (((((((E7 * r + E6) * r + E5) * r + E4) * r + E3) * &
                r + E2) * r + E1) * r + E0)
            den = (((((((F7 * r + F6) * r + F5) * r + F4) * r + F3) * &
                r + F2) * r + F1) * r + F0)
        endif
        if (x < 0._sp) then
            ans= real(-num/den,sp)
            return
        else 
            ans = real(num/den,sp)
            return
        endif    
        
    end subroutine erfinv
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    end module erfinv_m
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
