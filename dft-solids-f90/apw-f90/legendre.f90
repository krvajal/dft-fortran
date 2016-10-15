module legendre
use double

implicit none
private
public legendre0
contains

    function legendre0(x) result(fx)
        implicit none
        real(dp),intent(in) :: x
        real(dp) :: fx
        fx = 1
        

    end function legendre0

    function legendre1(x) result(fx)
        implicit none
        real(dp),intent(in) :: x
        real(dp) :: fx
        fx = x
        

    end function legendre1

    function legendre2(x) result(fx)
        implicit none
        real(dp),intent(in) :: x
        real(dp) :: fx
        fx = 0.5 * (3 *x**2 - 1)    

    end function legendre2
    

end module legendre