module utils

use double
use quad
implicit none
private 
public FillRadialFunc, solve_radial, norm, SaveMatrix
contains


! Fill F=F(E,l,V) to be used in the solution of
! the radial schroedinger equation
subroutine FillRadialFunc(E,V,l,r,RadialFunc)
    implicit none
    real(dp),intent(in) :: E
    real(dp), intent(in) :: V(:), r(:)
    integer, intent(in) :: l 
    real(dp), intent(out) :: RadialFunc(:)
    
    RadialFunc(2:) = -2* (E + V(2:)/r(2:)/2) + l*(l + 1)/(r(2:)**2) 
    
end subroutine FillRadialFunc


! solve the radial function from zero outwards
! up to an rmax
! we choose u(r= 0) = 0 so it is solved for any energy
 subroutine  solve_radial(r, RadialFunc,PhiStart, PhiNext,PhiRadialPoints)
        implicit none
    
        real(dp), intent(out)     :: PhiRadialPoints(:)
        real(dp), intent(in)       :: r(:), RadialFunc(:)
        real(dp),intent(in) :: PhiStart, PhiNext
        real(dp),dimension(size(r)):: w 
        integer  :: n , i
        real(dp) :: h
        real(dp) :: rmax 
        
        ! print *,'solving for E = ' , params(1)

        n= size(PhiRadialPoints)
        h = r(2) - r(1)
        rmax = r(n)
        
        ! solve using numerov stepper
        PhiRadialPoints(1) = PhiStart
        PhiRadialPoints(2) =  PhiNext
        
        w(1) = 0 
        w(2) = (1- h*h/12_dp * RadialFunc(2))* PhiRadialPoints(2)
        ! integrate 
        do i = 3,n
          w(i) = 2*w(i-1) - w(i-2) + h*h*RadialFunc(i - 1)*PhiRadialPoints(i-1)
          PhiRadialPoints(i) = w(i)/(1 - h*h/12_dp * RadialFunc(i))
        enddo

end subroutine solve_radial

  !==============================================
    ! compute norm for a reciprocal lattice point 
    ! with indices l,m,n for a bcc lattice
    ! Inputs:
    ! a real the cell 
    ! l,m,n integers
    !==============================================
    real(dp) function norm(a,l,m,n)
        real(dp) :: a 
        integer :: l,m,n

        norm = 2 * pi / a
        norm = norm  * sqrt((3*l**2 + 3* m**2 + 3*n**2 &
                - 2*l*m - 2* l*n - 2*m*n)*1.0_dp)


    end function norm

subroutine SaveMatrix(filename,N,A)
    implicit none
    character(len=*) :: filename
    integer,intent(in) :: N
    real(dp) :: A(N,N)
    
    integer :: u, i, j
    print *, filename
    open(newunit=u,file=filename,status = "replace")
    
    do i =1,N
        write (u,*) (A(i,j),j=1,N)
    enddo
    close(u)

end subroutine SaveMatrix

end module utils