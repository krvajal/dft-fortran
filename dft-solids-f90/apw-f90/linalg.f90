
module linalg
use double



private 

public det

contains

    !compute the determinant of a square
    !matrix 
    !Parameters: 
    !           M: integer > 0
    !           A: a real square matrix with dimension (M,M)
    !Return:    
    !       the determinant of the matrix
    real(dp) function det(M,A) 
        implicit none
        integer,intent(in) :: M
        real(dp),dimension(M,M) :: A
        real(dp),dimension(M,M) :: LU
        integer :: info,i
        integer,dimension(M) :: pivot


        LU = A
        call dgetrf(M,M,LU,M,pivot, info)

        if (info /= 0)then
            stop "det: lu decomposition failed"
        endif
        det = 1
        ! det  = (-1)^(S) * prod(diag(L))*prod(diag(U))
        ! but diag of L is always 1
        ! so compute prod(diag(U))
        ! S is the number of row exchanges
        do i =1,M

            det = det * LU(i,i)
            if(pivot(i) /= i) then
                ! row exchange here
                det = - det;
            endif
        enddo
    end function det
end module linalg