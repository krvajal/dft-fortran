

integer, parameter:: dp=kind(0.d0)                   ! double precision

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
    integer :: info,i,j
    real(dp) :: det


    call dgetrf(N,N,transpose(matrix),N,pivot)
    det = 1
    do i =1,N
        det = det * matrix(i,i);
        if(pivot(i) /= i) then
            det = - det;
        endif
    enddo
end function det