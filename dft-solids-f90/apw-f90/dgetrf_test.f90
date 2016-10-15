program getrf_test
    implicit none
    integer,parameter :: N = 3

    real(8),dimension(N,N) :: matrix, U,L
    integer,dimension(N) :: pivot
    integer :: info,i,j
    real(8) :: det

    L = 0
    U = 0
    matrix(1,:) = (/1, 2,4/)
    matrix(2,:) = (/3 ,8,14/)
    matrix(3,:) = (/2,6,13/)

    write (*,10) ((matrix(i,j), j =1,N) ,i=1,N)
    ! init matrix

    call dgetrf(N,N,transpose(matrix),N,pivot)
    do i = 1,N
        do j = 1, i-1
            L(i,j) = matrix(i,j)
        enddo
        L(i,i) = 1
    enddo


    do i = 1,N
        do j = i, N
            U(i,j) = matrix(i,j)
        enddo
    enddo

    det = 1
    do i =1,N
        det = det * matrix(i,i);
        if(pivot(i) /= i) then
            det = - det;
        endif
    enddo
    print *, "det=",det


    print *,"==========="
    write (*,10) ((matrix(i,j), j =1,N) ,i=1,N)
    print *,"==========="
    write (*,10) ((L(i,j), j =1,N) ,i=1,N)
    print *,"==========="
    write (*,10) ((U(i,j), j =1,N) ,i=1,N)
    print *, info

    10 format(1x,3f12.6)
    
end program getrf_test
