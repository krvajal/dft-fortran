program apw 
    use double
    implicit none
    
    integer,parameter :: N = 27 !basis dimension
    integer,parameter :: lmax = 4
    real(dp),parameter :: d = 6.822, R = 1 ! a.u.
    real(dp) :: K(N,3)
    real(dp) :: A(N,N), B(N,N), C(N,N,lmax)
    integer,parameter :: num_points = 1000
    real(dp) ::  p(num_points,3)
    real(dp) ::  rho_point(3), x_point(3)
    integer :: i
    real(dp),parameter :: omega = 3*d**3/4

    rho_point  = (/ 0.0_dp, 0.0_dp, 0.0_dp/); 
    x_point = (/ 2*pi/d, 0.0_dp, 0.0_dp/);

    call generate_K(N, K)
    call generate_p(num_points, p, rho_point, x_point )

    do i = 1,num_points
        call compute_A()
        call compute_B()
        call compute_C()

    enddo




    contains
    subroutine generate_K(N, K)
        implicit none
        integer,intent(in) :: N
        real(dp),intent(out) :: K(N,3)
        integer :: l,m, nn
        integer :: count   = 1
        real(dp) :: cutoff
        cutoff = 3.0_dp
        do l = -6,6
            do m =  -6,6
                  do nn = -6,6
                    if (norm(l,m,nn) < cutoff) then
                        print *,count
                        print *, l,m,nn
                        K(count,:) = 2*pi/d * (/-l + m  + nn, -m + l + nn, m + l  - nn  /)
                        count = count + 1
                        
                    endif 
                enddo
            enddo
        enddo
    end subroutine generate_K

    subroutine generate_p(N,p,pmin,pmax)
        implicit none
        integer, intent(in) :: N
        real(dp),intent(out) :: p(N,3)
        integer:: i
        real(dp) :: spacing
        real(dp),intent(in) :: pmax(3)
        real(dp),intent(in) :: pmin(3)
        real(dp) :: len
        real(dp) :: direct_vec(3)
        direct_vec  = pmax - pmin
        len = sqrt( sum( (direct_vec)**2))
        direct_vec = direct_vec/len !normalize vector

        spacing = len/(N-1)

        do i = 1,N-1
            p(i+1,:)= pmin + i * spacing*direct_vec
        enddo

    end subroutine generate_p


    subroutine compute_A(N,K,p,A)
        implicit none
        integer, intent(in) :: N
        real(dp), intent(in) :: K(N,3), p(3)
        real(dp),intent(out) :: A(N,N)
        integer :: i,j
        real(dp) ::  kij

        real(dp) :: jli, jlj

        do i =1,N
            do j=1,N
                A(i,j)  = - 4*pi*R**2/omega
                kij  = sqrt( sum( (K(i,:) - K(j,:))**2))
                if (i == j ) then
                    A(i,j) = A(i,j) * bessel_jn(1, kij*R)/kij
                else
                     A(i,j) = A(i,j)  * (1.0_dp/3) ! jl(x)->x/3 for small x
                endif
                call pm_polynomial_value()
            enddo
            A(i,i) = A(i,i) + 1
        enddo



    end subroutine compute_A

    subroutine compute_B(N,K,p,A,B)
        implicit none
        integer, intent(in) :: N
        real(dp), intent(in) :: K(N,3), p(3), A(N,N)
        real(dp),intent(out) :: B(N,N)
        real(dp) :: qi(3), qj(3)
        integer :: i,j
        real(dp) ::  kij
        
        do i =1,N
            qi = p + K(i,:)
            do j = 1,N
                qj= p + K(j,:)
                B(i,j) = 0.5_dp * A(i,j) * dot_product( qi,qj)
            enddo
        enddo


    end subroutine compute_B

    subroutine compute_C(N,lmax,K,p,C)
        implicit none
        integer, intent(in) :: N,lmax
        real(dp),intent(in) :: K(N,3), p(3)
        real(dp), intent(out) :: C(N,N, lmax + 1)
        integer :: i,j, l
        real(dp) :: qi(3), qj(3), lqi,lqj
        real(dp) :: x

        C = 2*pi * R*R/omega
        do i =1,N
            qi = p + K(i,:)
            lqi = sqrt(sum(qi**2))
            do j = 1, N
                qj = p + K(j,:)
                lqj = sqrt(sum(qj**2))
                do l = 0,lmax
                    x = dot_product(qi,qj)/(lqi*lqj)
                    C(i,j,l + 1) = C(i,j,l + 1) * (2*l + 1) 
                    C(i,j,l+1) = C(i,j,l+1) * pl( x,l)
                    C(i,j,l+1) = C(i,j,l+1) * bessel_jn(l, lqi*R) * bessel_jn(l, lqj*R)
                enddo
            enddo
        enddo

        
    end subroutine compute_C


    real(dp) function norm(l,m,n)
        integer :: l,m,n

        norm = 2 * pi / d
        norm = norm  * sqrt((3*l**2 + 3* m**2 + 3*n**2 &
                - 2*l*m - 2* l*n - 2*m*n)*1.0_dp)


    end function norm


    real(dp)    function pl(x,n)
    !======================================
    ! calculates Legendre polynomials Pn(x)
    ! using the recurrence relation
    ! if n > 100 the function retuns 0.0
    !======================================
        integer :: n, k
        real(dp) :: x
        real(dp),dimension(0:n) :: pln


        pln(0) = 1.0
        pln(1) = x

        if (n <= 1) then
            pl = pln(n)
        else
            do k=1,n-1
                pln(k+1) = ((2.0*k+1.0)*x*pln(k) - float(k)*pln(k-1))/(float(k+1))
            enddo
            pl = pln(n)
        end if
        return
    end function pl

end program apw