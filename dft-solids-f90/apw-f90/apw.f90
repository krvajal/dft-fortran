program apw 
    use double
    use utils
    use potential
    use linalg
    use special
    use legendre_pols
    implicit none
    
    integer,parameter ::NumRadialPoints =  1000
    integer,parameter :: N = 27 !basis dimension
    integer,parameter :: lmax = 4
    real(dp),parameter :: d = 6.822, R = 1 ! a.u.
    real(dp) :: K(N,3)
    real(dp) :: A(N,N), B(N,N), C(N,N,lmax + 1), Hamilt(N,N)
    integer,parameter :: NumPPoints = 1000
    real(dp) ::  PPoints(NumPPoints,3)
    real(dp) ::  RhoPoint(3), XPoint(3)
    integer :: i,j
    real(dp),parameter :: omega = 3*d**3/4
    real(dp) ::RadialFunc(NumRadialPoints),PhiRadialPoints(NumRadialPoints)
        
    real(dp) :: EnergyPoints(NumPPoints)
    real(dp) :: MinE, MaxE, DeltaE, HDet

    MinE = 0.02
    MaxE = 0.021
    DeltaE = (MaxE - MinE)/(NumPPoints-1)
    
    EnergyPoints = MinE + (/ (i*DeltaE,i=0,NumPPoints-1)/)

    !points in the reciprocal unit cell
    RhoPoint  = (/ 0.0_dp, 0.0_dp, 0.0_dp/); 
    XPoint = (/ 2*pi/d, 0.0_dp, 0.0_dp/);

    call LoadPotential
    
    call generate_K(N, K)
    call generate_p(NumPPoints, PPoints, RhoPoint, XPoint )
    ! they dont depend on the energy
    call compute_A(N,K,RhoPoint,A)

    call compute_B(N,K,RhoPoint,A,B)
    
    call compute_C(N,Lmax,K,RhoPoint,C)
    

    do i = 1, NumPPoints

        call SolveAtomic(EnergyPoints(i), RhoPoint, HDet)
        print *, EnergyPoints(i), HDet
    enddo


    contains
    

    subroutine SolveAtomic(E,p, HamiltDeterm)
        implicit none
        real(dp) :: E
        real(dp) :: PhiEnd(lmax + 1), PhiDerEnd(lmax + 1)
        integer :: l,i,j
        real(dp) :: Delta
        real(dp),intent(out) :: HamiltDeterm
        real(dp) :: p(3)
        real(dp), dimension(NumRadialPoints) :: PhiRadialPoints
        real(dp) :: PhiStart,PhiNext
        Delta = RadialPoints(NumRadialPoints) - RadialPoints(NumRadialPoints-1)

        !solve for each l

        do l = 0,lmax
            !generate the function points to solve the radial equation

            call FillRadialFunc(E, Pot, l, RadialPoints, RadialFunc)

            !solve the radial equation for a given L, V, y E
            !initial values for numerov
            
            if (l .eq. 0) then
                 PhiStart = 2*Delta*Delta/12._dp*29.0_dp
            !  Analytic expression; 29 is the nuclear charge of copper
            else
                PhiStart = 0_dp
            endif
            
            PhiNext = Delta ** (l+1)
            PhiRadialPoints = 0
            call solve_radial(RadialPoints,RadialFunc,PhiStart,PhiNext, PhiRadialPoints)
            ! print *,RadialPoints(NumPotPoints)

            PhiEnd(l + 1) = PhiRadialPoints(NumPotPoints)
            ! print *, PhiEnd
            PhiEnd(l + 1) = PhiEnd(l + 1)/RMuffinTin
          
            !compute derivative here
            PhiRadialPoints = PhiRadialPoints/RadialPoints

            PhiDerEnd(l+1) = +11.0_dp/6.0_dp*PhiRadialPoints(NumRadialPoints) &
                             - 3.0_dp*PhiRadialPoints(NumRadialPoints-1) & 
                             + 1.5_dp*PhiRadialPoints(NumRadialPoints-2) &
                             - PhiRadialPoints(NumRadialPoints-3)/3.0_dp
            PhiDerEnd(l+1) = PhiDerEnd(l+1)/Delta

            ! now calculate the determinant of |H - EI| for a mesh of E values 
            ! and find the roots of it


        enddo

        
        call evaluate_H(N,lmax,E,A,B,C,PhiEnd,PhiDerEnd,Hamilt)
        

        HamiltDeterm = det(N,Hamilt)
        ! HamiltDeterm = Hamilt(1,1)
    end subroutine SolveAtomic    

    subroutine evaluate_H(N,lmax,energy,MatA,MatB,MatC, PhiEnd, PhiDerEnd, Hamilt )
        implicit none
        integer, intent(in) :: N, lmax
        real(dp),intent(in) :: energy, MatA(N,N), MatB(N,N), MatC(N,N,lmax+1)
        real(dp),intent(in) :: PhiEnd(lmax+1),PhiDerEnd(lmax+1)
        real(dp), intent(out) :: Hamilt(N,N)
    
        integer :: i
        real(dp) :: CC(N,N)

        
        CC = 0.0_dp !init to zero

        do i = 0, lmax
            CC = CC + MatC(:,:,i+1) * PhiDerEnd(i+1)/PhiEnd(i+1)
        enddo
               
        Hamilt = - energy * MatA + MatB + CC
        do i = 1,N
            do j = 1,i
                Hamilt(j,i)= Hamilt(i,j)
            enddo
        enddo
        
    end subroutine evaluate_H
    
    subroutine generate_K(N, K)
        implicit none
        integer,intent(in)  :: N
        real(dp),intent(out)    :: K(N,3)
        integer     :: l,m, nn
        integer     :: count   = 1
        real(dp)    :: cutoff
        cutoff = 3.0_dp
        do l = -6,6
            do m =  -6,6
                  do nn = -6,6
                    if (norm(d,l,m,nn) < cutoff) then
                        ! print *,count
                        ! print *, l,m,nn
                        K(count,:) = 2*pi/d * (/-l + m  + nn, -m + l + nn, m + l  - nn  /)
                        count = count + 1
                        
                    endif 
                enddo
            enddo
        enddo
    end subroutine generate_K

    subroutine generate_p(N,p,pmin,pmax)
        implicit none
        integer, intent(in)     :: N
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


    !===================================
    ! Compute Aij for for the Hamiltonian
    ! Inputs:
    !  N dimension of the basis
    !  K array with the values of the reciprocal lattice points that
    !    are smaller than a cutoff energy
    !  p the value of k inside the first Bruillint zone
    !  A a real matrix with dimension NxN
    !===================================
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
                A(i,j)  = - 4*pi*RMuffinTin**2/omega
                kij  = sqrt( sum( (K(i,:) - K(j,:))**2))
                
                if (i /= j ) then
                    A(i,j) = A(i,j) * spherical_bessel_jn(1, kij*RMuffinTin)/kij
                else
                     A(i,j) = A(i,j)  /3.0_dp ! jl(x)->x/3 for small x
                endif
                
            enddo
            A(i,i) = A(i,i) + 1
        enddo



    end subroutine compute_A

    !===================================
    ! Compute Bij for for the Hamiltonian
    ! Inputs:
    !  N dimension of the basis
    !  K array with the values of the reciprocal lattice points that
    !    are smaller than a cutoff energy
    !  p the value of k inside the first Bruillint zone
    !  A a real matrix with dimension NxN
    !  B a real matrix with dimension NxN
    !===================================

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

    !===================================
    ! Compute Cijl for for the Hamiltonian
    ! Inputs:
    !  N dimension of the basis
    !  lmx integer the cutoff value for l 
    !  K array with the values of the reciprocal lattice points that
    !    are smaller than a cutoff energy
    !  p the value of k inside the first Bruillint zone
    !  C real cube with dimension NxNx(lmax + 1)
    !===================================

    subroutine compute_C(N,lmax,K,p,C)
        implicit none
        integer, intent(in) :: N,lmax
        real(dp),intent(in) :: K(N,3), p(3)
        real(dp), intent(out) :: C(N,N, lmax + 1)
        integer :: i,j, l
        real(dp) :: qi(3), qj(3), LengthQi,LengthQj
        real(dp) :: x

        C = 2*pi * RMuffinTin**2/omega
        do i =1,N
            qi = p + K(i,:)
            LengthQi = sqrt(sum(qi**2))
            do j = 1, N
                qj = p + K(j,:)
                LengthQj = sqrt(sum(qj**2))
                do l = 0,lmax
                    if(LengthQj > 1e-8 .and. LengthQi > 1e-8) then
                        x = dot_product(qi,qj)/(LengthQi*LengthQj)
                    else 
                        x = 1.0-1.D-10
                    endif

                    C(i,j,l + 1) = C(i,j,l + 1) * (2*l + 1) 
                    C(i,j,l+1) = C(i,j,l+1) * Legendre(l,x)
                    C(i,j,l+1) = C(i,j,l+1) * spherical_bessel_jn(l, LengthQi*RMuffinTin) &
                                 * spherical_bessel_jn(l, LengthQj*RMuffinTin)
                enddo
            enddo
        enddo

        
    end subroutine compute_C

  

end program apw
