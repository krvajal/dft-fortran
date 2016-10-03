! solve the Schroedinger equation using the Crank-Nicholson algorithm

program hw4
    implicit none
    real, parameter :: h = 0.001, L = 1.;
    real :: t = 0
    integer, parameter :: n = L/h
    integer :: iter, MAX = 1000;
    complex,dimension(n) ::  Wf
    real, dimension(n) :: V !the potential

    call init_pot1(V)
    do iter =  1, MAX
        call ckniterate(Wf,v,L,h);
    enddo


    contains
    ! use the potential for item 1
    subroutine init_pot1(V)
        implicit none
        real :: V(:)
        V = 0.0 !! just make it zero
    end subroutine init_pot1

    subroutine ckniterate(Wf, V, L, h)
        implicit none
        complex :: Wf(:)
        real :: V(:), L, h
        integer :: i, N
        complex,allocatable :: b(:)

        N = size(Wf)
        allocate(b(N))
        b = 1 +   2;




        if(allocated(b)) deallocate(b)
        
        
    end subroutine ckniterate

end program hw4
