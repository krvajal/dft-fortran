module potential
use double
implicit none
real(dp),allocatable :: Pot(:), RadialPoints(:)
integer  :: NumPotPoints = 1000
real(dp) :: RMuffinTin = 2.4151_dp
real(dp) :: Delta

contains
    subroutine LoadPotential()
        implicit none   
        integer :: u
        integer :: i
        open(newunit=u, file="linpot.txt", status="old")
        allocate(Pot(NumPotPoints))
        allocate(RadialPoints(NumPotPoints))

        print *, NumPotPoints

        do i =1 ,NumPotPoints
            read(u,*) RadialPoints(i), Pot(i)
        enddo
        
        
        close(u)

    end subroutine LoadPotential

end module potential