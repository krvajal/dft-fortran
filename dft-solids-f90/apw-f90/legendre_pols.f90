module legendre_pols

use double

private 
public Legendre
contains

  real(dp)  function Legendre (L, X)
  ! returns the Legendre polynomial P_l(x) as a function of l and x.
  ! Upward recursion is used. 
    implicit none

    integer L, helpL

    real(dp)X, PL, PlMin1, PlMin2

     if (L.EQ.0) then
       Legendre = 1_dp
     else if (L.EQ.1) then
       Legendre = X
     else
       PlMin1 = 1.D0
       Pl     = X
       DO HelpL =2, L
         PlMin2 = PlMin1
         PlMin1 = Pl
         Pl = (2*HelpL-1)*X*PlMin1 - (HelpL-1)*PlMin2
         Pl = Pl/(HelpL)
       ENDDO
       Legendre = Pl
     END if
  end function Legendre
end module legendre_pols