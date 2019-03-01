!************************************************************
! From Allen & Tildesley but originally from Knuth - `The
! art of computer programming', Addison-Wesley, Reading.
! Produces a normal (Gaussian) distribution with zero mean
! and unit variance.
!************************************************************

FUNCTION gauss(idum)

USE commondata, ONLY: dummy

IMPLICIT NONE

DOUBLE PRECISION :: gauss
!DOUBLE PRECISION, INTENT(IN) :: idum
real, INTENT(IN) :: idum
INTEGER :: j
DOUBLE PRECISION :: a1,a3,a5,a7,a9,x,r,rsq

a1=3.949846138
a3=0.252408784
a5=0.076542912
a7=0.008355968
a9=0.029899776
x=0.d0

do j=1,12
   call ransvu(dummy)
   x=x+dummy
enddo   

r=(x-6.0d0)/4.0d0
rsq=r*r
gauss=((((a9*rsq+a7)*rsq+a5)*rsq+a3)*rsq+a1)*r

return
END FUNCTION
