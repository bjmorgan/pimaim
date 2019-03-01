!************************************************************
!   From `Numerical Recipes'.
!   Produces a uniform random number in the 0-1 range.
!   Called with idum = -1 to seed or re-seed the sequence.
!************************************************************

FUNCTION ran1(idum)

IMPLICIT NONE

DOUBLE PRECISION :: ran1
DOUBLE PRECISION, DIMENSION(100) :: r
INTEGER, INTENT(IN) :: idum
INTEGER :: ix1,ix2,ix3,j
INTEGER, PARAMETER :: m1=259200,ia1=7141,ic1=54773
INTEGER, PARAMETER :: m2=134456,ia2=8121,ic2=28411
INTEGER, PARAMETER :: m3=243000,ia3=4561,ic3=51349
DOUBLE PRECISION, PARAMETER :: rm1=1.0d0/m1,rm2=1.0d0/m2
SAVE r,ix1,ix2,ix3

if(idum.lt.0.) then
   ix1=mod(ic1-idum,m1)
   ix1=mod(ia1*ix1+ic1,m1)
   ix2=mod(ix1,m2)
   ix1=mod(ia1*ix1+ic1,m1)
   ix3=mod(ix1,m3)
   do j=1,97
      ix1=mod(ia1*ix1+ic1,m1)
      ix2=mod(ia2*ix2+ic2,m2)
      r(j)=(dble(ix1)+dble(ix2)*rm2)*rm1
   enddo   
endif
 
ix1=mod(ia1*ix1+ic1,m1)
ix2=mod(ia2*ix2+ic2,m2)
ix3=mod(ia3*ix3+ic3,m3)

j=1+(97*ix3)/m3
r(j)=(dble(ix1)+dble(ix2)*rm2)*rm1
ran1=r(j)

return
END FUNCTION
