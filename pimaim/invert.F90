!************************************************************
!
!  Amazingly unsubtle matrix inversion algorithm
! removing the need for nag f01abf in the inversion
! of the short-range force matrix.
!
!  h(3,3) is the cell vector matrix.
!  hi(3,3) is the inverse of the above.
!  d is the determinant of h().
!
!  Modified from CCP5 code.
!
!************************************************************

SUBROUTINE invert(hloc,hiloc)

IMPLICIT NONE
DOUBLE PRECISION :: d,r
DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN) :: hloc
DOUBLE PRECISION, DIMENSION(3,3), INTENT(OUT) :: hiloc

!
! Calculate adjoint matrix.
!
hiloc(1,1)=hloc(2,2)*hloc(3,3)-hloc(3,2)*hloc(2,3)
hiloc(2,1)=hloc(3,1)*hloc(2,3)-hloc(2,1)*hloc(3,3)
hiloc(3,1)=hloc(2,1)*hloc(3,2)-hloc(3,1)*hloc(2,2)
hiloc(1,2)=hloc(3,2)*hloc(1,3)-hloc(1,2)*hloc(3,3)
hiloc(2,2)=hloc(1,1)*hloc(3,3)-hloc(3,1)*hloc(1,3)
hiloc(3,2)=hloc(3,1)*hloc(1,2)-hloc(1,1)*hloc(3,2)
hiloc(1,3)=hloc(1,2)*hloc(2,3)-hloc(2,2)*hloc(1,3)
hiloc(2,3)=hloc(2,1)*hloc(1,3)-hloc(1,1)*hloc(2,3)
hiloc(3,3)=hloc(1,1)*hloc(2,2)-hloc(2,1)*hloc(1,2)
!
! Calculate determinant.
!
d=hloc(1,1)*hiloc(1,1)+hloc(1,2)*hiloc(2,1)+hloc(1,3)*hiloc(3,1)
r=0.d0
if(abs(d).gt.0.d0)r=1.d0/d
!
! Complete inverse matrix.
!
hiloc=r*hiloc

return
END SUBROUTINE
