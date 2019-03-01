subroutine shortrange

!************************************************************
!
!  Adjusts the short range tensor and inverts the matrix via
! 
!
!************************************************************

USE lightdata, only: srxx,sryy,srzz,srxy,srxz,sryz,polarundist, &
                     xlipol
USE commondata, only: nanion
implicit none


! *** define our local variables
integer :: i,j,k,l,m
! *** end of local definitions

!=====================================================================
! Local double precision arrays.
double precision, dimension(3,3) :: inptarr,outarr
double precision pr

do i=1,3
   do j=1,3
      inptarr(i,j)=0.0d0
      outarr(i,j)=0.0d0
   enddo
enddo

!write(6,*) polarundist,xlipol(1)

!=====================================================================

do i=1,nanion

!=====================================================================
! Add the undistorted polarizability to the trace terms.
   inptarr(1,1)=srxx(i) + polarundist
   inptarr(1,2)=srxy(i)
   inptarr(1,3)=srxz(i)
   inptarr(2,2)=sryy(i) + polarundist
   inptarr(2,3)=sryz(i)
   inptarr(3,3)=srzz(i) + polarundist
   inptarr(2,1)=inptarr(1,2)
   inptarr(3,1)=inptarr(1,3)
   inptarr(3,2)=inptarr(2,3)

!write(10,*) i

!write(10,*) ((inptarr(l,m),m=1,3),l=1,3)


!=====================================================================
! Invert the matrix inptarr putting the result in outarr.
  
   call invert(inptarr,outarr)

!write(10,*) ' '
!=====================================================================
! Subtract the average polarizability (polarav) from the trace terms.

!write(10,*) ((outarr(l,m),m=1,3),l=1,3)


   srxx(i)=outarr(1,1)-xlipol(1)
   srxy(i)=outarr(1,2)
   srxz(i)=outarr(1,3)
   sryy(i)=outarr(2,2)-xlipol(1)
   sryz(i)=outarr(2,3)
   srzz(i)=outarr(3,3)-xlipol(1)



enddo

return

end subroutine
