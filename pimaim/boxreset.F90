!***********************************************************
!
! Sets/resets all box length dependent variables.
!
!***********************************************************
SUBROUTINE boxreset

USE commondata, ONLY: twopi,fourpi,onethird,eightpi,rsqmax
USE boxdata, ONLY: boxlenx,boxleny,boxlenz,h,cellvol,halfboxx,halfboxy,  &
                   halfboxz,halfboxxrec,halfboxyrec,halfboxzrec,twopiboxx, &
                   twopiboxy,twopiboxz,fourpicell,cellvol3rec,hi,fullhi, &
                   halfboxminsq,bh,fac,fullh

IMPLICIT NONE

INTEGER :: i,j
DOUBLE PRECISION, DIMENSION(3) :: lengths
DOUBLE PRECISION, DIMENSION(10) :: dcellinfo

lengths(1)=boxlenx
lengths(2)=boxleny
lengths(3)=boxlenz

call dcell(h,dcellinfo)

cellvol=dcellinfo(10)*boxlenx*boxleny*boxlenz
halfboxx=boxlenx/2.0d0
halfboxxrec=1.0d0/halfboxx
halfboxy=boxleny/2.0d0
halfboxyrec=1.0d0/halfboxy
halfboxz=boxlenz/2.0d0
halfboxzrec=1.0d0/halfboxz
twopiboxx=twopi/boxlenx
twopiboxy=twopi/boxleny
twopiboxz=twopi/boxlenz
fourpicell=fourpi/cellvol
cellvol3rec=onethird/cellvol
!
! Now invert the cell matrix for use with the force calculations
! and the Ewald summation.
!
call invert(h,hi)
do i=1,3
   fullh(i,:)=h(i,:)*lengths(:)
enddo
call invert(fullh,fullhi)
!
! Calculate factor for reciprocal space forces.
!
fac=eightpi/cellvol


return
END SUBROUTINE
