SUBROUTINE separations

use error_function 
!---> Memmory Reduction_S
!---> Parallelization_S
!USE commondata, ONLY: num,x,y,z,dxsav,dysav,dzsav
USE commondata, ONLY: num,x,y,z,dxsav,dysav,dzsav,num2,numx, &
                      cimlog,daimlog,quaimlog,ooaimlog,epplog, &
                      erfc,eta
!<--- Parallelization_E
!<--- Memmory Reduction_E
USE boxdata, ONLY: halfboxx,halfboxy,halfboxz,boxlenx,boxleny, &
                   boxlenz,hlab2

!---> Parallelization_S
use mpipara
!<--- Parallelization_E

IMPLICIT NONE

INTEGER :: i,j
DOUBLE PRECISION :: dxcf, dycf, dzcf, dr, drsq

!---> Parallelization_S
!do j=2,num
!   do i=1,j-1
do j=jst,jed
   do i=ist(j),ied(j)
!<--- Parallelization_E

! GWW x,y,z distance along vectors - therefore can add boxlen 
      dxcf=x(i)-x(j)
      dycf=y(i)-y(j)
      dzcf=z(i)-z(j)

! GWW therefore can add / subtract boxlen 
      if(abs(dxcf).gt.halfboxx) dxcf=dxcf-sign(boxlenx,dxcf)
      if(abs(dycf).gt.halfboxy) dycf=dycf-sign(boxleny,dycf)
      if(abs(dzcf).gt.halfboxz) dzcf=dzcf-sign(boxlenz,dzcf)

      numx = numadr(i,j)

! GWW convert to cartisain by mnultiplying by strain distortion matrix (det(h)=1) 
      dxsav(numx)=hlab2(1,1)*dxcf+hlab2(1,2)*dycf+hlab2(1,3)*dzcf
      dysav(numx)=hlab2(2,1)*dxcf+hlab2(2,2)*dycf+hlab2(2,3)*dzcf
      dzsav(numx)=hlab2(3,1)*dxcf+hlab2(3,2)*dycf+hlab2(3,3)*dzcf

! DT Adding erfc array evaluation here because it's notably ABSENT in RIM Ewald sum
      drsq=dxsav(numx)*dxsav(numx)+dysav(numx)*dysav(numx)+dzsav(numx)*dzsav(numx)
      dr=dsqrt(drsq)
      erfc(numx)=erfunc(eta*dr)

      write(911,'(i10,4(X,F12.6))') numx,dxsav(numx),dysav(numx),dzsav(numx), dr 

   enddo   
enddo    

!---> Parallelization_S
if( nprocs .ne. 1 ) then
   if( cimlog .or. daimlog .or. quaimlog ) then
      do j=jst2,jed2
         do i=ist2(j),ied2(j)
            numx = numadr(i,j)
            dxcf=x(i)-x(j)
            dycf=y(i)-y(j)
            dzcf=z(i)-z(j)
            if(abs(dxcf).gt.halfboxx) dxcf=dxcf-sign(boxlenx,dxcf)
            if(abs(dycf).gt.halfboxy) dycf=dycf-sign(boxleny,dycf)
            if(abs(dzcf).gt.halfboxz) dzcf=dzcf-sign(boxlenz,dzcf)
            dxsav(numx)=hlab2(1,1)*dxcf+hlab2(1,2)*dycf+hlab2(1,3)*dzcf
            dysav(numx)=hlab2(2,1)*dxcf+hlab2(2,2)*dycf+hlab2(2,3)*dzcf
            dzsav(numx)=hlab2(3,1)*dxcf+hlab2(3,2)*dycf+hlab2(3,3)*dzcf
         enddo
      enddo
   endif
   if( ooaimlog ) then
      if( cimlog .or. daimlog .or. quaimlog ) then
         do j=jst3,jed3
            do i=ist3(j),ied3(j)
               numx = numadr(i,j)
               dxcf=x(i)-x(j)
               dycf=y(i)-y(j)
               dzcf=z(i)-z(j)
               if(abs(dxcf).gt.halfboxx) dxcf=dxcf-sign(boxlenx,dxcf)
               if(abs(dycf).gt.halfboxy) dycf=dycf-sign(boxleny,dycf)
               if(abs(dzcf).gt.halfboxz) dzcf=dzcf-sign(boxlenz,dzcf)
               dxsav(numx)=hlab2(1,1)*dxcf+hlab2(1,2)*dycf+hlab2(1,3)*dzcf
               dysav(numx)=hlab2(2,1)*dxcf+hlab2(2,2)*dycf+hlab2(2,3)*dzcf
               dzsav(numx)=hlab2(3,1)*dxcf+hlab2(3,2)*dycf+hlab2(3,3)*dzcf
            enddo
         enddo
      endif
   elseif( .not. epplog ) then
      do j=jst3,jed3
         do i=ist3(j),ied3(j)
            numx = numadr(i,j)
            dxcf=x(i)-x(j)
            dycf=y(i)-y(j)
            dzcf=z(i)-z(j)
            if(abs(dxcf).gt.halfboxx) dxcf=dxcf-sign(boxlenx,dxcf)
            if(abs(dycf).gt.halfboxy) dycf=dycf-sign(boxleny,dycf)
            if(abs(dzcf).gt.halfboxz) dzcf=dzcf-sign(boxlenz,dzcf)
            dxsav(numx)=hlab2(1,1)*dxcf+hlab2(1,2)*dycf+hlab2(1,3)*dzcf
            dysav(numx)=hlab2(2,1)*dxcf+hlab2(2,2)*dycf+hlab2(2,3)*dzcf
            dzsav(numx)=hlab2(3,1)*dxcf+hlab2(3,2)*dycf+hlab2(3,3)*dzcf
         enddo
      enddo
   endif
   if( .not. epplog ) then
      do j=jst4,jed4
         do i=ist4(j),ied4(j)
            numx = numadr(i,j)
            dxcf=x(i)-x(j)
            dycf=y(i)-y(j)
            dzcf=z(i)-z(j)
            if(abs(dxcf).gt.halfboxx) dxcf=dxcf-sign(boxlenx,dxcf)
            if(abs(dycf).gt.halfboxy) dycf=dycf-sign(boxleny,dycf)
            if(abs(dzcf).gt.halfboxz) dzcf=dzcf-sign(boxlenz,dzcf)
            dxsav(numx)=hlab2(1,1)*dxcf+hlab2(1,2)*dycf+hlab2(1,3)*dzcf
            dysav(numx)=hlab2(2,1)*dxcf+hlab2(2,2)*dycf+hlab2(2,3)*dzcf
            dzsav(numx)=hlab2(3,1)*dxcf+hlab2(3,2)*dycf+hlab2(3,3)*dzcf
         enddo
      enddo
   endif
endif
!<--- Parallelization_E

return
END SUBROUTINE
