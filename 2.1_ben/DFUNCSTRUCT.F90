SUBROUTINE DFUNCSTRUCT(p,xi)

USE commondata, only: num,frrx,frry,frrz,pint2,stpxx,stpyy,stpzz,stcxx, &
                      stcyy,stczz,stsrxx,stsryy,stsrzz,stp2xx,stp2yy, &
                      stp2zz,stpsrxx,stpsryy,stpsrzz,stqsrxx,stqsryy, &
                      stqsrzz,stpqquadxx,stpqquadyy,stpqquadzz, &
                      stpdipquadxx,stpdipquadyy,stpdipquadzz, &
                      stpquadquadxx,stpquadquadyy,stpquadquadzz,stewzz, &
                      pextstruct, &
                      cellconstraints,stpxy,stpxz,stpyz,stcxy, &
                      stcxz,stcyz,stsrxy,stsrxz,stsryz,stp2xy,stp2xz, &
                      stp2yz,stpsrxy,stpsrxz,stpsryz,stqsrxy,stqsrxz, &
                      stqsryz,stpqquadxy,stpqquadxz,stpqquadyz, &
                      stpdipquadxy,stpdipquadxz,stpdipquadyz, &
                      stpquadquadxy,stpquadquadxz,stpquadquadyz
USE boxdata, ONLY: cellvol3rec

IMPLICIT NONE

INTEGER :: i
DOUBLE PRECISION, DIMENSION(num*3+6), INTENT(IN) :: p
DOUBLE PRECISION, DIMENSION(num*3+6), INTENT(OUT) :: xi
      
pint2(1,1)=3.0d0*(stpxx+ &
                 stcxx+stsrxx+stp2xx+stpsrxx+stqsrxx+ &
                 stpqquadxx+stpdipquadxx+stpquadquadxx) &
                 *cellvol3rec
pint2(2,2)=3.0d0*(stpyy+ &
                 stcyy+stsryy+stp2yy+stpsryy+stqsryy+ &
                 stpqquadyy+stpdipquadyy+stpquadquadyy) &
                 *cellvol3rec
pint2(3,3)=3.0d0*(stpzz+ &
                 stczz+stsrzz+stp2zz+stpsrzz+stqsrzz+ &
                 stpqquadzz+stpdipquadzz+stpquadquadzz) &
                 *cellvol3rec+stewzz
pint2(1,2)=3.0d0*(stpxy+ &
                 stcxy+stsrxy+stp2xy+stpsrxy+stqsrxy+ &
                 stpqquadxy+stpdipquadxy+stpquadquadxy) &
                 *cellvol3rec
pint2(1,3)=3.0d0*(stpxz+ &
                 stcxz+stsrxz+stp2xz+stpsrxz+stqsrxz+ &
                 stpqquadxz+stpdipquadxz+stpquadquadxz) &
                 *cellvol3rec
pint2(2,3)=3.0d0*(stpyz+ &
                 stcyz+stsryz+stp2yz+stpsryz+stqsryz+ &
                 stpqquadxz+stpdipquadxz+stpquadquadxz) &
                 *cellvol3rec

xi(1:num)=-frrx
xi(num+1:2*num)=-frry
xi(2*num+1:3*num)=-frrz
xi(3*num+1:3*num+6)=0.0d0
if(cellconstraints(1).GT.0) xi(3*num+1)=-(pint2(1,1)-pextstruct)
if(cellconstraints(2).GT.0) then 
   if(cellconstraints(2).EQ.cellconstraints(1)) then
      xi(3*num+1)=-((pint2(1,1)+pint2(2,2))/2.0d0-pextstruct)
   else
      xi(3*num+2)=-(pint2(2,2)-pextstruct)
   end if
end if
if(cellconstraints(3).GT.0) then
   if(cellconstraints(3).EQ.cellconstraints(1)) then
      if(cellconstraints(3).EQ.cellconstraints(2)) then
         xi(3*num+1)=-((pint2(1,1)+pint2(2,2)+pint2(3,3))/3.0d0-pextstruct)
      else
         xi(3*num+1)=-((pint2(1,1)+pint2(3,3))/2.0d0-pextstruct)
      end if
   else 
      if(cellconstraints(3).EQ.cellconstraints(2)) then
         xi(3*num+2)=-((pint2(2,2)+pint2(3,3))/2.0d0-pextstruct)
      else
         xi(3*num+3)=-(pint2(3,3)-pextstruct)
      end if
   end if
end if
if(cellconstraints(4).GT.0) xi(3*num+4)=-pint2(1,2)
if(cellconstraints(5).GT.0) then
   if(cellconstraints(5).EQ.cellconstraints(4)) then
      xi(3*num+4)=-(pint2(1,2)+pint2(1,3))/2.0d0
   else
      xi(3*num+5)=-pint2(1,3)
   end if
end if
if(cellconstraints(6).GT.0) then
   if(cellconstraints(6).EQ.cellconstraints(4)) then
      if(cellconstraints(6).EQ.cellconstraints(5)) then
         xi(3*num+4)=-(pint2(1,2)+pint2(1,3)+pint2(2,3))/3.0d0
      else
         xi(3*num+4)=-(pint2(1,2)+pint2(2,3))/2.0d0
      end if
   else
      if(cellconstraints(6).EQ.cellconstraints(5)) then
         xi(3*num+5)=-(pint2(1,3)+pint2(2,3))/2.0d0
      else
         xi(3*num+6)=-pint2(2,3)
      end if
   end if
end if

return
END SUBROUTINE
