!************************************************************
!  Performs an energy minimilisation on the startup configuration
! via a conjugant gradient method.
!************************************************************

SUBROUTINE conjstruct

USE commondata, ONLY: num,x,y,z,cellconstraints,relaxcell
USE boxdata, ONLY: boxlenx,boxleny,boxlenz,h,hlab2,hlab3

!---> Parallelization_S
use mpipara
!---> Parallelization_E
IMPLICIT NONE

INTEGER :: i
DOUBLE PRECISION, DIMENSION(3*num+6) :: p,xi
DOUBLE PRECISION, DIMENSION(10) :: dcellinfo
DOUBLE PRECISION gama,beta,alpha

call dcell(h,dcellinfo)
gama  = acos(dcellinfo(4))
beta  = acos(dcellinfo(5))
alpha = acos(dcellinfo(6))

p(1:num)=x
p(num+1:2*num)=y
p(2*num+1:3*num)=z
if(relaxcell) then
 if(cellconstraints(1).GT.0) p(3*num+1)=boxlenx
 if(cellconstraints(2).GT.0) p(3*num+2)=boxleny
 if(cellconstraints(3).GT.0) p(3*num+3)=boxlenz
 if(cellconstraints(4).GT.0) p(3*num+4)=gama
 if(cellconstraints(5).GT.0) p(3*num+5)=beta
 if(cellconstraints(6).GT.0) p(3*num+6)=alpha
end if

call FRPRMNSTRUCT(p,3*num+6)

x=p(1:num)
y=p(num+1:2*num)
z=p(2*num+1:3*num)
if(relaxcell) then
 if(cellconstraints(1).GT.0) boxlenx=p(3*num+1)
 if(cellconstraints(2).GT.0) boxleny=p(3*num+2)
 if(cellconstraints(3).GT.0) boxlenz=p(3*num+3)
 if((cellconstraints(2).GT.0).AND.(cellconstraints(2).EQ.cellconstraints(1))) &
    boxleny=boxlenx
 if((cellconstraints(3).GT.0).AND.(cellconstraints(3).EQ.cellconstraints(1))) &
    boxlenz=boxlenx
 if((cellconstraints(3).GT.0).AND.(cellconstraints(3).EQ.cellconstraints(2))) &
    boxlenz=boxleny
 if(cellconstraints(4).GT.0) gama   =p(3*num+4)
 if(cellconstraints(5).GT.0) beta   =p(3*num+5)
 if(cellconstraints(6).GT.0) alpha  =p(3*num+6)
 if((cellconstraints(5).GT.0).AND.(cellconstraints(5).EQ.cellconstraints(4))) &
    beta=gama
 if((cellconstraints(6).GT.0).AND.(cellconstraints(6).EQ.cellconstraints(4))) &
    alpha=gama
 if((cellconstraints(6).GT.0).AND.(cellconstraints(6).EQ.cellconstraints(5))) &
    alpha=beta
 h(1,1)=1.0d0
 h(2,1)=0.0d0
 h(3,1)=0.0d0
 h(1,2)=cos(gama)
 h(2,2)=sin(gama)
 h(3,2)=0.0d0
 h(1,3)=cos(beta)
 h(2,3)=(cos(alpha)-cos(beta)*cos(gama))/sin(gama)
 h(3,3)=sqrt(1.0d0-h(1,3)**2-h(2,3)**2)
 hlab2=h
 hlab3=h
 call boxreset
end if

do i=1,num
   if(x(i).lt.0.0d0) x(i)=x(i)+boxlenx
   if(y(i).lt.0.0d0) y(i)=y(i)+boxleny
   if(z(i).lt.0.0d0) z(i)=z(i)+boxlenz
   if(x(i).gt.boxlenx) x(i)=x(i)-boxlenx
   if(y(i).gt.boxleny) y(i)=y(i)-boxleny
   if(z(i).gt.boxlenz) z(i)=z(i)-boxlenz
enddo   

!---> Parallelization_S
if( iam .eq. 0 ) then

do i=1,num
   write(17,*)x(i),y(i),z(i)
   write(18,*)xi(i),xi(i+num),xi(i+2*num)
enddo   

write(*,*)
write(*,*)'**** Ion annealing completed ****'
write(*,*)

endif
!---> Parallelization_E
return
END SUBROUTINE
