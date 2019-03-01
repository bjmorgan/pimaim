FUNCTION F1DIMSTRUCT(XXX)

USE commondata, ONLY: num,x,y,z,conjgradlog,conjgradaimlog,environmentalpimlog, &
                      cellconstraints,relaxcell,frrx,frry,frrz
USE boxdata, ONLY: boxlenx,boxleny,boxlenz,h,hlab2,hlab3
USE cgdata, ONLY: PCOMSTRUCT,XICOMSTRUCT

IMPLICIT NONE

INTEGER :: i,J
DOUBLE PRECISION, INTENT(IN) :: XXX
DOUBLE PRECISION :: F1DIMSTRUCT
DOUBLE PRECISION :: FUNCSTRUCT,alpha,beta,gama
DOUBLE PRECISION, DIMENSION(num*3+6) :: XT
DOUBLE PRECISION, DIMENSION(10):: dcellinfo

call dcell(h,dcellinfo)
gama  = acos(dcellinfo(4))
beta  = acos(dcellinfo(5))
alpha = acos(dcellinfo(6))

XT=PCOMSTRUCT+XXX*XICOMSTRUCT

x=xt(1:num)
y=xt(num+1:2*num)
z=xt(2*num+1:3*num)
if(relaxcell) then
 if(cellconstraints(1).GT.0) boxlenx=xt(3*num+1)
 if(cellconstraints(2).GT.0) boxleny=xt(3*num+2)
 if(cellconstraints(3).GT.0) boxlenz=xt(3*num+3)
 if((cellconstraints(2).GT.0).AND.(cellconstraints(2).EQ.cellconstraints(1))) &
    boxleny=boxlenx
 if((cellconstraints(3).GT.0).AND.(cellconstraints(3).EQ.cellconstraints(1))) &
    boxlenz=boxlenx
 if((cellconstraints(3).GT.0).AND.(cellconstraints(3).EQ.cellconstraints(2))) &
    boxlenz=boxleny
 if(cellconstraints(4).GT.0) gama   =xt(3*num+4)
 if(cellconstraints(5).GT.0) beta   =xt(3*num+5)
 if(cellconstraints(6).GT.0) alpha  =xt(3*num+6)
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

if(environmentalpimlog) then
   call conjgradpimaim
else
   if(conjgradaimlog) call conjgradaim
   if(conjgradlog) call conjgrad
endif
call ener

F1DIMSTRUCT=FUNCSTRUCT(XT)
RETURN
END FUNCTION
