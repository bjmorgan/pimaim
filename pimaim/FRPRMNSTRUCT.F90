SUBROUTINE FRPRMNSTRUCT(P,N)

USE commondata, ONLY: num,x,y,z,conjgradaimlog,conjgradlog, &
                      environmentalpimlog,ftolstruc,cellconstraints, &
                      relaxcell,frrx,frry,frrz
USE boxdata, ONLY: boxlenx,boxleny,boxlenz,h,hlab2,hlab3

!---> Parallelization_S
use mpipara
!<--- Parallelization_E
IMPLICIT NONE

INTEGER :: N,ITER,i,J,ITS,icounter
INTEGER, PARAMETER :: ITMAX=1000
DOUBLE PRECISION, DIMENSION(num*3+6), INTENT(INOUT) :: P
DOUBLE PRECISION, DIMENSION(num*3+6) :: G,HHH,XI
DOUBLE PRECISION :: FP,GG,DGG,GAM,FUNCSTRUCT,FRET,alpha,beta,gama
DOUBLE PRECISION, PARAMETER :: EPS=1.E-10
DOUBLE PRECISION, DIMENSION(10):: dcellinfo
                                                                                             
call dcell(h,dcellinfo)
gama  = acos(dcellinfo(4))
beta  = acos(dcellinfo(5))
alpha = acos(dcellinfo(6))

icounter=0

call ener
FP=FUNCSTRUCT(P)
CALL DFUNCSTRUCT(P,XI)
!---> Parallelization_S
if( iam .eq. 0 ) then

   write(743,*) ' FRPRMNSTRUCT ' 
do i=1,num
   write(743,*) p(i),p(i+num),p(i+2*num)
   write(744,*) xi(i),xi(i+num),xi(i+2*num)
enddo

endif
!<--- Parallelization_E

G=-XI
HHH=G
XI=HHH

DO ITS=1,ITMAX
   ITER=ITS
   CALL LINMINSTRUCT(P,XI,FRET)

   if(iam.eq.0) write (*,*)' Linmin -', FRET 
   IF(2.*ABS(FRET-FP).LE.FTOLSTRUC*(ABS(FRET)+ABS(FP)+EPS)) THEN
      icounter=icounter+1
   ELSE
      icounter=0
   END IF
   IF(icounter.GE.3)RETURN
!---> Parallelization_S
   if( iam .eq. 0 ) then

   write(*,*) 'ITS,E,dE,tol:',ITS,real(FRET),real(FRET-FP),real(FTOLSTRUC*(ABS(FRET)+ABS(FP))/2.0d0)

     if(relaxcell) then
       write(745,*) ITS,real(boxlenx),real(boxleny),real(boxlenz)
       write(746,*) ITS,real(gama),real(beta),real(alpha)
     end if
   endif
!<--- Parallelization_E

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

   if(environmentalpimlog) then
      call conjgradpimaim
   else
      if(conjgradaimlog) call conjgradaim
      if(conjgradlog) call conjgrad
   endif
   call ener
!   FP=FUNCSTRUCT(P)
   FP=FRET
   CALL DFUNCSTRUCT(P,XI)
!---> Parallelization_S
   if( iam .eq. 0 ) then

   write(742,*) ITS,FP
   write(743,*) ITS,FP
   write(744,*) ITS,FP
   do i=1,num
      write(743,*) x(i),y(i),z(i)
      write(744,*) xi(i),xi(i+num),xi(i+2*num)
   enddo
 if(relaxcell) then
   write(745,*) ITS,real(boxlenx),real(boxleny),real(boxlenz)
   write(746,*) ITS,real(gama),real(beta),real(alpha)
 end if

   endif
!<--- Parallelization_E
   GG=SUM(G*G)
   DGG=SUM((XI+G)*XI)
   IF(GG.EQ.0.)RETURN
   GAM=DGG/GG
   G=-XI
   HHH=G+GAM*HHH
   ! steepest descent
!   if(mod(ITS,5).eq.0) HHH=G
   if(mod(ITS,20).eq.0) HHH=G
   XI=HHH
ENDDO    
!---> Parallelization_S
if( iam .eq. 0 ) then

write(*,*) 'FRPRSTRUCT maximum iterations exceeded'

endif
!<--- Parallelization_E
RETURN
END SUBROUTINE
