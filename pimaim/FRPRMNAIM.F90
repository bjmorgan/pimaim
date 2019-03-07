SUBROUTINE FRPRMNAIM(P,N)

USE commondata, ONLY: num,engpetot,environmentalpimlog,ftolaim,delta, &
                      epsilonx,epsilony,epsilonz,quaimxx,quaimyy,quaimzz, &
                      quaimxy,quaimxz,quaimyz,verbose

!---> Parallelization_S
use mpipara
!<--- Parallelization_E
IMPLICIT NONE

INTEGER :: N,ITER,i,J,ITS
INTEGER, PARAMETER :: ITMAX=200
DOUBLE PRECISION, DIMENSION(N), INTENT(INOUT) :: P
DOUBLE PRECISION, DIMENSION(N) :: G,HHH,XI
DOUBLE PRECISION :: FP,GG,DGG,GAM,FUNCAIM,FRET
DOUBLE PRECISION, PARAMETER :: EPS=1.E-10

engpetot=0.d0

if(.not.environmentalpimlog) call separations
!---> Parallelization_S
!!engtmp=0.0d0
!<--- Parallelization_E
call cgsr_energy

FP=FUNCAIM(P)

!---> Parallelization_S
if( iam .eq. 0 ) then

if(verbose) write(6,*)'iteration 0',FP

endif
!<--- Parallelization_E
!sgi
xi=0.0d0
!sgi
CALL DFUNCAIM(P,XI)

G=-XI
HHH=G
XI=HHH

DO ITS=1,ITMAX
   ITER=ITS
   CALL LINMINAIM(P,XI,FRET)

   IF(2.*ABS(FRET-FP).LE.FTOLAIM*(ABS(FRET)+ABS(FP)+EPS))RETURN

   delta=p(1:num)
   epsilonx=p(num+1:2*num)
   epsilony=p(2*num+1:3*num)
   epsilonz=p(3*num+1:4*num)
   quaimxx=p(4*num+1:5*num)
   quaimyy=p(5*num+1:6*num)
   quaimzz=p(6*num+1:7*num)
   quaimxy=p(7*num+1:8*num)
   quaimxz=p(8*num+1:9*num)
   quaimyz=p(9*num+1:10*num)

   engpetot=0.d0
!---> Parallelization_S
!! engtmp=0.0d0
!<--- Parallelization_E
   call cgsr_energy
   FP=FUNCAIM(P)
!---> Parallelization_S
   if( iam .eq. 0 ) then

   if(verbose) write(6,*)'iteration',ITS,FP

   endif
!<--- Parallelization_E
   CALL DFUNCAIM(P,XI)
   GG=SUM(G*G)
   DGG=SUM(XI*(G+XI))
   IF(GG.EQ.0.)RETURN
   GAM=DGG/GG
   G=-XI
   HHH=G+GAM*HHH
   XI=HHH
ENDDO    
RETURN
END SUBROUTINE
