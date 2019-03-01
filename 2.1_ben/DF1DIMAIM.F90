FUNCTION DF1DIMAIM(XXX)

USE commondata, ONLY: num,delta,epsilonx,epsilony,epsilonz,quaimxx,quaimyy, &
                      quaimzz,quaimxy,quaimxz,quaimyz
USE cgdata, ONLY: PCOMAIM,XICOMAIM

IMPLICIT NONE

INTEGER :: i,J
double precision, INTENT(IN) :: XXX
double precision :: DF1DIMAIM
double precision, DIMENSION(num*10) :: XT,DF

XT=PCOMAIM+XXX*XICOMAIM

!sgi
df=0.0d0
!sgi
CALL DFUNCAIM(XT,DF)
DF1DIMAIM=SUM(DF*XICOMAIM)
RETURN
END FUNCTION
