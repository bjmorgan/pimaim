FUNCTION F1DIMAIM(XXX)

USE commondata, ONLY: engpetot,delta,epsilonx,epsilony,epsilonz,quaimxx, &
                      quaimyy,quaimzz,quaimxy,quaimxz,quaimyz,num
USE cgdata, ONLY: PCOMAIM,XICOMAIM
!---> Parallelization_S
!!use mpipara
!<--- Parallelization_E
IMPLICIT NONE

INTEGER :: i,J
DOUBLE PRECISION, INTENT(IN) :: XXX
DOUBLE PRECISION :: F1DIMAIM
DOUBLE PRECISION :: FUNCAIM
DOUBLE PRECISION, DIMENSION(num*10) :: XT

XT=PCOMAIM+XXX*XICOMAIM

delta=xt(1:num)
epsilonx=xt(num+1:2*num)
epsilony=xt(2*num+1:3*num)
epsilonz=xt(3*num+1:4*num)
quaimxx=xt(4*num+1:5*num)
quaimyy=xt(5*num+1:6*num)
quaimzz=xt(6*num+1:7*num)
quaimxy=xt(7*num+1:8*num)
quaimxz=xt(8*num+1:9*num)
quaimyz=xt(9*num+1:10*num)

engpetot=0.d0 
!---> Parallelization_S
!!engtmp=0.0d0
!<--- Parallelization_E
call cgsr_energy

F1DIMAIM=FUNCAIM(XT)
RETURN
END FUNCTION
