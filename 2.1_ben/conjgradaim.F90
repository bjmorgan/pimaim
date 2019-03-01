SUBROUTINE conjgradaim

USE commondata, ONLY: num,delta,epsilonx,epsilony,epsilonz,quaimxx,quaimyy, &
                      quaimzz,quaimxy,quaimxz,quaimyz

!---> Parallelization_S
use mpipara
!<--- Parallelization_E
IMPLICIT NONE

INTEGER :: i
DOUBLE PRECISION, DIMENSION(num*10) :: p

p(1:num)=delta
p(num+1:2*num)=epsilonx
p(2*num+1:3*num)=epsilony
p(3*num+1:4*num)=epsilonz
p(4*num+1:5*num)=quaimxx
p(5*num+1:6*num)=quaimyy
p(6*num+1:7*num)=quaimzz
p(7*num+1:8*num)=quaimxy
p(8*num+1:9*num)=quaimxz
p(9*num+1:10*num)=quaimyz

call frprmnaim(p,num*10)

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

!---> Parallelization_S
if( iam .eq. 0 ) then

write(*,*)
write(*,*)'**** Anion annealing completed ****'
write(*,*)

endif
!<--- Parallelization_E

return
END SUBROUTINE
