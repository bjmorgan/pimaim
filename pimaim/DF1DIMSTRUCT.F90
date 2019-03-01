FUNCTION DF1DIMSTRUCT(XXX)

USE commondata, ONLY: num,x,y,z,cellconstraints,relaxcell, &
                      conjgradaimlog,conjgradlog,environmentalpimlog
USE boxdata, ONLY: boxlenx,boxleny,boxlenz,h,hlab2,hlab3
USE cgdata, ONLY: PCOMSTRUCT,XICOMSTRUCT

IMPLICIT NONE

INTEGER :: i,J
DOUBLE PRECISION, INTENT(IN) :: XXX
DOUBLE PRECISION :: DF1DIMSTRUCT,alpha,beta,gama
DOUBLE PRECISION, DIMENSION(num*3+6) :: XT,DF

XT=PCOMSTRUCT+XXX*XICOMSTRUCT
CALL DFUNCSTRUCT(XT,DF)
DF1DIMSTRUCT=SUM(DF*XICOMSTRUCT)
RETURN
END FUNCTION
