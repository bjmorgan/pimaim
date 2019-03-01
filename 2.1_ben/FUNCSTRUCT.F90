FUNCTION FUNCSTRUCT(p)

USE commondata, ONLY: engpetot,selfengtot,epsselfengtot,quaimselfengtot, &
                      reseng,dipsqeng,dipquadeng,quadeng,num,pextstruct

USE boxdata, ONLY: cellvol

IMPLICIT NONE

DOUBLE PRECISION :: FUNCSTRUCT
DOUBLE PRECISION, DIMENSION(num*3+6), INTENT(IN) :: p
INTEGER :: i

FUNCSTRUCT=engpetot+selfengtot+epsselfengtot+quaimselfengtot

FUNCSTRUCT=FUNCSTRUCT+SUM(reseng+dipsqeng+dipquadeng+quadeng)

FUNCSTRUCT=FUNCSTRUCT+pextstruct*cellvol

return
END FUNCTION
