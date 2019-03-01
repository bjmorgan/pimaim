subroutine pramp
!************************************************************
!  pressure ramping.
!************************************************************
USE commondata, ONLY: pext, deltapramp
!---> Parallelization_S
use mpipara
!<--- Parallelization_E

implicit none

pext=pext+deltapramp
if(iam.eq.0) write(6,*)'pext now ',pext

return
end
