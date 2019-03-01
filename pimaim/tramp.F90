SUBROUTINE tramp
!************************************************************
!  temperature ramping.
!************************************************************
USE commondata, ONLY: num,trantemp,trantkb,boltz,CUEp,CUEp2,  &
          CUEprec,CUEp2rec,CUEb,CUEbrec,CUEb2,CUEb2rec,  &
          W,Wgo,Wrec,Wgorec,nib,nth,nab,free,relax,relaxb,relaxb2, &
          dom,deltatramp
!---> Parallelization_S
use mpipara
!<--- Parallelization_E

IMPLICIT NONE

DOUBLE PRECISION :: relaxsq,relaxsqrec,gtran

trantemp=trantemp+deltatramp
if(iam.eq.0) write(*,*) 'T ramp: temperature set to ',trantemp,'K.'
trantkb=boltz*trantemp
!
! Thermostat parameters.
!
relaxsq=relax*relax
if(nth)relaxsqrec=1.0d0/relaxsq
gtran=3.0d0*dble(num-1)*trantkb

W=(free+3.0d0)*trantkb*(relaxb**2.0d0)
Wgo=(free+3.0d0)*trantkb*(relaxb2**2.0d0)/3.0d0
if(nib)then
   Wrec=1.0d0/W
   Wgorec=1.0d0/Wgo
endif
if(.not.nab) then
   Wgo=0.d0
   Wgorec=0.d0
endif

if(nth)CUEp=gtran/relaxsqrec
if(nth)CUEp2=CUEp/free
if(nth)CUEprec=1.d0/CUEp
if(nth)CUEp2rec=1.d0/CUEp2
if(nib)CUEb=dom*trantkb/relaxsqrec
if(nib)CUEb2=CUEb/dom
if(nib)CUEbrec=1.0d0/CUEb
if(nib)CUEb2rec=1.0d0/CUEb2
if(.not.nth) then
   CUEbrec=0.d0
   CUEb2rec=0.d0
endif

return
END SUBROUTINE
