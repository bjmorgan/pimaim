SUBROUTINE cgfirst1

USE commondata, ONLY: elecxq,elecyq,eleczq, &
                      engeff,exxq,eyyq,ezzq,exyq,exzq,eyzq, &
                      twopi,num,q,etaconst,eta
USE boxdata, ONLY: fullhi,fac,fourpicell
USE recipdata, ONLY: kmaxx,kmaxy,kmaxz,elcall,elsall,emcall,emsall,encall, &
                     ensall,rksqmax
!---> Parallelization_S
use mpipara
!<--- Parallelization_E

IMPLICIT NONE

integer :: i,l,m,n,ll,mm,nn,nmin,mmin
DOUBLE PRECISION :: rl,rkx1,rky1,rkz1,rm,rkx2,rky2, &
                 rkz2,rn,rkx3,rky3,rkz3,xkk,ckcs,ckss, &
                 akk,arl,arm,arn, &
                 arll,armm,arnn,arlm,arln,armn,egrad

! GWW - stack overflows whwn num is large so replace with allocateable 
!
! DOUBLE PRECISION, DIMENSION(num) :: clmall,slmall,ckcnoqall,cksnoqall,temp

DOUBLE PRECISION, allocatable :: clmall(:),slmall(:),ckcnoqall(:),cksnoqall(:),temp(:)

 allocate ( clmall(num),slmall(num),ckcnoqall(num),cksnoqall(num),temp(num))

!---> Parallelization_S
!do ll=0,kmaxx
 do ll=kmaxx_s,kmaxx_e
!<--- Parallelization_E
   l=iabs(ll)
   rl=dble(ll)*twopi
 
   rkx1=rl*fullhi(1,1)
   rky1=rl*fullhi(1,2)
   rkz1=rl*fullhi(1,3)

!---> Parallelization_S
!  do mm=mmin,kmaxy
   do mm=kmaxy_s(ll),kmaxy_e(ll)
!<--- Parallelization_E
      m=iabs(mm)
      rm=dble(mm)*twopi

      rkx2=rkx1+rm*fullhi(2,1)
      rky2=rky1+rm*fullhi(2,2)
      rkz2=rkz1+rm*fullhi(2,3)

      if(mm.ge.0) then
         clmall(:)=elcall(:,l)*emcall(:,m)-elsall(:,l)*emsall(:,m)
         slmall(:)=elsall(:,l)*emcall(:,m)+emsall(:,m)*elcall(:,l)
      else
         clmall(:)=elcall(:,l)*emcall(:,m)+elsall(:,l)*emsall(:,m)
         slmall(:)=elsall(:,l)*emcall(:,m)-emsall(:,m)*elcall(:,l)
      endif

!---> Parallelization_S
!     do nn=nmin,kmaxz
      do nn=kmaxz_s(mm,ll),kmaxz_e(mm,ll)
!<--- Parallelization_E
         n=iabs(nn)
         rn=dble(nn)*twopi

         rkx3=rkx2+rn*fullhi(3,1)
         rky3=rky2+rn*fullhi(3,2)
         rkz3=rkz2+rn*fullhi(3,3)

         xkk=rkx3*rkx3+rky3*rky3+rkz3*rkz3

         if(xkk.le.rksqmax) then
            if(nn.ge.0) then
               ckcnoqall(:)=clmall(:)*encall(:,n)-slmall(:)*ensall(:,n)
               cksnoqall(:)=slmall(:)*encall(:,n)+clmall(:)*ensall(:,n)
            else
               ckcnoqall(:)=clmall(:)*encall(:,n)+slmall(:)*ensall(:,n)
               cksnoqall(:)=slmall(:)*encall(:,n)-clmall(:)*ensall(:,n)
            endif

            ckcs=SUM(ckcnoqall*q)
            ckss=SUM(cksnoqall*q)

            akk=exp(etaconst*xkk)/xkk

            arl=akk*rkx3
            arm=akk*rky3
            arn=akk*rkz3
            arll=arl*rkx3
            armm=arm*rky3
            arnn=arn*rkz3
            arlm=arl*rky3
            arln=arl*rkz3
            armn=arm*rkz3
!
! force and electric field due to charges:
!
            temp=cksnoqall*ckcs-ckcnoqall*ckss

!---> Parallelization_S
!           elecxq=elecxq+arl*temp
!           elecyq=elecyq+arm*temp
!           eleczq=eleczq+arn*temp
            eltmp(1:num)=eltmp(1:num)+arl*temp
            eltmp(num+1:2*num)=eltmp(num+1:2*num)+arm*temp
            eltmp(2*num+1:3*num)=eltmp(2*num+1:3*num)+arn*temp
!<--- Parallelization_E
!
! This is the muTq term:
! Force and electric field gradient due to charges (AGUADO)
!
            temp=ckcnoqall*ckcs+cksnoqall*ckss

!---> Parallelization_S
!           exxq=exxq+arll*temp
!           eyyq=eyyq+armm*temp
!           ezzq=ezzq+arnn*temp
!           exyq=exyq+arlm*temp
!           exzq=exzq+arln*temp
!           eyzq=eyzq+armn*temp
            eltmp(3*num+1:4*num)=eltmp(3*num+1:4*num)+arll*temp
            eltmp(4*num+1:5*num)=eltmp(4*num+1:5*num)+armm*temp
            eltmp(5*num+1:6*num)=eltmp(5*num+1:6*num)+arnn*temp
            eltmp(6*num+1:7*num)=eltmp(6*num+1:7*num)+arlm*temp
            eltmp(7*num+1:8*num)=eltmp(7*num+1:8*num)+arln*temp
            eltmp(8*num+1:9*num)=eltmp(8*num+1:9*num)+armn*temp
!<--- Parallelization_E

            temp=ckcnoqall*ckcs+cksnoqall*ckss
!---> Parallelization_S
!           engeff=engeff+temp*akk*2.0d0
            eltmp(9*num+1:10*num)=eltmp(9*num+1:10*num)+temp*akk*2.0d0
!<--- Parallelization_E
         endif
      enddo   
!---> Parallelization_S
!     nmin=-kmaxz
!<--- Parallelization_E
   enddo   
!---> Parallelization_S
!  mmin=-kmaxy
!<--- Parallelization_E
enddo   

eltmp(1:9*num)=eltmp(1:9*num)*fac
eltmp(9*num+1:10*num)=eltmp(9*num+1:10*num)*fourpicell

return
END SUBROUTINE
