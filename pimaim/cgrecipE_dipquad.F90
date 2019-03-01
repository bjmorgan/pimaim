SUBROUTINE cgrecipE_dipquad

USE commondata, ONLY: xmu,ymu,zmu,quadxx,quadyy,quadzz,quadxy,quadxz, &
                      quadyz,twopi,num,q,etaconst,onethird
USE boxdata, ONLY: fullhi,fac
USE recipdata, ONLY: kmaxx,kmaxy,kmaxz,elcall,elsall,emcall,emsall,encall, &
                     ensall,rksqmax
!---> Parallelization_S
use mpipara
!<--- Parallelization_E

IMPLICIT NONE

integer :: i,l,m,n,ll,mm,nn,nmin,mmin
DOUBLE PRECISION :: rl,rkx1,rky1,rkz1,rm,rkx2,rky2, &
                 rkz2,rn,rkx3,rky3,rkz3,xkk,ckcs,ckss,cx,cy,cz,sx,sy,sz,cxx, &
                 cyy,czz,sxx,syy,szz,cxy,cxz,cyz,sxy,sxz,syz,akk,arl,arm,arn, &
                 arll,armm,arnn,arlm,arln,armn,arijdotsij,arijdotcij,clx,cly, &
                 clz,slx,sly,slz,rdotmus,rdotmuc,egrad
DOUBLE PRECISION, DIMENSION(num) :: clmall,slmall,ckcnoqall,cksnoqall,temp

!---> Parallelization_S
!mmin=0
!nmin=1
eltmp=0.d0
!<--- Parallelization_E

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

            cx=SUM(ckcnoqall*xmu)
            cy=SUM(ckcnoqall*ymu)
            cz=SUM(ckcnoqall*zmu)

            sx=SUM(cksnoqall*xmu)
            sy=SUM(cksnoqall*ymu)
            sz=SUM(cksnoqall*zmu)

            cxx=SUM(ckcnoqall*quadxx)
            cyy=SUM(ckcnoqall*quadyy)
            czz=SUM(ckcnoqall*quadzz)
            cxy=SUM(ckcnoqall*quadxy)
            cxz=SUM(ckcnoqall*quadxz)
            cyz=SUM(ckcnoqall*quadyz)

            sxx=SUM(cksnoqall*quadxx)
            syy=SUM(cksnoqall*quadyy)
            szz=SUM(cksnoqall*quadzz)
            sxy=SUM(cksnoqall*quadxy)
            sxz=SUM(cksnoqall*quadxz)
            syz=SUM(cksnoqall*quadyz)

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

            arijdotsij=arll*sxx+armm*syy+arnn*szz+ &
                       2.0d0*(arlm*sxy+arln*sxz+armn*syz)
            arijdotcij=arll*cxx+armm*cyy+arnn*czz+ &
                       2.0d0*(arlm*cxy+arln*cxz+armn*cyz)

            clx=cx*rkx3
            cly=cy*rky3
            clz=cz*rkz3
            slx=sx*rkx3
            sly=sy*rky3
            slz=sz*rkz3

            rdotmuc=clx+cly+clz
            rdotmus=slx+sly+slz
!
! this is the qTmu term:
! force and electric field due to dipoles:
!
            temp=-ckcnoqall*rdotmuc-cksnoqall*rdotmus

!---> Parallelization_S
!           elecx=elecx+arl*temp
!           elecy=elecy+arm*temp
!           elecz=elecz+arn*temp
            eltmp(1:num)=eltmp(1:num)+arl*temp
            eltmp(num+1:2*num)=eltmp(num+1:2*num)+arm*temp
            eltmp(2*num+1:3*num)=eltmp(2*num+1:3*num)+arn*temp
!<--- Parallelization_E
!
! This is the muTmu term
! Force and electric field gradient due to dipoles (AGUADO)
!
            temp=cksnoqall*rdotmuc-ckcnoqall*rdotmus

!---> Parallelization_S
!           exx=exx+arll*temp
!           eyy=eyy+armm*temp
!           ezz=ezz+arnn*temp
!           exy=exy+arlm*temp
!           exz=exz+arln*temp
!           eyz=eyz+armn*temp
            eltmp(3*num+1:4*num)=eltmp(3*num+1:4*num)+arll*temp
            eltmp(4*num+1:5*num)=eltmp(4*num+1:5*num)+armm*temp
            eltmp(5*num+1:6*num)=eltmp(5*num+1:6*num)+arnn*temp
            eltmp(6*num+1:7*num)=eltmp(6*num+1:7*num)+arlm*temp
            eltmp(7*num+1:8*num)=eltmp(7*num+1:8*num)+arln*temp
            eltmp(8*num+1:9*num)=eltmp(8*num+1:9*num)+armn*temp
!<--- Parallelization_E
!
! This is the qTquad term:
! Force and electric field due to quadrupoles: (AGUADO)
!
            temp=onethird*(ckcnoqall*arijdotsij-cksnoqall*arijdotcij)

!---> Parallelization_S
!           elecx=elecx+rkx3*temp
!           elecy=elecy+rky3*temp
!           elecz=elecz+rkz3*temp
            eltmp(1:num)=eltmp(1:num)+rkx3*temp
            eltmp(num+1:2*num)=eltmp(num+1:2*num)+rky3*temp
            eltmp(2*num+1:3*num)=eltmp(2*num+1:3*num)+rkz3*temp
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

eltmp=eltmp*fac

return
END SUBROUTINE
