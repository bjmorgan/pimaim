SUBROUTINE cgrealE_dippim

USE commondata, ONLY: x,y,z, &
                      q,xmu,ymu,zmu,dxsav,dysav,dzsav, &
!---> Memmory Reduction_S
!                     num,ewlog,twopi,ntype,rsqmax,etapi,etasq,erfc,nspec
                      num,ewlog,twopi,ntype,rsqmax,etapi,etasq,erfc,nspec,num2,numx
!<--- Memmory Reduction_E
USE boxdata, ONLY: hlab2,boxlenx,boxleny,boxlenz
!---> Parallelization_S
use mpipara
!<--- Parallelization_E

IMPLICIT NONE

INTEGER :: i,j,kk,ipoint,jpoint
DOUBLE PRECISION :: dxunitsq,dyunitsq,dzunitsq,frrxdipqij,frrxdipqji, &
                    frrydipqij,frrydipqji,frrzdipqij,frrzdipqji,xmuidotrij, &
                    xmujdotrij,egam,erfcr,expewld, &
                    totmz,totmzd,totmzu,facewvac,qq2api
DOUBLE PRECISION :: drsq,drsqrec,dr,dr3rec,dr4rec,dr5rec, &
                    chgtem,efact,efact2, &
                    sxx,syy,szz,sxy,sxz,syz,qirec,qjrec,drrec
DOUBLE PRECISION :: chgdipengij,chgdipengji
DOUBLE PRECISION :: qdotmuxij,qdotmuxji,qdotmuyij,qdotmuyji, &
                    qdotmuzij,qdotmuzji
 
!---> Parallelization_S
!do j=2,num
do j=jst,jed
!<--- Parallelization_E
   jpoint=ntype(j)
   qjrec=1.0d0/q(j)
!---> Parallelization_S
!  do i=1,j-1
   do i=ist(j),ied(j)
!<--- Parallelization_E
      ipoint=ntype(i)
      qirec=1.0d0/q(i)

!---> Memmory Reduction_S
      numx = numadr(i,j)
!     drsq=dxsav(i,j)*dxsav(i,j)+dysav(i,j)*dysav(i,j)+dzsav(i,j)*dzsav(i,j)
      drsq=dxsav(numx)*dxsav(numx)+dysav(numx)*dysav(numx)+dzsav(numx)*dzsav(numx)
!<--- Memmory Reduction_E

      if(drsq.ge.rsqmax) CYCLE

      drsqrec=1.0d0/drsq
      dr=dsqrt(drsq)
      drrec=1.0d0/dr
      dr3rec=drrec*drsqrec
      dr4rec=drsqrec*drsqrec
      dr5rec=dr4rec*drrec

      chgtem=q(i)*q(j)
      qq2api=etapi*chgtem
      expewld=exp(-etasq*drsq)
!---> Memmory Reduction_S
!     erfcr=chgtem*erfc(i,j)*drrec
      erfcr=chgtem*erfc(numx)*drrec
!<--- Memmory Reduction_E
      egam=(erfcr+qq2api*expewld)*drsqrec

!---> Memmory Reduction_S
!     xmuidotrij=xmu(i)*dxsav(i,j)+ymu(i)*dysav(i,j)+zmu(i)*dzsav(i,j)
!     xmujdotrij=xmu(j)*dxsav(i,j)+ymu(j)*dysav(i,j)+zmu(j)*dzsav(i,j)
      xmuidotrij=xmu(i)*dxsav(numx)+ymu(i)*dysav(numx)+zmu(i)*dzsav(numx)
      xmujdotrij=xmu(j)*dxsav(numx)+ymu(j)*dysav(numx)+zmu(j)*dzsav(numx)
!<--- Memmory Reduction_E
 
!---> Memmory Reduction_S
!     efact=etapi*expewld*dr+erfc(i,j)
      efact=etapi*expewld*dr+erfc(numx)
!<--- Memmory Reduction_E
      chgdipengij=(q(j)*xmuidotrij)*dr3rec
      chgdipengji=-(q(i)*xmujdotrij)*dr3rec

      qdotmuxij=q(j)*xmu(i)
      qdotmuxji=-q(i)*xmu(j)
      qdotmuyij=q(j)*ymu(i)
      qdotmuyji=-q(i)*ymu(j)
      qdotmuzij=q(j)*zmu(i)
      qdotmuzji=-q(i)*zmu(j)
      efact2=2.0d0*etasq*drsqrec*etapi*expewld

!---> Memmory Reduction_S
!     frrxdipqij=(qdotmuxij*dr3rec &
!                -3.0d0*dxsav(i,j)*drsqrec*chgdipengij)*efact &
!                -dxsav(i,j)*xmuidotrij*q(j)*efact2
!     frrxdipqji=(qdotmuxji*dr3rec &
!                -3.0d0*dxsav(i,j)*drsqrec*chgdipengji)*efact &
!                +dxsav(i,j)*xmujdotrij*q(i)*efact2
!     frrydipqij=(qdotmuyij*dr3rec &
!                -3.0d0*dysav(i,j)*drsqrec*chgdipengij)*efact &
!                -dysav(i,j)*xmuidotrij*q(j)*efact2
!     frrydipqji=(qdotmuyji*dr3rec &
!                -3.0d0*dysav(i,j)*drsqrec*chgdipengji)*efact &
!                +dysav(i,j)*xmujdotrij*q(i)*efact2
!     frrzdipqij=(qdotmuzij*dr3rec &
!                -3.0d0*dzsav(i,j)*drsqrec*chgdipengij)*efact &
!                -dzsav(i,j)*xmuidotrij*q(j)*efact2
!     frrzdipqji=(qdotmuzji*dr3rec &
!                -3.0d0*dzsav(i,j)*drsqrec*chgdipengji)*efact &
!                +dzsav(i,j)*xmujdotrij*q(i)*efact2
      frrxdipqij=(qdotmuxij*dr3rec &
                 -3.0d0*dxsav(numx)*drsqrec*chgdipengij)*efact &
                 -dxsav(numx)*xmuidotrij*q(j)*efact2
      frrxdipqji=(qdotmuxji*dr3rec &
                 -3.0d0*dxsav(numx)*drsqrec*chgdipengji)*efact &
                 +dxsav(numx)*xmujdotrij*q(i)*efact2
      frrydipqij=(qdotmuyij*dr3rec &
                 -3.0d0*dysav(numx)*drsqrec*chgdipengij)*efact &
                 -dysav(numx)*xmuidotrij*q(j)*efact2
      frrydipqji=(qdotmuyji*dr3rec &
                 -3.0d0*dysav(numx)*drsqrec*chgdipengji)*efact &
                 +dysav(numx)*xmujdotrij*q(i)*efact2
      frrzdipqij=(qdotmuzij*dr3rec &
                 -3.0d0*dzsav(numx)*drsqrec*chgdipengij)*efact &
                 -dzsav(numx)*xmuidotrij*q(j)*efact2
      frrzdipqji=(qdotmuzji*dr3rec &
                 -3.0d0*dzsav(numx)*drsqrec*chgdipengji)*efact &
                 +dzsav(numx)*xmujdotrij*q(i)*efact2
!<--- Memmory Reduction_E

!---> Parallelization_S
!     elecx(j)=elecx(j)-frrxdipqij*qjrec
!     elecy(j)=elecy(j)-frrydipqij*qjrec
!     elecz(j)=elecz(j)-frrzdipqij*qjrec
!     elecx(i)=elecx(i)+frrxdipqji*qirec
!     elecy(i)=elecy(i)+frrydipqji*qirec
!     elecz(i)=elecz(i)+frrzdipqji*qirec
      eltmp(j)=eltmp(j)-frrxdipqij*qjrec
      eltmp(num+j)=eltmp(num+j)-frrydipqij*qjrec
      eltmp(2*num+j)=eltmp(2*num+j)-frrzdipqij*qjrec
      eltmp(i)=eltmp(i)+frrxdipqji*qirec
      eltmp(num+i)=eltmp(num+i)+frrydipqji*qirec
      eltmp(2*num+i)=eltmp(2*num+i)+frrzdipqji*qirec
!<--- Parallelization_E

!---> Memmory Reduction_S
!     dxunitsq=dxsav(i,j)*dxsav(i,j)*drsqrec
!     dyunitsq=dysav(i,j)*dysav(i,j)*drsqrec
!     dzunitsq=dzsav(i,j)*dzsav(i,j)*drsqrec
      dxunitsq=dxsav(numx)*dxsav(numx)*drsqrec
      dyunitsq=dysav(numx)*dysav(numx)*drsqrec
      dzunitsq=dzsav(numx)*dzsav(numx)*drsqrec
!<--- Memmory Reduction_E
      sxx=1.0d0-3.0d0*dxunitsq
      syy=1.0d0-3.0d0*dyunitsq
      szz=1.0d0-3.0d0*dzunitsq
!---> Memmory Reduction_S
!     sxy=-3.0d0*dxsav(i,j)*dysav(i,j)*drsqrec
!     sxz=-3.0d0*dxsav(i,j)*dzsav(i,j)*drsqrec
!     syz=-3.0d0*dzsav(i,j)*dysav(i,j)*drsqrec
      sxy=-3.0d0*dxsav(numx)*dysav(numx)*drsqrec
      sxz=-3.0d0*dxsav(numx)*dzsav(numx)*drsqrec
      syz=-3.0d0*dzsav(numx)*dysav(numx)*drsqrec
!<--- Memmory Reduction_E
   enddo   
enddo   

totmz=0.d0
totmzd=0.d0
facewvac=0.0d0
if(ewlog) then
facewvac=twopi/(boxlenx*boxleny*boxlenz)
totmzu=SUM(zmu)
totmzd=SUM(q*(hlab2(3,1)*x+hlab2(3,2)*y+hlab2(3,3)*z))
totmz=totmzd+totmzu
endif

return
END SUBROUTINE
