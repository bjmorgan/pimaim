SUBROUTINE cgrealE_clus_quadpim

USE commondata, ONLY: elecxsr,elecysr,eleczsr,exxsr,eyysr,ezzsr,exysr,exzsr, &
                      eyzsr,elecx,elecy,elecz,exx,eyy,ezz,exy,exz,eyz,x,y,z, &
                      q,xmu,ymu,zmu,quadxx,quadyy,quadzz,quadxy,quadxz, &
                      quadyz,num,ntype,rsqmax,etapi,etasq,erfc, &
                      engeff,onethird,dxsav,dysav,dzsav,nspec, &
!---> Memmory Reduction_S
!                     nkfg,nkdamp,dampa,dampfac,fgb,fgc
                      nkfg,nkdamp,dampa,dampfac,fgb,fgc,num2,numx
!<--- Memmory Reduction_E
USE boxdata, ONLY: hlab2,boxlenx,boxleny,boxlenz
!---> Memmory Reduction_S
use mpipara
!<--- Memmory Reduction_E

IMPLICIT NONE

INTEGER :: i,j,kk,ipoint,jpoint
DOUBLE PRECISION :: dxunitsq,dyunitsq,dzunitsq,frrxdipqij,frrxdipqji, &
                    frrydipqij,frrydipqji,frrzdipqij,frrzdipqji,xmuidotrij, &
                    xmujdotrij,ettx,etty,ettz,egam,erfcr,expewld,eng, &
                    qq2api
DOUBLE PRECISION :: efact,efact2,efact3,efact4, &
                    chgquad2ij,chgquad2ji, &
                    chgquadengij,chgquadengji, &
                    chgdipengij,chgdipengji,chgtem 
DOUBLE PRECISION :: drsq,drsqrec,dr,dr3rec,dr4rec,dr5rec, &
                    dr7rec,dr9rec,dx,dy,dz,dxx,dyy,dzz,dxy,dxz,dyz, &
                    sxx,syy,szz,sxy,sxz,syz,threedr7,qirec,qjrec,drrec, &
                    fqquadxij,fqquadxji,fqquadyij,fqquadyji,fqquadzij,fqquadzji
DOUBLE PRECISION :: txxt,tyyt,tzzt,txyt,txzt,tyzt, &
                    elecsrxi,elecsryi,elecsrzi,elecsrxj,elecsryj,elecsrzj, &
                    r3dampi,r3dampj,factorial,xf, &
                    dampsumfi,dampsumfj,dampfunci,dampfuncj,dampaexpi,dampaexpj
DOUBLE PRECISION :: qdotmuxij,qdotmuxji,qdotmuyij,qdotmuyji, &
                    qdotmuzij,qdotmuzji
DOUBLE PRECISION :: txxxx,tyyyy,tzzzz,txxxy,txxxz,txxyy,txxzz, &
                    txyyy,txzzz,tyyyz,tyzzz,tyyzz,txxyz,txyyz,txyzz,txxx,tyyy, &
                 tzzz,txxy,txxz,txyy,txzz,txyz,tyyz,tyzz,txx,tyy,tzz,txy,txz,tyz
DOUBLE PRECISION :: exxsrg,eyysrg,ezzsrg,exysrg,exzsrg,eyzsrg, &
                    exxquadi,eyyquadi,ezzquadi,exyquadi,exzquadi,eyzquadi, &
                    exxquadj,eyyquadj,ezzquadj,exyquadj,exzquadj,eyzquadj, &
                    exxdipi,eyydipi,ezzdipi,exydipi,exzdipi,eyzdipi, &
                    exxdipj,eyydipj,ezzdipj,exydipj,exzdipj,eyzdipj
 
elecx=0.d0
elecy=0.d0
elecz=0.d0
exx=0.d0
eyy=0.d0
ezz=0.d0
exy=0.d0
exz=0.d0
eyz=0.d0
elecxsr=0.d0
elecysr=0.d0
eleczsr=0.d0
exxsr=0.d0
eyysr=0.d0
ezzsr=0.d0
exysr=0.d0
exzsr=0.d0
eyzsr=0.d0
engeff=0.d0

do j=2,num
   jpoint=ntype(j)
   qjrec=1.0d0/q(j)
   do i=1,j-1
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
      dr7rec=dr5rec*drsqrec
      dr9rec=dr7rec*drsqrec

      chgtem=q(i)*q(j)
      qq2api=etapi*chgtem
      expewld=exp(-etasq*drsq)
!---> Memmory Reduction_S
!     erfcr=chgtem*erfc(i,j)*drrec
      erfcr=chgtem*erfc(numx)*drrec
!<--- Memmory Reduction_E
      egam=(erfcr+qq2api*expewld)*drsqrec

!---> Memmory Reduction_S
!     ettx=egam*dxsav(i,j)
!     etty=egam*dysav(i,j)
!     ettz=egam*dzsav(i,j)
      ettx=egam*dxsav(numx)
      etty=egam*dysav(numx)
      ettz=egam*dzsav(numx)
!<--- Memmory Reduction_E

      elecx(j)=elecx(j)-(ettx*qjrec)
      elecy(j)=elecy(j)-(etty*qjrec)
      elecz(j)=elecz(j)-(ettz*qjrec)
      elecx(i)=elecx(i)+(ettx*qirec)
      elecy(i)=elecy(i)+(etty*qirec)
      elecz(i)=elecz(i)+(ettz*qirec)

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

      elecx(j)=elecx(j)-frrxdipqij*qjrec
      elecy(j)=elecy(j)-frrydipqij*qjrec
      elecz(j)=elecz(j)-frrzdipqij*qjrec
      elecx(i)=elecx(i)+frrxdipqji*qirec
      elecy(i)=elecy(i)+frrydipqji*qirec
      elecz(i)=elecz(i)+frrzdipqji*qirec

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
      txx=-sxx*dr3rec
      tyy=-syy*dr3rec
      tzz=-szz*dr3rec
      txy=-sxy*dr3rec
      txz=-sxz*dr3rec
      tyz=-syz*dr3rec

      threedr7=3.0d0*dr7rec
!---> Memmory Reduction_S
!     txxx=(5.0d0*dxsav(i,j)*dxsav(i,j)*dxsav(i,j) &
!         -3.0d0*dxsav(i,j)*drsq)*threedr7
!     tyyy=(5.0d0*dysav(i,j)*dysav(i,j)*dysav(i,j) &
!         -3.0d0*dysav(i,j)*drsq)*threedr7
!     tzzz=(5.0d0*dzsav(i,j)*dzsav(i,j)*dzsav(i,j) &
!         -3.0d0*dzsav(i,j)*drsq)*threedr7
!     txxy=(5.0d0*dxsav(i,j)*dxsav(i,j)*dysav(i,j) &
!         -drsq*dysav(i,j))*threedr7
!     txyy=(5.0d0*dxsav(i,j)*dysav(i,j)*dysav(i,j) &
!         -drsq*dxsav(i,j))*threedr7
!     txxz=(5.0d0*dxsav(i,j)*dxsav(i,j)*dzsav(i,j) &
!         -drsq*dzsav(i,j))*threedr7
!     txzz=(5.0d0*dxsav(i,j)*dzsav(i,j)*dzsav(i,j) &
!         -drsq*dxsav(i,j))*threedr7
!     tyyz=(5.0d0*dysav(i,j)*dysav(i,j)*dzsav(i,j) &
!         -drsq*dzsav(i,j))*threedr7
!     tyzz=(5.0d0*dysav(i,j)*dzsav(i,j)*dzsav(i,j) &
!         -drsq*dysav(i,j))*threedr7
!     txyz=15.0d0*dr7rec*dxsav(i,j)*dysav(i,j)*dzsav(i,j)
      txxx=(5.0d0*dxsav(numx)*dxsav(numx)*dxsav(numx) &
          -3.0d0*dxsav(numx)*drsq)*threedr7
      tyyy=(5.0d0*dysav(numx)*dysav(numx)*dysav(numx) &
          -3.0d0*dysav(numx)*drsq)*threedr7
      tzzz=(5.0d0*dzsav(numx)*dzsav(numx)*dzsav(numx) &
          -3.0d0*dzsav(numx)*drsq)*threedr7
      txxy=(5.0d0*dxsav(numx)*dxsav(numx)*dysav(numx) &
          -drsq*dysav(numx))*threedr7
      txyy=(5.0d0*dxsav(numx)*dysav(numx)*dysav(numx) &
          -drsq*dxsav(numx))*threedr7
      txxz=(5.0d0*dxsav(numx)*dxsav(numx)*dzsav(numx) &
          -drsq*dzsav(numx))*threedr7
      txzz=(5.0d0*dxsav(numx)*dzsav(numx)*dzsav(numx) &
          -drsq*dxsav(numx))*threedr7
      tyyz=(5.0d0*dysav(numx)*dysav(numx)*dzsav(numx) &
          -drsq*dzsav(numx))*threedr7
      tyzz=(5.0d0*dysav(numx)*dzsav(numx)*dzsav(numx) &
          -drsq*dysav(numx))*threedr7
      txyz=15.0d0*dr7rec*dxsav(numx)*dysav(numx)*dzsav(numx)
!<--- Memmory Reduction_E
      efact3= 2.0d0*efact2*(drsqrec+etasq)

      engeff(i)=engeff(i)+erfcr*qirec
      engeff(j)=engeff(j)+erfcr*qjrec

!---> Memmory Reduction_S
!     dx=dxsav(i,j)
!     dy=dysav(i,j)
!     dz=dzsav(i,j)
      dx=dxsav(numx)
      dy=dysav(numx)
      dz=dzsav(numx)
!<--- Memmory Reduction_E
      dxx=dx*dx
      dyy=dy*dy
      dzz=dz*dz
      dxy=dx*dy
      dxz=dx*dz
      dyz=dy*dz

      chgquadengij=((txx*quadxx(i)+txy*quadxy(i)+txy*quadxy(i) &
                    +txz*quadxz(i)+txz*quadxz(i)+tyy*quadyy(i) &
                    +tyz*quadyz(i)+tzz*quadzz(i)+tyz*quadyz(i))*q(j))
      chgquadengji=((txx*quadxx(j)+txy*quadxy(j)+txy*quadxy(j) &
                    +txz*quadxz(j)+txz*quadxz(j)+tyy*quadyy(j) &
                    +tyz*quadyz(j)+tzz*quadzz(j)+tyz*quadyz(j))*q(i))

      chgquad2ij=((dxx*quadxx(i)+dxy*quadxy(i)+dxy*quadxy(i) &
                  +dxz*quadxz(i)+dxz*quadxz(i)+dyy*quadyy(i) &
                  +dyz*quadyz(i)+dzz*quadzz(i)+dyz*quadyz(i))*q(j))
      chgquad2ji=((dxx*quadxx(j)+dxy*quadxy(j)+dxy*quadxy(j) &
                  +dxz*quadxz(j)+dxz*quadxz(j)+dyy*quadyy(j) &
                  +dyz*quadyz(j)+dzz*quadzz(j)+dyz*quadyz(j))*q(i))

      fqquadxij=(q(j)*(quadxx(i)*txxx+quadxy(i)*txxy+quadxz(i)*txxz &
                    +quadxy(i)*txxy+quadyy(i)*txyy+quadyz(i)*txyz &
                    +quadxz(i)*txxz+quadyz(i)*txyz+quadzz(i)*txzz))
      fqquadxji=(q(i)*(quadxx(j)*txxx+quadxy(j)*txxy+quadxz(j)*txxz &
                    +quadxy(j)*txxy+quadyy(j)*txyy+quadyz(j)*txyz &
                    +quadxz(j)*txxz+quadyz(j)*txyz+quadzz(j)*txzz))
      fqquadyij=(q(j)*(quadxx(i)*txxy+quadxy(i)*txyy+quadxz(i)*txyz &
                    +quadxy(i)*txyy+quadyy(i)*tyyy+quadyz(i)*tyyz &
                    +quadxz(i)*txyz+quadyz(i)*tyyz+quadzz(i)*tyzz))
      fqquadyji=(q(i)*(quadxx(j)*txxy+quadxy(j)*txyy+quadxz(j)*txyz &
                    +quadxy(j)*txyy+quadyy(j)*tyyy+quadyz(j)*tyyz &
                    +quadxz(j)*txyz+quadyz(j)*tyyz+quadzz(j)*tyzz))
      fqquadzij=(q(j)*(quadxx(i)*txxz+quadxy(i)*txyz+quadxz(i)*txzz &
                    +quadxy(i)*txyz+quadyy(i)*tyyz+quadyz(i)*tyzz &
                    +quadxz(i)*txzz+quadyz(i)*tyzz+quadzz(i)*tzzz))
      fqquadzji=(q(i)*(quadxx(j)*txxz+quadxy(j)*txyz+quadxz(j)*txzz &
                    +quadxy(j)*txyz+quadyy(j)*tyyz+quadyz(j)*tyzz &
                    +quadxz(j)*txzz+quadyz(j)*tyzz+quadzz(j)*tzzz))

      efact4=efact2*dr*drsq

!---> Memmory Reduction_S
!     fqquadxij=fqquadxij*efact+chgquadengij*efact4*dxsav(i,j) &
!              +chgquad2ij*efact3*dxsav(i,j)         
!     fqquadyij=fqquadyij*efact+chgquadengij*efact4*dysav(i,j) &
!              +chgquad2ij*efact3*dysav(i,j)
!     fqquadzij=fqquadzij*efact+chgquadengij*efact4*dzsav(i,j) &
!              +chgquad2ij*efact3*dzsav(i,j)
!     fqquadxji=fqquadxji*efact+chgquadengji*efact4*dxsav(i,j) &
!              +chgquad2ji*efact3*dxsav(i,j)
!     fqquadyji=fqquadyji*efact+chgquadengji*efact4*dysav(i,j) &
!              +chgquad2ji*efact3*dysav(i,j)
!     fqquadzji=fqquadzji*efact+chgquadengji*efact4*dzsav(i,j) &
!              +chgquad2ji*efact3*dzsav(i,j)
      fqquadxij=fqquadxij*efact+chgquadengij*efact4*dxsav(numx) &
               +chgquad2ij*efact3*dxsav(numx)         
      fqquadyij=fqquadyij*efact+chgquadengij*efact4*dysav(numx) &
               +chgquad2ij*efact3*dysav(numx)
      fqquadzij=fqquadzij*efact+chgquadengij*efact4*dzsav(numx) &
               +chgquad2ij*efact3*dzsav(numx)
      fqquadxji=fqquadxji*efact+chgquadengji*efact4*dxsav(numx) &
               +chgquad2ji*efact3*dxsav(numx)
      fqquadyji=fqquadyji*efact+chgquadengji*efact4*dysav(numx) &
               +chgquad2ji*efact3*dysav(numx)
      fqquadzji=fqquadzji*efact+chgquadengji*efact4*dzsav(numx) &
               +chgquad2ji*efact3*dzsav(numx)
!<--- Memmory Reduction_E

!---> Memmory Reduction_S
!     fqquadxij=fqquadxij-2.0d0*efact2*q(j)* &
!               (quadxx(i)*dxsav(i,j)+quadxy(i)*dysav(i,j)+ &
!                quadxz(i)*dzsav(i,j))
!     fqquadyij=fqquadyij-2.0d0*efact2*q(j)* &
!               (quadxy(i)*dxsav(i,j)+quadyy(i)*dysav(i,j)+ &
!                quadyz(i)*dzsav(i,j))
!     fqquadzij=fqquadzij-2.0d0*efact2*q(j)* &
!               (quadxz(i)*dxsav(i,j)+quadyz(i)*dysav(i,j)+ &
!                quadzz(i)*dzsav(i,j))
!     fqquadxji=fqquadxji-2.0d0*efact2*q(i)* &
!               (quadxx(j)*dxsav(i,j)+quadxy(j)*dysav(i,j)+ &
!                quadxz(j)*dzsav(i,j))
!     fqquadyji=fqquadyji-2.0d0*efact2*q(i)* &
!               (quadxy(j)*dxsav(i,j)+quadyy(j)*dysav(i,j)+ &
!                quadyz(j)*dzsav(i,j))
!     fqquadzji=fqquadzji-2.0d0*efact2*q(i)* &
!               (quadxz(j)*dxsav(i,j)+quadyz(j)*dysav(i,j)+ &
!                quadzz(j)*dzsav(i,j))
      fqquadxij=fqquadxij-2.0d0*efact2*q(j)* &
                (quadxx(i)*dxsav(numx)+quadxy(i)*dysav(numx)+ &
                 quadxz(i)*dzsav(numx))
      fqquadyij=fqquadyij-2.0d0*efact2*q(j)* &
                (quadxy(i)*dxsav(numx)+quadyy(i)*dysav(numx)+ &
                 quadyz(i)*dzsav(numx))
      fqquadzij=fqquadzij-2.0d0*efact2*q(j)* &
                (quadxz(i)*dxsav(numx)+quadyz(i)*dysav(numx)+ &
                 quadzz(i)*dzsav(numx))
      fqquadxji=fqquadxji-2.0d0*efact2*q(i)* &
                (quadxx(j)*dxsav(numx)+quadxy(j)*dysav(numx)+ &
                 quadxz(j)*dzsav(numx))
      fqquadyji=fqquadyji-2.0d0*efact2*q(i)* &
                (quadxy(j)*dxsav(numx)+quadyy(j)*dysav(numx)+ &
                 quadyz(j)*dzsav(numx))
      fqquadzji=fqquadzji-2.0d0*efact2*q(i)* &
                (quadxz(j)*dxsav(numx)+quadyz(j)*dysav(numx)+ &
                 quadzz(j)*dzsav(numx))
!<--- Memmory Reduction_E

      fqquadxij=fqquadxij*onethird
      fqquadyij=fqquadyij*onethird
      fqquadzij=fqquadzij*onethird
      fqquadxji=fqquadxji*onethird
      fqquadyji=fqquadyji*onethird
      fqquadzji=fqquadzji*onethird

      elecx(j)=elecx(j)-fqquadxij*qjrec
      elecy(j)=elecy(j)-fqquadyij*qjrec
      elecz(j)=elecz(j)-fqquadzij*qjrec
      elecx(i)=elecx(i)+fqquadxji*qirec
      elecy(i)=elecy(i)+fqquadyji*qirec
      elecz(i)=elecz(i)+fqquadzji*qirec

      txxxx=(105.0d0*dx*dx*dx*dx*dr9rec)+(9.0d0*dr5rec) &
           -(90.0d0*dx*dx*dr7rec)
      tyyyy=(105.0d0*dy*dy*dy*dy*dr9rec)+(9.0d0*dr5rec) &
           -(90.0d0*dy*dy*dr7rec)
      tzzzz=(105.0d0*dz*dz*dz*dz*dr9rec)+(9.0d0*dr5rec) &
           -(90.0d0*dz*dz*dr7rec)

      txxxy=(105.0d0*dx*dx*dx*dy*dr9rec)-(45.0d0*dx*dy*dr7rec)
      txxxz=(105.0d0*dx*dx*dx*dz*dr9rec)-(45.0d0*dx*dz*dr7rec)
      txyyy=(105.0d0*dx*dy*dy*dy*dr9rec)-(45.0d0*dx*dy*dr7rec)
      txzzz=(105.0d0*dx*dz*dz*dz*dr9rec)-(45.0d0*dx*dz*dr7rec)
      tyyyz=(105.0d0*dz*dy*dy*dy*dr9rec)-(45.0d0*dz*dy*dr7rec)
      tyzzz=(105.0d0*dz*dz*dz*dy*dr9rec)-(45.0d0*dz*dy*dr7rec)

      txxyy=(105.0d0*dx*dx*dy*dy*dr9rec)+(3.0d0*dr5rec) &
           -(15.0d0*((dx*dx)+(dy*dy))*dr7rec)
      txxzz=(105.0d0*dx*dx*dz*dz*dr9rec)+(3.0d0*dr5rec) &
           -(15.0d0*((dx*dx)+(dz*dz))*dr7rec)
      tyyzz=(105.0d0*dz*dz*dy*dy*dr9rec)+(3.0d0*dr5rec) &
           -(15.0d0*((dz*dz)+(dy*dy))*dr7rec)

      txxyz=(105.0d0*dx*dx*dy*dz*dr9rec)-(15.0d0*dy*dz*dr7rec)
      txyyz=(105.0d0*dx*dy*dy*dz*dr9rec)-(15.0d0*dx*dz*dr7rec)
      txyzz=(105.0d0*dx*dz*dy*dz*dr9rec)-(15.0d0*dy*dx*dr7rec)

      dampsumfi=1.0d0
      xf=1.0d0
      factorial=1.0d0
      do kk=1,nkdamp(jpoint,ipoint)
         xf=(xf*dampa(jpoint,ipoint)*dr)
         factorial=factorial*float(kk)
         dampsumfi=dampsumfi+(xf/factorial)
      enddo   

      dampaexpi=dexp(-dampa(jpoint,ipoint)*dr)
      dampfunci=-dampsumfi*dampaexpi*dampfac(jpoint,ipoint)
      r3dampi=dr3rec*dampfunci*q(j)
!---> Memmory Reduction_S
!     elecsrxi=dxsav(i,j)*r3dampi
!     elecsryi=dysav(i,j)*r3dampi
!     elecsrzi=dzsav(i,j)*r3dampi
      elecsrxi=dxsav(numx)*r3dampi
      elecsryi=dysav(numx)*r3dampi
      elecsrzi=dzsav(numx)*r3dampi
!<--- Memmory Reduction_E

      elecxsr(i)=elecxsr(i)+elecsrxi
      elecysr(i)=elecysr(i)+elecsryi
      eleczsr(i)=eleczsr(i)+elecsrzi

      dampsumfj=1.0d0
      xf=1.0d0
      factorial=1.0d0
      do kk=1,nkdamp(ipoint,jpoint)
         xf=(xf*dampa(ipoint,jpoint)*dr)
         factorial=factorial*float(kk)
         dampsumfj=dampsumfj+(xf/factorial)
      enddo  

      dampaexpj=dexp(-dampa(ipoint,jpoint)*dr)
      dampfuncj=-dampsumfj*dampaexpj*dampfac(ipoint,jpoint)
      r3dampj=dr3rec*dampfuncj*q(i)
!---> Memmory Reduction_S
!     elecsrxj=-dxsav(i,j)*r3dampj
!     elecsryj=-dysav(i,j)*r3dampj
!     elecsrzj=-dzsav(i,j)*r3dampj
      elecsrxj=-dxsav(numx)*r3dampj
      elecsryj=-dysav(numx)*r3dampj
      elecsrzj=-dzsav(numx)*r3dampj
!<--- Memmory Reduction_E

      elecxsr(j)=elecxsr(j)+elecsrxj
      elecysr(j)=elecysr(j)+elecsryj
      eleczsr(j)=eleczsr(j)+elecsrzj

      txxt=txx*efact+efact2*dxx
      tyyt=tyy*efact+efact2*dyy
      tzzt=tzz*efact+efact2*dzz
      txyt=txy*efact+efact2*dxy
      txzt=txz*efact+efact2*dxz
      tyzt=tyz*efact+efact2*dyz

      exx(j)=exx(j)-txxt*q(i)
      eyy(j)=eyy(j)-tyyt*q(i)
      ezz(j)=ezz(j)-tzzt*q(i)
      exy(j)=exy(j)-txyt*q(i)
      exz(j)=exz(j)-txzt*q(i)
      eyz(j)=eyz(j)-tyzt*q(i)
      exx(i)=exx(i)-txxt*q(j)
      eyy(i)=eyy(i)-tyyt*q(j)
      ezz(i)=ezz(i)-tzzt*q(j)
      exy(i)=exy(i)-txyt*q(j)
      exz(i)=exz(i)-txzt*q(j)
      eyz(i)=eyz(i)-tyzt*q(j)

      exxdipi=txxx*xmu(i)+txxy*ymu(i)+txxz*zmu(i)
      eyydipi=tyyy*ymu(i)+txyy*xmu(i)+tyyz*zmu(i)
      ezzdipi=tzzz*zmu(i)+tyzz*ymu(i)+txzz*xmu(i)
      exydipi=txxy*xmu(i)+txyy*ymu(i)+txyz*zmu(i)
      exzdipi=txxz*xmu(i)+txyz*ymu(i)+txzz*zmu(i)
      eyzdipi=txyz*xmu(i)+tyyz*ymu(i)+tyzz*zmu(i)

      exx(j)=exx(j)+exxdipi
      eyy(j)=eyy(j)+eyydipi
      ezz(j)=ezz(j)+ezzdipi
      exy(j)=exy(j)+exydipi
      exz(j)=exz(j)+exzdipi
      eyz(j)=eyz(j)+eyzdipi

      exxdipj=-txxx*xmu(j)-txxy*ymu(j)-txxz*zmu(j)
      eyydipj=-tyyy*ymu(j)-txyy*xmu(j)-tyyz*zmu(j)
      ezzdipj=-tzzz*zmu(j)-tyzz*ymu(j)-txzz*xmu(j)
      exydipj=-txxy*xmu(j)-txyy*ymu(j)-txyz*zmu(j)
      exzdipj=-txxz*xmu(j)-txyz*ymu(j)-txzz*zmu(j)
      eyzdipj=-txyz*xmu(j)-tyyz*ymu(j)-tyzz*zmu(j)

      exx(i)=exx(i)+exxdipj
      eyy(i)=eyy(i)+eyydipj
      ezz(i)=ezz(i)+ezzdipj
      exy(i)=exy(i)+exydipj
      exz(i)=exz(i)+exzdipj
      eyz(i)=eyz(i)+eyzdipj

      exxquadi=-txxxx*quadxx(i) &
               -2.0d0*txxxy*quadxy(i)-2.0d0*txxxz*quadxz(i) &
               -txxyy*quadyy(i)-txxzz*quadzz(i)-2.0d0*txxyz*quadyz(i)
      eyyquadi=-txxyy*quadxx(i) &
               -2.0d0*txyyy*quadxy(i)-2.0d0*txyyz*quadxz(i) &
               -tyyyy*quadyy(i)-tyyzz*quadzz(i)-2.0d0*tyyyz*quadyz(i)
      ezzquadi=-txxzz*quadxx(i) &
               -2.0d0*txyzz*quadxy(i)-2.0d0*txzzz*quadxz(i) &
               -tyyzz*quadyy(i)-tzzzz*quadzz(i)-2.0d0*tyzzz*quadyz(i)
      exyquadi=-txxxy*quadxx(i) &
               -2.0d0*txxyy*quadxy(i)-2.0d0*txxyz*quadxz(i) &
               -txyyy*quadyy(i)-txyzz*quadzz(i)-2.0d0*txyyz*quadyz(i)
      exzquadi=-txxxz*quadxx(i) &
               -2.0d0*txxyz*quadxy(i)-2.0d0*txxzz*quadxz(i) &
               -txyyz*quadyy(i)-txzzz*quadzz(i)-2.0d0*txyzz*quadyz(i)
      eyzquadi=-txxyz*quadxx(i) &
               -2.0d0*txyyz*quadxy(i)-2.0d0*txyzz*quadxz(i) &
               -tyyyz*quadyy(i)-tyzzz*quadzz(i)-2.0d0*tyyzz*quadyz(i)

      exx(j)=exx(j)+exxquadi*onethird
      eyy(j)=eyy(j)+eyyquadi*onethird
      ezz(j)=ezz(j)+ezzquadi*onethird
      exy(j)=exy(j)+exyquadi*onethird
      exz(j)=exz(j)+exzquadi*onethird
      eyz(j)=eyz(j)+eyzquadi*onethird

      exxquadj=-txxxx*quadxx(j) &
               -2.0d0*txxxy*quadxy(j)-2.0d0*txxxz*quadxz(j) &
               -txxyy*quadyy(j)-txxzz*quadzz(j)-2.0d0*txxyz*quadyz(j)
      eyyquadj=-txxyy*quadxx(j) &
               -2.0d0*txyyy*quadxy(j)-2.0d0*txyyz*quadxz(j) &
               -tyyyy*quadyy(j)-tyyzz*quadzz(j)-2.0d0*tyyyz*quadyz(j)
      ezzquadj=-txxzz*quadxx(j) &
               -2.0d0*txyzz*quadxy(j)-2.0d0*txzzz*quadxz(j) &
               -tyyzz*quadyy(j)-tzzzz*quadzz(j)-2.0d0*tyzzz*quadyz(j)
      exyquadj=-txxxy*quadxx(j) &
               -2.0d0*txxyy*quadxy(j)-2.0d0*txxyz*quadxz(j) &
               -txyyy*quadyy(j)-txyzz*quadzz(j)-2.0d0*txyyz*quadyz(j)
      exzquadj=-txxxz*quadxx(j) &
               -2.0d0*txxyz*quadxy(j)-2.0d0*txxzz*quadxz(j) &
               -txyyz*quadyy(j)-txzzz*quadzz(j)-2.0d0*txyzz*quadyz(j)
      eyzquadj=-txxyz*quadxx(j) &
               -2.0d0*txyyz*quadxy(j)-2.0d0*txyzz*quadxz(j) &
               -tyyyz*quadyy(j)-tyzzz*quadzz(j)-2.0d0*tyyzz*quadyz(j)

      exx(i)=exx(i)+exxquadj*onethird
      eyy(i)=eyy(i)+eyyquadj*onethird
      ezz(i)=ezz(i)+ezzquadj*onethird
      exy(i)=exy(i)+exyquadj*onethird
      exz(i)=exz(i)+exzquadj*onethird
      eyz(i)=eyz(i)+eyzquadj*onethird

      dampsumfi=1.0d0
      xf=1.0d0
      factorial=1.0d0
      do kk=1,nkfg(jpoint,ipoint)
         xf=(xf*fgb(jpoint,ipoint)*dr)
         factorial=factorial*float(kk)
         dampsumfi=dampsumfi+(xf/factorial)
      enddo   

      dampaexpi=dexp(-fgb(jpoint,ipoint)*dr)
      dampfunci=-dampsumfi*dampaexpi*fgc(jpoint,ipoint)
      exxsrg=txx*dampfunci*q(j)
      eyysrg=tyy*dampfunci*q(j)
      ezzsrg=tzz*dampfunci*q(j)
      exysrg=txy*dampfunci*q(j)
      exzsrg=txz*dampfunci*q(j)
      eyzsrg=tyz*dampfunci*q(j)

      exxsr(i)=exxsr(i)-exxsrg
      eyysr(i)=eyysr(i)-eyysrg
      ezzsr(i)=ezzsr(i)-ezzsrg
      exysr(i)=exysr(i)-exysrg
      exzsr(i)=exzsr(i)-exzsrg
      eyzsr(i)=eyzsr(i)-eyzsrg

      dampsumfj=1.0d0
      xf=1.0d0
      factorial=1.0d0
      do kk=1,nkfg(ipoint,jpoint)
         xf=(xf*fgb(ipoint,jpoint)*dr)
         factorial=factorial*float(kk)
         dampsumfj=dampsumfj+(xf/factorial)
      enddo   

      dampaexpj=dexp(-fgb(ipoint,jpoint)*dr)
      dampfuncj=-dampsumfj*dampaexpj*fgc(ipoint,jpoint)
      exxsrg=txx*dampfuncj*q(i)
      eyysrg=tyy*dampfuncj*q(i)
      ezzsrg=tzz*dampfuncj*q(i)
      exysrg=txy*dampfuncj*q(i)
      exzsrg=txz*dampfuncj*q(i)
      eyzsrg=tyz*dampfuncj*q(i)

      exxsr(j)=exxsr(j)-exxsrg
      eyysr(j)=eyysr(j)-eyysrg
      ezzsr(j)=ezzsr(j)-ezzsrg
      exysr(j)=exysr(j)-exysrg
      exzsr(j)=exzsr(j)-exzsrg
      eyzsr(j)=eyzsr(j)-eyzsrg

   enddo   
enddo   

elecx=elecx+elecxsr
elecy=elecy+elecysr
elecz=elecz+eleczsr
exx=exx+exxsr
eyy=eyy+eyysr
ezz=ezz+ezzsr
exy=exy+exysr
exz=exz+exzsr
eyz=eyz+eyzsr

return
END SUBROUTINE
