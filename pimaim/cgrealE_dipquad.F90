!************************************************************
!   Calculates just that part of the real energy needed by the
!   conjgrad routine (AGUADO 2002)
!************************************************************

SUBROUTINE cgrealE_dipquad

USE commondata, ONLY: x,y,z, &
                      q,xmu,ymu,zmu,quadxx,quadyy,quadzz,quadxy,quadxz, &
                      quadyz,num,ewlog,twopi,ntype,rsqmax,etapi,etasq,erfc, &
!---> Memmory Reduction_S
!                     onethird,dxsav,dysav,dzsav,nspec
                      onethird,dxsav,dysav,dzsav,nspec,num2,numx
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
DOUBLE PRECISION :: efact,efact2,efact3,efact4, &
                    chgquad2ij,chgquad2ji, &
                    chgquadengij,chgquadengji, &
                    chgdipengij,chgdipengji,chgtem
DOUBLE PRECISION :: drsq,drsqrec,dr,dr3rec,dr4rec,dr5rec, &
                    dr7rec,dr9rec,dx,dy,dz,dxx,dyy,dzz,dxy,dxz,dyz, &
                    sxx,syy,szz,sxy,sxz,syz,threedr7,qirec,qjrec,drrec
DOUBLE PRECISION :: fqquadxij,fqquadxji,fqquadyij,fqquadyji,fqquadzij,fqquadzji
DOUBLE PRECISION :: qdotmuxij,qdotmuxji,qdotmuyij,qdotmuyji, &
                    qdotmuzij,qdotmuzji,txijdotthetaij,txijdotthetaji, &
                  tyijdotthetaij,tyijdotthetaji,tzijdotthetaij,tzijdotthetaji, &
                    thetaxidotrij,thetaxjdotrij,thetayidotrij,thetayjdotrij, &
                    thetazidotrij,thetazjdotrij
DOUBLE PRECISION :: txxxx,tyyyy,tzzzz,txxxy,txxxz,txxyy,txxzz, &
                    txyyy,txzzz,tyyyz,tyzzz,tyyzz,txxyz,txyyz,txyzz,txxx,tyyy, &
                 tzzz,txxy,txxz,txyy,txzz,txyz,tyyz,tyzz,txx,tyy,tzz,txy,txz,tyz
DOUBLE PRECISION :: exxsrg,eyysrg,ezzsrg,exysrg,exzsrg,eyzsrg, &
                    exxquadi,eyyquadi,ezzquadi,exyquadi,exzquadi,eyzquadi, &
                    exxquadj,eyyquadj,ezzquadj,exyquadj,exzquadj,eyzquadj, &
                    exxdipi,eyydipi,ezzdipi,exydipi,exzdipi,eyzdipi, &
                    exxdipj,eyydipj,ezzdipj,exydipj,exzdipj,eyzdipj
 
!Calculate real space energy, forces, stress tensor, electric fields
! and field gradients...
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
      if (drsq.ge.rsqmax) CYCLE

      drsqrec=1.d0/drsq
      dr=dsqrt(drsq)
      drrec=1.d0/dr
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
!     xmuidotrij=xmu(i)*dxsav(i,j)+ymu(i)*dysav(i,j)+zmu(i)*dzsav(i,j)
!     xmujdotrij=xmu(j)*dxsav(i,j)+ymu(j)*dysav(i,j)+zmu(j)*dzsav(i,j)
      xmuidotrij=xmu(i)*dxsav(numx)+ymu(i)*dysav(numx)+zmu(i)*dzsav(numx)
      xmujdotrij=xmu(j)*dxsav(numx)+ymu(j)*dysav(numx)+zmu(j)*dzsav(numx)
!<--- Memmory Reduction_E
!
! Permanent charge - dipole energy:
!
!---> Memmory Reduction_S
!     efact=etapi*expewld*dr+erfc(i,j)
      efact=etapi*expewld*dr+erfc(numx)
!<--- Memmory Reduction_E
      chgdipengij=(q(j)*xmuidotrij)*dr3rec
      chgdipengji=-(q(i)*xmujdotrij)*dr3rec
!
! Permanent charge - dipole forces:
!
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
!
! Dipole-dipole interaction tensor:
!
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
!
! Dipole-quadrupole interaction tensor...
!
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
!
! End of dipole only interaction calcs.
! Start of quadrupole handling section.
!
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
!
! quadrupole-charge interaction:
!   theta_alpha_beta(i)*T_alpha_beta*q(j)+
!   q(i)*T_alpha_beta*theta_alpha_beta(j)
!
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
!
! Quadrupole-charge forces:
!  -theta_beta_gamma(i)*T_alpha_beta_gamma*q(j)
!
      fqquadxij=(q(j)*(quadxx(i)*txxx+quadxy(i)*txxy+quadxz(i)*txxz &
                    +quadxy(i)*txxy+quadyy(i)*txyy+quadyz(i)*txyz &
                    +quadxz(i)*txxz+quadyz(i)*txyz+quadzz(i)*txzz))
!
!  -q(i)*T_alpha_beta_gamma*theta_beta_gamma(j)
!
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
!
! Save fqquad's for use in dipole-quadrupole Ewald sum (AGUADO)
!
      txijdotthetaij=fqquadxji*qirec
      txijdotthetaji=fqquadxij*qjrec
      tyijdotthetaij=fqquadyji*qirec
      tyijdotthetaji=fqquadyij*qjrec
      tzijdotthetaij=fqquadzji*qirec
      tzijdotthetaji=fqquadzij*qjrec            
!
! Correct the forces for PBCs (AGUADO)
!
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

!     thetaxidotrij=quadxx(i)*dxsav(i,j)+quadxy(i)*dysav(i,j)+ &
!                   quadxz(i)*dzsav(i,j)
!     thetaxjdotrij=quadxx(j)*dxsav(i,j)+quadxy(j)*dysav(i,j)+ &
!                   quadxz(j)*dzsav(i,j)
!     thetayidotrij=quadxy(i)*dxsav(i,j)+quadyy(i)*dysav(i,j)+ &
!                   quadyz(i)*dzsav(i,j)
!     thetayjdotrij=quadxy(j)*dxsav(i,j)+quadyy(j)*dysav(i,j)+ &
!                   quadyz(j)*dzsav(i,j)
!     thetazidotrij=quadxz(i)*dxsav(i,j)+quadyz(i)*dysav(i,j)+ &
!                   quadzz(i)*dzsav(i,j)
!     thetazjdotrij=quadxz(j)*dxsav(i,j)+quadyz(j)*dysav(i,j)+ &
!                   quadzz(j)*dzsav(i,j)
      thetaxidotrij=quadxx(i)*dxsav(numx)+quadxy(i)*dysav(numx)+ &
                    quadxz(i)*dzsav(numx)
      thetaxjdotrij=quadxx(j)*dxsav(numx)+quadxy(j)*dysav(numx)+ &
                    quadxz(j)*dzsav(numx)
      thetayidotrij=quadxy(i)*dxsav(numx)+quadyy(i)*dysav(numx)+ &
                    quadyz(i)*dzsav(numx)
      thetayjdotrij=quadxy(j)*dxsav(numx)+quadyy(j)*dysav(numx)+ &
                    quadyz(j)*dzsav(numx)
      thetazidotrij=quadxz(i)*dxsav(numx)+quadyz(i)*dysav(numx)+ &
                    quadzz(i)*dzsav(numx)
      thetazjdotrij=quadxz(j)*dxsav(numx)+quadyz(j)*dysav(numx)+ &
                    quadzz(j)*dzsav(numx)
!<--- Memmory Reduction_E

      fqquadxij=fqquadxij-2.0d0*efact2*q(j)*thetaxidotrij
      fqquadyij=fqquadyij-2.0d0*efact2*q(j)*thetayidotrij
      fqquadzij=fqquadzij-2.0d0*efact2*q(j)*thetazidotrij
      fqquadxji=fqquadxji-2.0d0*efact2*q(i)*thetaxjdotrij
      fqquadyji=fqquadyji-2.0d0*efact2*q(i)*thetayjdotrij
      fqquadzji=fqquadzji-2.0d0*efact2*q(i)*thetazjdotrij

!    Factor of one-third from Buckingham

      fqquadxij=fqquadxij*onethird
      fqquadyij=fqquadyij*onethird
      fqquadzij=fqquadzij*onethird
      fqquadxji=fqquadxji*onethird
      fqquadyji=fqquadyji*onethird
      fqquadzji=fqquadzji*onethird

!---> Parallelization_S
!     elecx(j)=elecx(j)-fqquadxij*qjrec
!     elecy(j)=elecy(j)-fqquadyij*qjrec
!     elecz(j)=elecz(j)-fqquadzij*qjrec
!     elecx(i)=elecx(i)+fqquadxji*qirec
!     elecy(i)=elecy(i)+fqquadyji*qirec
!     elecz(i)=elecz(i)+fqquadzji*qirec
      eltmp(j)=eltmp(j)-fqquadxij*qjrec
      eltmp(num+j)=eltmp(num+j)-fqquadyij*qjrec
      eltmp(2*num+j)=eltmp(2*num+j)-fqquadzij*qjrec
      eltmp(i)=eltmp(i)+fqquadxji*qirec
      eltmp(num+i)=eltmp(num+i)+fqquadyji*qirec
      eltmp(2*num+i)=eltmp(2*num+i)+fqquadzji*qirec
!<--- Parallelization_E
!
! Quadrupole-quadrupole interaction tensor 
!
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
!
! 2) Dipole contribution +T_3 mu.
!    With PBC's corrections (AGUADO).
!
      exxdipi=(txxx*xmu(i)+txxy*ymu(i)+txxz*zmu(i))*efact+ &
               dxx*xmuidotrij*efact3-2.0d0*dx*xmu(i)*efact2+ &
               txx*xmuidotrij*efact4
      eyydipi=(tyyy*ymu(i)+txyy*xmu(i)+tyyz*zmu(i))*efact+ &
               dyy*xmuidotrij*efact3-2.0d0*dy*ymu(i)*efact2+ &
               tyy*xmuidotrij*efact4
      ezzdipi=(tzzz*zmu(i)+tyzz*ymu(i)+txzz*xmu(i))*efact+ &
               dzz*xmuidotrij*efact3-2.0d0*dz*zmu(i)*efact2+ &
               tzz*xmuidotrij*efact4
      exydipi=(txxy*xmu(i)+txyy*ymu(i)+txyz*zmu(i))*efact+ &
              dxy*xmuidotrij*efact3-(dy*xmu(i)+dx*ymu(i))*efact2+ &
               txy*xmuidotrij*efact4
      exzdipi=(txxz*xmu(i)+txyz*ymu(i)+txzz*zmu(i))*efact+ &
              dxz*xmuidotrij*efact3-(dz*xmu(i)+dx*zmu(i))*efact2+ &
               txz*xmuidotrij*efact4
      eyzdipi=(txyz*xmu(i)+tyyz*ymu(i)+tyzz*zmu(i))*efact+ &
              dyz*xmuidotrij*efact3-(dz*ymu(i)+dy*zmu(i))*efact2+ &
               tyz*xmuidotrij*efact4

!---> Parallelization_S
!     exx(j)=exx(j)+exxdipi
!     eyy(j)=eyy(j)+eyydipi
!     ezz(j)=ezz(j)+ezzdipi
!     exy(j)=exy(j)+exydipi
!     exz(j)=exz(j)+exzdipi
!     eyz(j)=eyz(j)+eyzdipi
      eltmp(3*num+j)=eltmp(3*num+j)+exxdipi
      eltmp(4*num+j)=eltmp(4*num+j)+eyydipi
      eltmp(5*num+j)=eltmp(5*num+j)+ezzdipi
      eltmp(6*num+j)=eltmp(6*num+j)+exydipi
      eltmp(7*num+j)=eltmp(7*num+j)+exzdipi
      eltmp(8*num+j)=eltmp(8*num+j)+eyzdipi
!<--- Parallelization_E

      exxdipj=(-txxx*xmu(j)-txxy*ymu(j)-txxz*zmu(j))*efact &
              -dxx*xmujdotrij*efact3+2.0d0*dx*xmu(j)*efact2- &
              txx*xmujdotrij*efact4
      eyydipj=(-tyyy*ymu(j)-txyy*xmu(j)-tyyz*zmu(j))*efact &
              -dyy*xmujdotrij*efact3+2.0d0*dy*ymu(j)*efact2- &
              tyy*xmujdotrij*efact4
      ezzdipj=(-tzzz*zmu(j)-tyzz*ymu(j)-txzz*xmu(j))*efact &
              -dzz*xmujdotrij*efact3+2.0d0*dz*zmu(j)*efact2- &
              tzz*xmujdotrij*efact4
      exydipj=(-txxy*xmu(j)-txyy*ymu(j)-txyz*zmu(j))*efact &
             -dxy*xmujdotrij*efact3+(dy*xmu(j)+dx*ymu(j))*efact2- &
              txy*xmujdotrij*efact4
      exzdipj=(-txxz*xmu(j)-txyz*ymu(j)-txzz*zmu(j))*efact &
             -dxz*xmujdotrij*efact3+(dz*xmu(j)+dx*zmu(j))*efact2- &
              txz*xmujdotrij*efact4
      eyzdipj=(-txyz*xmu(j)-tyyz*ymu(j)-tyzz*zmu(j))*efact &
             -dyz*xmujdotrij*efact3+(dz*ymu(j)+dy*zmu(j))*efact2- &
              tyz*xmujdotrij*efact4

!---> Parallelization_S
!     exx(i)=exx(i)+exxdipj
!     eyy(i)=eyy(i)+eyydipj
!     ezz(i)=ezz(i)+ezzdipj
!     exy(i)=exy(i)+exydipj
!     exz(i)=exz(i)+exzdipj
!     eyz(i)=eyz(i)+eyzdipj
      eltmp(3*num+i)=eltmp(3*num+i)+exxdipj
      eltmp(4*num+i)=eltmp(4*num+i)+eyydipj
      eltmp(5*num+i)=eltmp(5*num+i)+ezzdipj
      eltmp(6*num+i)=eltmp(6*num+i)+exydipj
      eltmp(7*num+i)=eltmp(7*num+i)+exzdipj
      eltmp(8*num+i)=eltmp(8*num+i)+eyzdipj
!<--- Parallelization_E

!
! 3) Quadrupole contribution - T_4 quad.
!    With PBC's corrections (AGUADO)
!
        exxquadi=-txxxx*quadxx(i)-2.0d0*txxxy*quadxy(i) &
                 -2.0d0*txxxz*quadxz(i)-txxyy*quadyy(i) &
                 -txxzz*quadzz(i)-2.0d0*txxyz*quadyz(i)
        eyyquadi=-txxyy*quadxx(i)-2.0d0*txyyy*quadxy(i) &
                 -2.0d0*txyyz*quadxz(i)-tyyyy*quadyy(i) &
                 -tyyzz*quadzz(i)-2.0d0*tyyyz*quadyz(i)
        ezzquadi=-txxzz*quadxx(i)-2.0d0*txyzz*quadxy(i) &
                 -2.0d0*txzzz*quadxz(i)-tyyzz*quadyy(i) &
                 -tzzzz*quadzz(i)-2.0d0*tyzzz*quadyz(i)
        exyquadi=-txxxy*quadxx(i)-2.0d0*txxyy*quadxy(i) &
                 -2.0d0*txxyz*quadxz(i)-txyyy*quadyy(i) &
                 -txyzz*quadzz(i)-2.0d0*txyyz*quadyz(i)
        exzquadi=-txxxz*quadxx(i)-2.0d0*txxyz*quadxy(i) &
                 -2.0d0*txxzz*quadxz(i)-txyyz*quadyy(i) &
                 -txzzz*quadzz(i)-2.0d0*txyzz*quadyz(i)
        eyzquadi=-txxyz*quadxx(i)-2.0d0*txyyz*quadxy(i) &
                 -2.0d0*txyzz*quadxz(i)-tyyyz*quadyy(i) &
                 -tyzzz*quadzz(i)-2.0d0*tyyzz*quadyz(i)

!---> Parallelization_S
!     exx(j)=exx(j)+exxquadi*onethird
!     eyy(j)=eyy(j)+eyyquadi*onethird
!     ezz(j)=ezz(j)+ezzquadi*onethird
!     exy(j)=exy(j)+exyquadi*onethird
!     exz(j)=exz(j)+exzquadi*onethird
!     eyz(j)=eyz(j)+eyzquadi*onethird
      eltmp(3*num+j)=eltmp(3*num+j)+exxquadi*onethird
      eltmp(4*num+j)=eltmp(4*num+j)+eyyquadi*onethird
      eltmp(5*num+j)=eltmp(5*num+j)+ezzquadi*onethird
      eltmp(6*num+j)=eltmp(6*num+j)+exyquadi*onethird
      eltmp(7*num+j)=eltmp(7*num+j)+exzquadi*onethird
      eltmp(8*num+j)=eltmp(8*num+j)+eyzquadi*onethird
!<--- Parallelization_E

       exxquadj=-txxxx*quadxx(j)-2.0d0*txxxy*quadxy(j) &
                -2.0d0*txxxz*quadxz(j)-txxyy*quadyy(j) &
                -txxzz*quadzz(j)-2.0d0*txxyz*quadyz(j)
       eyyquadj=-txxyy*quadxx(j)-2.0d0*txyyy*quadxy(j) &
                -2.0d0*txyyz*quadxz(j)-tyyyy*quadyy(j) &
                -tyyzz*quadzz(j)-2.0d0*tyyyz*quadyz(j)
       ezzquadj=-txxzz*quadxx(j)-2.0d0*txyzz*quadxy(j) &
                -2.0d0*txzzz*quadxz(j)-tyyzz*quadyy(j) &
                -tzzzz*quadzz(j)-2.0d0*tyzzz*quadyz(j)
       exyquadj=-txxxy*quadxx(j)-2.0d0*txxyy*quadxy(j) &
                -2.0d0*txxyz*quadxz(j)-txyyy*quadyy(j) &
                -txyzz*quadzz(j)-2.0d0*txyyz*quadyz(j)
       exzquadj=-txxxz*quadxx(j)-2.0d0*txxyz*quadxy(j) &
                -2.0d0*txxzz*quadxz(j)-txyyz*quadyy(j) &
                -txzzz*quadzz(j)-2.0d0*txyzz*quadyz(j)
       eyzquadj=-txxyz*quadxx(j)-2.0d0*txyyz*quadxy(j) &
                -2.0d0*txyzz*quadxz(j)-tyyyz*quadyy(j) &
                -tyzzz*quadzz(j)-2.0d0*tyyzz*quadyz(j)

!---> Parallelization_S
!     exx(i)=exx(i)+exxquadj*onethird
!     eyy(i)=eyy(i)+eyyquadj*onethird
!     ezz(i)=ezz(i)+ezzquadj*onethird
!     exy(i)=exy(i)+exyquadj*onethird
!     exz(i)=exz(i)+exzquadj*onethird
!     eyz(i)=eyz(i)+eyzquadj*onethird
      eltmp(3*num+i)=eltmp(3*num+i)+exxquadj*onethird
      eltmp(4*num+i)=eltmp(4*num+i)+eyyquadj*onethird
      eltmp(5*num+i)=eltmp(5*num+i)+ezzquadj*onethird
      eltmp(6*num+i)=eltmp(6*num+i)+exyquadj*onethird
      eltmp(7*num+i)=eltmp(7*num+i)+exzquadj*onethird
      eltmp(8*num+i)=eltmp(8*num+i)+eyzquadj*onethird
!<--- Parallelization_E
   enddo   
enddo   

totmz=0.d0
totmzd=0.d0
if(ewlog) then
facewvac=twopi/(boxlenx*boxleny*boxlenz)
totmzu=SUM(zmu)
totmzd=SUM(q*(hlab2(3,1)*x+hlab2(3,2)*y+hlab2(3,3)*z))
totmz=totmzd+totmzu
else
facewvac=0.d0
endif

return
END SUBROUTINE
