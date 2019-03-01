SUBROUTINE cgfirst2

USE commondata, ONLY: elecxq,elecyq,eleczq,exxq,eyyq,ezzq,exyq,exzq,eyzq,x,y,z, &
                      q,zmu, &
                      num,ewlog,twopi,ntype,rsqmax,etapi,etasq,erfc, &
                      engeff,dxsav,dysav,dzsav,nspec, &
!---> Memmory Reduction_S
!                     nkdamp,dampa,dampfac,nkfg,fgb,fgc
                      nkdamp,dampa,dampfac,nkfg,fgb,fgc,num2,numx, &
                      ewcorr
!<--- Memmory Reduction_E
USE boxdata, ONLY: hlab2,boxlenx,boxleny,boxlenz
!---> Parallelization_S
use mpipara
!<--- Parallelization_E

IMPLICIT NONE

INTEGER :: i,j,kk,ipoint,jpoint
DOUBLE PRECISION :: dxunitsq,dyunitsq,dzunitsq, &
                    ettx,etty,ettz,egam,erfcr,expewld, &
                    totmz,totmzd,totmzu,facewvac,qq2api
DOUBLE PRECISION :: efact,efact2,chgtem
DOUBLE PRECISION :: drsq,drsqrec,dr,dr3rec, &
                    dx,dy,dz,dxx,dyy,dzz,dxy,dxz,dyz, &
                    sxx,syy,szz,sxy,sxz,syz,qirec,qjrec,drrec
DOUBLE PRECISION :: txxt,tyyt,tzzt,txyt,txzt,tyzt, &
                    elecsrxi,elecsryi,elecsrzi,elecsrxj,elecsryj,elecsrzj, &
                    r3dampi,r3dampj,factorial,xf, &
                    dampsumfi,dampsumfj,dampfunci,dampfuncj,dampaexpi,dampaexpj
DOUBLE PRECISION :: txx,tyy,tzz,txy,txz,tyz
DOUBLE PRECISION :: exxsrg,eyysrg,ezzsrg,exysrg,exzsrg,eyzsrg
 
totmz=0.d0
totmzd=0.d0

!Calculate terms for Vacuum boundary conditions in Ewald summation
if(ewlog) then
   facewvac=twopi/(boxlenx*boxleny*boxlenz)
   totmzu=SUM(zmu)
   totmzd=SUM(q*(hlab2(3,1)*x+hlab2(3,2)*y+hlab2(3,3)*z))
   totmz=totmzd+totmzu
else
   facewvac=0.d0
endif
!Calculate real space energy, forces, stress tensor, electric fields
! and field gradients...
!---> Memmory Reduction_S
!<--- Memmory Reduction_E
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
!MPIAGUADO-->Here we could have severe load inbalance problems. The distribution
!            of ion pairs into nodes does not care about the spatial proximity of
!            those pairs, so some nodes may just contain distant particles and
!            finish immediately!
!            There should be a better way to distribute the work in terms of
!            connected clusters of particles, at least for a solid. For a liquid,
!            this problem does not apply (or is impossible to solve!).
      if (drsq.ge.rsqmax) CYCLE
!<--MPIAGUADO

      drsqrec=1.d0/drsq
      dr=dsqrt(drsq)
      drrec=1.d0/dr
      dr3rec=drrec*drsqrec

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

!---> Parallelization_S
!     elecxq(j)=elecxq(j)-(ettx*qjrec)
!     elecyq(j)=elecyq(j)-(etty*qjrec)
!     eleczq(j)=eleczq(j)-(ettz*qjrec)

!     elecxq(i)=elecxq(i)+(ettx*qirec)
!     elecyq(i)=elecyq(i)+(etty*qirec)
!     eleczq(i)=eleczq(i)+(ettz*qirec)
      eltmp(j)=eltmp(j)-(ettx*qjrec)
      eltmp(num+j)=eltmp(num+j)-(etty*qjrec)
      eltmp(2*num+j)=eltmp(2*num+j)-(ettz*qjrec)

      eltmp(i)=eltmp(i)+(ettx*qirec)
      eltmp(num+i)=eltmp(num+i)+(etty*qirec)
      eltmp(2*num+i)=eltmp(2*num+i)+(ettz*qirec)
!<--- Parallelization_E
!
! Permanent charge - dipole energy:
!
!---> Memmory Reduction_S
!     efact=etapi*expewld*dr+erfc(i,j)
      efact=etapi*expewld*dr+erfc(numx)
!<--- Memmory Reduction_E
!
! Permanent charge - dipole forces:
!
      efact2=2.0d0*etasq*drsqrec*etapi*expewld
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
      sxx=1.0d0-3.0d0*dxunitsq
      syy=1.0d0-3.0d0*dyunitsq
      szz=1.0d0-3.0d0*dzunitsq
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

!---> Parallelization_S
!     engeff(i)=engeff(i)+erfcr*qirec
!     engeff(j)=engeff(j)+erfcr*qjrec
      eltmp(9*num+i)=eltmp(9*num+i)+erfcr*qirec
      eltmp(9*num+j)=eltmp(9*num+j)+erfcr*qjrec
!<--- Parallelization_E
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
! Polarization damping terms
! TMP i-j term - i is defined as at the centre of the cluster.
!
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
!---> Parallelization_S
!     elecxsr(i)=elecxsr(i)+elecsrxi
!     elecysr(i)=elecysr(i)+elecsryi
!     eleczsr(i)=eleczsr(i)+elecsrzi
      eltmp(10*num+i)=eltmp(10*num+i)+elecsrxi
      eltmp(11*num+i)=eltmp(11*num+i)+elecsryi
      eltmp(12*num+i)=eltmp(12*num+i)+elecsrzi
!<--- Parallelization_E

! TMP j-i term - j now at the centre of the cluster.

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
!---> Parallelization_S
!     elecxsr(j)=elecxsr(j)+elecsrxj
!     elecysr(j)=elecysr(j)+elecsryj
!     eleczsr(j)=eleczsr(j)+elecsrzj
      eltmp(10*num+j)=eltmp(10*num+j)+elecsrxj
      eltmp(11*num+j)=eltmp(11*num+j)+elecsryj
      eltmp(12*num+j)=eltmp(12*num+j)+elecsrzj
!<--- Parallelization_E
!
! Field gradients. 
! 1) Permanent charge contribution -T_2q.
!
      txxt=txx*efact+efact2*dxx
      tyyt=tyy*efact+efact2*dyy
      tzzt=tzz*efact+efact2*dzz
      txyt=txy*efact+efact2*dxy
      txzt=txz*efact+efact2*dxz
      tyzt=tyz*efact+efact2*dyz

!---> Parallelization_S
!     exxq(j)=exxq(j)-txxt*q(i)
!     eyyq(j)=eyyq(j)-tyyt*q(i)
!     ezzq(j)=ezzq(j)-tzzt*q(i)
!     exyq(j)=exyq(j)-txyt*q(i)
!     exzq(j)=exzq(j)-txzt*q(i)
!     eyzq(j)=eyzq(j)-tyzt*q(i)
!     exxq(i)=exxq(i)-txxt*q(j)
!     eyyq(i)=eyyq(i)-tyyt*q(j)
!     ezzq(i)=ezzq(i)-tzzt*q(j)
!     exyq(i)=exyq(i)-txyt*q(j)
!     exzq(i)=exzq(i)-txzt*q(j)
!     eyzq(i)=eyzq(i)-tyzt*q(j)
      eltmp(3*num+j)=eltmp(3*num+j)-txxt*q(i)
      eltmp(4*num+j)=eltmp(4*num+j)-tyyt*q(i)
      eltmp(5*num+j)=eltmp(5*num+j)-tzzt*q(i)
      eltmp(6*num+j)=eltmp(6*num+j)-txyt*q(i)
      eltmp(7*num+j)=eltmp(7*num+j)-txzt*q(i)
      eltmp(8*num+j)=eltmp(8*num+j)-tyzt*q(i)
      eltmp(3*num+i)=eltmp(3*num+i)-txxt*q(j)
      eltmp(4*num+i)=eltmp(4*num+i)-tyyt*q(j)
      eltmp(5*num+i)=eltmp(5*num+i)-tzzt*q(j)
      eltmp(6*num+i)=eltmp(6*num+i)-txyt*q(j)
      eltmp(7*num+i)=eltmp(7*num+i)-txzt*q(j)
      eltmp(8*num+i)=eltmp(8*num+i)-tyzt*q(j)
!<--- Parallelization_E
! Field grad. SR contribution.

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

!---> Parallelization_S
!     exxsr(i)=exxsr(i)-exxsrg
!     eyysr(i)=eyysr(i)-eyysrg
!     ezzsr(i)=ezzsr(i)-ezzsrg
!     exysr(i)=exysr(i)-exysrg
!     exzsr(i)=exzsr(i)-exzsrg
!     eyzsr(i)=eyzsr(i)-eyzsrg
      eltmp(13*num+i)=eltmp(13*num+i)-exxsrg
      eltmp(14*num+i)=eltmp(14*num+i)-eyysrg
      eltmp(15*num+i)=eltmp(15*num+i)-ezzsrg
      eltmp(16*num+i)=eltmp(16*num+i)-exysrg
      eltmp(17*num+i)=eltmp(17*num+i)-exzsrg
      eltmp(18*num+i)=eltmp(18*num+i)-eyzsrg
!<--- Parallelization_E

! TMP j-i term - j now at the centre of the cluster.

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

!---> Parallelization_S
!     exxsr(j)=exxsr(j)-exxsrg
!     eyysr(j)=eyysr(j)-eyysrg
!     ezzsr(j)=ezzsr(j)-ezzsrg
!     exysr(j)=exysr(j)-exysrg
!     exzsr(j)=exzsr(j)-exzsrg
!     eyzsr(j)=eyzsr(j)-eyzsrg
      eltmp(13*num+j)=eltmp(13*num+j)-exxsrg
      eltmp(14*num+j)=eltmp(14*num+j)-eyysrg
      eltmp(15*num+j)=eltmp(15*num+j)-ezzsrg
      eltmp(16*num+j)=eltmp(16*num+j)-exysrg
      eltmp(17*num+j)=eltmp(17*num+j)-exzsrg
      eltmp(18*num+j)=eltmp(18*num+j)-eyzsrg
!<--- Parallelization_E

   enddo   

enddo   

ewcorr=facewvac*2.0d0*totmz

return
END SUBROUTINE
