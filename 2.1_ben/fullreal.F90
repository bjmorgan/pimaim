SUBROUTINE fullreal

USE commondata, ONLY: num,ntype, &
                      q,x,y,z,xmu,ymu,zmu,quadxx,quadyy,quadzz,quadxy, &
                      quadxz,quadyz,elecx,elecy,elecz,exx,eyy,ezz,exy,exz,eyz, &
                      ewlog,twopi,rsqmax,dxsav,dysav,dzsav,etapi,etasq,erfc, &
                      dddamp,dddamp2,dddamp3,dddamp4,dddamp5,dddamp6,dqdamp, &
                      dqdamp2,dqdamp3,dqdamp4,dqdamp5,dqdamp6,dqdamp7,dqdamp8, &
                      dispalpsq,dispalp4,dispalp6,ftc,ftd,onethird, &
                      srdipx,srdipy,srdipz,asdipx,asdipy,asdipz,srquadxx, &
                      srquadyy,srquadzz,srquadxy,srquadxz,srquadyz,asquadxx, &
                      asquadyy,asquadzz,asquadxy,asquadxz,asquadyz,alppolar, &
                      Cpolar,stewzz,elecxsr,elecysr,eleczsr,exxsr, &
                      eyysr,ezzsr,exysr,exzsr,eyzsr,nkdamp,dampa,dampfac,nkfg, &
! force check
                      gw_force, & 
!---> Memmory Reduction_S
!                     fgb,fgc
                      fgb,fgc,num2,numx
!<--- Memmory Reduction_E
USE boxdata, ONLY: hlab2,boxlenx,boxleny,boxlenz
!---> Optimize_S
use mpipara
!<--- Optimize_E

IMPLICIT NONE

INTEGER :: i,j,kk,ipoint,jpoint,m,n
DOUBLE PRECISION :: dxunitsq,dyunitsq,dzunitsq,frrxdipqij,frrxdipqji, &
  frrydipqij,frrydipqji,frrzdipqij,frrzdipqji,xmuidotrij, &
  xmujdotrij,egam,erfcr,expewld, &
  totmz,totmzd,totmzu,facewvac,facewvac2,qq2api
DOUBLE PRECISION :: efact,efact2,efact3,efact4,efact5,efact6, &
  xmuiTthetaj,xmujTthetai,chgquad2,chgquad2ij,chgquad2ji, &
  chgquadeng,chgquadengij,chgquadengji,dipdipeng,chgdipeng, &
  chgdipengij,chgdipengji,chgtem
DOUBLE PRECISION :: drsq,drsqrec,dr,dr3rec,dr4rec,dr5rec, &
  dr6rec,dr7rec,dr9rec,dr11rec,dx,dy,dz,dxx,dyy,dzz,dxy,dxz,dyz, &
  sxx,syy,szz,sxy,sxz,syz,threedr7,qirec,qjrec,drrec
DOUBLE PRECISION :: quadxxitquadj,quadyyitquadj,quadzzitquadj, &
  quadxyitquadj,quadxzitquadj,quadyzitquadj,txxijquadj, &
  tyyijquadj,tzzijquadj,txyijquadj,txzijquadj,tyzijquadj, &
  fqquadxij,fqquadxji,fqquadyij,fqquadyji,fqquadzij,fqquadzji
DOUBLE PRECISION :: elecsrxi,elecsryi,elecsrzi,elecsrxj,elecsryj,elecsrzj, &
  r3dampi,r3dampj,factorial,xf, &
  dampsumfi,dampsumfj,dampfunci,dampfuncj,dampaexpi,dampaexpj
DOUBLE PRECISION :: quadquadeng,quadquadeng1,quadquadeng2, &
  quadquadeng3,quadquadeng4,quadquadeng5,quadquadeng6, &
  quadquadeng7,quadquadeng6ij,quadquadeng6ji,quadquadengxx2, &
  quadquadengyy2,quadquadengzz2,quadquadengxy2,quadquadengxz2, &
  quadquadengyz2,quaddipengij,quaddipengji,quaddipeng
DOUBLE PRECISION :: qdotmuxij,qdotmuxji,qdotmuyij,qdotmuyji, &
  qdotmuzij,qdotmuzji,txijdotthetaij,txijdotthetaji, &
  tyijdotthetaij,tyijdotthetaji,tzijdotthetaij,tzijdotthetaji, &
  thetaxidotrij,thetaxjdotrij,thetayidotrij,thetayjdotrij, &
  thetazidotrij,thetazjdotrij,xmuidotT,ymuidotT,zmuidotT
DOUBLE PRECISION :: txxxx,tyyyy,tzzzz,txxxy,txxxz,txxyy,txxzz, &
  txyyy,txzzz,tyyyz,tyzzz,tyyzz,txxyz,txyyz,txyzz,txxx,tyyy, &
  tzzz,txxy,txxz,txyy,txzz,txyz,tyyz,tyzz,txx,tyy,tzz,txy,txz,tyz
DOUBLE PRECISION ::  &
  exxsrg,eyysrg,ezzsrg,exysrg,exzsrg,eyzsrg
DOUBLE PRECISION :: fqquadxsr,fqquadysr,fqquadzsr,frrxsrqij, &
  frrxsrqji,frrysrqij,frrysrqji,frrzsrqij,frrzsrqji,frrxsri, &
  frrxsrj,frrysri,frrysrj,frrzsri,frrzsrj
DOUBLE PRECISION :: txxxxx,tyyyyy,tzzzzz,txxxxy,txxxxz,txzzzz,txyyyy, &
  tyyyyz,tyzzzz,txxxyy,txxxzz,txxzzz,txxyyy,tyyyzz,tyyzzz,txxxyz, &
  txyyyz,txyzzz,txxyyz,txyyzz,txxyzz,fttx,ftty,fttz
DOUBLE PRECISION :: fquadquad1,fquadquad2a,fquadquad2b,fquadquad3, &
  fquadquad4xxx,fquadquad4yyx,fquadquad4zzx,fquadquad4xyx, &
  fquadquad4xzx,fquadquad4yzx,fquadquad4x,fquadquad4y,fquadquad4z, &
  fquadquad4xxy,fquadquad4yyy,fquadquad4zzy,fquadquad4xyy, &
  fquadquad4xzy,fquadquad4yzy,fquadquad4xzz,fquadquad4yzz, &
  fquadquad4xxz,fquadquad4yyz,fquadquad4zzz,fquadquad4xyz, &
  fquadquad5,fquadquad6x,fquadquad6y,fquadquad6z,fquadquad7, &
  fquadquad8,fquadquad9a1,fquadquad9a2,fquadquad9b1,fquadquad9b2
DOUBLE PRECISION :: fquadquad10,fquadquad11xij,fquadquad11xji, &
  fquadquad11x,fquadquad11yij,fquadquad11yji,fquadquad11y, &
  fquadquad11zij,fquadquad11zji,fquadquad11z,fquadquad12, &
  fquadquad13x,fquadquad13y,fquadquad13z,fquadquadx,fquadquady, &
  fquadquadz,fquadquadxxbyrx,fquadquadxxbyry,fquadquadxxbyrz, &
  fquadquadyybyrx,fquadquadyybyry,fquadquadyybyrz,fquadquadzzbyrx, &
  fquadquadzzbyry,fquadquadzzbyrz,fquadquadxybyrx,fquadquadxybyry, &
  fquadquadxybyrz,fquadquadxzbyrx,fquadquadxzbyry,fquadquadxzbyrz, &
  fquadquadyzbyrx,fquadquadyzbyry,fquadquadyzbyrz

DOUBLE PRECISION :: dipi,dipj,diprfunci,diprfuncj,T2dotquadi, &
  T2dotquadj,dampfuncdiffi,dampfuncdiffj,efact7,efact8, &
  fdipdipdamp,gamtot,sgam,sgam1,cpe,cpe1,dispexp,dqdampforce, &
  dqdampexp,dampsumdq,dddampforce,dddampexp,dampsum,f6,f8

double precision :: cpe2, cpe3, cpe4, sgam2, sgam3, sgam4 

DOUBLE PRECISION :: fdipquadij1,fdipquadji1,fdipquadij2,fdipquadji2, &
  fdipquadij3,fdipquadji3,fdipquadij4,fdipquadji4,fdipquadij5, &
  fdipquadji5,fdipquadij6,fdipquadji6,fdipquadij7,fdipquadji7, &
  fdipquadij8,fdipquadji8,fdipquadijx,fdipquadjix,fdipquadijy, &
  fdipquadjiy,fdipquadijz,fdipquadjiz,fdipquadx,fdipquady, &
  fdipquadz
DOUBLE PRECISION :: txxijquadi,tyyijquadi,tzzijquadi,txyijquadi, &
  txzijquadi,tyzijquadi,thetaxidotmuj,thetaxjdotmui,thetayidotmuj, &
  thetayjdotmui,thetazidotmuj,thetazjdotmui,fmumux,fmumuy,fmumuz, &
  fqquadx,fqquady,fqquadz,frrxdipq,frrydipq,frrzdipq


#ifdef dipole
dimtmp=0.0d0
#endif 

#ifdef ewald_surface
!  totmz=0.d0
!  totmzd=0.d0
!  facewvac=0.d0
!  facewvac2=0.d0
  if(ewlog) then
    facewvac=twopi/(boxlenx*boxleny*boxlenz)
    facewvac2=facewvac/(boxlenx*boxleny*boxlenz)
    totmzu=SUM(zmu)
    totmzd=SUM(q*(hlab2(3,1)*x+hlab2(3,2)*y+hlab2(3,3)*z))
    totmz=totmzd+totmzu
  endif
#endif 

!do j=2,num
do j=jst,jed
   jpoint=ntype(j)
! GWW - why not just store 1/q(ntype(j)) - i.e. 1/q for each type 
   qjrec=1.0d0/q(j)
!  do i=1,j-1
   do i=ist(j),ied(j)
      ipoint=ntype(i)
      qirec=1.0d0/q(i)

      numx = numadr(i,j)
!     drsq=dxsav(i,j)*dxsav(i,j)+dysav(i,j)*dysav(i,j)+dzsav(i,j)*dzsav(i,j)
      drsq=dxsav(numx)*dxsav(numx)+dysav(numx)*dysav(numx)+dzsav(numx)*dzsav(numx)

      if(drsq.ge.rsqmax) CYCLE

      drsqrec=1.0d0/drsq
      dr=dsqrt(drsq)
      drrec=1.0d0/dr
      dr3rec=drrec*drsqrec
      dr4rec=drsqrec*drsqrec
      dr5rec=dr4rec*drrec
      dr6rec=dr5rec*drrec
      dr7rec=dr5rec*drsqrec
      dr9rec=dr7rec*drsqrec
      dr11rec=dr9rec*drsqrec


! GWW - terms for real space charge interaction 
      chgtem=q(i)*q(j)
      qq2api=etapi*chgtem
      expewld=exp(-etasq*drsq)
!     erfcr=chgtem*erfc(i,j)*drrec
      erfcr=chgtem*erfc(numx)*drrec

! force terms for stress f/r 
      egam=(erfcr+qq2api*expewld)*drsqrec



! damping function for C6 
      dampsum=1.0d0+(dddamp(ipoint,jpoint)*dr) &
               +(dddamp2(ipoint,jpoint)*drsq) &
               +(dddamp3(ipoint,jpoint)*drsq*dr) &
               +(dddamp4(ipoint,jpoint)*drsq*drsq) &
               +(dddamp5(ipoint,jpoint)*drsq*drsq*dr) &
               +(dddamp6(ipoint,jpoint)*drsq*drsq*drsq)
      dddampexp=dexp(-dddamp(ipoint,jpoint)*dr)
      f6=1.0d0-(dddampexp*dampsum)
      dddampforce=dddamp(ipoint,jpoint)*dddamp6(ipoint,jpoint) &
                 *drsq*drsq*drsq*dddampexp

! damping function for C8
      dampsumdq=1.0d0+(dqdamp(ipoint,jpoint)*dr) &
               +(dqdamp2(ipoint,jpoint)*drsq) &
               +(dqdamp3(ipoint,jpoint)*drsq*dr) &
               +(dqdamp4(ipoint,jpoint)*drsq*drsq) &
               +(dqdamp5(ipoint,jpoint)*drsq*drsq*dr) &
               +(dqdamp6(ipoint,jpoint)*drsq*drsq*drsq) &
               +(dqdamp7(ipoint,jpoint)*drsq*drsq*drsq*dr) &
               +(dqdamp8(ipoint,jpoint)*drsq*drsq*drsq*drsq)
      dqdampexp=dexp(-dqdamp(ipoint,jpoint)*dr)
      f8=1.0d0-(dqdampexp*dampsumdq)
      dqdampforce=dqdamp(ipoint,jpoint)*dqdamp8(ipoint,jpoint) &
                 *drsq*drsq*drsq*drsq*dqdampexp


#ifdef vdw_ewald 
! GWW - Ewald real space C6 vdw
! unfortunately in setup definies the reciprical space inidividual atoms
! C6 terms using rules hene C22 is not done the same in real and reciprical.
! In addition there is not sign of damping in the reciprical space term
! rather shite really
! Will have to come back and see if I can fix the reciprical space ewald vdw
! at some poitn - but for now provide the possibility of a real space using the
! input parameters and actually using the damping !
      dispexp=exp(-dispalpsq*drsq)
      cpe1  = -ftc(ipoint,jpoint)*dispexp* &
               (dr4rec+dispalpsq*drsqrec+dispalp4*0.5d0)*drsqrec

      sgam1 = -(6.0d0*(dr6rec+dispalpsq*dr4rec+0.5d0*dispalp4* &
               drsqrec)+dispalp6)*drsqrec*dispexp*ftc(ipoint,jpoint)
! add C6 and C8 terms together, energy and force term (F/r) 
! for Ewald need to calculte difference between full and damped and take off


! full realspace C6 energy 
      cpe2  = -ftc(ipoint,jpoint)*dr6rec
! full realspace - damped energy 
      cpe4   = (1-f6)*cpe2 

! damped  energy   difference between Ewald and (realspace - damped) 
      cpe    = cpe1 - cpe4 

! and now damped C8 term 
      cpe3  = -ftd(ipoint,jpoint)*drsqrec*dr6rec
      cpe    = cpe + f8*cpe3 

! forces 

! C6 undamped realspace force  
      sgam2 = -6.0d0*cpe2*drsqrec 

! full real space - damped C6
      sgam4  =  -sgam2 +(f6 * sgam2 - cpe2*dddampforce*drrec)      

! Final C6 force term, Ewald - (realspace - damped) 
      sgam =  sgam1- sgam4 

! and now damped C8 forces 
      sgam3 = -8.0d0*cpe3*drsqrec 
      sgam = sgam -(f8 * sgam3 - cpe3*dqdampforce*drrec) 

#else 

! real space only vdw 
! add C6 and C8 terms together, energy and force term (F/r) 

! full real space C6  vdw - either use as vdw or use to calculate for 
      cpe2  = -ftc(ipoint,jpoint)*dr6rec
      cpe3  = -ftd(ipoint,jpoint)*drsqrec*dr6rec
      cpe   = f6*cpe2 + f8*cpe3

      sgam2 = -6.0d0*cpe2*drsqrec 
      sgam3 = -8.0d0*cpe3*drsqrec 

      sgam  = -(f6 * sgam2 - cpe2*dddampforce*drrec)     & 
              -(f8 * sgam3 - cpe3*dqdampforce*drrec) 
#endif 

! GWW - vdw forces - used in stress matrix
!     fttx=sgam*dxsav(i,j)
!     ftty=sgam*dysav(i,j)
!     fttz=sgam*dzsav(i,j)
      fttx=sgam*dxsav(numx)
      ftty=sgam*dysav(numx)
      fttz=sgam*dzsav(numx)

! GWW - seperate vdw force so can get total vdw - as reciprical space (when
! used) is stored in same place  
      eltmp(3*num+j)=eltmp(3*num+j)-fttx
      eltmp(4*num+j)=eltmp(4*num+j)-ftty
      eltmp(5*num+j)=eltmp(5*num+j)-fttz
      eltmp(3*num+i)=eltmp(3*num+i)+fttx
      eltmp(4*num+i)=eltmp(4*num+i)+ftty
      eltmp(5*num+i)=eltmp(5*num+i)+fttz

! GWW - vdw seperately above 
! f/r - charges + vdw 
      gamtot=egam  !  +sgam

! GWW - energy for charge, vdw and dipole (c-d, d-d) 
! moved above dipole section for c and v 
!     eng=eng+erfcr+cpe-chgdipeng+dipdipeng
      sctmp(55)=sctmp(55)+erfcr
      sctmp(56)=sctmp(56)+cpe

! GWW vdw stress 
!     stsrxx=stsrxx+fttx*dxsav(i,j)
!     stsrxy=stsrxy+fttx*dysav(i,j)
!     stsrxz=stsrxz+fttx*dzsav(i,j)
!     stsryy=stsryy+ftty*dysav(i,j)
!     stsryz=stsryz+ftty*dzsav(i,j)
!     stsrzz=stsrzz+fttz*dzsav(i,j)
      sctmp(1)=sctmp(1)+fttx*dxsav(numx)
      sctmp(2)=sctmp(2)+fttx*dysav(numx)
      sctmp(3)=sctmp(3)+fttx*dzsav(numx)
      sctmp(4)=sctmp(4)+ftty*dysav(numx)
      sctmp(5)=sctmp(5)+ftty*dzsav(numx)
      sctmp(6)=sctmp(6)+fttz*dzsav(numx)

! GWW - charge -charge stress
!     stcxx=stcxx+egam*dxsav(i,j)*dxsav(i,j)
!     stcxy=stcxy+egam*dxsav(i,j)*dysav(i,j)
!     stcxz=stcxz+egam*dxsav(i,j)*dzsav(i,j)
!     stcyy=stcyy+egam*dysav(i,j)*dysav(i,j)
!     stcyz=stcyz+egam*dysav(i,j)*dzsav(i,j)
!     stczz=stczz+egam*dzsav(i,j)*dzsav(i,j)
      sctmp(7) =sctmp(7) +egam*dxsav(numx)*dxsav(numx)
      sctmp(8) =sctmp(8) +egam*dxsav(numx)*dysav(numx)
      sctmp(9) =sctmp(9) +egam*dxsav(numx)*dzsav(numx)
      sctmp(10)=sctmp(10)+egam*dysav(numx)*dysav(numx)
      sctmp(11)=sctmp(11)+egam*dysav(numx)*dzsav(numx)
      sctmp(12)=sctmp(12)+egam*dzsav(numx)*dzsav(numx)


#ifdef dipole 
!     xmuidotrij=xmu(i)*dxsav(i,j)+ymu(i)*dysav(i,j)+zmu(i)*dzsav(i,j)
!     xmujdotrij=xmu(j)*dxsav(i,j)+ymu(j)*dysav(i,j)+zmu(j)*dzsav(i,j)
      xmuidotrij=xmu(i)*dxsav(numx)+ymu(i)*dysav(numx)+zmu(i)*dzsav(numx)
      xmujdotrij=xmu(j)*dxsav(numx)+ymu(j)*dysav(numx)+zmu(j)*dzsav(numx)

!     efact=etapi*expewld*dr+erfc(i,j)
      efact=etapi*expewld*dr+erfc(numx)

      chgdipengij=(q(j)*xmuidotrij)*dr3rec
      chgdipengji=-(q(i)*xmujdotrij)*dr3rec
      chgdipeng=(chgdipengij+chgdipengji)*efact

      qdotmuxij=q(j)*xmu(i)
      qdotmuxji=-q(i)*xmu(j)
      qdotmuyij=q(j)*ymu(i)
      qdotmuyji=-q(i)*ymu(j)
      qdotmuzij=q(j)*zmu(i)
      qdotmuzji=-q(i)*zmu(j)
      efact2=2.0d0*etasq*drsqrec*etapi*expewld

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

      frrxdipq=frrxdipqij+frrxdipqji
      frrydipq=frrydipqij+frrydipqji
      frrzdipq=frrzdipqij+frrzdipqji

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


      xmuidotT=xmu(i)*sxx+ymu(i)*sxy+zmu(i)*sxz
      ymuidotT=xmu(i)*sxy+ymu(i)*syy+zmu(i)*syz
      zmuidotT=xmu(i)*sxz+ymu(i)*syz+zmu(i)*szz

      dipdipeng=xmuidotT*xmu(j)+ymuidotT*ymu(j)+zmuidotT*zmu(j)

      fdipdipdamp=dipdipeng
      dipdipeng=dipdipeng*efact*dr3rec-xmuidotrij*xmujdotrij*efact2

      threedr7=3.0d0*dr7rec

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
 
      fmumux=(xmu(j)*(xmu(i)*txxx+ymu(i)*txxy+zmu(i)*txxz)) &
            +(ymu(j)*(xmu(i)*txxy+ymu(i)*txyy+zmu(i)*txyz)) &
            +(zmu(j)*(xmu(i)*txxz+ymu(i)*txyz+zmu(i)*txzz))
      fmumuy=(xmu(j)*(xmu(i)*txxy+ymu(i)*txyy+zmu(i)*txyz)) &
            +(ymu(j)*(xmu(i)*txyy+ymu(i)*tyyy+zmu(i)*tyyz)) &
            +(zmu(j)*(xmu(i)*txyz+ymu(i)*tyyz+zmu(i)*tyzz))
      fmumuz=(xmu(j)*(xmu(i)*txxz+ymu(i)*txyz+zmu(i)*txzz)) &
            +(ymu(j)*(xmu(i)*txyz+ymu(i)*tyyz+zmu(i)*tyzz)) &
            +(zmu(j)*(xmu(i)*txzz+ymu(i)*tyzz+zmu(i)*tzzz))

      efact3= 2.0d0*efact2*(drsqrec+etasq)

!     fmumux = -fmumux*efact+fdipdipdamp*efact2*dxsav(i,j) &
!              +xmu(i)*xmujdotrij*efact2+xmu(j)*xmuidotrij*efact2 &
!              -xmuidotrij*xmujdotrij*dxsav(i,j)*efact3
!     fmumuy = -fmumuy*efact+fdipdipdamp*efact2*dysav(i,j) &
!              +ymu(i)*xmujdotrij*efact2+ymu(j)*xmuidotrij*efact2 &
!              -xmuidotrij*xmujdotrij*dysav(i,j)*efact3
!     fmumuz = -fmumuz*efact+fdipdipdamp*efact2*dzsav(i,j) &
!              +zmu(i)*xmujdotrij*efact2+zmu(j)*xmuidotrij*efact2 &
!              -xmuidotrij*xmujdotrij*dzsav(i,j)*efact3
      fmumux = -fmumux*efact+fdipdipdamp*efact2*dxsav(numx) &
               +xmu(i)*xmujdotrij*efact2+xmu(j)*xmuidotrij*efact2 &
               -xmuidotrij*xmujdotrij*dxsav(numx)*efact3
      fmumuy = -fmumuy*efact+fdipdipdamp*efact2*dysav(numx) &
               +ymu(i)*xmujdotrij*efact2+ymu(j)*xmuidotrij*efact2 &
               -xmuidotrij*xmujdotrij*dysav(numx)*efact3
      fmumuz = -fmumuz*efact+fdipdipdamp*efact2*dzsav(numx) &
               +zmu(i)*xmujdotrij*efact2+zmu(j)*xmuidotrij*efact2 &
               -xmuidotrij*xmujdotrij*dzsav(numx)*efact3

! GWW - energy for charge, vdw and dipole (c-d, d-d) 
! charge and dipole move to before dipole section 
!     eng=eng+erfcr+cpe-chgdipeng+dipdipeng
!      sctmp(55)=sctmp(55)+erfcr
!      sctmp(56)=sctmp(56)+cpe
      sctmp(57)=sctmp(57)-chgdipeng+dipdipeng

!     stpxx=stpxx+frrxdipq*dxsav(i,j)
!     stpxy=stpxy+0.5d0*(frrxdipq*dysav(i,j)+frrydipq*dxsav(i,j))
!     stpxz=stpxz+0.5d0*(frrxdipq*dzsav(i,j)+frrzdipq*dxsav(i,j))
!     stpyy=stpyy+frrydipq*dysav(i,j)
!     stpyz=stpyz+0.5d0*(frrydipq*dzsav(i,j)+frrzdipq*dysav(i,j))
!     stpzz=stpzz+frrzdipq*dzsav(i,j)
      sctmp(13)=sctmp(13)+frrxdipq*dxsav(numx)
      sctmp(14)=sctmp(14)+0.5d0*(frrxdipq*dysav(numx)+frrydipq*dxsav(numx))
      sctmp(15)=sctmp(15)+0.5d0*(frrxdipq*dzsav(numx)+frrzdipq*dxsav(numx))
      sctmp(16)=sctmp(16)+frrydipq*dysav(numx)
      sctmp(17)=sctmp(17)+0.5d0*(frrydipq*dzsav(numx)+frrzdipq*dysav(numx))
      sctmp(18)=sctmp(18)+frrzdipq*dzsav(numx)

!     stp2xx=stp2xx+fmumux*dxsav(i,j)
!     stp2xy=stp2xy+0.5d0*(fmumux*dysav(i,j)+fmumuy*dxsav(i,j))
!     stp2xz=stp2xz+0.5d0*(fmumux*dzsav(i,j)+fmumuz*dxsav(i,j))
!     stp2yy=stp2yy+fmumuy*dysav(i,j)
!     stp2yz=stp2yz+0.5d0*(fmumuy*dzsav(i,j)+fmumuz*dysav(i,j))
!     stp2zz=stp2zz+fmumuz*dzsav(i,j)
      sctmp(19)=sctmp(19)+fmumux*dxsav(numx)
      sctmp(20)=sctmp(20)+0.5d0*(fmumux*dysav(numx)+fmumuy*dxsav(numx))
      sctmp(21)=sctmp(21)+0.5d0*(fmumux*dzsav(numx)+fmumuz*dxsav(numx))
      sctmp(22)=sctmp(22)+fmumuy*dysav(numx)
      sctmp(23)=sctmp(23)+0.5d0*(fmumuy*dzsav(numx)+fmumuz*dysav(numx))
      sctmp(24)=sctmp(24)+fmumuz*dzsav(numx)
#endif 

! GWW - changed to not contain vdw 
! forces for charge and vdw  
!     fttx=gamtot*dxsav(i,j)+frrxdipq+fmumux
!     ftty=gamtot*dysav(i,j)+frrydipq+fmumuy
!     fttz=gamtot*dzsav(i,j)+frrzdipq+fmumuz
      fttx=gamtot*dxsav(numx)
      ftty=gamtot*dysav(numx)
      fttz=gamtot*dzsav(numx)

! check forces 
      gw_force(j)=gw_force(j)-fttx
      gw_force(num+j)=gw_force(num+j)-ftty
      gw_force(2*num+j)=gw_force(2*num+j)-fttz
      gw_force(i)=gw_force(i)+fttx
      gw_force(num+i)=gw_force(num+i)+ftty
      gw_force(2*num+i)=gw_force(2*num+i)+fttz

#ifdef dipole
! check q-mu forces 
      gw_force(3*num+j)=gw_force(3*num+j)-frrxdipq
      gw_force(4*num+j)=gw_force(4*num+j)-frrydipq
      gw_force(5*num+j)=gw_force(5*num+j)-frrzdipq
      gw_force(3*num+i)=gw_force(3*num+i)+frrxdipq
      gw_force(4*num+i)=gw_force(4*num+i)+frrydipq
      gw_force(5*num+i)=gw_force(5*num+i)+frrzdipq

! check mu-mu forces 
      gw_force(6*num+j)=gw_force(6*num+j)-fmumux       
      gw_force(7*num+j)=gw_force(7*num+j)-fmumuy       
      gw_force(8*num+j)=gw_force(8*num+j)-fmumuz       
      gw_force(6*num+i)=gw_force(6*num+i)+fmumux       
      gw_force(7*num+i)=gw_force(7*num+i)+fmumuy       
      gw_force(8*num+i)=gw_force(8*num+i)+fmumuz       

! Gww - add dipole to forces 
      fttx=fttx+frrxdipq+fmumux
      ftty=ftty+frrydipq+fmumuy
      fttz=fttz+frrzdipq+fmumuz
#endif 

! GWW add forces to atomic running totals 
!     frrx(j)=frrx(j)-fttx
!     frry(j)=frry(j)-ftty
!     frrz(j)=frrz(j)-fttz
!     frrx(i)=frrx(i)+fttx
!     frry(i)=frry(i)+ftty
!     frrz(i)=frrz(i)+fttz
      eltmp(j)=eltmp(j)-fttx
      eltmp(num+j)=eltmp(num+j)-ftty
      eltmp(2*num+j)=eltmp(2*num+j)-fttz
      eltmp(i)=eltmp(i)+fttx
      eltmp(num+i)=eltmp(num+i)+ftty
      eltmp(2*num+i)=eltmp(2*num+i)+fttz


#ifdef quadrupole 
!     dx=dxsav(i,j)
!     dy=dysav(i,j)
!     dz=dzsav(i,j)
      dx=dxsav(numx)
      dy=dysav(numx)
      dz=dzsav(numx)
      dxx=dx*dx
      dyy=dy*dy
      dzz=dz*dz
      dxy=dx*dy
      dxz=dx*dz
      dyz=dy*dz

      txx=-sxx*dr3rec
      tyy=-syy*dr3rec
      tzz=-szz*dr3rec
      txy=-sxy*dr3rec
      txz=-sxz*dr3rec
      tyz=-syz*dr3rec

      chgquadengij=((txx*quadxx(i)+txy*quadxy(i)+txy*quadxy(i) &
                    +txz*quadxz(i)+txz*quadxz(i)+tyy*quadyy(i) &
                    +tyz*quadyz(i)+tzz*quadzz(i)+tyz*quadyz(i))*q(j))
      chgquadengji=((txx*quadxx(j)+txy*quadxy(j)+txy*quadxy(j) &
                    +txz*quadxz(j)+txz*quadxz(j)+tyy*quadyy(j) &
                    +tyz*quadyz(j)+tzz*quadzz(j)+tyz*quadyz(j))*q(i))

      chgquadeng=chgquadengij+chgquadengji

      chgquad2ij=((dxx*quadxx(i)+dxy*quadxy(i)+dxy*quadxy(i) &
                  +dxz*quadxz(i)+dxz*quadxz(i)+dyy*quadyy(i) &
                  +dyz*quadyz(i)+dzz*quadzz(i)+dyz*quadyz(i))*q(j))
      chgquad2ji=((dxx*quadxx(j)+dxy*quadxy(j)+dxy*quadxy(j) &
                  +dxz*quadxz(j)+dxz*quadxz(j)+dyy*quadyy(j) &
                  +dyz*quadyz(j)+dzz*quadzz(j)+dyz*quadyz(j))*q(i))

      chgquad2=chgquad2ij+chgquad2ji

      chgquadeng=(chgquadeng*efact)+chgquad2*efact2
      chgquadeng=chgquadeng*onethird

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

      txijdotthetaij=fqquadxji*qirec
      txijdotthetaji=fqquadxij*qjrec
      tyijdotthetaij=fqquadyji*qirec
      tyijdotthetaji=fqquadyij*qjrec
      tzijdotthetaij=fqquadzji*qirec
      tzijdotthetaji=fqquadzij*qjrec            

      efact4=efact2*dr*drsq

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

      fqquadxij=fqquadxij-2.0d0*efact2*q(j)*thetaxidotrij
      fqquadyij=fqquadyij-2.0d0*efact2*q(j)*thetayidotrij
      fqquadzij=fqquadzij-2.0d0*efact2*q(j)*thetazidotrij
      fqquadxji=fqquadxji-2.0d0*efact2*q(i)*thetaxjdotrij
      fqquadyji=fqquadyji-2.0d0*efact2*q(i)*thetayjdotrij
      fqquadzji=fqquadzji-2.0d0*efact2*q(i)*thetazjdotrij

      fqquadxij=fqquadxij*onethird
      fqquadyij=fqquadyij*onethird
      fqquadzij=fqquadzij*onethird
      fqquadxji=fqquadxji*onethird
      fqquadyji=fqquadyji*onethird
      fqquadzji=fqquadzji*onethird

      fqquadx=fqquadxij+fqquadxji
      fqquady=fqquadyij+fqquadyji
      fqquadz=fqquadzij+fqquadzji

!     stpqquadxx=stpqquadxx+fqquadx*dxsav(i,j)
!     stpqquadxy=stpqquadxy+0.5d0*(fqquadx*dysav(i,j)+fqquady*dxsav(i,j))
!     stpqquadxz=stpqquadxz+0.5d0*(fqquadx*dzsav(i,j)+fqquadz*dxsav(i,j))
!     stpqquadyy=stpqquadyy+fqquady*dysav(i,j)
!     stpqquadyz=stpqquadyz+0.5d0*(fqquadz*dysav(i,j)+fqquady*dzsav(i,j))
!     stpqquadzz=stpqquadzz+fqquadz*dzsav(i,j)
      sctmp(25)=sctmp(25)+fqquadx*dxsav(numx)
      sctmp(26)=sctmp(26)+0.5d0*(fqquadx*dysav(numx)+fqquady*dxsav(numx))
      sctmp(27)=sctmp(27)+0.5d0*(fqquadx*dzsav(numx)+fqquadz*dxsav(numx))
      sctmp(28)=sctmp(28)+fqquady*dysav(numx)
      sctmp(29)=sctmp(29)+0.5d0*(fqquadz*dysav(numx)+fqquady*dzsav(numx))
      sctmp(30)=sctmp(30)+fqquadz*dzsav(numx)

      xmuiTthetaj=xmu(i)*(txxx*quadxx(j)+(2.0d0*txxy*quadxy(j)) &
                +(2.0d0*txxz*quadxz(j))+(txyy*quadyy(j))+(txzz*quadzz(j)) &
                +(2.0d0*txyz*quadyz(j)))+ymu(i)*(txxy*quadxx(j) &
                +(2.0d0*txyy*quadxy(j))+(2.0d0*txyz*quadxz(j)) &
                +(tyyy*quadyy(j))+(tyzz*quadzz(j))+(2.0d0*tyyz*quadyz(j))) &
                +zmu(i)*(txxz*quadxx(j)+(2.0d0*txyz*quadxy(j)) &
                +(2.0d0*txzz*quadxz(j))+(tyyz*quadyy(j))+(tzzz*quadzz(j)) &
                +(2.0d0*tyzz*quadyz(j)))

      xmujTthetai=xmu(j)*(txxx*quadxx(i)+(2.0d0*txxy*quadxy(i)) &
                +(2.0d0*txxz*quadxz(i))+(txyy*quadyy(i))+(txzz*quadzz(i)) &
                +(2.0d0*txyz*quadyz(i)))+ymu(j)*(txxy*quadxx(i) &
                +(2.0d0*txyy*quadxy(i))+(2.0d0*txyz*quadxz(i)) &
                +(tyyy*quadyy(i))+(tyzz*quadzz(i))+(2.0d0*tyyz*quadyz(i))) &
                +zmu(j)*(txxz*quadxx(i)+(2.0d0*txyz*quadxy(i)) &
                +(2.0d0*txzz*quadxz(i))+(tyyz*quadyy(i))+(tzzz*quadzz(i)) &
                +(2.0d0*tyzz*quadyz(i)))

      quaddipengij=xmuiTthetaj*efact &
                  +efact3*xmuidotrij*(chgquad2ji*qirec) &
                  -2.0d0*efact2*(xmu(i)*thetaxjdotrij+ &
                                 ymu(i)*thetayjdotrij+ &
                                 zmu(i)*thetazjdotrij) &
                  +efact4*xmuidotrij*(chgquadengji*qirec)
      quaddipengji=xmujTthetai*efact &
                  +efact3*xmujdotrij*(chgquad2ij*qjrec) &
                  -2.0d0*efact2*(xmu(j)*thetaxidotrij+ &
                                 ymu(j)*thetayidotrij+ &
                                 zmu(j)*thetazidotrij) &
                  +efact4*xmujdotrij*(chgquadengij*qjrec)

      quaddipeng=(quaddipengij-quaddipengji)*onethird

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

      fdipquadijx=xmu(i)*(txxxx*quadxx(j)+(2.0d0*txxxy*quadxy(j)) &
                  +(2.0d0*txxxz*quadxz(j))+(txxyy*quadyy(j)) &
                  +(txxzz*quadzz(j))+(2.0d0*txxyz*quadyz(j))) &
          +ymu(i)*((txxxy*quadxx(j))+(2.0d0*txxyy*quadxy(j)) &
                  +(2.0d0*txxyz*quadxz(j))+txyyy*quadyy(j) &
                  +(txyzz*quadzz(j))+(2.0d0*txyyz*quadyz(j))) &
          +zmu(i)*((txxxz*quadxx(j))+(2.0d0*txxyz*quadxy(j)) &
                  +(2.0d0*txxzz*quadxz(j))+(txyyz*quadyy(j)) &
                         +txzzz*quadzz(j)+(2.0d0*txyzz*quadyz(j)))
      fdipquadjix=xmu(j)*(txxxx*quadxx(i)+(2.0d0*txxxy*quadxy(i)) &
                  +(2.0d0*txxxz*quadxz(i))+(txxyy*quadyy(i)) &
                  +(txxzz*quadzz(i))+(2.0d0*txxyz*quadyz(i))) &
          +ymu(j)*((txxxy*quadxx(i))+(2.0d0*txxyy*quadxy(i)) &
                  +(2.0d0*txxyz*quadxz(i))+txyyy*quadyy(i) &
                  +(txyzz*quadzz(i))+(2.0d0*txyyz*quadyz(i))) &
          +zmu(j)*((txxxz*quadxx(i))+(2.0d0*txxyz*quadxy(i)) &
                  +(2.0d0*txxzz*quadxz(i))+(txyyz*quadyy(i)) &
                         +txzzz*quadzz(i)+(2.0d0*txyzz*quadyz(i)))
      fdipquadijy=xmu(i)*(txxxy*quadxx(j)+(2.0d0*txxyy*quadxy(j)) &
                  +(2.0d0*txxyz*quadxz(j))+(txyyy*quadyy(j)) &
                  +(txyzz*quadzz(j))+(2.0d0*txyyz*quadyz(j))) &
          +ymu(i)*((txxyy*quadxx(j))+(2.0d0*txyyy*quadxy(j)) &
                  +(2.0d0*txyyz*quadxz(j))+tyyyy*quadyy(j) &
                  +(tyyzz*quadzz(j))+(2.0d0*tyyyz*quadyz(j))) &
          +zmu(i)*((txxyz*quadxx(j))+(2.0d0*txyyz*quadxy(j)) &
                  +(2.0d0*txyzz*quadxz(j))+(tyyyz*quadyy(j)) &
                         +tyzzz*quadzz(j)+(2.0d0*tyyzz*quadyz(j)))
      fdipquadjiy=xmu(j)*(txxxy*quadxx(i)+(2.0d0*txxyy*quadxy(i)) &
                  +(2.0d0*txxyz*quadxz(i))+(txyyy*quadyy(i)) &
                  +(txyzz*quadzz(i))+(2.0d0*txyyz*quadyz(i))) &
          +ymu(j)*((txxyy*quadxx(i))+(2.0d0*txyyy*quadxy(i)) &
                  +(2.0d0*txyyz*quadxz(i))+tyyyy*quadyy(i) &
                  +(tyyzz*quadzz(i))+(2.0d0*tyyyz*quadyz(i))) &
          +zmu(j)*((txxyz*quadxx(i))+(2.0d0*txyyz*quadxy(i)) &
                  +(2.0d0*txyzz*quadxz(i))+(tyyyz*quadyy(i)) &
                         +tyzzz*quadzz(i)+(2.0d0*tyyzz*quadyz(i)))
      fdipquadijz=xmu(i)*(txxxz*quadxx(j)+(2.0d0*txxyz*quadxy(j)) &
                  +(2.0d0*txxzz*quadxz(j))+(txyyz*quadyy(j)) &
                  +(txzzz*quadzz(j))+(2.0d0*txyzz*quadyz(j))) &
          +ymu(i)*((txxyz*quadxx(j))+(2.0d0*txyyz*quadxy(j)) &
                  +(2.0d0*txyzz*quadxz(j))+tyyyz*quadyy(j) &
                  +(tyzzz*quadzz(j))+(2.0d0*tyyzz*quadyz(j))) &
          +zmu(i)*((txxzz*quadxx(j))+(2.0d0*txyzz*quadxy(j)) &
                  +(2.0d0*txzzz*quadxz(j))+(tyyzz*quadyy(j)) &
                         +tzzzz*quadzz(j)+(2.0d0*tyzzz*quadyz(j)))
      fdipquadjiz=xmu(j)*(txxxz*quadxx(i)+(2.0d0*txxyz*quadxy(i)) &
                  +(2.0d0*txxzz*quadxz(i))+(txyyz*quadyy(i)) &
                  +(txzzz*quadzz(i))+(2.0d0*txyzz*quadyz(i))) &
          +ymu(j)*((txxyz*quadxx(i))+(2.0d0*txyyz*quadxy(i)) &
                  +(2.0d0*txyzz*quadxz(i))+tyyyz*quadyy(i) &
                  +(tyzzz*quadzz(i))+(2.0d0*tyyzz*quadyz(i))) &
          +zmu(j)*((txxzz*quadxx(i))+(2.0d0*txyzz*quadxy(i)) &
                  +(2.0d0*txzzz*quadxz(i))+(tyyzz*quadyy(i)) &
                         +tzzzz*quadzz(i)+(2.0d0*tyzzz*quadyz(i)))

      efact5=2.0d0*efact3*(etasq+drsqrec)+4.0d0*efact2*dr4rec
      efact6=dr*(3.0d0*efact2-drsq*efact3)
      thetaxidotmuj=quadxx(i)*xmu(j)+quadxy(i)*ymu(j)+quadxz(i)*zmu(j)
      thetaxjdotmui=quadxx(j)*xmu(i)+quadxy(j)*ymu(i)+quadxz(j)*zmu(i)
      thetayidotmuj=quadxy(i)*xmu(j)+quadyy(i)*ymu(j)+quadyz(i)*zmu(j)
      thetayjdotmui=quadxy(j)*xmu(i)+quadyy(j)*ymu(i)+quadyz(j)*zmu(i)
      thetazidotmuj=quadxz(i)*xmu(j)+quadyz(i)*ymu(j)+quadzz(i)*zmu(j)
      thetazjdotmui=quadxz(j)*xmu(i)+quadyz(j)*ymu(i)+quadzz(j)*zmu(i)

      fdipquadij1=efact5*xmuidotrij*(chgquad2ji*qirec)
      fdipquadji1=efact5*xmujdotrij*(chgquad2ij*qjrec)
      fdipquadij2=efact3*chgquad2ji*qirec
      fdipquadji2=efact3*chgquad2ij*qjrec
      fdipquadij3=2.0d0*efact3*xmuidotrij
      fdipquadji3=2.0d0*efact3*xmujdotrij
      fdipquadij4=2.0d0*efact3*(xmu(i)*thetaxjdotrij+ &
                                ymu(i)*thetayjdotrij+zmu(i)*thetazjdotrij)
      fdipquadji4=2.0d0*efact3*(xmu(j)*thetaxidotrij+ &
                                ymu(j)*thetayidotrij+zmu(j)*thetazidotrij)
      fdipquadij5=efact6*xmuidotrij*(chgquadengji*qirec)
      fdipquadji5=efact6*xmujdotrij*(chgquadengij*qjrec)
      fdipquadij6=efact4*(chgquadengji*qirec)
      fdipquadji6=efact4*(chgquadengij*qjrec)
      fdipquadij7=efact4*xmuidotrij
      fdipquadji7=efact4*xmujdotrij
      fdipquadij8=efact4*xmuiTthetaj
      fdipquadji8=efact4*xmujTthetai

!     fdipquadijx=fdipquadijx*efact+fdipquadij1*dxsav(i,j) &
!                -fdipquadij2*xmu(i)-fdipquadij3*thetaxjdotrij &
!                -fdipquadij4*dxsav(i,j)+2.0d0*efact2*thetaxjdotmui &
!                -fdipquadij5*dxsav(i,j)-fdipquadij6*xmu(i) &
!                +fdipquadij7*txijdotthetaij+fdipquadij8*dxsav(i,j)
!     fdipquadjix=fdipquadjix*efact+fdipquadji1*dxsav(i,j) &
!                -fdipquadji2*xmu(j)-fdipquadji3*thetaxidotrij &
!                -fdipquadji4*dxsav(i,j)+2.0d0*efact2*thetaxidotmuj &
!                -fdipquadji5*dxsav(i,j)-fdipquadji6*xmu(j) &
!                +fdipquadji7*txijdotthetaji+fdipquadji8*dxsav(i,j)
!     fdipquadijy=fdipquadijy*efact+fdipquadij1*dysav(i,j) &
!                -fdipquadij2*ymu(i)-fdipquadij3*thetayjdotrij &
!                -fdipquadij4*dysav(i,j)+2.0d0*efact2*thetayjdotmui &
!                -fdipquadij5*dysav(i,j)-fdipquadij6*ymu(i) &
!                +fdipquadij7*tyijdotthetaij+fdipquadij8*dysav(i,j)
!     fdipquadjiy=fdipquadjiy*efact+fdipquadji1*dysav(i,j) &
!                -fdipquadji2*ymu(j)-fdipquadji3*thetayidotrij &
!                -fdipquadji4*dysav(i,j)+2.0d0*efact2*thetayidotmuj &
!                -fdipquadji5*dysav(i,j)-fdipquadji6*ymu(j) &
!                +fdipquadji7*tyijdotthetaji+fdipquadji8*dysav(i,j)
!     fdipquadijz=fdipquadijz*efact+fdipquadij1*dzsav(i,j) &
!                -fdipquadij2*zmu(i)-fdipquadij3*thetazjdotrij &
!                -fdipquadij4*dzsav(i,j)+2.0d0*efact2*thetazjdotmui &
!                -fdipquadij5*dzsav(i,j)-fdipquadij6*zmu(i) &
!                +fdipquadij7*tzijdotthetaij+fdipquadij8*dzsav(i,j)
!     fdipquadjiz=fdipquadjiz*efact+fdipquadji1*dzsav(i,j) &
!                -fdipquadji2*zmu(j)-fdipquadji3*thetazidotrij &
!                -fdipquadji4*dzsav(i,j)+2.0d0*efact2*thetazidotmuj &
!                -fdipquadji5*dzsav(i,j)-fdipquadji6*zmu(j) &
!                +fdipquadji7*tzijdotthetaji+fdipquadji8*dzsav(i,j)
      fdipquadijx=fdipquadijx*efact+fdipquadij1*dxsav(numx) &
                 -fdipquadij2*xmu(i)-fdipquadij3*thetaxjdotrij &
                 -fdipquadij4*dxsav(numx)+2.0d0*efact2*thetaxjdotmui &
                 -fdipquadij5*dxsav(numx)-fdipquadij6*xmu(i) &
                 +fdipquadij7*txijdotthetaij+fdipquadij8*dxsav(numx)
      fdipquadjix=fdipquadjix*efact+fdipquadji1*dxsav(numx) &
                 -fdipquadji2*xmu(j)-fdipquadji3*thetaxidotrij &
                 -fdipquadji4*dxsav(numx)+2.0d0*efact2*thetaxidotmuj &
                 -fdipquadji5*dxsav(numx)-fdipquadji6*xmu(j) &
                 +fdipquadji7*txijdotthetaji+fdipquadji8*dxsav(numx)
      fdipquadijy=fdipquadijy*efact+fdipquadij1*dysav(numx) &
                 -fdipquadij2*ymu(i)-fdipquadij3*thetayjdotrij &
                 -fdipquadij4*dysav(numx)+2.0d0*efact2*thetayjdotmui &
                 -fdipquadij5*dysav(numx)-fdipquadij6*ymu(i) &
                 +fdipquadij7*tyijdotthetaij+fdipquadij8*dysav(numx)
      fdipquadjiy=fdipquadjiy*efact+fdipquadji1*dysav(numx) &
                 -fdipquadji2*ymu(j)-fdipquadji3*thetayidotrij &
                 -fdipquadji4*dysav(numx)+2.0d0*efact2*thetayidotmuj &
                 -fdipquadji5*dysav(numx)-fdipquadji6*ymu(j) &
                 +fdipquadji7*tyijdotthetaji+fdipquadji8*dysav(numx)
      fdipquadijz=fdipquadijz*efact+fdipquadij1*dzsav(numx) &
                 -fdipquadij2*zmu(i)-fdipquadij3*thetazjdotrij &
                 -fdipquadij4*dzsav(numx)+2.0d0*efact2*thetazjdotmui &
                 -fdipquadij5*dzsav(numx)-fdipquadij6*zmu(i) &
                 +fdipquadij7*tzijdotthetaij+fdipquadij8*dzsav(numx)
      fdipquadjiz=fdipquadjiz*efact+fdipquadji1*dzsav(numx) &
                 -fdipquadji2*zmu(j)-fdipquadji3*thetazidotrij &
                 -fdipquadji4*dzsav(numx)+2.0d0*efact2*thetazidotmuj &
                 -fdipquadji5*dzsav(numx)-fdipquadji6*zmu(j) &
                 +fdipquadji7*tzijdotthetaji+fdipquadji8*dzsav(numx)

      fdipquadx=onethird*(fdipquadijx-fdipquadjix)
      fdipquady=onethird*(fdipquadijy-fdipquadjiy)
      fdipquadz=onethird*(fdipquadijz-fdipquadjiz)

!     stpdipquadxx=stpdipquadxx-fdipquadx*dxsav(i,j)
!     stpdipquadxy=stpdipquadxy-0.5d0*(fdipquadx*dysav(i,j) &
!                       +fdipquady*dxsav(i,j))
!     stpdipquadxz=stpdipquadxz-0.5d0*(fdipquadx*dzsav(i,j) &
!                       +fdipquadz*dxsav(i,j))
!     stpdipquadyy=stpdipquadyy-fdipquady*dysav(i,j)
!     stpdipquadyz=stpdipquadyz-0.5d0*(fdipquadz*dysav(i,j) &
!                       +fdipquady*dzsav(i,j))
!     stpdipquadzz=stpdipquadzz-fdipquadz*dzsav(i,j)
      sctmp(31)=sctmp(31)-fdipquadx*dxsav(numx)
      sctmp(32)=sctmp(32)-0.5d0*(fdipquadx*dysav(numx) &
                        +fdipquady*dxsav(numx))
      sctmp(33)=sctmp(33)-0.5d0*(fdipquadx*dzsav(numx) &
                        +fdipquadz*dxsav(numx))
      sctmp(34)=sctmp(34)-fdipquady*dysav(numx)
      sctmp(35)=sctmp(35)-0.5d0*(fdipquadz*dysav(numx) &
                        +fdipquady*dzsav(numx))
      sctmp(36)=sctmp(36)-fdipquadz*dzsav(numx)

      txxijquadj=txxxx*quadxx(j)+(2.0d0*txxxy*quadxy(j)) &
                +(2.0d0*txxxz*quadxz(j))+(2.0d0*txxyz*quadyz(j)) &
                +txxyy*quadyy(j)+txxzz*quadzz(j)
      tyyijquadj=txxyy*quadxx(j)+(2.0d0*txyyy*quadxy(j)) &
                +(2.0d0*txyyz*quadxz(j))+(2.0d0*tyyyz*quadyz(j)) &
                +tyyyy*quadyy(j)+tyyzz*quadzz(j)
      tzzijquadj=txxzz*quadxx(j)+(2.0d0*txyzz*quadxy(j)) &
                +(2.0d0*txzzz*quadxz(j))+(2.0d0*tyzzz*quadyz(j)) &
                +tyyzz*quadyy(j)+tzzzz*quadzz(j)
      txyijquadj=txxxy*quadxx(j)+(2.0d0*txxyy*quadxy(j)) &
                +(2.0d0*txxyz*quadxz(j))+txyyy*quadyy(j) &
                +(2.0d0*txyyz*quadyz(j))+txyzz*quadzz(j) 
      txzijquadj=txxxz*quadxx(j)+(2.0d0*txxyz*quadxy(j)) &
                +(2.0d0*txxzz*quadxz(j))+txyyz*quadyy(j) &
                +(2.0d0*txyzz*quadyz(j))+txzzz*quadzz(j)
      tyzijquadj=txxyz*quadxx(j)+(2.0d0*txyyz*quadxy(j)) &
                +(2.0d0*txyzz*quadxz(j))+tyyyz*quadyy(j) &
                +(2.0d0*tyyzz*quadyz(j))+tyzzz*quadzz(j)

      quadxxiTquadj=quadxx(i)*txxijquadj
      quadyyiTquadj=quadyy(i)*tyyijquadj
      quadzziTquadj=quadzz(i)*tzzijquadj
      quadxyiTquadj=quadxy(i)*txyijquadj
      quadxziTquadj=quadxz(i)*txzijquadj
      quadyziTquadj=quadyz(i)*tyzijquadj

      txxijquadi=txxxx*quadxx(i)+(2.0d0*txxxy*quadxy(i)) &
                 +(2.0d0*txxxz*quadxz(i))+(2.0d0*txxyz*quadyz(i)) &
                 +txxyy*quadyy(i)+txxzz*quadzz(i)
      tyyijquadi=txxyy*quadxx(i)+(2.0d0*txyyy*quadxy(i)) &
                 +(2.0d0*txyyz*quadxz(i))+(2.0d0*tyyyz*quadyz(i)) &
                 +tyyyy*quadyy(i)+tyyzz*quadzz(i)
      tzzijquadi=txxzz*quadxx(i)+(2.0d0*txyzz*quadxy(i)) &
                 +(2.0d0*txzzz*quadxz(i))+(2.0d0*tyzzz*quadyz(i)) &
                 +tyyzz*quadyy(i)+tzzzz*quadzz(i)
      txyijquadi=txxxy*quadxx(i)+(2.0d0*txxyy*quadxy(i)) &
                 +(2.0d0*txxyz*quadxz(i))+txyyy*quadyy(i) &
                 +(2.0d0*txyyz*quadyz(i))+txyzz*quadzz(i)
      txzijquadi=txxxz*quadxx(i)+(2.0d0*txxyz*quadxy(i)) &
                 +(2.0d0*txxzz*quadxz(i))+txyyz*quadyy(i) &
                 +(2.0d0*txyzz*quadyz(i))+txzzz*quadzz(i)
      tyzijquadi=txxyz*quadxx(i)+(2.0d0*txyyz*quadxy(i)) &
                 +(2.0d0*txyzz*quadxz(i))+tyyyz*quadyy(i) &
                 +(2.0d0*tyyzz*quadyz(i))+tyzzz*quadzz(i)

      quadquadeng1=efact5*(chgquad2ij*qjrec)*(chgquad2ji*qirec)
      quadquadengxx2=3.0d0*dxx*quadxx(j)+((dxx+dyy)/2.0d0)*quadyy(j)+ &
                     ((dxx+dzz)/2.0d0)*quadzz(j)+3.0d0*dxy*quadxy(j)+ &
                     3.0d0*dxz*quadxz(j)+dyz*quadyz(j)
      quadquadengyy2=3.0d0*dyy*quadyy(j)+((dxx+dyy)/2.0d0)*quadxx(j)+ &
                     ((dyy+dzz)/2.0d0)*quadzz(j)+3.0d0*dxy*quadxy(j)+ &
                     3.0d0*dyz*quadyz(j)+dxz*quadxz(j)
      quadquadengzz2=3.0d0*dzz*quadzz(j)+((dxx+dzz)/2.0d0)*quadxx(j)+ &
                     ((dyy+dzz)/2.0d0)*quadyy(j)+3.0d0*dxz*quadxz(j)+ &
                     3.0d0*dyz*quadyz(j)+dxy*quadxy(j)
      quadquadengxy2=(3.0d0*dxy/2.0d0)*(quadxx(j)+quadyy(j))+ &
                     (dxy/2.0d0)*quadzz(j)+(dxx+dyy)*quadxy(j)+ &
                     dxz*quadyz(j)+dyz*quadxz(j)
      quadquadengxz2=(3.0d0*dxz/2.0d0)*(quadxx(j)+quadzz(j))+ &
                     (dxz/2.0d0)*quadyy(j)+(dxx+dzz)*quadxz(j)+ &
                     dxy*quadyz(j)+dyz*quadxy(j)
      quadquadengyz2=(3.0d0*dyz/2.0d0)*(quadyy(j)+quadzz(j))+ &
                     (dyz/2.0d0)*quadxx(j)+(dyy+dzz)*quadyz(j)+ &
                     dxy*quadxz(j)+dxz*quadxy(j)
      quadquadeng2=-efact3*(quadxx(i)*quadquadengxx2+ &
                            quadyy(i)*quadquadengyy2+ &
                            quadzz(i)*quadquadengzz2+ &
                      2.0d0*quadxy(i)*quadquadengxy2+ &
                      2.0d0*quadxz(i)*quadquadengxz2+ &
                      2.0d0*quadyz(i)*quadquadengyz2)
      quadquadeng3=-2.0d0*efact3*( &
                   dxx*(quadxx(i)*quadxx(j)+quadxy(i)*quadxy(j)+ &
                        quadxz(i)*quadxz(j))+ &
                   dyy*(quadxy(i)*quadxy(j)+quadyy(i)*quadyy(j)+ &
                        quadyz(i)*quadyz(j))+ &
                   dzz*(quadxz(i)*quadxz(j)+quadyz(i)*quadyz(j)+ &
                        quadzz(i)*quadzz(j))+ &
                   dxy*(quadxx(i)*quadxy(j)+quadxx(j)*quadxy(i)+ &
                        quadyy(i)*quadxy(j)+quadyy(j)*quadxy(i)+ &
                        quadxz(i)*quadyz(j)+quadxz(j)*quadyz(i))+ &
                   dxz*(quadxx(i)*quadxz(j)+quadxx(j)*quadxz(i)+ &
                        quadzz(i)*quadxz(j)+quadzz(j)*quadxz(i)+ &
                        quadxy(i)*quadyz(j)+quadxy(j)*quadyz(i))+ &
                   dyz*(quadyy(i)*quadyz(j)+quadyy(j)*quadyz(i)+ &
                        quadzz(i)*quadyz(j)+quadzz(j)*quadyz(i)+ &
                        quadxy(i)*quadxz(j)+quadxy(j)*quadxz(i)))
      quadquadeng4=2.0d0*efact2*(quadxx(i)*quadxx(j)+ &
                                 quadyy(i)*quadyy(j)+ &
                                 quadzz(i)*quadzz(j)+ &
                          2.0d0*(quadxy(i)*quadxy(j)+ &
                                 quadxz(i)*quadxz(j)+ &
                                 quadyz(i)*quadyz(j)))
      quadquadeng5=-efact6* &
                   ((chgquad2ij*qjrec)*(chgquadengji*qirec)+ &
                   (chgquad2ji*qirec)*(chgquadengij*qjrec))/2.0d0
      quadquadeng6ij=efact4*( &
                     quadxx(i)*(dx*txijdotthetaij)+ &
                     quadyy(i)*(dy*tyijdotthetaij)+ &
                     quadzz(i)*(dz*tzijdotthetaij)+ &
               quadxy(i)*(dx*tyijdotthetaij+dy*txijdotthetaij)+ &
               quadxz(i)*(dx*tzijdotthetaij+dz*txijdotthetaij)+ &
               quadyz(i)*(dy*tzijdotthetaij+dz*tyijdotthetaij))
      quadquadeng6ji=efact4*( &
                     quadxx(j)*(dx*txijdotthetaji)+ &
                     quadyy(j)*(dy*tyijdotthetaji)+ &
                     quadzz(j)*(dz*tzijdotthetaji)+ &
               quadxy(j)*(dx*tyijdotthetaji+dy*txijdotthetaji)+ &
               quadxz(j)*(dx*tzijdotthetaji+dz*txijdotthetaji)+ &
               quadyz(j)*(dy*tzijdotthetaji+dz*tyijdotthetaji))
      quadquadeng6=quadquadeng6ij+quadquadeng6ji
      quadquadeng7=efact*(quadxxiTquadj+quadyyiTquadj+ &
                          quadzziTquadj+2.0d0*(quadxyiTquadj+ &
                          quadxziTquadj+quadyziTquadj))

      quadquadeng=quadquadeng1+quadquadeng2+quadquadeng3+quadquadeng4 &
                 +quadquadeng5+quadquadeng6+quadquadeng7
      quadquadeng=quadquadeng*onethird*onethird

!     eng=eng+chgquadeng+quadquadeng-quaddipeng 
      sctmp(58)=sctmp(58)+chgquadeng+quadquadeng-quaddipeng 

      txxxxx=-945.0d0*dx*dx*dx*dx*dx*dr11rec &
            +1050d0*dx*dx*dx*dr9rec-225.0d0*dx*dr7rec
      tyyyyy=-945.0d0*dy*dy*dy*dy*dy*dr11rec &
            +1050d0*dy*dy*dy*dr9rec-225.0d0*dy*dr7rec
      tzzzzz=-945.0d0*dz*dz*dz*dz*dz*dr11rec &
            +1050d0*dz*dz*dz*dr9rec-225.0d0*dz*dr7rec
      txxxxy=-945.0d0*dx*dx*dx*dx*dy*dr11rec &
            -45.0d0*dy*dr7rec+630.0d0*dx*dx*dy*dr9rec
      txxxxz=-945.0d0*dx*dx*dx*dx*dz*dr11rec &
            -45.0d0*dz*dr7rec+630.0d0*dx*dx*dz*dr9rec
      tyyyyz=-945.0d0*dy*dy*dy*dy*dz*dr11rec &
            -45.0d0*dz*dr7rec+630.0d0*dy*dy*dz*dr9rec
      txyyyy=-945.0d0*dx*dy*dy*dy*dy*dr11rec &
            -45.0d0*dx*dr7rec+630.0d0*dx*dy*dy*dr9rec
      txzzzz=-945.0d0*dx*dz*dz*dz*dz*dr11rec &
            -45.0d0*dx*dr7rec+630.0d0*dx*dz*dz*dr9rec
      tyzzzz=-945.0d0*dy*dz*dz*dz*dz*dr11rec &
            -45.0d0*dy*dr7rec+630.0d0*dy*dz*dz*dr9rec
      txxxyy=-945.0d0*dx*dx*dx*dy*dy*dr11rec+105.0d0*dx*dx*dx*dr9rec &
            +315.0d0*dx*dy*dy*dr9rec-45.0d0*dx*dr7rec
      txxxzz=-945.0d0*dx*dx*dx*dz*dz*dr11rec+105.0d0*dx*dx*dx*dr9rec &
            +315.0d0*dx*dz*dz*dr9rec-45.0d0*dx*dr7rec
      tyyyzz=-945.0d0*dy*dy*dy*dz*dz*dr11rec+105.0d0*dy*dy*dy*dr9rec &
            +315.0d0*dy*dz*dz*dr9rec-45.0d0*dy*dr7rec
      txxyyy=-945.0d0*dx*dx*dy*dy*dy*dr11rec+105.0d0*dy*dy*dy*dr9rec &
            +315.0d0*dx*dx*dy*dr9rec-45.0d0*dy*dr7rec
      txxzzz=-945.0d0*dx*dx*dz*dz*dz*dr11rec+105.0d0*dz*dz*dz*dr9rec &
            +315.0d0*dx*dx*dz*dr9rec-45.0d0*dz*dr7rec
      tyyzzz=-945.0d0*dy*dy*dz*dz*dz*dr11rec+105.0d0*dz*dz*dz*dr9rec &
            +315.0d0*dy*dy*dz*dr9rec-45.0d0*dz*dr7rec
      txxxyz=-945.0d0*dx*dx*dx*dy*dz*dr11rec+315.0d0*dx*dy*dz*dr9rec
      txyyyz=-945.0d0*dy*dy*dy*dx*dz*dr11rec+315.0d0*dx*dy*dz*dr9rec
      txyzzz=-945.0d0*dz*dz*dz*dx*dy*dr11rec+315.0d0*dx*dy*dz*dr9rec
      txxyyz=-945.0d0*dx*dx*dy*dy*dz*dr11rec &
            -15.0d0*dz*dr7rec+105.0d0*dz*dr9rec*((dx*dx)+(dy*dy))
      txxyzz=-945.0d0*dx*dx*dy*dz*dz*dr11rec &
            -15.0d0*dy*dr7rec+105.0d0*dy*dr9rec*((dx*dx)+(dz*dz))
      txyyzz=-945.0d0*dx*dy*dy*dz*dz*dr11rec &
            -15.0d0*dx*dr7rec+105.0d0*dx*dr9rec*((dy*dy)+(dz*dz))

      fquadquadxxbyrx=quadxx(i)*(txxxxx*quadxx(j)+(2.0d0*txxxxy*quadxy(j)) &
                     +(2.0d0*txxxxz*quadxz(j))+(2.0d0*txxxyz*quadyz(j)) &
                     +txxxyy*quadyy(j)+txxxzz*quadzz(j))
      fquadquadxxbyry=quadxx(i)*(txxxxy*quadxx(j)+(2.0d0*txxxyy*quadxy(j)) &
                     +(2.0d0*txxxyz*quadxz(j))+(2.0d0*txxyyz*quadyz(j)) &
                     +txxyyy*quadyy(j)+txxyzz*quadzz(j))
      fquadquadxxbyrz=quadxx(i)*(txxxxz*quadxx(j)+(2.0d0*txxxyz*quadxy(j)) &
                     +(2.0d0*txxxzz*quadxz(j))+(2.0d0*txxyzz*quadyz(j)) &
                     +txxyyz*quadyy(j)+txxzzz*quadzz(j))
      fquadquadyybyrx=quadyy(i)*(txxxyy*quadxx(j)+(2.0d0*txxyyy*quadxy(j)) &
                     +(2.0d0*txxyyz*quadxz(j))+(2.0d0*txyyyz*quadyz(j)) &
                     +txyyyy*quadyy(j)+txyyzz*quadzz(j))
      fquadquadyybyry=quadyy(i)*(txxyyy*quadxx(j)+(2.0d0*txyyyy*quadxy(j)) &
                     +(2.0d0*txyyyz*quadxz(j))+(2.0d0*tyyyyz*quadyz(j)) &
                     +tyyyyy*quadyy(j)+tyyyzz*quadzz(j))
      fquadquadyybyrz=quadyy(i)*(txxyyz*quadxx(j)+(2.0d0*txyyyz*quadxy(j)) &
                     +(2.0d0*txyyzz*quadxz(j))+(2.0d0*tyyyzz*quadyz(j)) &
                     +tyyyyz*quadyy(j)+tyyzzz*quadzz(j))
      fquadquadzzbyrx=quadzz(i)*(txxxzz*quadxx(j)+(2.0d0*txxyzz*quadxy(j)) &
                     +(2.0d0*txxzzz*quadxz(j))+(2.0d0*txyzzz*quadyz(j)) &
                     +txyyzz*quadyy(j)+txzzzz*quadzz(j))
      fquadquadzzbyry=quadzz(i)*(txxyzz*quadxx(j)+(2.0d0*txyyzz*quadxy(j)) &
                     +(2.0d0*txyzzz*quadxz(j))+(2.0d0*tyyzzz*quadyz(j)) &
                     +tyyyzz*quadyy(j)+tyzzzz*quadzz(j))
      fquadquadzzbyrz=quadzz(i)*(txxzzz*quadxx(j)+(2.0d0*txyzzz*quadxy(j)) &
                     +(2.0d0*txzzzz*quadxz(j))+(2.0d0*tyzzzz*quadyz(j)) &
                     +tyyzzz*quadyy(j)+tzzzzz*quadzz(j))
      fquadquadxybyrx=quadxy(i)*((2.0d0*txxxxy*quadxx(j)) &
                    +(4.0d0*txxxyy*quadxy(j))+(4.0d0*txxxyz*quadxz(j)) &
                    +(2.0d0*txxyyy*quadyy(j))+(4.0d0*txxyyz*quadyz(j)) &
                    +(2.0d0*txxyzz*quadzz(j)))
      fquadquadxybyry=quadxy(i)*((2.0d0*txxxyy*quadxx(j)) &
                    +(4.0d0*txxyyy*quadxy(j))+(4.0d0*txxyyz*quadxz(j)) &
                    +(2.0d0*txyyyy*quadyy(j))+(4.0d0*txyyyz*quadyz(j)) &
                    +(2.0d0*txyyzz*quadzz(j)))
      fquadquadxybyrz=quadxy(i)*((2.0d0*txxxyz*quadxx(j)) &
                    +(4.0d0*txxyyz*quadxy(j))+(4.0d0*txxyzz*quadxz(j)) &
                    +(2.0d0*txyyyz*quadyy(j))+(4.0d0*txyyzz*quadyz(j)) &
                    +(2.0d0*txyzzz*quadzz(j)))
      fquadquadxzbyrx=quadxz(i)*((2.0d0*txxxxz*quadxx(j)) &
                    +(4.0d0*txxxyz*quadxy(j))+(4.0d0*txxxzz*quadxz(j)) &
                    +(2.0d0*txxyyz*quadyy(j))+(4.0d0*txxyzz*quadyz(j)) &
                    +(2.0d0*txxzzz*quadzz(j)))
      fquadquadxzbyry=quadxz(i)*((2.0d0*txxxyz*quadxx(j)) &
                    +(4.0d0*txxyyz*quadxy(j))+(4.0d0*txxyzz*quadxz(j)) &
                    +(2.0d0*txyyyz*quadyy(j))+(4.0d0*txyyzz*quadyz(j)) &
                    +(2.0d0*txyzzz*quadzz(j)))
      fquadquadxzbyrz=quadxz(i)*((2.0d0*txxxzz*quadxx(j)) &
                    +(4.0d0*txxyzz*quadxy(j))+(4.0d0*txxzzz*quadxz(j)) &
                    +(2.0d0*txyyzz*quadyy(j))+(4.0d0*txyzzz*quadyz(j)) &
                    +(2.0d0*txzzzz*quadzz(j)))
      fquadquadyzbyrx=quadyz(i)*((2.0d0*txxxyz*quadxx(j)) &
                    +(4.0d0*txxyyz*quadxy(j))+(4.0d0*txxyzz*quadxz(j)) &
                    +(2.0d0*txyyyz*quadyy(j))+(4.0d0*txyyzz*quadyz(j)) &
                    +(2.0d0*txyzzz*quadzz(j)))
      fquadquadyzbyry=quadyz(i)*((2.0d0*txxyyz*quadxx(j)) &
                    +(4.0d0*txyyyz*quadxy(j))+(4.0d0*txyyzz*quadxz(j)) &
                    +(2.0d0*tyyyyz*quadyy(j))+(4.0d0*tyyyzz*quadyz(j)) &
                    +(2.0d0*tyyzzz*quadzz(j)))
      fquadquadyzbyrz=quadyz(i)*((2.0d0*txxyzz*quadxx(j)) &
                    +(4.0d0*txyyzz*quadxy(j))+(4.0d0*txyzzz*quadxz(j)) &
                    +(2.0d0*tyyyzz*quadyy(j))+(4.0d0*tyyzzz*quadyz(j)) &
                    +(2.0d0*tyzzzz*quadzz(j)))

      efact7=2.0d0*efact5*(etasq+drsqrec)+8.0d0*efact3*dr4rec+ &
             16.0d0*efact2*dr6rec
      efact8=7.0d0*efact2*drrec-2.0d0*dr*efact3*(2.0d0-drsq*etasq)

      fquadquad1=-efact7*(chgquad2ij*qjrec)*(chgquad2ji*qirec)
      fquadquad2a=2.0d0*efact5*(chgquad2ij*qjrec)
      fquadquad2b=2.0d0*efact5*(chgquad2ji*qirec)
      fquadquad3=-(quadquadeng2/efact3)*efact5
      fquadquad4xxx=6.0d0*dx*quadxx(j)+dx*(quadyy(j)+quadzz(j))+ &
                    3.0d0*dy*quadxy(j)+3.0d0*dz*quadxz(j)
      fquadquad4yyx=dx*quadxx(j)+3.0d0*dy*quadxy(j)+dz*quadxz(j)
      fquadquad4zzx=dx*quadxx(j)+dy*quadxy(j)+3.0d0*dz*quadxz(j)
      fquadquad4xyx=3.0d0*dy*(quadxx(j)+quadyy(j))+4.0d0*dx*quadxy(j)+ &
                    dz*quadyz(j)+dy*quadzz(j)
      fquadquad4xzx=3.0d0*dz*(quadxx(j)+quadzz(j))+4.0d0*dx*quadxz(j)+ &
                    dy*quadyz(j)+dz*quadyy(j)
      fquadquad4yzx=dy*quadxz(j)+dz*quadxy(j)
      fquadquad4x=-(quadxx(i)*fquadquad4xxx+quadyy(i)*fquadquad4yyx+ &
                    quadzz(i)*fquadquad4zzx+quadxy(i)*fquadquad4xyx+ &
                    quadxz(i)*fquadquad4xzx+quadyz(i)*fquadquad4yzx)
      fquadquad4xxy=dy*quadyy(j)+3.0d0*dx*quadxy(j)+dz*quadyz(j)
      fquadquad4yyy=6.0d0*dy*quadyy(j)+dy*(quadxx(j)+quadzz(j))+ &
                    3.0d0*dx*quadxy(j)+3.0d0*dz*quadyz(j)
      fquadquad4zzy=dy*quadyy(j)+dx*quadxy(j)+3.0d0*dz*quadyz(j)
      fquadquad4xyy=3.0d0*dx*(quadxx(j)+quadyy(j))+4.0d0*dy*quadxy(j)+ &
                    dz*quadxz(j)+dx*quadzz(j)
      fquadquad4xzy=dx*quadyz(j)+dz*quadxy(j)
      fquadquad4yzy=3.0d0*dz*(quadyy(j)+quadzz(j))+4.0d0*dy*quadyz(j)+ &
                    dx*quadxz(j)+dz*quadxx(j)
      fquadquad4y=-(quadxx(i)*fquadquad4xxy+quadyy(i)*fquadquad4yyy+ &
                    quadzz(i)*fquadquad4zzy+quadxy(i)*fquadquad4xyy+ &
                    quadxz(i)*fquadquad4xzy+quadyz(i)*fquadquad4yzy)
      fquadquad4xxz=dz*quadzz(j)+3.0d0*dx*quadxz(j)+dy*quadyz(j)
      fquadquad4yyz=dz*quadzz(j)+dx*quadxz(j)+3.0d0*dy*quadyz(j)
      fquadquad4zzz=6.0d0*dz*quadzz(j)+dz*(quadxx(j)+quadyy(j))+ &
                    3.0d0*dx*quadxz(j)+3.0d0*dy*quadyz(j)
      fquadquad4xyz=dx*quadyz(j)+dy*quadxz(j)
      fquadquad4xzz=3.0d0*dx*(quadxx(j)+quadzz(j))+4.0d0*dz*quadxz(j)+ &
                    dy*quadxy(j)+dx*quadyy(j)
      fquadquad4yzz=3.0d0*dy*(quadyy(j)+quadzz(j))+4.0d0*dz*quadyz(j)+ &
                    dx*quadxy(j)+dy*quadxx(j)
      fquadquad4z=-(quadxx(i)*fquadquad4xxz+quadyy(i)*fquadquad4yyz+ &
                    quadzz(i)*fquadquad4zzz+quadxy(i)*fquadquad4xyz+ &
                    quadxz(i)*fquadquad4xzz+quadyz(i)*fquadquad4yzz)
      fquadquad5=-(quadquadeng3/efact3)*efact5
      fquadquad6x=-2.0d0*(2.0d0*dx*(quadxx(i)*quadxx(j)+ &
                  quadxy(i)*quadxy(j)+quadxz(i)*quadxz(j))+ &
                  dy*(quadxx(i)*quadxy(j)+quadxx(j)*quadxy(i)+ &
                      quadyy(i)*quadxy(j)+quadyy(j)*quadxy(i)+ &
                      quadxz(i)*quadyz(j)+quadxz(j)*quadyz(i))+ &
                  dz*(quadxx(i)*quadxz(j)+quadxx(j)*quadxz(i)+ &
                      quadzz(i)*quadxz(j)+quadzz(j)*quadxz(i)+ &
                      quadxy(i)*quadyz(j)+quadxy(j)*quadyz(i)))
      fquadquad6y=-2.0d0*(2.0d0*dy*(quadxy(i)*quadxy(j)+ &
                  quadyy(i)*quadyy(j)+quadyz(i)*quadyz(j))+ &
                  dx*(quadxx(i)*quadxy(j)+quadxx(j)*quadxy(i)+ &
                      quadyy(i)*quadxy(j)+quadyy(j)*quadxy(i)+ &
                      quadxz(i)*quadyz(j)+quadxz(j)*quadyz(i))+ &
                  dz*(quadyy(i)*quadyz(j)+quadyy(j)*quadyz(i)+ &
                      quadzz(i)*quadyz(j)+quadzz(j)*quadyz(i)+ &
                      quadxz(i)*quadxy(j)+quadxz(j)*quadxy(i)))
      fquadquad6z=-2.0d0*(2.0d0*dz*(quadxz(i)*quadxz(j)+ &
                  quadyz(i)*quadyz(j)+quadzz(i)*quadzz(j))+ &
                  dx*(quadxx(i)*quadxz(j)+quadxx(j)*quadxz(i)+ &
                      quadzz(i)*quadxz(j)+quadzz(j)*quadxz(i)+ &
                      quadxy(i)*quadyz(j)+quadxy(j)*quadyz(i))+ &
                  dy*(quadyy(i)*quadyz(j)+quadyy(j)*quadyz(i)+ &
                      quadzz(i)*quadyz(j)+quadzz(j)*quadyz(i)+ &
                      quadxz(i)*quadxy(j)+quadxz(j)*quadxy(i)))
      fquadquad7=-quadquadeng4*(efact3/efact2)
      fquadquad8=(efact8/efact6)*quadquadeng5
      fquadquad9a1=-efact6*(chgquadengji*qirec)
      fquadquad9a2=-efact6*(chgquadengij*qjrec)
      fquadquad9b1=efact6*(chgquad2ij*qjrec)
      fquadquad9b2=efact6*(chgquad2ji*qirec)
      fquadquad10=-(efact6/efact4)*quadquadeng6
      fquadquad11xij=-(quadxx(i)*dx*txxijquadj+quadyy(i)*dy*txyijquadj+ &
                      quadzz(i)*dz*txzijquadj+quadxy(i)*(dx*txyijquadj+ &
                                 dy*txxijquadj)+quadxz(i)*(dx*txzijquadj+ &
                                 dz*txxijquadj)+quadyz(i)*(dy*txzijquadj+ &
                                 dz*txyijquadj))
      fquadquad11xji=-(quadxx(j)*dx*txxijquadi+quadyy(j)*dy*txyijquadi+ &
                      quadzz(j)*dz*txzijquadi+quadxy(j)*(dx*txyijquadi+ &
                                 dy*txxijquadi)+quadxz(j)*(dx*txzijquadi+ &
                                 dz*txxijquadi)+quadyz(j)*(dy*txzijquadi+ &
                                 dz*txyijquadi))
      fquadquad11x=fquadquad11xij+fquadquad11xji
      fquadquad11yij=-(quadxx(i)*dx*txyijquadj+quadyy(i)*dy*tyyijquadj+ &
                      quadzz(i)*dz*tyzijquadj+quadxy(i)*(dx*tyyijquadj+ &
                                 dy*txyijquadj)+quadxz(i)*(dx*tyzijquadj+ &
                                 dz*txyijquadj)+quadyz(i)*(dy*tyzijquadj+ &
                                 dz*tyyijquadj))
      fquadquad11yji=-(quadxx(j)*dx*txyijquadi+quadyy(j)*dy*tyyijquadi+ &
                      quadzz(j)*dz*tyzijquadi+quadxy(j)*(dx*tyyijquadi+ &
                                 dy*txyijquadi)+quadxz(j)*(dx*tyzijquadi+ &
                                 dz*txyijquadi)+quadyz(j)*(dy*tyzijquadi+ &
                                 dz*tyyijquadi))
      fquadquad11y=fquadquad11yij+fquadquad11yji
      fquadquad11zij=-(quadxx(i)*dx*txzijquadj+quadyy(i)*dy*tyzijquadj+ &
                      quadzz(i)*dz*tzzijquadj+quadxy(i)*(dx*tyzijquadj+ &
                                 dy*txzijquadj)+quadxz(i)*(dx*tzzijquadj+ &
                                 dz*txzijquadj)+quadyz(i)*(dy*tzzijquadj+ &
                                 dz*tyzijquadj))
      fquadquad11zji=-(quadxx(j)*dx*txzijquadi+quadyy(j)*dy*tyzijquadi+ &
                      quadzz(j)*dz*tzzijquadi+quadxy(j)*(dx*tyzijquadi+ &
                                 dy*txzijquadi)+quadxz(j)*(dx*tzzijquadi+ &
                                 dz*txzijquadi)+quadyz(j)*(dy*tzzijquadi+ &
                                 dz*tyzijquadi))
      fquadquad11z=fquadquad11zij+fquadquad11zji
      fquadquad12=-(efact4/efact)*quadquadeng7
      fquadquad13x=fquadquadxxbyrx+fquadquadyybyrx &
                +fquadquadzzbyrx+fquadquadxybyrx+fquadquadxzbyrx+fquadquadyzbyrx
      fquadquad13y=fquadquadxxbyry+fquadquadyybyry &
                +fquadquadzzbyry+fquadquadxybyry+fquadquadxzbyry+fquadquadyzbyry
      fquadquad13z=fquadquadxxbyrz+fquadquadyybyrz &
                +fquadquadzzbyrz+fquadquadxybyrz+fquadquadxzbyrz+fquadquadyzbyrz

      fquadquadx=fquadquad1*dx+ &
                 fquadquad2a*thetaxjdotrij+fquadquad2b*thetaxidotrij+ &
                 fquadquad3*dx+fquadquad4x*efact3+ &
                 fquadquad5*dx+fquadquad6x*efact3+ &
                 fquadquad7*dx+fquadquad8*dx+ &
                 fquadquad9a1*thetaxidotrij+fquadquad9a2*thetaxjdotrij+ &
                 fquadquad9b1*txijdotthetaij+fquadquad9b2*txijdotthetaji+ &
                 fquadquad10*dx+fquadquad11x*efact4+ &
                 fquadquad12*dx+fquadquad13x*efact
      fquadquady=fquadquad1*dy+ &
                 fquadquad2a*thetayjdotrij+fquadquad2b*thetayidotrij+ &
                 fquadquad3*dy+fquadquad4y*efact3+ &
                 fquadquad5*dy+fquadquad6y*efact3+ &
                 fquadquad7*dy+fquadquad8*dy+ &
                 fquadquad9a1*thetayidotrij+fquadquad9a2*thetayjdotrij+ &
                 fquadquad9b1*tyijdotthetaij+fquadquad9b2*tyijdotthetaji+ &
                 fquadquad10*dy+fquadquad11y*efact4+ &
                 fquadquad12*dy+fquadquad13y*efact
      fquadquadz=fquadquad1*dz+ &
                 fquadquad2a*thetazjdotrij+fquadquad2b*thetazidotrij+ &
                 fquadquad3*dz+fquadquad4z*efact3+ &
                 fquadquad5*dz+fquadquad6z*efact3+ &
                 fquadquad7*dz+fquadquad8*dz+ &
                 fquadquad9a1*thetazidotrij+fquadquad9a2*thetazjdotrij+ &
                 fquadquad9b1*tzijdotthetaij+fquadquad9b2*tzijdotthetaji+ &
                 fquadquad10*dz+fquadquad11z*efact4+ &
                 fquadquad12*dz+fquadquad13z*efact

      fquadquadx=fquadquadx*onethird*onethird
      fquadquady=fquadquady*onethird*onethird
      fquadquadz=fquadquadz*onethird*onethird

!     stpquadquadxx=stpquadquadxx-fquadquadx*dxsav(i,j)
!     stpquadquadxy=stpquadquadxy-0.5d0*(fquadquadx*dysav(i,j) &
!                       +fquadquady*dxsav(i,j))
!     stpquadquadxz=stpquadquadxz-0.5d0*(fquadquadx*dzsav(i,j) &
!                       +fquadquadz*dxsav(i,j))
!     stpquadquadyy=stpquadquadyy-fquadquady*dysav(i,j)
!     stpquadquadyz=stpquadquadyz-0.5d0*(fquadquadz*dysav(i,j) &
!                       +fquadquady*dzsav(i,j))
!     stpquadquadzz=stpquadquadzz-fquadquadz*dzsav(i,j)
      sctmp(37)=sctmp(37)-fquadquadx*dxsav(numx)
      sctmp(38)=sctmp(38)-0.5d0*(fquadquadx*dysav(numx) &
                        +fquadquady*dxsav(numx))
      sctmp(39)=sctmp(39)-0.5d0*(fquadquadx*dzsav(numx) &
                        +fquadquadz*dxsav(numx))
      sctmp(40)=sctmp(40)-fquadquady*dysav(numx)
      sctmp(41)=sctmp(41)-0.5d0*(fquadquadz*dysav(numx) &
                        +fquadquady*dzsav(numx))
      sctmp(42)=sctmp(42)-fquadquadz*dzsav(numx)

      fttx=fqquadx-fdipquadx-fquadquadx
      ftty=fqquady-fdipquady-fquadquady
      fttz=fqquadz-fdipquadz-fquadquadz

!     frrx(j)=frrx(j)-fttx
!     frry(j)=frry(j)-ftty
!     frrz(j)=frrz(j)-fttz
!     frrx(i)=frrx(i)+fttx
!     frry(i)=frry(i)+ftty
!     frrz(i)=frrz(i)+fttz
      eltmp(j)=eltmp(j)-fttx
      eltmp(num+j)=eltmp(num+j)-ftty
      eltmp(2*num+j)=eltmp(2*num+j)-fttz
      eltmp(i)=eltmp(i)+fttx
      eltmp(num+i)=eltmp(num+i)+ftty
      eltmp(2*num+i)=eltmp(2*num+i)+fttz

#endif 

#ifdef dipole 
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
      dampfuncdiffi=dampa(jpoint,ipoint)*dampaexpi*dampfac(jpoint,ipoint) &
             *(((dampa(jpoint,ipoint)*dr)**nkdamp(jpoint,ipoint))/(factorial))

      r3dampi=dr3rec*dampfunci*q(j)

      elecsrxi=dxsav(numx)*r3dampi
      elecsryi=dysav(numx)*r3dampi
      elecsrzi=dzsav(numx)*r3dampi

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
      dampfuncdiffj=dampa(ipoint,jpoint)*dampaexpj*dampfac(ipoint,jpoint) &
             *(((dampa(ipoint,jpoint)*dr)**nkdamp(ipoint,jpoint))/(factorial))

      r3dampj=dr3rec*dampfuncj*q(i)

      elecsrxj=-dxsav(numx)*r3dampj
      elecsryj=-dysav(numx)*r3dampj
      elecsrzj=-dzsav(numx)*r3dampj

!     srdipeng(ipoint)=srdipeng(ipoint)-(elecsrxi*xmu(i) &
!                         +elecsryi*ymu(i)+elecsrzi*zmu(i))
!     srdipeng(jpoint)=srdipeng(jpoint)-(elecsrxj*xmu(j) &
!                         +elecsryj*ymu(j)+elecsrzj*zmu(j))
      dimtmp(ipoint,1)=dimtmp(ipoint,1)-(elecsrxi*xmu(i) &
                          +elecsryi*ymu(i)+elecsrzi*zmu(i))
      dimtmp(jpoint,1)=dimtmp(jpoint,1)-(elecsrxj*xmu(j) &
                          +elecsryj*ymu(j)+elecsrzj*zmu(j))

!      dipi=xmu(i)*dxsav(numx)+ymu(i)*dysav(numx)+zmu(i)*dzsav(numx)
!      dipj=xmu(j)*dxsav(numx)+ymu(j)*dysav(numx)+zmu(j)*dzsav(numx)
! GWW - not sure why the same calc is save in two different variables 

      diprfunci=xmuidotrij*drrec*q(j)*((-3.0d0*dampfunci*dr3rec*drrec) &
              +(dampfuncdiffi*dr3rec))
      diprfuncj=xmujdotrij*drrec*q(i)*((-3.0d0*dampfuncj*dr3rec*drrec) &
              +(dampfuncdiffj*dr3rec))

      frrxsri=(xmu(i)*r3dampi)+(dxsav(numx)*diprfunci)
      frrysri=(ymu(i)*r3dampi)+(dysav(numx)*diprfunci)
      frrzsri=(zmu(i)*r3dampi)+(dzsav(numx)*diprfunci)
      frrxsrj=(xmu(j)*r3dampj)+(dxsav(numx)*diprfuncj)
      frrysrj=(ymu(j)*r3dampj)+(dysav(numx)*diprfuncj)
      frrzsrj=(zmu(j)*r3dampj)+(dzsav(numx)*diprfuncj)

!     stpsrxx=stpsrxx+frrxsri*dxsav(i,j)-frrxsrj*dxsav(i,j)
!     stpsrxy=stpsrxy+0.5d0*(frrxsri*dysav(i,j)+frrysri*dxsav(i,j)) &
!                    -0.5d0*(frrxsrj*dysav(i,j)+frrysrj*dxsav(i,j))
!     stpsrxz=stpsrxz+0.5d0*(frrxsri*dzsav(i,j)+frrzsri*dxsav(i,j)) &
!                    -0.5d0*(frrxsrj*dzsav(i,j)+frrzsrj*dxsav(i,j))
!     stpsryy=stpsryy+frrysri*dysav(i,j)-frrysrj*dysav(i,j)
!     stpsryz=stpsryz+0.5d0*(frrysri*dzsav(i,j)+frrzsri*dysav(i,j)) &
!                    -0.5d0*(frrysrj*dzsav(i,j)+frrzsrj*dysav(i,j))
!     stpsrzz=stpsrzz+frrzsri*dzsav(i,j)-frrzsrj*dzsav(i,j)
      sctmp(43)=sctmp(43)+frrxsri*dxsav(numx)-frrxsrj*dxsav(numx)
      sctmp(44)=sctmp(44)+0.5d0*(frrxsri*dysav(numx)+frrysri*dxsav(numx)) &
                     -0.5d0*(frrxsrj*dysav(numx)+frrysrj*dxsav(numx))
      sctmp(45)=sctmp(45)+0.5d0*(frrxsri*dzsav(numx)+frrzsri*dxsav(numx)) &
                     -0.5d0*(frrxsrj*dzsav(numx)+frrzsrj*dxsav(numx))
      sctmp(46)=sctmp(46)+frrysri*dysav(numx)-frrysrj*dysav(numx)
      sctmp(47)=sctmp(47)+0.5d0*(frrysri*dzsav(numx)+frrzsri*dysav(numx)) &
                     -0.5d0*(frrysrj*dzsav(numx)+frrzsrj*dysav(numx))
      sctmp(48)=sctmp(48)+frrzsri*dzsav(numx)-frrzsrj*dzsav(numx)

!     frrx(i)=frrx(i)+frrxsri-frrxsrj
!     frry(i)=frry(i)+frrysri-frrysrj
!     frrz(i)=frrz(i)+frrzsri-frrzsrj
!     frrx(j)=frrx(j)-frrxsri+frrxsrj
!     frry(j)=frry(j)-frrysri+frrysrj
!     frrz(j)=frrz(j)-frrzsri+frrzsrj
      eltmp(j)=eltmp(j)-frrxsri+frrxsrj
      eltmp(num+j)=eltmp(num+j)-frrysri+frrysrj
      eltmp(2*num+j)=eltmp(2*num+j)-frrzsri+frrzsrj
      eltmp(i)=eltmp(i)+frrxsri-frrxsrj
      eltmp(num+i)=eltmp(num+i)+frrysri-frrysrj
      eltmp(2*num+i)=eltmp(2*num+i)+frrzsri-frrzsrj

! check q-mu damping forces 
      gw_force( 9*num+j)=gw_force( 9*num+j)-frrxsri+frrxsrj
      gw_force(10*num+j)=gw_force(10*num+j)-frrysri+frrysrj
      gw_force(11*num+j)=gw_force(11*num+j)-frrzsri+frrzsrj
      gw_force( 9*num+i)=gw_force( 9*num+i)+frrxsri-frrxsrj
      gw_force(10*num+i)=gw_force(10*num+i)+frrysri-frrysrj
      gw_force(11*num+i)=gw_force(11*num+i)+frrzsri-frrzsrj
#endif

#ifdef quadrupole  
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
      dampfuncdiffi=fgb(jpoint,ipoint)*dampaexpi*fgc(jpoint,ipoint) &
             *(((fgb(jpoint,ipoint)*dr)**nkfg(jpoint,ipoint))/(factorial))

      exxsrg=txx*dampfunci*q(j)
      eyysrg=tyy*dampfunci*q(j)
      ezzsrg=tzz*dampfunci*q(j)
      exysrg=txy*dampfunci*q(j)
      exzsrg=txz*dampfunci*q(j)
      eyzsrg=tyz*dampfunci*q(j)

!---> Parallelization_S
!     srquadeng(ipoint)=srquadeng(ipoint)+onethird*(exxsrg*quadxx(i) &
      dimtmp(ipoint,2)=dimtmp(ipoint,2)+onethird*(exxsrg*quadxx(i) &
                                  +eyysrg*quadyy(i)+ezzsrg*quadzz(i) &
                           +2.0d0*(exysrg*quadxy(i)+exzsrg*quadxz(i) &
                                        +eyzsrg*quadyz(i)))
!<--- Parallelization_E

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
      dampfuncdiffj=fgb(ipoint,jpoint)*dampaexpj*fgc(ipoint,jpoint) &
             *(((fgb(ipoint,jpoint)*dr)**nkfg(ipoint,jpoint))/(factorial))

      exxsrg=txx*dampfuncj*q(i)
      eyysrg=tyy*dampfuncj*q(i)
      ezzsrg=tzz*dampfuncj*q(i)
      exysrg=txy*dampfuncj*q(i)
      exzsrg=txz*dampfuncj*q(i)
      eyzsrg=tyz*dampfuncj*q(i)

!---> Parallelization_S
!     srquadeng(jpoint)=srquadeng(jpoint)+onethird*(exxsrg*quadxx(j) &
      dimtmp(jpoint,2)=dimtmp(jpoint,2)+onethird*(exxsrg*quadxx(j) &
                                  +eyysrg*quadyy(j)+ezzsrg*quadzz(j) &
                           +2.0d0*(exysrg*quadxy(j)+exzsrg*quadxz(j) &
                                  +eyzsrg*quadyz(j)))
!<--- Parallelization_E

      T2dotquadj=quadxx(j)*txx+quadyy(j)*tyy &
               +quadzz(j)*tzz+2.0d0*(quadxy(j)*txy &
               +quadxz(j)*txz+quadyz(j)*tyz)
      T2dotquadi=quadxx(i)*txx+quadyy(i)*tyy &
                +quadzz(i)*tzz+2.0d0*(quadxy(i)*txy &
                +quadxz(i)*txz+quadyz(i)*tyz)

      frrxsrqji=-((-txxx*quadxx(j)-txyy*quadyy(j)-txzz*quadzz(j) &
               -2.0d0*txxy*quadxy(j)-2.0d0*txxz*quadxz(j) &
               -2.0d0*txyz*quadyz(j))*dampfuncj &
               +(T2dotquadj*dampfuncdiffj*dx*drrec))*onethird*q(i)
      frrxsrqij=-((-txxx*quadxx(i)-txyy*quadyy(i)-txzz*quadzz(i) &
               -2.0d0*txxy*quadxy(i)-2.0d0*txxz*quadxz(i) &
               -2.0d0*txyz*quadyz(i))*dampfunci &
               +(T2dotquadi*dampfuncdiffi*dx*drrec))*onethird*q(j)
      frrysrqji=-((-txxy*quadxx(j)-tyyy*quadyy(j)-tyzz*quadzz(j) &
               -2.0d0*txyy*quadxy(j)-2.0d0*txyz*quadxz(j) &
               -2.0d0*tyyz*quadyz(j))*dampfuncj &
               +(T2dotquadj*dampfuncdiffj*dy*drrec))*onethird*q(i)
      frrysrqij=-((-txxy*quadxx(i)-tyyy*quadyy(i)-tyzz*quadzz(i) &
               -2.0d0*txyy*quadxy(i)-2.0d0*txyz*quadxz(i) &
               -2.0d0*tyyz*quadyz(i))*dampfunci &
               +(T2dotquadi*dampfuncdiffi*dy*drrec))*onethird*q(j)
      frrzsrqji=-((-txxz*quadxx(j)-tyyz*quadyy(j)-tzzz*quadzz(j) &
               -2.0d0*txyz*quadxy(j)-2.0d0*txzz*quadxz(j) &
               -2.0d0*tyzz*quadyz(j))*dampfuncj &
               +(T2dotquadj*dampfuncdiffj*dz*drrec))*onethird*q(i)
      frrzsrqij=-((-txxz*quadxx(i)-tyyz*quadyy(i)-tzzz*quadzz(i) &
               -2.0d0*txyz*quadxy(i)-2.0d0*txzz*quadxz(i) &
               -2.0d0*tyzz*quadyz(i))*dampfunci &
               +(T2dotquadi*dampfuncdiffi*dz*drrec))*onethird*q(j)

      fqquadxsr=frrxsrqji+frrxsrqij
      fqquadysr=frrysrqji+frrysrqij
      fqquadzsr=frrzsrqji+frrzsrqij

!---> Memmory Reduction_S
!---> Parallelization_S
!     stqsrxx=stqsrxx+fqquadxsr*dxsav(i,j)
!     stqsrxy=stqsrxy+0.5d0*(fqquadxsr*dysav(i,j) &
!                    +fqquadysr*dxsav(i,j))
!     stqsrxz=stqsrxz+0.5d0*(fqquadxsr*dzsav(i,j) &
!                    +fqquadzsr*dxsav(i,j))
!     stqsryy=stqsryy+fqquadysr*dysav(i,j)
!     stqsryz=stqsryz+0.5d0*(fqquadysr*dzsav(i,j) &
!                    +fqquadzsr*dysav(i,j))
!     stqsrzz=stqsrzz+fqquadzsr*dzsav(i,j)
      sctmp(49)=sctmp(49)+fqquadxsr*dxsav(numx)
      sctmp(50)=sctmp(50)+0.5d0*(fqquadxsr*dysav(numx) &
                     +fqquadysr*dxsav(numx))
      sctmp(51)=sctmp(51)+0.5d0*(fqquadxsr*dzsav(numx) &
                     +fqquadzsr*dxsav(numx))
      sctmp(52)=sctmp(52)+fqquadysr*dysav(numx)
      sctmp(53)=sctmp(53)+0.5d0*(fqquadysr*dzsav(numx) &
                     +fqquadzsr*dysav(numx))
      sctmp(54)=sctmp(54)+fqquadzsr*dzsav(numx)
!<--- Parallelization_E
!<--- Memmory Reduction_E

!---> Parallelization_S
!     frrx(i)=frrx(i)+fqquadxsr
!     frry(i)=frry(i)+fqquadysr
!     frrz(i)=frrz(i)+fqquadzsr
!     frrx(j)=frrx(j)-fqquadxsr
!     frry(j)=frry(j)-fqquadysr
!     frrz(j)=frrz(j)-fqquadzsr
      eltmp(j)=eltmp(j)-fqquadxsr
      eltmp(num+j)=eltmp(num+j)-fqquadysr
      eltmp(2*num+j)=eltmp(2*num+j)-fqquadzsr
      eltmp(i)=eltmp(i)+fqquadxsr
      eltmp(num+i)=eltmp(num+i)+fqquadysr
      eltmp(2*num+i)=eltmp(2*num+i)+fqquadzsr
#endif  

   enddo   
enddo   


asdipx=alppolar*elecx
asdipy=alppolar*elecy
asdipz=alppolar*elecz
srdipx=alppolar*elecxsr
srdipy=alppolar*elecysr
srdipz=alppolar*eleczsr
asquadxx=Cpolar*exx
asquadyy=Cpolar*eyy
asquadzz=Cpolar*ezz
asquadxy=Cpolar*exy
asquadxz=Cpolar*exz
asquadyz=Cpolar*eyz
srquadxx=Cpolar*exxsr
srquadyy=Cpolar*eyysr
srquadzz=Cpolar*ezzsr
srquadxy=Cpolar*exysr
srquadxz=Cpolar*exzsr
srquadyz=Cpolar*eyzsr

#ifdef ewald_surface
  eltmp(2*num+1:3*num)=eltmp(2*num+1:3*num)-facewvac*2.0d0*q*totmz
  sctmp(55)=sctmp(55)+facewvac*(totmz*totmz)
  stewzz= facewvac2*totmz*(totmz-2.0d0*totmzd)
#else
  stewzz = 0.0 
#endif 


!   print*,'Total charge interaction = ',sctmp(55) 
return
END SUBROUTINE
