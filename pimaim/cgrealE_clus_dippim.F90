SUBROUTINE cgrealE_clus_dippim

USE commondata, ONLY: elecxsr,elecysr,eleczsr,elecx,elecy,elecz,x,y,z, &
                      q,xmu,ymu,zmu,engeff,dxsav,dysav,dzsav, &
                      num,twopi,ntype,rsqmax,etapi,etasq,erfc,nspec, &
!---> Memmory Reduction_S
                      nkdamp,dampa,dampfac,num2,numx
!<--- Memmory Reduction_E
USE boxdata, ONLY: hlab2,boxlenx,boxleny,boxlenz
!---> Memmory Reduction_S
use mpipara
!<--- Memmory Reduction_E

implicit none

integer :: i,j,kk,ipoint,jpoint
double precision :: dxunitsq,dyunitsq,dzunitsq,frrxdipqij,frrxdipqji, &
                    frrydipqij,frrydipqji,frrzdipqij,frrzdipqji,xmuidotrij, &
                    xmujdotrij,ettx,etty,ettz,egam,erfcr,expewld,qq2api
double precision :: drsq,drsqrec,dr,dr3rec, &
                    chgtem,efact,efact2, &
                    sxx,syy,szz,sxy,sxz,syz,qirec,qjrec,drrec
double precision :: chgdipengij,chgdipengji, &
                    elecsrxi,elecsryi,elecsrzi,elecsrxj,elecsryj,elecsrzj, &
                    r3dampi,r3dampj,factorial,xf,dampa1,dampa2,dampa3,dampa4, &
                    dampsumfi,dampsumfj,dampfunci,dampfuncj,dampaexpi,dampaexpj
double precision :: qdotmuxij,qdotmuxji,qdotmuyij,qdotmuyji, &
                    qdotmuzij,qdotmuzji
 
elecx=0.d0
elecy=0.d0
elecz=0.d0
elecxsr=0.d0
elecysr=0.d0
eleczsr=0.d0
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

      engeff(i)=engeff(i)+erfcr*qirec
      engeff(j)=engeff(j)+erfcr*qjrec

      dampa1=dampa(jpoint,ipoint)
      dampa2=dampa1*dampa1/2.0d0
      dampa3=dampa1*dampa1*dampa1/6.0d0
      dampa4=dampa1*dampa1*dampa1*dampa1/24.0d0

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

      dampa1=dampa(ipoint,jpoint)
      dampa2=dampa1*dampa1/2.0d0
      dampa3=dampa1*dampa1*dampa1/6.0d0
      dampa4=dampa1*dampa1*dampa1*dampa1/24.0d0

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

   enddo   
enddo   

elecx=elecx+elecxsr
elecy=elecy+elecysr
elecz=elecz+eleczsr

return
END SUBROUTINE
