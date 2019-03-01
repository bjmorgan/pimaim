SUBROUTINE sr_energy

USE commondata, ONLY: num,nspec,epplog,xftlog,rvplog,cimlog,daimlog, &
                      quaimlog,ntype,dxsav,dysav,dzsav,rsqmaxsr, &
                      ftalp, &
                      ftb,rvph,rvpn,rvpr4,ftalpx,ftbx,nsp,ftb2,ftb3,ftbeta, &
                      ftgamma,delta,epsilonx,epsilony,epsilonz,quaimxx, &
                      quaimyy,quaimzz,quaimxy,quaimxz,quaimyz,ooaimlog, &
!---> Memmory Reduction_S
!                     engpetot,nanion,nrpower,extraalpha,extrab
                      engpetot,nanion,nrpower,extraalpha,extrab,num2,numx,nrpower2
!<--- Memmory Reduction_E
!---> Parallelization_S
use mpipara
!<--- Parallelization_E

IMPLICIT NONE

INTEGER :: i,j,ipoint,jpoint
DOUBLE PRECISION :: dr,drrec,drsq,drsqrec,expft,sgam, &
                 fttx,ftty,fttz,expftx,expft2,expft3,expft4,dx,dy,dz, &
                 dxx,dyy,dzz,dxy,dxz,dyz,dr4rec,txx,tyy,tzz,txy,txz,tyz, &
                 dxunit,dyunit,dzunit,sgam2x,sgam2y,sgam2z, &
                 sgam2xx,sgam2yy,sgam2zz,sgam3,drnrec,repsilon,rrtheta, &
                 repsiloni,rrthetai,repsilonj,rrthetaj,repsiloncat,rrthetacat

! Select SR energy to be calculated.
if(epplog) then
   if((.not.(xftlog)).and.(.not.(rvplog))) then
!---> Parallelization_S
!     do j=2,num
      do j=jst,jed
!<--- Parallelization_E
         jpoint=ntype(j)
!---> Parallelization_S
!        do i=1,j-1
         do i=ist(j),ied(j)
!<--- Parallelization_E
            ipoint=ntype(i)

!---> Memmory Reduction_S
            numx = numadr(i,j)
!           drsq=dxsav(i,j)*dxsav(i,j)+dysav(i,j)*dysav(i,j) &
!               +dzsav(i,j)*dzsav(i,j)
            drsq=dxsav(numx)*dxsav(numx)+dysav(numx)*dysav(numx) &
                +dzsav(numx)*dzsav(numx)
!<--- Memmory Reduction_E

            if(drsq.ge.rsqmaxsr) CYCLE

            drsqrec=1.0d0/drsq
            dr=dsqrt(drsq)
            drrec=1.0d0/dr

            expft=ftb(ipoint,jpoint)*exp(-ftalp(ipoint,jpoint)*dr)

            sgam=(expft*ftalp(ipoint,jpoint)*drrec)

!---> Memmory Reduction_S
!           fttx=sgam*dxsav(i,j)
!           ftty=sgam*dysav(i,j)
!           fttz=sgam*dzsav(i,j)
            fttx=sgam*dxsav(numx)
            ftty=sgam*dysav(numx)
            fttz=sgam*dzsav(numx)
!<--- Memmory Reduction_E

!---> Parallelization_S
!           engsr=engsr+expft
            sctmp(59)=sctmp(59)+expft
!<--- Parallelization_E

!---> Memmory Reduction_S
!---> Parallelization_S
!           stsrxx=stsrxx+fttx*dxsav(i,j)
!           stsrxy=stsrxy+fttx*dysav(i,j)
!           stsrxz=stsrxz+fttx*dzsav(i,j)
!           stsryy=stsryy+ftty*dysav(i,j)
!           stsryz=stsryz+ftty*dzsav(i,j)
!           stsrzz=stsrzz+fttz*dzsav(i,j)
            sctmp(1)=sctmp(1)+fttx*dxsav(numx)
            sctmp(2)=sctmp(2)+fttx*dysav(numx)
            sctmp(3)=sctmp(3)+fttx*dzsav(numx)
            sctmp(4)=sctmp(4)+ftty*dysav(numx)
            sctmp(5)=sctmp(5)+ftty*dzsav(numx)
            sctmp(6)=sctmp(6)+fttz*dzsav(numx)
!<--- Parallelization_E
!<--- Memmory Reduction_E

!---> Parallelization_S
!           frrx(j)=frrx(j)-fttx
!           frry(j)=frry(j)-ftty
!           frrz(j)=frrz(j)-fttz
!           frrx(i)=frrx(i)+fttx
!           frry(i)=frry(i)+ftty
!           frrz(i)=frrz(i)+fttz
            eltmp(j)=eltmp(j)-fttx
            eltmp(num+j)=eltmp(num+j)-ftty
            eltmp(2*num+j)=eltmp(2*num+j)-fttz
            eltmp(i)=eltmp(i)+fttx
            eltmp(num+i)=eltmp(num+i)+ftty
            eltmp(2*num+i)=eltmp(2*num+i)+fttz
!<--- Parallelization_E
         enddo   
      enddo   
   endif

   if(rvplog) then
!---> Parallelization_S
!     do j=2,num
      do j=jst,jed
!<--- Parallelization_E
         jpoint=ntype(j)
!---> Parallelization_S
!        do i=1,j-1
         do i=ist(j),ied(j)
!<--- Parallelization_E
            ipoint=ntype(i)

!---> Memmory Reduction_S
            numx = numadr(i,j)
!           drsq=dxsav(i,j)*dxsav(i,j)+dysav(i,j)*dysav(i,j) &
!               +dzsav(i,j)*dzsav(i,j)
            drsq=dxsav(numx)*dxsav(numx)+dysav(numx)*dysav(numx) &
                +dzsav(numx)*dzsav(numx)
!<--- Memmory Reduction_E

            if(drsq.ge.rsqmaxsr) CYCLE

            drsqrec=1.0d0/drsq
            dr=dsqrt(drsq)
            drrec=1.0d0/dr

            expft=(rvph(ipoint,jpoint)/(dr**rvpn(ipoint,jpoint))) &
                 -(rvpr4(ipoint,jpoint)/(dr**4.0d0))

            sgam=((rvpn(ipoint,jpoint)*rvph(ipoint,jpoint) &
             /(dr**(rvpn(ipoint,jpoint)+1.0d0))) &
             -(4.0d0*rvpr4(ipoint,jpoint)/(dr**5.0d0)))*drrec

!---> Memmory Reduction_S
!           fttx=sgam*dxsav(i,j)
!           ftty=sgam*dysav(i,j)
!           fttz=sgam*dzsav(i,j)
            fttx=sgam*dxsav(numx)
            ftty=sgam*dysav(numx)
            fttz=sgam*dzsav(numx)
!<--- Memmory Reduction_E

!---> Parallelization_S
!           engsr=engsr+expft
            sctmp(59)=sctmp(59)+expft
!<--- Parallelization_E

!---> Memmory Reduction_S
!---> Parallelization_S
!           stsrxx=stsrxx+fttx*dxsav(i,j)
!           stsrxy=stsrxy+fttx*dysav(i,j)
!           stsrxz=stsrxz+fttx*dzsav(i,j)
!           stsryy=stsryy+ftty*dysav(i,j)
!           stsryz=stsryz+ftty*dzsav(i,j)
!           stsrzz=stsrzz+fttz*dzsav(i,j)
            sctmp(1)=sctmp(1)+fttx*dxsav(numx)
            sctmp(2)=sctmp(2)+fttx*dysav(numx)
            sctmp(3)=sctmp(3)+fttx*dzsav(numx)
            sctmp(4)=sctmp(4)+ftty*dysav(numx)
            sctmp(5)=sctmp(5)+ftty*dzsav(numx)
            sctmp(6)=sctmp(6)+fttz*dzsav(numx)
!<--- Parallelization_E
!<--- Memmory Reduction_E

!---> Parallelization_S
!           frrx(j)=frrx(j)-fttx
!           frry(j)=frry(j)-ftty
!           frrz(j)=frrz(j)-fttz
!           frrx(i)=frrx(i)+fttx
!           frry(i)=frry(i)+ftty
!           frrz(i)=frrz(i)+fttz
            eltmp(j)=eltmp(j)-fttx
            eltmp(num+j)=eltmp(num+j)-ftty
            eltmp(2*num+j)=eltmp(2*num+j)-fttz
            eltmp(i)=eltmp(i)+fttx
            eltmp(num+i)=eltmp(num+i)+ftty
            eltmp(2*num+i)=eltmp(2*num+i)+fttz
!<--- Parallelization_E
         enddo   
      enddo   
   endif

   if(xftlog) then
!---> Parallelization_S
!     do j=2,num
      do j=jst,jed
!<--- Parallelization_E
         jpoint=ntype(j)
!---> Parallelization_S
!        do i=1,j-1
         do i=ist(j),ied(j)
!<--- Parallelization_E
            ipoint=ntype(i)

!---> Memmory Reduction_S
            numx = numadr(i,j)
!           drsq=dxsav(i,j)*dxsav(i,j)+dysav(i,j)*dysav(i,j) &
!               +dzsav(i,j)*dzsav(i,j)
            drsq=dxsav(numx)*dxsav(numx)+dysav(numx)*dysav(numx) &
                +dzsav(numx)*dzsav(numx)
!<--- Memmory Reduction_E

            if(drsq.ge.rsqmaxsr) CYCLE

            drsqrec=1.0d0/drsq
            dr=dsqrt(drsq)
            drrec=1.0d0/dr

            drnrec=drrec**nrpower(ipoint,jpoint)

            expft=ftb(ipoint,jpoint)*exp(-ftalp(ipoint,jpoint)*dr)*drnrec

            if(ftbx(ipoint,jpoint).gt.0.0d0) then 
               expftx=ftbx(ipoint,jpoint)*exp(-ftalpx(ipoint,jpoint)*drsq)
            else
               expftx=0.d0
            endif

            sgam=expft*ftalp(ipoint,jpoint)*drrec &
               +(expft*nrpower(ipoint,jpoint)*drsqrec) &
               +expftx*2.0d0*ftalpx(ipoint,jpoint)

!---> Memmory Reduction_S
!           fttx=sgam*dxsav(i,j)
!           ftty=sgam*dysav(i,j)
!           fttz=sgam*dzsav(i,j)
            fttx=sgam*dxsav(numx)
            ftty=sgam*dysav(numx)
            fttz=sgam*dzsav(numx)
!<--- Memmory Reduction_E

!---> Parallelization_S
!           engsr=engsr+expftx+expft
            sctmp(59)=sctmp(59)+expftx+expft
!<--- Parallelization_E

!---> Memmory Reduction_S
!---> Parallelization_S
!           stsrxx=stsrxx+fttx*dxsav(i,j)
!           stsrxy=stsrxy+fttx*dysav(i,j)
!           stsrxz=stsrxz+fttx*dzsav(i,j)
!           stsryy=stsryy+ftty*dysav(i,j)
!           stsryz=stsryz+ftty*dzsav(i,j)
!           stsrzz=stsrzz+fttz*dzsav(i,j)
            sctmp(1)=sctmp(1)+fttx*dxsav(numx)
            sctmp(2)=sctmp(2)+fttx*dysav(numx)
            sctmp(3)=sctmp(3)+fttx*dzsav(numx)
            sctmp(4)=sctmp(4)+ftty*dysav(numx)
            sctmp(5)=sctmp(5)+ftty*dzsav(numx)
            sctmp(6)=sctmp(6)+fttz*dzsav(numx)
!<--- Parallelization_E
!<--- Memmory Reduction_E

!---> Parallelization_S
!           frrx(j)=frrx(j)-fttx
!           frry(j)=frry(j)-ftty
!           frrz(j)=frrz(j)-fttz
!           frrx(i)=frrx(i)+fttx
!           frry(i)=frry(i)+ftty
!           frrz(i)=frrz(i)+fttz
            eltmp(j)=eltmp(j)-fttx
            eltmp(num+j)=eltmp(num+j)-ftty
            eltmp(2*num+j)=eltmp(2*num+j)-fttz
            eltmp(i)=eltmp(i)+fttx
            eltmp(num+i)=eltmp(num+i)+ftty
            eltmp(2*num+i)=eltmp(2*num+i)+fttz
!<--- Parallelization_E
         enddo 
      enddo 
   endif
endif

! CIM energy + force calc.
if(cimlog) then
!---> Parallelization_S
!  do j=nsp(1)+1,num
   do j=jst2,jed2
!<--- Parallelization_E
      jpoint=ntype(j)
!---> Parallelization_S
!     do i=1,nsp(1)
      do i=ist2(j),ied2(j)
!<--- Parallelization_E
         ipoint=ntype(i)

!---> Memmory Reduction_S
         numx = numadr(i,j)
!        drsq=dxsav(i,j)*dxsav(i,j)+dysav(i,j)*dysav(i,j)+dzsav(i,j)*dzsav(i,j)
         drsq=dxsav(numx)*dxsav(numx)+dysav(numx)*dysav(numx)+dzsav(numx)*dzsav(numx)
!<--- Memmory Reduction_E

         if(drsq.ge.rsqmaxsr) CYCLE

         drsqrec=1.0d0/drsq
         dr=dsqrt(drsq)
         drrec=1.0d0/dr

         expft=(ftb(ipoint,jpoint)) &
               *exp(-ftalp(ipoint,jpoint)*(dr-delta(i)-delta(j)))
         expft2=(ftb2(ipoint,jpoint)) &
               *exp(-ftbeta(ipoint,jpoint)*(dr-delta(i)-delta(j)))
         expft3=(ftb3(ipoint,jpoint)) &
               *exp(-ftgamma(ipoint,jpoint)*(dr-delta(i)-delta(j)))
         expft4=extrab*exp(-extraalpha*dr)

         if(ftbx(ipoint,jpoint).gt.0.0d0) then 
            expftx=ftbx(ipoint,jpoint)*exp(-ftalpx(ipoint,jpoint)*drsq)
         else
            expftx=0.d0
         endif

         sgam=(ftalp(ipoint,jpoint)*expft) &
             +(ftbeta(ipoint,jpoint)*expft2) &
             +(ftgamma(ipoint,jpoint)*expft3) &
             +(extraalpha*expft4) &
             +expftx*2.0d0*ftalpx(ipoint,jpoint)

!---> Memmory Reduction_S
!        fttx=sgam*(dxsav(i,j)*drrec)
!        ftty=sgam*(dysav(i,j)*drrec)
!        fttz=sgam*(dzsav(i,j)*drrec)
         fttx=sgam*(dxsav(numx)*drrec)
         ftty=sgam*(dysav(numx)*drrec)
         fttz=sgam*(dzsav(numx)*drrec)
!<--- Memmory Reduction_E

!---> Parallelization_S
!        engsr=engsr+expft+expft2+expft3+expft4
         sctmp(59)=sctmp(59)+expft+expft2+expft3+expft4+expftx
!<--- Parallelization_E

!---> Memmory Reduction_S
!---> Parallelization_S
!        stsrxx=stsrxx+fttx*dxsav(i,j)
!        stsrxy=stsrxy+fttx*dysav(i,j)
!        stsrxz=stsrxz+fttx*dzsav(i,j)
!        stsryy=stsryy+ftty*dysav(i,j)
!        stsryz=stsryz+ftty*dzsav(i,j)
!        stsrzz=stsrzz+fttz*dzsav(i,j)
!        stsrxx=stsrxx+fttx*dxsav(numx)
!        stsrxy=stsrxy+fttx*dysav(numx)
!        stsrxz=stsrxz+fttx*dzsav(numx)
!        stsryy=stsryy+ftty*dysav(numx)
!        stsryz=stsryz+ftty*dzsav(numx)
!        stsrzz=stsrzz+fttz*dzsav(numx)
         sctmp(1)=sctmp(1)+fttx*dxsav(numx)
         sctmp(2)=sctmp(2)+fttx*dysav(numx)
         sctmp(3)=sctmp(3)+fttx*dzsav(numx)
         sctmp(4)=sctmp(4)+ftty*dysav(numx)
         sctmp(5)=sctmp(5)+ftty*dzsav(numx)
         sctmp(6)=sctmp(6)+fttz*dzsav(numx)
!<--- Parallelization_E
!<--- Memmory Reduction_E

!---> Parallelization_S
!        frrx(j)=frrx(j)-fttx
!        frry(j)=frry(j)-ftty
!        frrz(j)=frrz(j)-fttz
!        frrx(i)=frrx(i)+fttx
!        frry(i)=frry(i)+ftty
!        frrz(i)=frrz(i)+fttz
            eltmp(j)=eltmp(j)-fttx
            eltmp(num+j)=eltmp(num+j)-ftty
            eltmp(2*num+j)=eltmp(2*num+j)-fttz
            eltmp(i)=eltmp(i)+fttx
            eltmp(num+i)=eltmp(num+i)+ftty
            eltmp(2*num+i)=eltmp(2*num+i)+fttz
!<--- Parallelization_E
      enddo   
   enddo   
endif
! DAIM energy + force calc.
if(daimlog) then
!---> Parallelization_S
!  do j=nsp(1)+1,num
   do j=jst2,jed2
!<--- Parallelization_E
      jpoint=ntype(j)
!---> Parallelization_S
!     do i=1,nsp(1)
      do i=ist2(j),ied2(j)
!<--- Parallelization_E
         ipoint=ntype(i)

!---> Memmory Reduction_S
         numx = numadr(i,j)
!        drsq=dxsav(i,j)*dxsav(i,j)+dysav(i,j)*dysav(i,j)+dzsav(i,j)*dzsav(i,j)
         drsq=dxsav(numx)*dxsav(numx)+dysav(numx)*dysav(numx)+dzsav(numx)*dzsav(numx)
!<--- Memmory Reduction_E

         if(drsq.ge.rsqmaxsr) CYCLE

         drsqrec=1.0d0/drsq
         dr=dsqrt(drsq)
         drrec=1.0d0/dr
!---> Memmory Reduction_S
!        dx=dxsav(i,j)
!        dy=dysav(i,j)
!        dz=dzsav(i,j)
         dx=dxsav(numx)
         dy=dysav(numx)
         dz=dzsav(numx)
!<--- Memmory Reduction_E

         repsilon=epsilonx(i)*dx+epsilony(i)*dy+epsilonz(i)*dz
         repsiloncat=epsilonx(j)*dx+epsilony(j)*dy+epsilonz(j)*dz
         repsilon=repsilon*drrec
         repsiloncat=repsiloncat*drrec

         expft=ftb(ipoint,jpoint) &
               *exp(-ftalp(ipoint,jpoint)*(dr-delta(i)-delta(j) &
                 -repsilon+repsiloncat))
         expft2=ftb2(ipoint,jpoint) &
               *exp(-ftbeta(ipoint,jpoint)*(dr-delta(i)-delta(j) &
                 -repsilon+repsiloncat))
         expft3=ftb3(ipoint,jpoint) &
               *exp(-ftgamma(ipoint,jpoint)*(dr-delta(i)-delta(j) &
                 -repsilon+repsiloncat))
         expft4=extrab*exp(-extraalpha*dr)


         if(ftbx(ipoint,jpoint).gt.0.0d0) then 
            expftx=ftbx(ipoint,jpoint)*exp(-ftalpx(ipoint,jpoint)*drsq)
         else
            expftx=0.d0
         endif

         dxunit=dx*drrec
         dyunit=dy*drrec
         dzunit=dz*drrec

         sgam=ftalp(ipoint,jpoint)*expft &
               +ftbeta(ipoint,jpoint)*expft2 &
               +ftgamma(ipoint,jpoint)*expft3 &
               +expftx*2.0d0*ftalpx(ipoint,jpoint)
         sgam2x=drrec*(epsilonx(i)-epsilonx(j))- &
                dx*drsqrec*(repsilon-repsiloncat)
         sgam2y=drrec*(epsilony(i)-epsilony(j))- &
                dy*drsqrec*(repsilon-repsiloncat)
         sgam2z=drrec*(epsilonz(i)-epsilonz(j))- &
                dz*drsqrec*(repsilon-repsiloncat)
         sgam3=extraalpha*drrec*expft4



         fttx=sgam*(dx*drrec-sgam2x)+sgam3*dx
         ftty=sgam*(dy*drrec-sgam2y)+sgam3*dy
         fttz=sgam*(dz*drrec-sgam2z)+sgam3*dz

!---> Parallelization_S
!        engsr=engsr+expft+expft2+expft3+expft4
         sctmp(59)=sctmp(59)+expft+expft2+expft3+expft4+expftx
!<--- Parallelization_E

!---> Memmory Reduction_S
!---> Parallelization_S
!        stsrxx=stsrxx+fttx*dxsav(i,j)
!        stsrxy=stsrxy+fttx*dysav(i,j)
!        stsrxz=stsrxz+fttx*dzsav(i,j)
!        stsryy=stsryy+ftty*dysav(i,j)
!        stsryz=stsryz+ftty*dzsav(i,j)
!        stsrzz=stsrzz+fttz*dzsav(i,j)
!        stsrxx=stsrxx+fttx*dxsav(numx)
!        stsrxy=stsrxy+fttx*dysav(numx)
!        stsrxz=stsrxz+fttx*dzsav(numx)
!        stsryy=stsryy+ftty*dysav(numx)
!        stsryz=stsryz+ftty*dzsav(numx)
!        stsrzz=stsrzz+fttz*dzsav(numx)
         sctmp(1)=sctmp(1)+fttx*dxsav(numx)
         sctmp(2)=sctmp(2)+fttx*dysav(numx)
         sctmp(3)=sctmp(3)+fttx*dzsav(numx)
         sctmp(4)=sctmp(4)+ftty*dysav(numx)
         sctmp(5)=sctmp(5)+ftty*dzsav(numx)
         sctmp(6)=sctmp(6)+fttz*dzsav(numx)
!<--- Parallelization_E
!<--- Memmory Reduction_E

!---> Parallelization_S
!        frrx(j)=frrx(j)-fttx
!        frry(j)=frry(j)-ftty
!        frrz(j)=frrz(j)-fttz
!        frrx(i)=frrx(i)+fttx
!        frry(i)=frry(i)+ftty
!        frrz(i)=frrz(i)+fttz
            eltmp(j)=eltmp(j)-fttx
            eltmp(num+j)=eltmp(num+j)-ftty
            eltmp(2*num+j)=eltmp(2*num+j)-fttz
            eltmp(i)=eltmp(i)+fttx
            eltmp(num+i)=eltmp(num+i)+ftty
            eltmp(2*num+i)=eltmp(2*num+i)+fttz
!<--- Parallelization_E
      enddo   
   enddo   
endif
! QAIM energy + force calc.
if(quaimlog) then
!---> Parallelization_S
!  do j=nsp(1)+1,num
   do j=jst2,jed2
!<--- Parallelization_E
      jpoint=ntype(j)
!---> Parallelization_S
!     do i=1,nsp(1)
      do i=ist2(j),ied2(j)
!<--- Parallelization_E
         ipoint=ntype(i)

!---> Memmory Reduction_S
         numx = numadr(i,j)
!        drsq=dxsav(i,j)*dxsav(i,j)+dysav(i,j)*dysav(i,j)+dzsav(i,j)*dzsav(i,j)
         drsq=dxsav(numx)*dxsav(numx)+dysav(numx)*dysav(numx)+dzsav(numx)*dzsav(numx)
!<--- Memmory Reduction_E

         if(drsq.ge.rsqmaxsr) CYCLE

         drsqrec=1.0d0/drsq
         dr=dsqrt(drsq)
         drrec=1.0d0/dr
         dr4rec=drsqrec*drsqrec
!---> Memmory Reduction_S
!        dx=dxsav(i,j)
!        dy=dysav(i,j)
!        dz=dzsav(i,j)
         dx=dxsav(numx)
         dy=dysav(numx)
         dz=dzsav(numx)
!<--- Memmory Reduction_E

         repsilon=epsilonx(i)*dx+epsilony(i)*dy+epsilonz(i)*dz
         repsiloncat=epsilonx(j)*dx+epsilony(j)*dy+epsilonz(j)*dz
         repsilon=repsilon*drrec
         repsiloncat=repsiloncat*drrec

         dxx=dx*dx
         dyy=dy*dy
         dzz=dz*dz
         dxy=dx*dy
         dxz=dx*dz
         dyz=dy*dz
         txx=3.0d0*dxx*drsqrec-1.0d0
         tyy=3.0d0*dyy*drsqrec-1.0d0
         tzz=3.0d0*dzz*drsqrec-1.0d0
         txy=3.0d0*dxy*drsqrec
         txz=3.0d0*dxz*drsqrec
         tyz=3.0d0*dyz*drsqrec

         rrtheta=quaimxx(i)*txx+quaimyy(i)*tyy+quaimzz(i)*tzz+ &
              2.0d0*(quaimxy(i)*txy+quaimxz(i)*txz+quaimyz(i)*tyz)
         rrthetacat=quaimxx(j)*txx+quaimyy(j)*tyy+ &
                    quaimzz(j)*tzz+2.0d0*(quaimxy(j)*txy+ &
                    quaimxz(j)*txz+quaimyz(j)*tyz)

         expft=ftb(ipoint,jpoint) &
              *exp(-ftalp(ipoint,jpoint)*(dr-delta(i)-delta(j)-repsilon &
           +repsiloncat-rrtheta-rrthetacat))
         expft2=ftb2(ipoint,jpoint) &
              *exp(-ftbeta(ipoint,jpoint)*(dr-delta(i)-delta(j)-repsilon &
           +repsiloncat-rrtheta-rrthetacat))
         expft3=ftb3(ipoint,jpoint) &
              *exp(-ftgamma(ipoint,jpoint)*(dr-delta(i)-delta(j)-repsilon &
           +repsiloncat-rrtheta-rrthetacat))
         expft4=extrab*exp(-extraalpha*dr)

         if(ftbx(ipoint,jpoint).gt.0.0d0) then 
            expftx=ftbx(ipoint,jpoint)*exp(-ftalpx(ipoint,jpoint)*drsq)
         else
            expftx=0.d0
         endif

         dxunit=dx*drrec
         dyunit=dy*drrec
         dzunit=dz*drrec

         sgam=ftalp(ipoint,jpoint)*expft &
               +ftbeta(ipoint,jpoint)*expft2 &
               +ftgamma(ipoint,jpoint)*expft3 &
               +expftx*2.0d0*ftalpx(ipoint,jpoint)
         sgam2x=drrec*(epsilonx(i)-epsilonx(j))- &
                 dx*drsqrec*(repsilon-repsiloncat)
         sgam2y=drrec*(epsilony(i)-epsilony(j))- &
                 dy*drsqrec*(repsilon-repsiloncat)
         sgam2z=drrec*(epsilonz(i)-epsilonz(j))- &
                 dz*drsqrec*(repsilon-repsiloncat)

          sgam2xx=(quaimxx(j)+quaimxx(i)) &
             *(6.0d0*dx*drsqrec-6.0d0*dxx*dx*dr4rec)- &
              (quaimyy(j)+quaimyy(i))*6.0d0*dx*dyy*dr4rec- &
              (quaimzz(j)+quaimzz(i))*6.0d0*dx*dzz*dr4rec+ &
              (quaimxy(j)+quaimxy(i))* &
              (6.0d0*dy*drsqrec-12.0d0*dy*dxx*dr4rec)+ &
              (quaimxz(j)+quaimxz(i))* &
              (6.0d0*dz*drsqrec-12.0d0*dz*dxx*dr4rec)- &
              (quaimyz(j)+quaimyz(i))*12.0d0*dx*dy*dz*dr4rec

          sgam2yy=(quaimyy(j)+quaimyy(i))*(6.0d0*dy*drsqrec &
             -6.0d0*dyy*dy*dr4rec)- &
              (quaimxx(j)+quaimxx(i))*6.0d0*dy*dxx*dr4rec- &
              (quaimzz(j)+quaimzz(i))*6.0d0*dy*dzz*dr4rec+ &
              (quaimxy(j)+quaimxy(i))* &
              (6.0d0*dx*drsqrec-12.0d0*dx*dyy*dr4rec)+ &
              (quaimyz(j)+quaimyz(i))* &
              (6.0d0*dz*drsqrec-12.0d0*dz*dyy*dr4rec)- &
              (quaimxz(j)+quaimxz(i))*12.0d0*dx*dy*dz*dr4rec

          sgam2zz=(quaimzz(j)+quaimzz(i))*(6.0d0*dz*drsqrec &
             -6.0d0*dzz*dz*dr4rec)- &
              (quaimxx(j)+quaimxx(i))*6.0d0*dz*dxx*dr4rec- &
              (quaimyy(j)+quaimyy(i))*6.0d0*dz*dyy*dr4rec+ &
              (quaimxz(j)+quaimxz(i))* &
              (6.0d0*dx*drsqrec-12.0d0*dx*dzz*dr4rec)+ &
              (quaimyz(j)+quaimyz(i))* &
              (6.0d0*dy*drsqrec-12.0d0*dy*dzz*dr4rec)- &
              (quaimxy(j)+quaimxy(i))*12.0d0*dx*dy*dz*dr4rec

         sgam3=extraalpha*drrec*expft4

         fttx=sgam*(dx*drrec-sgam2x-sgam2xx)+sgam3*dx
         ftty=sgam*(dy*drrec-sgam2y-sgam2yy)+sgam3*dy
         fttz=sgam*(dz*drrec-sgam2z-sgam2zz)+sgam3*dz

!---> Parallelization_S
!        engsr=engsr+expft+expft2+expft3+expft4
         sctmp(59)=sctmp(59)+expft+expft2+expft3+expft4+expftx
!<--- Parallelization_E

!---> Parallelization_S
!        frrx(j)=frrx(j)-fttx
!        frry(j)=frry(j)-ftty
!        frrz(j)=frrz(j)-fttz
!        frrx(i)=frrx(i)+fttx
!        frry(i)=frry(i)+ftty
!        frrz(i)=frrz(i)+fttz
            eltmp(j)=eltmp(j)-fttx
            eltmp(num+j)=eltmp(num+j)-ftty
            eltmp(2*num+j)=eltmp(2*num+j)-fttz
            eltmp(i)=eltmp(i)+fttx
            eltmp(num+i)=eltmp(num+i)+ftty
            eltmp(2*num+i)=eltmp(2*num+i)+fttz
!<--- Parallelization_E

!---> Memmory Reduction_S
!---> Parallelization_S
!        stsrxx=stsrxx+fttx*dxsav(i,j)
!        stsrxy=stsrxy+fttx*dysav(i,j)
!        stsrxz=stsrxz+fttx*dzsav(i,j)
!        stsryy=stsryy+ftty*dysav(i,j)
!        stsryz=stsryz+ftty*dzsav(i,j)
!        stsrzz=stsrzz+fttz*dzsav(i,j)
!        stsrxx=stsrxx+fttx*dxsav(numx)
!        stsrxy=stsrxy+fttx*dysav(numx)
!        stsrxz=stsrxz+fttx*dzsav(numx)
!        stsryy=stsryy+ftty*dysav(numx)
!        stsryz=stsryz+ftty*dzsav(numx)
!        stsrzz=stsrzz+fttz*dzsav(numx)
         sctmp(1)=sctmp(1)+fttx*dxsav(numx)
         sctmp(2)=sctmp(2)+fttx*dysav(numx)
         sctmp(3)=sctmp(3)+fttx*dzsav(numx)
         sctmp(4)=sctmp(4)+ftty*dysav(numx)
         sctmp(5)=sctmp(5)+ftty*dzsav(numx)
         sctmp(6)=sctmp(6)+fttz*dzsav(numx)
!<--- Parallelization_E
!<--- Memmory Reduction_E
      enddo   
   enddo   
endif
! anion-anion EPP
! TMP - CIM, AIM.
if(ooaimlog) then
   if(cimlog) then
!---> Parallelization_S
!     do j=2,nanion
      do j=jst3,jed3
!<--- Parallelization_E
         jpoint=ntype(j)
!---> Parallelization_S
!        do i=1,j-1
         do i=ist3(j),ied3(j)
!<--- Parallelization_E
            ipoint=ntype(i)

!---> Memmory Reduction_S
            numx = numadr(i,j)
!           drsq=dxsav(i,j)*dxsav(i,j)+dysav(i,j)*dysav(i,j) &
!               +dzsav(i,j)*dzsav(i,j)
            drsq=dxsav(numx)*dxsav(numx)+dysav(numx)*dysav(numx) &
                +dzsav(numx)*dzsav(numx)
!<--- Memmory Reduction_E

            if(drsq.ge.rsqmaxsr) CYCLE

            drsqrec=1.0d0/drsq
            dr=dsqrt(drsq)
            drrec=1.0d0/dr

            expft=ftb(ipoint,jpoint)*exp(-ftalp(ipoint,jpoint) &
                 *(dr-delta(i)-delta(j)))
            expft2=ftb2(ipoint,jpoint)*exp(-ftbeta(ipoint,jpoint) &
                 *(dr-delta(i)-delta(j)))
            expft3=ftb3(ipoint,jpoint)*exp(-ftgamma(ipoint,jpoint) &
                 *(dr-delta(i)-delta(j)))

            if(ftbx(ipoint,jpoint).gt.0.0d0) then 
               expftx=ftbx(ipoint,jpoint)*exp(-ftalpx(ipoint,jpoint)*drsq)
            else
               expftx=0.d0
            endif

            sgam=expft*ftalp(ipoint,jpoint)+ &
                 expft2*ftbeta(ipoint,jpoint)+ &
                 expft3*ftgamma(ipoint,jpoint) &
                 +expftx*2.0d0*ftalpx(ipoint,jpoint)

            fttx=sgam*(dx*drrec)
            ftty=sgam*(dy*drrec)
            fttz=sgam*(dz*drrec)

!---> Parallelization_S
!           engsr=engsr+expft+expft2+expft3
            sctmp(59)=sctmp(59)+expft+expft2+expft3+expftx
!<--- Parallelization_E

!---> Memmory Reduction_S
!---> Parallelization_S
!           stsrxx=stsrxx+fttx*dxsav(i,j)
!           stsrxy=stsrxy+fttx*dysav(i,j)
!           stsrxz=stsrxz+fttx*dzsav(i,j)
!           stsryy=stsryy+ftty*dysav(i,j)
!           stsryz=stsryz+ftty*dzsav(i,j)
!           stsrzz=stsrzz+fttz*dzsav(i,j)
!        stsrxx=stsrxx+fttx*dxsav(numx)
!        stsrxy=stsrxy+fttx*dysav(numx)
!        stsrxz=stsrxz+fttx*dzsav(numx)
!        stsryy=stsryy+ftty*dysav(numx)
!        stsryz=stsryz+ftty*dzsav(numx)
!        stsrzz=stsrzz+fttz*dzsav(numx)
            sctmp(1)=sctmp(1)+fttx*dxsav(numx)
            sctmp(2)=sctmp(2)+fttx*dysav(numx)
            sctmp(3)=sctmp(3)+fttx*dzsav(numx)
            sctmp(4)=sctmp(4)+ftty*dysav(numx)
            sctmp(5)=sctmp(5)+ftty*dzsav(numx)
            sctmp(6)=sctmp(6)+fttz*dzsav(numx)
!<--- Parallelization_E
!<--- Memmory Reduction_E

!---> Parallelization_S
!           frrx(j)=frrx(j)-fttx
!           frry(j)=frry(j)-ftty
!           frrz(j)=frrz(j)-fttz
!           frrx(i)=frrx(i)+fttx
!           frry(i)=frry(i)+ftty
!           frrz(i)=frrz(i)+fttz
            eltmp(j)=eltmp(j)-fttx
            eltmp(num+j)=eltmp(num+j)-ftty
            eltmp(2*num+j)=eltmp(2*num+j)-fttz
            eltmp(i)=eltmp(i)+fttx
            eltmp(num+i)=eltmp(num+i)+ftty
            eltmp(2*num+i)=eltmp(2*num+i)+fttz
!<--- Parallelization_E
         enddo   
      enddo   
   else if(daimlog) then
!---> Parallelization_S
!     do j=2,nanion
      do j=jst3,jed3
!<--- Parallelization_E
         jpoint=ntype(j)
!---> Parallelization_S
!        do i=1,j-1
         do i=ist3(j),ied3(j)
!<--- Parallelization_E
            ipoint=ntype(i)

!---> Memmory Reduction_S
            numx = numadr(i,j)
!           drsq=dxsav(i,j)*dxsav(i,j)+dysav(i,j)*dysav(i,j) &
!               +dzsav(i,j)*dzsav(i,j)
            drsq=dxsav(numx)*dxsav(numx)+dysav(numx)*dysav(numx) &
                +dzsav(numx)*dzsav(numx)
!<--- Memmory Reduction_E

            if(drsq.ge.rsqmaxsr) CYCLE

            drsqrec=1.0d0/drsq
            dr=dsqrt(drsq)
            drrec=1.0d0/dr
!---> Memmory Reduction_S
!           dx=dxsav(i,j)
!           dy=dysav(i,j)
!           dz=dzsav(i,j)
            dx=dxsav(numx)
            dy=dysav(numx)
            dz=dzsav(numx)
!<--- Memmory Reduction_E

            repsiloni=epsilonx(i)*dx+epsilony(i)*dy+epsilonz(i)*dz
            repsilonj=epsilonx(j)*dx+epsilony(j)*dy+epsilonz(j)*dz
            repsiloni=repsiloni*drrec
            repsilonj=repsilonj*drrec

            expft=ftb(ipoint,jpoint)*exp(-ftalp(ipoint,jpoint) &
                *(dr-delta(i)-delta(j)-repsiloni+repsilonj))
            expft2=ftb2(ipoint,jpoint)*exp(-ftbeta(ipoint,jpoint) &
                *(dr-delta(i)-delta(j)-repsiloni+repsilonj))
            expft3=ftb3(ipoint,jpoint)*exp(-ftgamma(ipoint,jpoint) &
                *(dr-delta(i)-delta(j)-repsiloni+repsilonj))

            if(ftbx(ipoint,jpoint).gt.0.0d0) then 
               expftx=ftbx(ipoint,jpoint)*exp(-ftalpx(ipoint,jpoint)*drsq)
            else
               expftx=0.d0
            endif

            dxunit=dx*drrec
            dyunit=dy*drrec
            dzunit=dz*drrec

            sgam=expft*ftalp(ipoint,jpoint) &
                +expft2*ftbeta(ipoint,jpoint) &
                +expft3*ftgamma(ipoint,jpoint) &
                +expftx*2.0d0*ftalpx(ipoint,jpoint)
            sgam2x=drrec*(epsilonx(i)-epsilonx(j)) &
                  -dx*drsqrec*(repsiloni-repsilonj)
            sgam2y=drrec*(epsilony(i)-epsilony(j)) &
                  -dy*drsqrec*(repsiloni-repsilonj)
            sgam2z=drrec*(epsilonz(i)-epsilonz(j)) &
                  -dz*drsqrec*(repsiloni-repsilonj)

            fttx=sgam*(dx*drrec-sgam2x)
            ftty=sgam*(dy*drrec-sgam2y)
            fttz=sgam*(dz*drrec-sgam2z)

!---> Parallelization_S
!           engsr=engsr+expft+expft2+expft3
            sctmp(59)=sctmp(59)+expft+expft2+expft3+expftx
!<--- Parallelization_E

!---> Memmory Reduction_S
!---> Parallelization_S
!           stsrxx=stsrxx+fttx*dxsav(i,j)
!           stsrxy=stsrxy+fttx*dysav(i,j)
!           stsrxz=stsrxz+fttx*dzsav(i,j)
!           stsryy=stsryy+ftty*dysav(i,j)
!           stsryz=stsryz+ftty*dzsav(i,j)
!           stsrzz=stsrzz+fttz*dzsav(i,j)
!        stsrxx=stsrxx+fttx*dxsav(numx)
!        stsrxy=stsrxy+fttx*dysav(numx)
!        stsrxz=stsrxz+fttx*dzsav(numx)
!        stsryy=stsryy+ftty*dysav(numx)
!        stsryz=stsryz+ftty*dzsav(numx)
!        stsrzz=stsrzz+fttz*dzsav(numx)
            sctmp(1)=sctmp(1)+fttx*dxsav(numx)
            sctmp(2)=sctmp(2)+fttx*dysav(numx)
            sctmp(3)=sctmp(3)+fttx*dzsav(numx)
            sctmp(4)=sctmp(4)+ftty*dysav(numx)
            sctmp(5)=sctmp(5)+ftty*dzsav(numx)
            sctmp(6)=sctmp(6)+fttz*dzsav(numx)
!<--- Parallelization_E
!<--- Memmory Reduction_E

!---> Parallelization_S
!           frrx(j)=frrx(j)-fttx
!           frry(j)=frry(j)-ftty
!           frrz(j)=frrz(j)-fttz
!           frrx(i)=frrx(i)+fttx
!           frry(i)=frry(i)+ftty
!           frrz(i)=frrz(i)+fttz
            eltmp(j)=eltmp(j)-fttx
            eltmp(num+j)=eltmp(num+j)-ftty
            eltmp(2*num+j)=eltmp(2*num+j)-fttz
            eltmp(i)=eltmp(i)+fttx
            eltmp(num+i)=eltmp(num+i)+ftty
            eltmp(2*num+i)=eltmp(2*num+i)+fttz
!<--- Parallelization_E
         enddo   
      enddo   
   else if(quaimlog) then
!---> Parallelization_S
!     do j=2,nanion
      do j=jst3,jed3
!<--- Parallelization_E
         jpoint=ntype(j)
!---> Parallelization_S
!        do i=1,j-1
         do i=ist3(j),ied3(j)
!<--- Parallelization_E
            ipoint=ntype(i)

!---> Memmory Reduction_S
            numx = numadr(i,j)
!           drsq=dxsav(i,j)*dxsav(i,j)+dysav(i,j)*dysav(i,j) &
!               +dzsav(i,j)*dzsav(i,j)
            drsq=dxsav(numx)*dxsav(numx)+dysav(numx)*dysav(numx) &
                +dzsav(numx)*dzsav(numx)
!<--- Memmory Reduction_E

            if(drsq.ge.rsqmaxsr) CYCLE

            drsqrec=1.0d0/drsq
            dr=dsqrt(drsq)
            drrec=1.0d0/dr
            dr4rec=drsqrec*drsqrec
!---> Memmory Reduction_S
!           dx=dxsav(i,j)
!           dy=dysav(i,j)
!           dz=dzsav(i,j)
            dx=dxsav(numx)
            dy=dysav(numx)
            dz=dzsav(numx)
!<--- Memmory Reduction_E

            repsiloni=epsilonx(i)*dx+epsilony(i)*dy+epsilonz(i)*dz
            repsilonj=epsilonx(j)*dx+epsilony(j)*dy+epsilonz(j)*dz
            repsiloni=repsiloni*drrec
            repsilonj=repsilonj*drrec

            dxx=dx*dx
            dyy=dy*dy
            dzz=dz*dz
            dxy=dx*dy
            dxz=dx*dz
            dyz=dy*dz
            txx=3.0d0*dxx*drsqrec-1.0d0
            tyy=3.0d0*dyy*drsqrec-1.0d0
            tzz=3.0d0*dzz*drsqrec-1.0d0
            txy=3.0d0*dxy*drsqrec
            txz=3.0d0*dxz*drsqrec
            tyz=3.0d0*dyz*drsqrec

            rrthetai=quaimxx(i)*txx+quaimyy(i)*tyy+quaimzz(i)*tzz+ &
                 2.0d0*(quaimxy(i)*txy+quaimxz(i)*txz+quaimyz(i)*tyz)
            rrthetaj=quaimxx(j)*txx+quaimyy(j)*tyy+quaimzz(j)*tzz+ &
                 2.0d0*(quaimxy(j)*txy+quaimxz(j)*txz+quaimyz(j)*tyz)

            expft=ftb(ipoint,jpoint)*exp(-ftalp(ipoint,jpoint) &
                *(dr-delta(i)-delta(j)-repsiloni+repsilonj &
                  -rrthetaj-rrthetai))
            expft2=ftb2(ipoint,jpoint)*exp(-ftbeta(ipoint,jpoint) &
                *(dr-delta(i)-delta(j)-repsiloni+repsilonj &
                  -rrthetaj-rrthetai))
            expft3=ftb3(ipoint,jpoint)*exp(-ftgamma(ipoint,jpoint) &
                *(dr-delta(i)-delta(j)-repsiloni+repsilonj &
                  -rrthetaj-rrthetai))

            if(ftbx(ipoint,jpoint).gt.0.0d0) then 
               expftx=ftbx(ipoint,jpoint)*exp(-ftalpx(ipoint,jpoint)*drsq)
            else
               expftx=0.d0
            endif

            dxunit=dx*drrec
            dyunit=dy*drrec
            dzunit=dz*drrec

            sgam=expft*ftalp(ipoint,jpoint) &
                +expft2*ftbeta(ipoint,jpoint) &
                +expft3*ftgamma(ipoint,jpoint) &
                +expftx*2.0d0*ftalpx(ipoint,jpoint)
            sgam2x=drrec*(epsilonx(i)-epsilonx(j)) &
                  -dx*drsqrec*(repsiloni-repsilonj)
            sgam2y=drrec*(epsilony(i)-epsilony(j)) &
                  -dy*drsqrec*(repsiloni-repsilonj)
            sgam2z=drrec*(epsilonz(i)-epsilonz(j)) &
                  -dz*drsqrec*(repsiloni-repsilonj)

             sgam2xx=(quaimxx(j)+quaimxx(i)) &
                *(6.0d0*dx*drsqrec-6.0d0*dxx*dx*dr4rec)- &
                 (quaimyy(j)+quaimyy(i))*6.0d0*dx*dyy*dr4rec- &
                 (quaimzz(j)+quaimzz(i))*6.0d0*dx*dzz*dr4rec+ &
                 (quaimxy(j)+quaimxy(i))* &
                     (6.0d0*dy*drsqrec-12.0d0*dy*dxx*dr4rec)+ &
                 (quaimxz(j)+quaimxz(i))* &
                     (6.0d0*dz*drsqrec-12.0d0*dz*dxx*dr4rec)- &
                 (quaimyz(j)+quaimyz(i))*12.0d0*dx*dy*dz*dr4rec

             sgam2yy=(quaimyy(j)+quaimyy(i))*(6.0d0*dy*drsqrec &
                -6.0d0*dyy*dy*dr4rec)- &
                 (quaimxx(j)+quaimxx(i))*6.0d0*dy*dxx*dr4rec- &
                 (quaimzz(j)+quaimzz(i))*6.0d0*dy*dzz*dr4rec+ &
                 (quaimxy(j)+quaimxy(i))* &
                     (6.0d0*dx*drsqrec-12.0d0*dx*dyy*dr4rec)+ &
                 (quaimyz(j)+quaimyz(i))* &
                     (6.0d0*dz*drsqrec-12.0d0*dz*dyy*dr4rec)- &
                 (quaimxz(j)+quaimxz(i))*12.0d0*dx*dy*dz*dr4rec

             sgam2zz=(quaimzz(j)+quaimzz(i))*(6.0d0*dz*drsqrec &
                -6.0d0*dzz*dz*dr4rec)- &
                 (quaimxx(j)+quaimxx(i))*6.0d0*dz*dxx*dr4rec- &
                 (quaimyy(j)+quaimyy(i))*6.0d0*dz*dyy*dr4rec+ &
                 (quaimxz(j)+quaimxz(i))* &
                     (6.0d0*dx*drsqrec-12.0d0*dx*dzz*dr4rec)+ &
                 (quaimyz(j)+quaimyz(i))* &
                     (6.0d0*dy*drsqrec-12.0d0*dy*dzz*dr4rec)- &
                 (quaimxy(j)+quaimxy(i))*12.0d0*dx*dy*dz*dr4rec

            fttx=sgam*(dx*drrec-sgam2x-sgam2xx)
            ftty=sgam*(dy*drrec-sgam2y-sgam2yy)
            fttz=sgam*(dz*drrec-sgam2z-sgam2zz)

!---> Parallelization_S
!           engsr=engsr+expft+expft2+expft3
            sctmp(59)=sctmp(59)+expft+expft2+expft3+expftx
!<--- Parallelization_E

!---> Memmory Reduction_S
!---> Parallelization_S
!           stsrxx=stsrxx+fttx*dxsav(i,j)
!           stsrxy=stsrxy+fttx*dysav(i,j)
!           stsrxz=stsrxz+fttx*dzsav(i,j)
!           stsryy=stsryy+ftty*dysav(i,j)
!           stsryz=stsryz+ftty*dzsav(i,j)
!           stsrzz=stsrzz+fttz*dzsav(i,j)
!        stsrxx=stsrxx+fttx*dxsav(numx)
!        stsrxy=stsrxy+fttx*dysav(numx)
!        stsrxz=stsrxz+fttx*dzsav(numx)
!        stsryy=stsryy+ftty*dysav(numx)
!        stsryz=stsryz+ftty*dzsav(numx)
!        stsrzz=stsrzz+fttz*dzsav(numx)
            sctmp(1)=sctmp(1)+fttx*dxsav(numx)
            sctmp(2)=sctmp(2)+fttx*dysav(numx)
            sctmp(3)=sctmp(3)+fttx*dzsav(numx)
            sctmp(4)=sctmp(4)+ftty*dysav(numx)
            sctmp(5)=sctmp(5)+ftty*dzsav(numx)
            sctmp(6)=sctmp(6)+fttz*dzsav(numx)
!<--- Parallelization_E
!<--- Memmory Reduction_E

!---> Parallelization_S
!           frrx(j)=frrx(j)-fttx
!           frry(j)=frry(j)-ftty
!           frrz(j)=frrz(j)-fttz
!           frrx(i)=frrx(i)+fttx
!           frry(i)=frry(i)+ftty
!           frrz(i)=frrz(i)+fttz
            eltmp(j)=eltmp(j)-fttx
            eltmp(num+j)=eltmp(num+j)-ftty
            eltmp(2*num+j)=eltmp(2*num+j)-fttz
            eltmp(i)=eltmp(i)+fttx
            eltmp(num+i)=eltmp(num+i)+ftty
            eltmp(2*num+i)=eltmp(2*num+i)+fttz
!<--- Parallelization_E
         enddo   
      enddo   
   endif
else if(.not.epplog)then
!---> Parallelization_S
!  do j=2,nanion
   do j=jst3,jed3
!<--- Parallelization_E
      jpoint=ntype(j)
!---> Parallelization_S
!     do i=1,j-1
      do i=ist3(j),ied3(j)
!<--- Parallelization_E
         ipoint=ntype(i)

!---> Memmory Reduction_S
         numx = numadr(i,j)
!        drsq=dxsav(i,j)*dxsav(i,j)+dysav(i,j)*dysav(i,j)+dzsav(i,j)*dzsav(i,j)
         drsq=dxsav(numx)*dxsav(numx)+dysav(numx)*dysav(numx)+dzsav(numx)*dzsav(numx)
!<--- Memmory Reduction_E

         if(drsq.ge.rsqmaxsr) CYCLE

         drsqrec=1.0d0/drsq
         dr=dsqrt(drsq)
         drrec=1.0d0/dr

         expft=ftb(ipoint,jpoint)*exp(-ftalp(ipoint,jpoint)*dr)

         if(ftbx(ipoint,jpoint).gt.0.0d0) then 
            expftx=ftbx(ipoint,jpoint)*exp(-ftalpx(ipoint,jpoint)*drsq)
         else
            expftx=0.d0
         endif

         sgam=(expft*ftalp(ipoint,jpoint)*drrec) &
              +expftx*2.0d0*ftalpx(ipoint,jpoint)

!---> Memmory Reduction_S
!        fttx=sgam*dxsav(i,j)
!        ftty=sgam*dysav(i,j)
!        fttz=sgam*dzsav(i,j)
         fttx=sgam*dxsav(numx)
         ftty=sgam*dysav(numx)
         fttz=sgam*dzsav(numx)
!<--- Memmory Reduction_E

!---> Parallelization_S
!        engsr=engsr+expft
         sctmp(59)=sctmp(59)+expft+expftx
!<--- Parallelization_E

!---> Memmory Reduction_S
!---> Parallelization_S
!        stsrxx=stsrxx+fttx*dxsav(i,j)
!        stsrxy=stsrxy+fttx*dysav(i,j)
!        stsrxz=stsrxz+fttx*dzsav(i,j)
!        stsryy=stsryy+ftty*dysav(i,j)
!        stsryz=stsryz+ftty*dzsav(i,j)
!        stsrzz=stsrzz+fttz*dzsav(i,j)
!        stsrxx=stsrxx+fttx*dxsav(numx)
!        stsrxy=stsrxy+fttx*dysav(numx)
!        stsrxz=stsrxz+fttx*dzsav(numx)
!        stsryy=stsryy+ftty*dysav(numx)
!        stsryz=stsryz+ftty*dzsav(numx)
!        stsrzz=stsrzz+fttz*dzsav(numx)
         sctmp(1)=sctmp(1)+fttx*dxsav(numx)
         sctmp(2)=sctmp(2)+fttx*dysav(numx)
         sctmp(3)=sctmp(3)+fttx*dzsav(numx)
         sctmp(4)=sctmp(4)+ftty*dysav(numx)
         sctmp(5)=sctmp(5)+ftty*dzsav(numx)
         sctmp(6)=sctmp(6)+fttz*dzsav(numx)
!<--- Parallelization_E
!<--- Memmory Reduction_E

!---> Parallelization_S
!        frrx(j)=frrx(j)-fttx
!        frry(j)=frry(j)-ftty
!        frrz(j)=frrz(j)-fttz
!        frrx(i)=frrx(i)+fttx
!        frry(i)=frry(i)+ftty
!        frrz(i)=frrz(i)+fttz
            eltmp(j)=eltmp(j)-fttx
            eltmp(num+j)=eltmp(num+j)-ftty
            eltmp(2*num+j)=eltmp(2*num+j)-fttz
            eltmp(i)=eltmp(i)+fttx
            eltmp(num+i)=eltmp(num+i)+ftty
            eltmp(2*num+i)=eltmp(2*num+i)+fttz
!<--- Parallelization_E
      enddo   
   enddo   
endif
! cation-cation EPP
if(.not.epplog)then
!---> Parallelization_S
!do j=nanion+2,num
do j=jst4,jed4
   jpoint=ntype(j)
!  do i=1+nanion,j-1
   do i=ist4(j),ied4(j)
!<--- Parallelization_E
      ipoint=ntype(i)

!---> Memmory Reduction_S
      numx = numadr(i,j)
!     drsq=dxsav(i,j)*dxsav(i,j)+dysav(i,j)*dysav(i,j)+dzsav(i,j)*dzsav(i,j)
      drsq=dxsav(numx)*dxsav(numx)+dysav(numx)*dysav(numx)+dzsav(numx)*dzsav(numx)
!<--- Memmory Reduction_E

      if(drsq.ge.rsqmaxsr) CYCLE

      drsqrec=1.0d0/drsq
      dr=dsqrt(drsq)
      drrec=1.0d0/dr

      expft=ftb(ipoint,jpoint)*exp(-ftalp(ipoint,jpoint)*dr)
            
      if(ftbx(ipoint,jpoint).gt.0.0d0) then 
         expftx=ftbx(ipoint,jpoint)*exp(-ftalpx(ipoint,jpoint)*drsq)
      else
         expftx=0.d0
      endif
      
      sgam=(expft*ftalp(ipoint,jpoint)*drrec) &
           +expftx*2.0d0*ftalpx(ipoint,jpoint)

!---> Memmory Reduction_S
!     fttx=sgam*dxsav(i,j)
!     ftty=sgam*dysav(i,j)
!     fttz=sgam*dzsav(i,j)
      fttx=sgam*dxsav(numx)
      ftty=sgam*dysav(numx)
      fttz=sgam*dzsav(numx)
!<--- Memmory Reduction_E

!---> Parallelization_S
!     engsr=engsr+expft
      sctmp(59)=sctmp(59)+expft+expftx
!<--- Parallelization_E

!---> Memmory Reduction_S
!---> Parallelization_S
!     stsrxx=stsrxx+fttx*dxsav(i,j)
!     stsrxy=stsrxy+fttx*dysav(i,j)
!     stsrxz=stsrxz+fttx*dzsav(i,j)
!     stsryy=stsryy+ftty*dysav(i,j)
!     stsryz=stsryz+ftty*dzsav(i,j)
!     stsrzz=stsrzz+fttz*dzsav(i,j)
!        stsrxx=stsrxx+fttx*dxsav(numx)
!        stsrxy=stsrxy+fttx*dysav(numx)
!        stsrxz=stsrxz+fttx*dzsav(numx)
!        stsryy=stsryy+ftty*dysav(numx)
!        stsryz=stsryz+ftty*dzsav(numx)
!        stsrzz=stsrzz+fttz*dzsav(numx)
      sctmp(1)=sctmp(1)+fttx*dxsav(numx)
      sctmp(2)=sctmp(2)+fttx*dysav(numx)
      sctmp(3)=sctmp(3)+fttx*dzsav(numx)
      sctmp(4)=sctmp(4)+ftty*dysav(numx)
      sctmp(5)=sctmp(5)+ftty*dzsav(numx)
      sctmp(6)=sctmp(6)+fttz*dzsav(numx)
!<--- Parallelization_E
!<--- Memmory Reduction_E

!---> Parallelization_S
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
!<--- Parallelization_E
   enddo   
enddo   
endif

return
END SUBROUTINE
