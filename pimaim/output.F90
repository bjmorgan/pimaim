SUBROUTINE output

USE commondata, ONLY: pereng,dippimlog,num,reseng,quadpimlog,dipsqeng, &
                      dipquadeng,quadeng,engpetot,tranke,tcell,tvol,tpzeta, &
                      tbzeta,PeeVee,pzeta,bzeta,tranke,cimlog,selfengtot, &
                      daimlog,epsselfengtot,quaimlog,quaimselfengtot,Pcomx, &
                      Pcomy,Pcomz,stsrxx,stsryy,stsrzz,stsrxy,stsrxz,stsryz, &
                      stcxx,stcyy,stczz,stcxy,stcxz,stcyz,stpxx,stpyy,stpzz, &
                      stpxy,stpxz,stpyz,stp2xx,stp2yy,stp2zz,stp2xy,stp2xz, &
                      stp2yz,stpsrxx,stpsryy,stpsrzz,stpsrxy,stpsrxz,stpsryz, &
                      stqsrxx,stqsryy,stqsrzz,stqsrxy,stqsrxz,stqsryz, &
                      stpqquadxx,stpqquadyy,stpqquadzz,stpqquadxy,stpqquadxz, &
                      stpqquadyz,stpdipquadxx,stpdipquadyy,stpdipquadzz, &
                      stpdipquadxy,stpdipquadxz,stpdipquadyz,stpquadquadxx, &
                      stpquadquadyy,stpquadquadzz,stpquadquadxy,stpquadquadxz, &
                      stpquadquadyz,stewzz,percell,pint2,vol3,nstep,x,y,z,vx, &
                      vy,vz,tke,nummon,elecx,elecy,elecz,exx,eyy,ezz,exy,exz, &
                      eyz,xmu,ymu,zmu,quadxx,quadyy,quadzz,quadxy,quadxz, &
                      quadyz,srdipx,srdipy,srdipz,alppolar,asdipx,asdipy, &
                      asdipz,Bpolar,Cpolar,gammapolar,srquadxx,srquadyy, &
                      srquadzz,srquadxy,srquadxz,srquadyz,asquadxx,asquadyy, &
                      asquadzz,asquadxy,asquadxz,asquadyz,delta,epsilonx, &
                      epsilony,epsilonz,quaimxx,quaimyy,quaimzz,quaimxy, &
                      quaimxz,quaimyz,selfeps,selfquaim,pervel,perfric,nmon, &
                      engeff,dtime, amass, chg, atmnam, weight,chge,frrx,frry,&
		      frrz,nunitcellx,nunitcelly,nunitcellz,istrj,jstrj,keytrj,&
                      jobname 
USE boxdata, ONLY: cellvol3rec,cellvol,boxlenx,boxleny,boxlenz,bee,hlab3,h,a0,b0,c0

IMPLICIT NONE

INTEGER :: i,nchan,nnchan,nnn,nt,nnnchan,nnnnchan,nnnn
DOUBLE PRECISION :: polengtot

    polengtot=0.0 
if(pereng) then
   pereng=.false.

   if(dippimlog) polengtot=SUM(reseng)
   if(quadpimlog) polengtot=SUM(reseng+dipsqeng+dipquadeng+quadeng)
! energies
!   write(21,*)nstep,real(engpetot),real(tranke)
   write(21,*)nstep,engpetot,tranke
   write(32,*)nstep,real(tcell),real(tvol)
   write(34,*)nstep,real(tpzeta),real(pzeta)
   write(44,*)nstep,real(tbzeta),real(bzeta)
   write(33,*)nstep,real(PeeVee),(engpetot+tranke+polengtot &
             +selfengtot+epsselfengtot+quaimselfengtot &
             +tcell+tvol+PeeVee+pzeta+tpzeta+tbzeta+bzeta)

   if(cimlog) then
      write(77,*)nstep,real(selfengtot)
      write(25,*)nstep,real(polengtot) &
                      ,real(engpetot+tranke+polengtot+selfengtot)
   else if(daimlog) then
      write(77,*)nstep,real(selfengtot)
      write(78,*)nstep,real(epsselfengtot)
      write(25,*)nstep,real(polengtot) &
                      ,real(engpetot+tranke+polengtot+selfengtot+epsselfengtot)
   else if(quaimlog) then
      write(77,*)nstep,real(selfengtot)
      write(78,*)nstep,real(epsselfengtot)
      write(79,*)nstep,real(quaimselfengtot)
      write(25,*)nstep,real(polengtot) &
                      ,real((engpetot+tranke+polengtot+selfengtot &
                           +epsselfengtot+quaimselfengtot))
   else
      write(25,*)nstep,real(polengtot),real(engpetot+tranke+polengtot)
   endif
! cell momentum
   write(36,*)nstep,real(Pcomx),real(Pcomy),real(Pcomz)

! stress tensor
!
! total term:
   write (56,*) nstep,(stpxx+stcxx+stsrxx+stpyy+stcyy+ &
                 stsryy+stpzz+stczz+stsrzz+stp2xx+ &
                 stp2yy+stp2zz+stpsrxx+stpsryy+stpsrzz+ &
                 stpqquadxx+stpqquadyy+stpqquadzz+ &
                 stpdipquadxx+stpdipquadyy+stpdipquadzz+ &
                 stpquadquadxx+stpquadquadyy+stpquadquadzz+ &
                 stqsrxx+stqsryy+stqsrzz)* &
                 cellvol3rec+stewzz/3.0d0
!
! components: (xx, yy, zz)
!
   write (57,*) nstep,real((stpxx+stcxx+stsrxx+ &
                 stp2xx+stpsrxx+stpqquadxx+stpdipquadxx+ &
                 stpquadquadxx+stqsrxx)/cellvol),real(( &
                 stpyy+stcyy+stsryy+stp2yy+stpsryy+ &
                 stpqquadyy+stpdipquadyy+stpquadquadyy+ &
                 stqsryy)/cellvol),real((stpzz+stczz+stsrzz+ &
                 stp2zz+stpsrzz+stpqquadzz+stpdipquadzz+ &
                 stpquadquadzz+stqsrzz)/cellvol+stewzz)
!
! components: (xy, xz, yz)
!
   write (58,*) nstep,real((stpxy+stcxy+stsrxy+ &
                 stp2xy+stpsrxy+stpqquadxy+stpdipquadxy+ &
                 stpquadquadxy+stqsrxy)/cellvol) &
                ,real((stpxz+stcxz+stsrxz+stp2xz+stpsrxz+ &
                 stpqquadxz+stpdipquadxz+stpquadquadxz+ &
                 stqsrxz)/cellvol),real((stpyz+stcyz+stsryz+ &
                 stp2yz+stpsryz+stpqquadyz+stpdipquadyz+ &
                 stpquadquadyz+stqsryz)/cellvol)
!
! charge-dipole + s-r dipole + dipole-dipole term:
!
   write (62,*) nstep,real((stpxx+stpyy+stpzz)*cellvol3rec) &
                ,real((stpsrxx+stpsryy+stpsrzz)*cellvol3rec) &
                ,real((stp2xx+stp2yy+stp2zz)*cellvol3rec)
!
! charge-charge + short-range term + VCB term:
!
   write (63,*) nstep,real((stcxx+stcyy+stczz)*cellvol3rec) &
               ,real((stsrxx+stsryy+stsrzz)*cellvol3rec) &
               ,real(stewzz)
!
! charge-quadrupole + dipole-quadrupole + quadrupole-quadrupole
!
   write (64,*) nstep,real((stpqquadxx+stpqquadyy+ &
                            stpqquadzz+stqsrxx+ &
                            stqsryy+stqsrzz)*cellvol3rec) &
               ,real((stpdipquadxx+stpdipquadyy+ &
                      stpdipquadzz)*cellvol3rec) &
               ,real((stpquadquadxx+stpquadquadyy+ &
                     stpquadquadxx)*cellvol3rec)
!
! Different contributions to the surface tension.
!
   write (66,*) nstep,real((stczz-(stcxx+stcyy)/2.0d0)/cellvol) &
               ,real((stsrzz-(stsrxx+stsryy)/2.0d0)/cellvol)

   write (67,*) nstep,real((stpzz-(stpxx+stpyy)/2.0d0)/cellvol) &
               ,real((stpsrzz-(stpsrxx+stpsryy)/2.0d0)/cellvol) &
               ,real((stp2zz-(stp2xx+stp2yy)/2.0d0)/cellvol)

   write (68,*) nstep,real((stpqquadzz-(stpqquadxx+stpqquadyy)/2.0d0)/cellvol) &
         ,real((stpdipquadzz-(stpdipquadxx+stpdipquadyy)/2.0d0)/cellvol) &
         ,real((stpquadquadzz-(stpquadquadxx+stpquadquadyy)/2.0d0)/cellvol) &
         ,real((stqsrzz-(stqsrxx+stqsryy)/2.0d0)/cellvol)
endif
! from trans_chains
if(percell) then
   percell=.false.
   write(51,*) nstep,real(pint2(1,1)),real(pint2(2,2)),real(pint2(3,3))
   write(52,*) nstep,((pint2(1,1)+pint2(2,2)+pint2(3,3))/3.0d0)
   write(55,*) nstep,real(bee(4)),real(bee(5)),real(bee(6))
   write(53,'(I12,1x,3(F30.16,1x))') nstep,boxlenx,boxleny,boxlenz
   write(54,*) nstep,vol3,(float(num)/vol3)
endif

if(pervel) then
   pervel=.false.
   write(61,*)hlab3(1,1),hlab3(1,2),hlab3(1,3)
   write(61,*)hlab3(2,1),hlab3(2,2),hlab3(2,3)
   write(61,*)hlab3(3,1),hlab3(3,2),hlab3(3,3)
   write(61,*)boxlenx,boxleny,boxlenz


   do i=1,num
      write(22,*)vx(i),vy(i),vz(i)
      write(23,*)x(i),y(i),z(i)
      write(60,*) hlab3(1,1)*x(i)+hlab3(1,2)*y(i)+hlab3(1,3)*z(i), &
                  hlab3(2,1)*x(i)+hlab3(2,2)*y(i)+hlab3(2,3)*z(i), &
                  hlab3(3,1)*x(i)+hlab3(3,2)*y(i)+hlab3(3,3)*z(i)
   enddo  
      if (dippimlog) then
         do i=1, num
            write(24,*) xmu(i), ymu(i), zmu(i)
         enddo
      endif

endif

if(perfric) then
   perfric=.false.
!sgi   write(30,*)nstep,real(tke)
   write(30,*)nstep,tke

   if(dippimlog) then
      nchan=81
      nnchan=84
      do i=1,nummon
         nnn=0
         write(nchan+nnn,*)nstep,real(xmu(nmon(i))) &
          ,real(elecx(nmon(i))*alppolar(i)) &
          ,real(elecx(nmon(i)))
         write(nchan+1+nnn,*)nstep,real(ymu(nmon(i))) &
          ,real(elecy(nmon(i))*alppolar(i)) &
          ,real(elecy(nmon(i)))
         write(nchan+2+nnn,*)nstep,real(zmu(nmon(i))) &
          ,real(elecz(nmon(i))*alppolar(i)) &
          ,real(elecz(nmon(i)))
         write(nnchan+nnn,*)nstep,real(asdipx(nmon(i))) &
          ,real(srdipx(nmon(i)))
         write(nnchan+1+nnn,*)nstep,real(asdipy(nmon(i))) &
          ,real(srdipy(nmon(i)))
         write(nnchan+2+nnn,*)nstep,real(asdipz(nmon(i))) &
          ,real(srdipz(nmon(i)))
      enddo   
   endif
   if(quadpimlog) then
      nchan=81
      nnchan=71
      nnnchan=84
      nnnnchan=91
      do i=1,nummon
         nnn=0
         nnnn=0
         write(nchan+nnn,*)nstep,real(xmu(nmon(i))) &
          ,real(elecx(nmon(i))*alppolar(i)) &
          ,real(elecx(nmon(i))*alppolar(i) &
       +(1.0d0/3.0d0)*Bpolar(i)*(elecx(nmon(i))*exx(nmon(i)) &
       -(0.5d0*elecx(nmon(i))*eyy(nmon(i))) &
       -(0.5d0*elecx(nmon(i))*ezz(nmon(i))) &
       +(1.5d0*elecy(nmon(i))*exy(nmon(i))) &
       +(1.5d0*elecz(nmon(i))*exz(nmon(i)))) &
       +(1.0d0/6.0d0)*gammapolar(i)* &
               (elecx(nmon(i))*elecx(nmon(i))*elecx(nmon(i)) &
               +elecx(nmon(i))*elecy(nmon(i))*elecy(nmon(i)) &
               +elecx(nmon(i))*elecz(nmon(i))*elecz(nmon(i))))
         write(nchan+1+nnn,*)nstep,real(ymu(nmon(i))) &
          ,real(elecy(nmon(i))*alppolar(i)) &
        ,real(elecy(nmon(i))*alppolar(i) &
       +(1.0d0/3.0d0)*Bpolar(i)*(elecy(nmon(i))*eyy(nmon(i)) &
       -(0.5d0*elecy(nmon(i))*exx(nmon(i))) &
       -(0.5d0*elecy(nmon(i))*ezz(nmon(i))) &
       +(1.5d0*elecx(nmon(i))*exy(nmon(i))) &
       +(1.5d0*elecz(nmon(i))*eyz(nmon(i)))) &
       +(1.0d0/6.0d0)*gammapolar(i)* &
               (elecy(nmon(i))*elecx(nmon(i))*elecx(nmon(i)) &
               +elecy(nmon(i))*elecy(nmon(i))*elecy(nmon(i)) &
               +elecy(nmon(i))*elecz(nmon(i))*elecz(nmon(i))))
         write(nchan+2+nnn,*)nstep,real(zmu(nmon(i))) &
          ,real(elecz(nmon(i))*alppolar(i)) &
        ,real(elecz(nmon(i))*alppolar(i) &
       +(1.0d0/3.0d0)*Bpolar(i)*(elecz(nmon(i))*ezz(nmon(i)) &
       -(0.5d0*elecz(nmon(i))*exx(nmon(i))) &
       -(0.5d0*elecz(nmon(i))*eyy(nmon(i))) &
       +(1.5d0*elecy(nmon(i))*eyz(nmon(i))) &
       +(1.5d0*elecx(nmon(i))*exz(nmon(i)))) &
       +(1.0d0/6.0d0)*gammapolar(i)* &
               (elecz(nmon(i))*elecx(nmon(i))*elecx(nmon(i)) &
               +elecz(nmon(i))*elecy(nmon(i))*elecy(nmon(i)) &
               +elecz(nmon(i))*elecz(nmon(i))*elecz(nmon(i))))

         write(nnchan+nnnn,*)nstep,real(quadxx(nmon(i))) &
        ,real((2.0d0/3.0d0)*Cpolar(i)*exx(nmon(i)) &
        -(1.0d0/3.0d0)*Cpolar(i)*eyy(nmon(i)) &
        -(1.0d0/3.0d0)*Cpolar(i)*ezz(nmon(i))) &
        ,real((2.0d0/3.0d0)*Cpolar(i)*exx(nmon(i)) &
        -(1.0d0/3.0d0)*Cpolar(i)*eyy(nmon(i)) &
        -(1.0d0/3.0d0)*Cpolar(i)*ezz(nmon(i)) &
        +0.5d0*Bpolar(i)*elecx(nmon(i))*elecx(nmon(i)) &
        -0.25d0*Bpolar(i)*elecy(nmon(i))*elecy(nmon(i)) &
        -0.25d0*Bpolar(i)*elecz(nmon(i))*elecz(nmon(i)))
         write(nnchan+1+nnnn,*)nstep,real(quadyy(nmon(i))) &
        ,real((2.0d0/3.0d0)*Cpolar(i)*eyy(nmon(i)) &
        -(1.0d0/3.0d0)*Cpolar(i)*exx(nmon(i)) &
        -(1.0d0/3.0d0)*Cpolar(i)*ezz(nmon(i))) &
        ,real((2.0d0/3.0d0)*Cpolar(i)*eyy(nmon(i)) &
        -(1.0d0/3.0d0)*Cpolar(i)*exx(nmon(i)) &
        -(1.0d0/3.0d0)*Cpolar(i)*ezz(nmon(i)) &
        +0.5d0*Bpolar(i)*elecy(nmon(i))*elecy(nmon(i)) &
        -0.25d0*Bpolar(i)*elecx(nmon(i))*elecx(nmon(i)) &
        -0.25d0*Bpolar(i)*elecz(nmon(i))*elecz(nmon(i)))
         write(nnchan+2+nnnn,*)nstep,real(quadzz(nmon(i))) &
        ,real((2.0d0/3.0d0)*Cpolar(i)*ezz(nmon(i)) &
        -(1.0d0/3.0d0)*Cpolar(i)*eyy(nmon(i)) &
        -(1.0d0/3.0d0)*Cpolar(i)*exx(nmon(i))) &
        ,real((2.0d0/3.0d0)*Cpolar(i)*ezz(nmon(i)) &
        -(1.0d0/3.0d0)*Cpolar(i)*eyy(nmon(i)) &
        -(1.0d0/3.0d0)*Cpolar(i)*exx(nmon(i)) &
        +0.5d0*Bpolar(i)*elecz(nmon(i))*elecz(nmon(i)) &
        -0.25d0*Bpolar(i)*elecy(nmon(i))*elecy(nmon(i)) &
        -0.25d0*Bpolar(i)*elecx(nmon(i))*elecx(nmon(i)))
         write(nnchan+3+nnnn,*)nstep,real(quadxy(nmon(i))) &
        ,real(Cpolar(i)*exy(nmon(i))) &
        ,real((0.75d0*Bpolar(i)*elecx(nmon(i))*elecy(nmon(i))) &
        +(Cpolar(i)*exy(nmon(i))))
         write(nnchan+4+nnnn,*)nstep,real(quadxz(nmon(i))) &
        ,real(Cpolar(i)*exz(nmon(i))) &
        ,real((0.75d0*Bpolar(i)*elecx(nmon(i))*elecz(nmon(i))) &
        +(Cpolar(i)*exz(nmon(i))))
         write(nnchan+5+nnnn,*)nstep,real(quadyz(nmon(i))) &
        ,real(Cpolar(i)*eyz(nmon(i))) &
        ,real((0.75d0*Bpolar(i)*elecy(nmon(i))*elecz(nmon(i))) &
        +(Cpolar(i)*eyz(nmon(i))))

         write(nnnchan+nnn,*)nstep,real(asdipx(nmon(i))) &
          ,real(srdipx(nmon(i)))
         write(nnnchan+1+nnn,*)nstep,real(asdipy(nmon(i))) &
          ,real(srdipy(nmon(i)))
         write(nnnchan+2+nnn,*)nstep,real(asdipz(nmon(i))) &
          ,real(srdipz(nmon(i)))

         write(nnnnchan+nnnn,*)nstep,real(asquadxx(nmon(i))) &
          ,real(srquadxx(nmon(i)))
         write(nnnnchan+1+nnnn,*)nstep,real(asquadyy(nmon(i))) &
          ,real(srquadyy(nmon(i)))
         write(nnnnchan+2+nnnn,*)nstep,real(asquadzz(nmon(i))) &
          ,real(srquadzz(nmon(i)))
         write(nnnnchan+3+nnnn,*)nstep,real(asquadxy(nmon(i))) &
          ,real(srquadxy(nmon(i)))
         write(nnnnchan+4+nnnn,*)nstep,real(asquadxz(nmon(i))) &
          ,real(srquadxz(nmon(i)))
         write(nnnnchan+5+nnnn,*)nstep,real(asquadyz(nmon(i))) &
          ,real(srquadyz(nmon(i)))
      enddo   
   endif
   if((cimlog).or.(daimlog).or.(quaimlog)) then
      nchan=87
      do i=1,nummon
         nnn=0
         write(nchan+nnn,*)nstep,real(delta(nmon(i))) &
                                ,real(engeff(nmon(i))) &
                                ,real(alppolar(nmon(i)))
      enddo   
   endif
   if((daimlog).or.(quaimlog)) then
      nchan=88
      do i=1,nummon
         nnn=0
         write(nchan+nnn,*)nstep,real(epsilonx(nmon(i))) &
           ,real(epsilony(nmon(i))),real(epsilonz(nmon(i))) &
           ,real(selfeps(nmon(i)))
      enddo   
   endif 
   if(quaimlog) then
      nchan=89
      do i=1,nummon
         nnn=0
         write(nchan+nnn,*)nstep,real(quaimxx(nmon(i))) &
                                ,real(quaimyy(nmon(i))) &
                                ,real(quaimzz(nmon(i))) &
                                ,real(selfquaim(nmon(i)))
         write(nchan+nnn+1,*)nstep,real(quaimxy(nmon(i))) &
                                ,real(quaimxz(nmon(i))) &
                                ,real(quaimyz(nmon(i)))
      enddo   
   endif
endif

return
END SUBROUTINE
