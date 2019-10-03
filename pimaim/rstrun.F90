SUBROUTINE rstrun

USE commondata, ONLY: num,x,y,z,vx,vy,vz,rescalelog,ntotstp,nsofar,nstpbrk, &
                      nstep,nrun,pzeta3,bzeta3,vpzeta3,vbzeta3,vpzeta2, &
                      vbzeta2,eps2,eps3,veps2,veps3,vol3,pi,xmu,ymu,zmu, &
                      quadxx,quadyy,quadzz,quadxy,quadxz,quadyz,delta, &
                      epsilonx,epsilony,epsilonz, &
                      quaimxx,quaimyy,quaimzz,quaimxy,quaimxz,quaimyz, &
                      quadpimlog,cimlog,daimlog,quaimlog,rimlog, &
                      conjgradlog,rsqmax,rcut,exit_code
USE boxdata, ONLY: vg3,h,boxlenx,boxleny,boxlenz,h3,hlab2,hlab3,halfboxminsq, &
                   b,h3,hi3,fullh
!---> Parallelization_S
use mpipara
use mpi_spawn
!<--- Parallelization_E

IMPLICIT NONE

INTEGER :: i,j,k
LOGICAL :: crdlog,vellog,chglog,fullrunlog

open(35,file='restart.dat',status='old')
   read(35,*)crdlog
   read(35,*)vellog
   read(35,*)chglog
   read(35,*)fullrunlog

   if(crdlog) then
      do i=1,num
         read(35,*)x(i),y(i),z(i)
      enddo   
   else
!---> Parallelization_S
      if( iam .eq. 0 ) then

      write(6,*)'No coordinate information in restart file!'

      endif
      exit_code = 1
      call close_down_mpi() 
!<--- Parallelization_E
      stop
   endif

   if(vellog) then
      do i=1,num
         read(35,*)vx(i),vy(i),vz(i)
      enddo   
   else
      if(rescalelog) call velset
   endif

   if(chglog) then
      do i=1,num
         read(35,*)xmu(i),ymu(i),zmu(i)
      enddo   
      if(quadpimlog) then
         do i=1,num
            read(35,*)quadxx(i),quadyy(i),quadzz(i)
            read(35,*)quadxy(i),quadxz(i),quadyz(i)
         enddo   
      endif
      if(cimlog.or.daimlog.or.quaimlog) then
         do i=1,num
            read(35,*)delta(i)
         enddo   
         if(daimlog.or.quaimlog) then
            do i=1,num
               read(35,*)epsilonx(i),epsilony(i),epsilonz(i)
            enddo    
            if(quaimlog) then
               do i=1,num
                  read(35,*)quaimxx(i),quaimyy(i),quaimzz(i)
                  read(35,*)quaimxy(i),quaimxz(i),quaimyz(i)
               enddo   
            endif
         endif
      endif
   else
      if(.not. conjgradlog .and. .not. rimlog) then
!---> Parallelization_S
         if( iam .eq. 0 ) then

         write(6,*)'Warning: No dipole variables read in,'
         write(6,*)'         and no annealing is scheduled.'

         endif
!<--- Parallelization_E
      endif
   endif

   if(fullrunlog) then
      ntotstp=0
      read(35,*)nsofar
      do i=1,nsofar
         read(35,*)nstpbrk(i)
         ntotstp=ntotstp+nstpbrk(i)
      enddo   
      nstep=nstep+ntotstp
      nrun=nrun+ntotstp

      do i=1,5
         read(35,*)vpzeta2(i)
         read(35,*)vpzeta3(i)
         read(35,*)pzeta3(i)
         read(35,*)vbzeta2(i)
         read(35,*)vbzeta3(i)
         read(35,*)bzeta3(i)
      enddo   
      read(35,*)eps2,eps3
      read(35,*)veps2,veps3
      read(35,*)vg3(1,1),vg3(2,1),vg3(3,1)
      read(35,*)vg3(1,2),vg3(2,2),vg3(3,2)
      read(35,*)vg3(1,3),vg3(2,3),vg3(3,3)
   endif

   read(35,*) h(1,1),h(1,2),h(1,3)
   read(35,*) h(2,1),h(2,2),h(2,3)
   read(35,*) h(3,1),h(3,2),h(3,3)
   read(35,*)boxlenx
   read(35,*)boxleny
   read(35,*)boxlenz

   call dcell(h,b)
   vol3=boxlenx*boxleny*boxlenz*b(10)
   hlab3=h
   hlab2=h
   h3(:,1)=h(:,1)*boxlenx/(vol3**(1.0d0/3.0d0))
   h3(:,2)=h(:,2)*boxleny/(vol3**(1.0d0/3.0d0))
   h3(:,3)=h(:,3)*boxlenz/(vol3**(1.0d0/3.0d0))
   call invert(h3,hi3)
   call boxreset
   call dcell(fullh,b)
   halfboxminsq=(0.5d0*dmin1(b(7),b(8),b(9)))**2
   rsqmax=rcut**2.0d0

   call kset

close(35)

!---> Parallelization_S
if( iam .eq. 0 ) then
#ifndef ppfit
open(50,file='startuppos.out',status='new')
   do i=1,num
      write(50,*) h(1,1)*x(i)+h(1,2)*y(i)+h(1,3)*z(i), &
                  h(2,1)*x(i)+h(2,2)*y(i)+h(2,3)*z(i), &
                  h(3,1)*x(i)+h(3,2)*y(i)+h(3,3)*z(i)
   enddo
close(50)
#endif

write(6,*)
write(6,*)'**** Run restarted ****'
write(6,*)

endif
!<--- Parallelization_E

return
END SUBROUTINE
