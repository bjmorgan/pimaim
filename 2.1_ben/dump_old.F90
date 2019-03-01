!************************************************************
! Dumps ion positions, velocities and information on
! the extra degrees of freedom to disk, for later use
! as a restart file.
!************************************************************

SUBROUTINE dump

use commondata
use boxdata
!USE commondata, ONLY: fileout,crddumplog,veldumplog,chgdumplog,fulldumplog, &
!                      num,x,y,z,vx,vy,vz,xmu,ymu,zmu,quadxx,quadyy,quadzz, &
!                      quadxy,quadxz,quadyz,delta,epsilonx,epsilony,epsilonz, &
!                      quaimxx,quaimyy,quaimzz,quaimxy,quaimxz,quaimyz,nsofar, &
!                      nstpbrk,nstep,ntotstp,vpzeta2,vbzeta2,vpzeta3,vbzeta3, &
!                      pzeta3,bzeta3,eps2,eps3,veps2,veps3,cimlog,daimlog, &
!                      quaimlog,quadpimlog
!USE boxdata, ONLY: vg3,boxlenx,boxleny,boxlenz,hlab3

IMPLICIT NONE

INTEGER :: i
double precision :: au_to_ang,au_to_ang3,au_to_ps,conv_v,au_to_force
!! The constants below are used to convert from A.U. to A and ps, which are used in 
!! DL_POLY
au_to_ang= 0.5291772108
au_to_ps= 2.418884325505e-05
conv_v=au_to_ang/au_to_ps
au_to_ang3=au_to_ang*au_to_ang*au_to_ang 
au_to_force=9.1093821545e-31/(1000*1.6605402e-27)    !!MABC: Not sure about this conversion factor
!!au_to_force=1.6605402e-27/(9.1093821545e-31)


! File HISTORY resembles that from DL_POLY and allows data manipulation
open (527,file='CONFIG',status='unknown')
rewind (527)
open (528,file='POSCAR',status='unknown')
rewind (528)

 write(527,'(a150)')jobname
 write(528,'(a150)')jobname
 write(528,'(f12.6)') 1.000
 write(527,'(3i10)') keytrj,3,num

  write(527,1001) au_to_ang*boxlenx*h(1,1),au_to_ang*boxleny*h(2,1),au_to_ang*boxlenz*h(3,1)  !MABC: write unit cell to CONFIG 
  write(527,1001) au_to_ang*boxlenx*h(1,2),au_to_ang*boxleny*h(2,2),au_to_ang*boxlenz*h(3,2)  !MABC: write unit cell to CONFIG 
  write(527,1001) au_to_ang*boxlenx*h(1,3),au_to_ang*boxleny*h(2,3),au_to_ang*boxlenz*h(3,3)  !MABC: write unit cell to CONFIG 

  write(528,1001) au_to_ang*boxlenx*h(1,1),au_to_ang*boxleny*h(2,1),au_to_ang*boxlenz*h(3,1)  !MABC: write unit cell to POSCAR 
  write(528,1001) au_to_ang*boxlenx*h(1,2),au_to_ang*boxleny*h(2,2),au_to_ang*boxlenz*h(3,2)  !MABC: write unit cell to POSCAR 
  write(528,1001) au_to_ang*boxlenx*h(1,3),au_to_ang*boxleny*h(2,3),au_to_ang*boxlenz*h(3,3)  !MABC: write unit cell to POSCAR 
  write(528,*) (spectype(i),i=1,nspec)
  write(528,*) (nsp(i),i=1,nspec)
  write(528,*) 'Direct'

          do i = 1,num

	 write(527,'(a8,i10,2f12.6)') atmnam(i),i,weight(i),chge(i)   !MABC: write atm label, number, mass, chg
	 write(527,'(1p,3e12.4)')au_to_ang*x(i),au_to_ang*y(i),au_to_ang*z(i)    !MABC: write positions in cell frame to HISTORY
         write(528,'(1p,3e12.5)')(au_to_ang*x(i))/(au_to_ang*boxlenx*h(1,1)),(au_to_ang*y(i))/(au_to_ang*boxlenx*h(2,2)),   & 
	      &    (au_to_ang*z(i))/(au_to_ang*boxlenx*h(3,3))   
            if(keytrj.ge.1)then
             write(527,'(1p,3e12.4)')vx(i)*conv_v,vy(i)*conv_v,vz(i)*conv_v !MABC: write to HISTORY
            endif
            if(keytrj.ge.2)then
             write(527,'(1p,3e12.4)')frrx(i)*au_to_force,frry(i)*au_to_force,frrz(i)*au_to_force !MABC: write forces to HISTORY
            endif

          enddo

close(527)
close(528)



open(37,file=fileout,status='unknown')

   write(37,*) crddumplog
   write(37,*) veldumplog
   write(37,*) chgdumplog
   write(37,*) fulldumplog

   if(crddumplog) then
      do i=1,num
         write(37,*)x(i),y(i),z(i)
      enddo   
   endif

   if(veldumplog) then
      do i=1,num
         write(37,*)vx(i),vy(i),vz(i)
      enddo   
   endif

   if(chgdumplog) then
      do i=1,num
         write(37,*)xmu(i),ymu(i),zmu(i)
      enddo    

      if (quadpimlog) then
         do i=1,num
            write(37,*)quadxx(i),quadyy(i),quadzz(i)
            write(37,*)quadxy(i),quadxz(i),quadyz(i)
         enddo   
      endif

      if(cimlog.or.daimlog.or.quaimlog) then
         do i=1,num
            write(37,*)delta(i)
         enddo   
      endif

      if(daimlog.or.quaimlog) then
         do i=1,num
            write(37,*)epsilonx(i),epsilony(i),epsilonz(i)
         enddo   
      endif

      if(quaimlog) then
         do i=1,num
            write(37,*)quaimxx(i),quaimyy(i),quaimzz(i)
            write(37,*)quaimxy(i),quaimxz(i),quaimyz(i)
         enddo   
      endif
   endif

   if(fulldumplog) then
      write(37,*)nsofar+1
      do i=1,nsofar
         write(37,*)nstpbrk(i)
      enddo   
      write(37,*)nstep-ntotstp
!write out barostat data
      do i=1,5
         write(37,*)vpzeta2(i)
         write(37,*)vpzeta3(i)
         write(37,*)pzeta3(i)
         write(37,*)vbzeta2(i)
         write(37,*)vbzeta3(i)
         write(37,*)bzeta3(i)
      enddo   
      write(37,*)eps2,eps3
      write(37,*)veps2,veps3

      write(37,*)vg3(1,1),vg3(2,1),vg3(3,1)
      write(37,*)vg3(1,2),vg3(2,2),vg3(3,2)
      write(37,*)vg3(1,3),vg3(2,3),vg3(3,3)
   endif
!write out new cell matrix
   write(37,*) hlab3(1,1),hlab3(1,2),hlab3(1,3)
   write(37,*) hlab3(2,1),hlab3(2,2),hlab3(2,3)
   write(37,*) hlab3(3,1),hlab3(3,2),hlab3(3,3)
   write (37,*) boxlenx
   write (37,*) boxleny
   write (37,*) boxlenz

close(37)

1001 format (g12.6,2X,g12.6,2X,g12.6,2X) !MABC

return
END SUBROUTINE
