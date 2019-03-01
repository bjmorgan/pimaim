SUBROUTINE output_MABC

USE commondata
USE boxdata

IMPLICIT NONE

integer :: i
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
open(666,file='HISTORY',position='append')

! File fulloutput.txt contains information such as runtype, energies, pressure,etc
open(667,file='fulloutput.dat',position='append')
if(nstep.eq.istrj) then
write(667,'(a80)') jobname
write(667,'(3a10)') 'Runtype= ', runtype,runtype2
write(667,'(a10,5a5)') 'Species: ', (spectype(i),i=1,nspec)
write(667,'(2(a15,l3,2X))') 'Isotropic: ',nib,'Anisotropic: ',nab
write(667,1000) 'Step No.','PE  ','KE','Tot. Eng.','Vol.(A^3)','Press.(GPa)', 'Temp.(K)'
endif
1000 format (10(a,10X))

if((nstep.ge.istrj).and.(mod(nstep-istrj-1,jstrj).eq.0))then
  write(666,'(a8,4i10,a2,g12.6)') 'timestep',nstep,num,keytrj,3,'    ', dtime*au_to_ps !!MABC: write to HISTORY
  write(666,'(3(g12.6,2X))') au_to_ang*boxlenx*h(1,1),au_to_ang*boxleny*h(2,1),au_to_ang*boxlenz*h(3,1)  !MABC: write unit cell to HISTORY
  write(666,'(3(g12.6,2X))') au_to_ang*boxlenx*h(1,2),au_to_ang*boxleny*h(2,2),au_to_ang*boxlenz*h(3,2)  !MABC: write unit cell to HISTORY
  write(666,'(3(g12.6,2X))') au_to_ang*boxlenx*h(1,3),au_to_ang*boxleny*h(2,3),au_to_ang*boxlenz*h(3,3)  !MABC: write unit cell to HISTORY
  write(667,'(i8,10(2X,g16.10))') nstep,engpetot,tranke,(engpetot+tranke),(vol3*au_to_ang3),&
                                 ((pint2(1,1)+pint2(2,2)+pint2(3,3))/3.0d0)*2.9421912E04,tke

          do i = 1,num

           write(666,'(a8,i10,2f12.6)') atmnam(i),i,weight(i),chge(i)   !MABC: write atm label, number, mass, chg
           write(666,'(1p,3e12.4)')au_to_ang*x(i),au_to_ang*y(i),au_to_ang*z(i)    !MABC: write positions in cell frame to HISTORY
            if(keytrj.ge.1)then
             write(666,'(1p,3e12.4)')vx(i)*conv_v,vy(i)*conv_v,vz(i)*conv_v !MABC: write to HISTORY
            endif
            if(keytrj.ge.2)then
             write(666,'(1p,3e12.4)')frrx(i)*au_to_force,frry(i)*au_to_force,frrz(i)*au_to_force !MABC: write forces to HISTORY
            endif

          enddo
endif

close(666)

END SUBROUTINE
