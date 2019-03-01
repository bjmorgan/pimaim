!
! Calls correct energy routine.
!

SUBROUTINE ener

USE commondata, ONLY: dippimlog,quadpimlog,cluslog,chargequadlog, &
                      dipquadlog,fullewaldlog,rimlog,num,frrx, &
                      frry,frrz,engpetot,stsrxx,stsryy,stsrzz, &
                      stsrxy,stsrxz,stsryz,stcxx,stcyy,stczz,stcxy,stcxz, &
                      stcyz,stpxx,stpyy,stpzz,stpxy,stpxz,stpyz,stp2xx, &
                      stp2yy,stp2zz,stp2xy,stp2xz,stp2yz,stpqquadxx, &
                      stpqquadyy,stpqquadzz,stpqquadxy,stpqquadxz,stpqquadyz, &
                      stpdipquadxx,stpdipquadyy,stpdipquadzz,stpdipquadxy, &
                      stpdipquadxz,stpdipquadyz,stpquadquadxx,stpquadquadyy, &
                      stpquadquadzz,stpquadquadxy,stpquadquadxz,stpquadquadyz, &
                      stpsrxx,stpsryy,stpsrzz,stpsrxy,stpsrxz,stpsryz, &
                      stqsrxx,stqsryy,stqsrzz,stqsrxy,stqsrxz,stqsryz,nspec, &
! force check     
                      gw_force, & 
! MS modif 31/03
                      dispcorrec1
! MS endmodif 31/03
! MS modif 31/03
USE boxdata, ONLY: cellvol
! MS endmodif 31/03
USE mpipara

IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(nspec,2) :: dimmtmp

DOUBLE PRECISION ::  time0, time1, time2 

Integer :: i 

#ifndef dipole
! GWW - currently need to call seperation here as for RIM - since done in
! multipole relaxation - must change in the future 
  call separations !Have to call this otherwise distances (d?sav) aren't set up!
#endif 

!--------------------------------------------------------------------------------
#ifdef time2
  time20 = mpi_wtime()
#endif 

    call fullrecip

#ifdef time2
  time21 = mpi_wtime()
  time22 = time21 - time20
  if( iam .eq. 0 ) then
     write(*,'(a,F20.8,a)') ' Real     Time : ',time22,' (Sec.)'
  endif
#endif 

!--------------------------------------------------------------------------------

#ifdef time2
  time20 = mpi_wtime()
#endif 

    call fullreal

#ifdef time2
  time21 = mpi_wtime()
  time22 = time21 - time20
  if( iam .eq. 0 ) then
     write(*,'(a,F20.8,a)') ' Recip E  Time : ',time22,' (Sec.)'
  endif
#endif 

!--------------------------------------------------------------------------------
#ifdef time2
  time20 = mpi_wtime()
#endif 
    call sr_energy

#ifdef time2
  time21 = mpi_wtime()
  time22 = time21 - time20
  if( iam .eq. 0 ) then
     write(*,'(a,F20.8,a)') ' Short R  Time : ',time22,' (Sec.)'
  endif
#endif 

!--------------------------------------------------------------------------------

#ifndef dipole 
! GWW strange thing to do for RIM - must check if this is actually required
  dimtmp=0.d0
#endif 


!--------------------------------------------------------------------------------
!if(quadpimlog) then
!   if(cluslog) then
!      call realE_clus_quadpim
!   else
!      if(chargequadlog) then 
!         call recipE_chargequad
!         call realE_chargequad
!      endif
!      if(dipquadlog) then 
!         call recipE_dipquad
!         call realE_dipquad
!      endif
!      if(fullewaldlog) then 
!         call fullrecip
!         call fullreal
!      endif
!   endif
!   call sr_energy
!endif
!
!if(dippimlog) then
!   if(cluslog) then
!      call realE_clus_dippim
!   else
!
!time0 = mpi_wtime()
!!      call recipE_dippim
!         call fullrecip
!time1 = mpi_wtime()
!time2 = time1 - time0
!if( iam .eq. 0 ) then
!   write(*,'(a,F20.8,a)') ' RecipE   Time : ',time2,' (Sec.)'
!endif
!
!time0 = mpi_wtime()
!      call realE_dippim
!time1 = mpi_wtime()
!time2 = time1 - time0
!if( iam .eq. 0 ) then
!   write(*,'(a,F20.8,a)') ' Real     Time : ',time2,' (Sec.)'
!endif
!
!   endif
!time0 = mpi_wtime()
!   call sr_energy
!time1 = mpi_wtime()
!time2 = time1 - time0
!if( iam .eq. 0 ) then
!   write(*,'(a,F20.8,a)') ' Short R  Time : ',time2,' (Sec.)'
!endif
!
!endif
!
!if(rimlog) then
!   call separations !Have to call this otherwise distances (d?sav) aren't set up!
!   if(cluslog) then
!      call realE_clus_rgd
!   else
!      call rgdrecipE
!!      call fullrecip 
!      call rgdrealE
!   endif
!   call sr_energy
!   dimtmp=0.d0
!endif

CALL MPI_ALLREDUCE(eltmp,elltmp,6*num,mpi_double_precision,MPI_SUM, &
                   MPI_COMM_WORLD,ierr)

frrx=elltmp(1:num)+elltmp(3*num+1:4*num)
frry=elltmp(num+1:2*num)+elltmp(4*num+1:5*num)
frrz=elltmp(2*num+1:3*num)+elltmp(5*num+1:6*num)

CALL MPI_ALLREDUCE(sctmp,scctmp,59,mpi_double_precision,MPI_SUM, &
                   MPI_COMM_WORLD,ierr)

CALL MPI_ALLREDUCE(dimtmp,dimmtmp,nspec*2,mpi_double_precision,MPI_SUM, &
                   MPI_COMM_WORLD,ierr)

! GWW - some stuff in engpetot  - need to distribute them to appropriate bits 
scctmp(57) = scctmp(57) + SUM(dimmtmp)

! Total, coulomb(q), VDW,  dipole ,  quadrupole, short range  
! GWW - sum energy
! 55 - charge- charge
! 56 - VDW
! 57 - dipole
! 58 - quadrupole 
! 59 - short range 
engpetot = scctmp(55) + scctmp(56) + scctmp(57) + scctmp(58) + scctmp(59)

!print*,'Total Coulombic E = ', scctmp(55) 

!#ifdef debug
  if( iam .eq. 0 ) then
    write(*,'(a,6(F15.8,X))') ' Energy = : ',engpetot, (scctmp(i),i=55,59) 
  endif
!#endif 

#ifdef debug
  if (iam.eq.0) then
    write(900,*)'# vdw forces' 
    write(901,*)'# total forces' 
    write(902,*)'# charge forces' 
    write(903,*)'# q-mu + mu-q  forces' 
    write(904,*)'# mu-mu forces' 
    write(905,*)'# mu-mu forces' 

    do i=1, num
    write(900,'(i6,3(2X,F20.12))') i, elltmp(3*num+i), elltmp(4*num+i), elltmp(5*num+i) 
    write(901,'(i6,3(2X,F20.12))') i, frrx(i), frry(i), frrz(i) 
    write(902,'(i6,3(2X,F20.12))') i, gw_force(i), gw_force(1*num+i), gw_force(2*num+i) 
    write(903,'(i6,3(2X,F20.12))') i, gw_force(3*num+i), gw_force(4*num+i), gw_force(5*num+i) 
    write(904,'(i6,3(2X,F20.12))') i, gw_force(6*num+i), gw_force(7*num+i), gw_force(8*num+i) 
    write(905,'(i6,3(2X,F20.12))') i, gw_force( 9*num+i), gw_force(10*num+i), gw_force(11*num+i) 
    enddo 

    write(900,'("Sum ",3(2X,F20.12))') sum(elltmp(3*num+1:4*num)), sum(elltmp(4*num+1:5*num)), sum(elltmp(5*num+1:6*num)) 
    write(901,'("Sum ",3(2X,F20.12))') sum(frrx), sum(frry), sum(frrz) 
    write(902,'("Sum ",3(2X,F20.12))') sum(gw_force(1:num)), sum(gw_force(num+1:2*num)), sum(gw_force(2*num+1:3*num)) 
    write(903,'("Sum ",3(2X,F20.12))') sum(gw_force(3*num+1:4*num)), sum(gw_force(4*num+1:5*num)), sum(gw_force(5*num+1:6*num)) 
    write(904,'("Sum ",3(2X,F20.12))') sum(gw_force(6*num+1:7*num)), sum(gw_force(7*num+1:8*num)), sum(gw_force(8*num+1:9*num)) 
    write(905,'("Sum ",3(2X,F20.12))') sum(gw_force(9*num+1:10*num)), sum(gw_force(10*num+1:11*num)), sum(gw_force(11*num+1:12*num)) 
  endif 

#endif 




! GWW correct for drift in forces 

!   frrx = frrx - sum(frrx) / num
!   frry = frry - sum(frry) / num
!   frrz = frrz - sum(frrz) / num


! GWW - vdw stress
  stsrxx   = scctmp(1)
  stsrxy   = scctmp(2)
  stsrxz   = scctmp(3)
  stsryy   = scctmp(4)
  stsryz   = scctmp(5)
  stsrzz   = scctmp(6)

#ifdef ewald_vdw
! GWW - correction of Ewald vdw stress tensor 
  stsrxx=stsrxx+dispcorrec1/cellvol
  stsryy=stsryy+dispcorrec1/cellvol
  stsrzz=stsrzz+dispcorrec1/cellvol
#endif 

! GWW - charge - charge stress tensor
  stcxx    = scctmp(7)
  stcxy    = scctmp(8)
  stcxz    = scctmp(9)
  stcyy    = scctmp(10)
  stcyz    = scctmp(11)
  stczz    = scctmp(12)

#ifdef dipole 
!GWW - charge - dipole stress 
  stpxx    = scctmp(13)
  stpxy    = scctmp(14)
  stpxz    = scctmp(15)
  stpyy    = scctmp(16)
  stpyz    = scctmp(17)
  stpzz    = scctmp(18)

!GWW - dipole - dipole stress 
  stp2xx   = scctmp(19)
  stp2xy   = scctmp(20)
  stp2xz   = scctmp(21)
  stp2yy   = scctmp(22)
  stp2yz   = scctmp(23)
  stp2zz   = scctmp(24)
#endif 

#ifdef quadrupole 
  stpqquadxx = scctmp(25)
  stpqquadxy = scctmp(26)
  stpqquadxz = scctmp(27)
  stpqquadyy = scctmp(28)
  stpqquadyz = scctmp(29)
  stpqquadzz = scctmp(30)

  stpdipquadxx = scctmp(31)
  stpdipquadxy = scctmp(32)
  stpdipquadxz = scctmp(33)
  stpdipquadyy = scctmp(34)
  stpdipquadyz = scctmp(35)
  stpdipquadzz = scctmp(36)

  stpquadquadxx = scctmp(37)
  stpquadquadxy = scctmp(38)
  stpquadquadxz = scctmp(39)
  stpquadquadyy = scctmp(40)
  stpquadquadyz = scctmp(41)
  stpquadquadzz = scctmp(42)
#endif 

#ifdef dipole 
! GWW - something to do with dipole damping - only in real space ?
  stpsrxx = scctmp(43)
  stpsrxy = scctmp(44)
  stpsrxz = scctmp(45)
  stpsryy = scctmp(46)
  stpsryz = scctmp(47)
  stpsrzz = scctmp(48)
#endif 

#ifdef quadrupole
! GWW - something to do with quadrupole damping - only in real space ?
  stqsrxx = scctmp(49)
  stqsrxy = scctmp(50)
  stqsrxz = scctmp(51)
  stqsryy = scctmp(52)
  stqsryz = scctmp(53)
  stqsrzz = scctmp(54)
#endif 

return
END SUBROUTINE
