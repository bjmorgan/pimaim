SUBROUTINE rstrun

  USE commondata, ONLY: num,x,y,z,vx,vy,vz,rescalelog,ntotstp,nsofar,nstpbrk, &
                      nstep,nrun,pzeta3,bzeta3,vpzeta3,vbzeta3,vpzeta2, &
                      vbzeta2,eps2,eps3,veps2,veps3,vol3,pi,xmu,ymu,zmu, &
                      quadxx,quadyy,quadzz,quadxy,quadxz,quadyz,delta, &
                      epsilonx,epsilony,epsilonz, &
                      quaimxx,quaimyy,quaimzz,quaimxy,quaimxz,quaimyz, &
                      dippimlog, quadpimlog,cimlog,daimlog,quaimlog,rimlog, &
                      conjgradlog,rsqmax,rcut, nsp, nspec 
  USE boxdata, ONLY: vg3,h,boxlenx,boxleny,boxlenz,h3,hlab2,hlab3,halfboxminsq, &
                   b,h3,hi3,fullh

  use mpipara

  IMPLICIT NONE

  INTEGER :: i,j,k
  LOGICAL :: coord_input_log,vel_input_log,dip_input_log, quad_input_log 
  LOGICAL :: cim_input_log, daim_input_log, quaim_input_log, fullrun_input_log
  CHARACTER(len=50) :: dummy_string 
  integer :: num_for_each_type(nspec)

   open(35,file='restart.dat',status='old')
   read(35,*)coord_input_log 
   read(35,*)vel_input_log
   read(35,*)dip_input_log, quad_input_log 
   read(35,*)cim_input_log, daim_input_log, quaim_input_log 
   read(35,*)fullrun_input_log




! read in lattice vectors - put in sensible order for Version 2.0 
    read (35,'(A50)') dummy_string 
    read(35,*)boxlenx, boxleny, boxlenz
    read(35,*) h(1,1),h(2,1),h(3,1)
    read(35,*) h(1,2),h(2,2),h(3,2)
    read(35,*) h(1,3),h(2,3),h(3,3)

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

    if (iam.eq.0.and.rsqmax.gt.halfboxminsq)  &
            write(6,*)'WARNING - cut off large than half smallest box width'

    call kset


! GWW - not reall sure the point of this  
!   if( iam .eq. 0 ) then
!
!   open(50,file='startuppos.out',status='new')
!      do i=1,num
!         write(50,*) h(1,1)*x(i)+h(2,1)*y(i)+h(3,1)*z(i), &
!                     h(1,2)*x(i)+h(2,2)*y(i)+h(3,2)*z(i), &
!                     h(1,3)*x(i)+h(2,3)*y(i)+h(3,3)*z(i)
!      enddo
!   close(50)
!   endif

    if(coord_input_log) then
! read in number of each type - check against runtime.inpt values  
      read (35,'(A50)') dummy_string 
      read (35,*) (num_for_each_type(i), i=1,nspec )

      do i=1,nspec
        if (nsp(i).ne.num_for_each_type(i)) then
           if (iam.eq.0) write(6,*) 'WARNING - number of atoms of type ',i,' does not match runtime.inpt'
       endif
      enddo 

! read in coordinates 
      read (35,'(A50)') dummy_string 
      do i=1,num
        read(35,*)x(i),y(i),z(i)
      enddo   
    else
      if( iam .eq. 0 ) then
        write(6,*)'No coordinate information in restart file!'
        write(6,*)'STOPPING !'
      endif
      call mpi_finalize(ierr)
      stop
    endif

! read in velocities 
    if(vel_input_log) then
      read (35,'(A50)') dummy_string 
      do i=1,num
        read(35,*)vx(i),vy(i),vz(i)
      enddo   
    else
      if(rescalelog) call velset
    endif


! read in dipole / quadrupole data 
    if(dip_input_log) then
      read (35,'(A50)') dummy_string 

! read dipoles if this is actually a dipole run 
      if(dippimlog.or.quadpimlog) then 
        do i=1,num
          read(35,*)xmu(i),ymu(i),zmu(i)
        enddo   
! or read the line but do not store - as at some point will only allocate
! appropriate storage for the run 
      else 
        if( iam .eq. 0 ) write(6,*) 'WARNING - dipole data present but not a dipole run'
        do i=1,num
          read (35,'(A50)') dummy_string 
        enddo   
      endif 
    endif 


   if(quad_input_log) then
      read (35,'(A50)') dummy_string 
! read quadrupoles if present and a quadrupole run 
      if(dippimlog.or.quadpimlog) then 
        do i=1,num
          read(35,*)quadxx(i),quadyy(i),quadzz(i)
          read(35,*)quadxy(i),quadxz(i),quadyz(i)
        enddo   
! or read the line but do not store - as at some point will only allocate
! appropriate storage for the run 
      else 
        if( iam .eq. 0 ) write(6,*) 'WARNING - quadrupole data present but not a quadrupole run'
        do i=1,num
          read (35,'(A50)') dummy_string 
          read (35,'(A50)') dummy_string 
        enddo   
      endif 
    endif 



! cim, diam and quaim input 

! cim input 
    if (cim_input_log.or.daim_input_log.or.quaim_input_log) then 
      if(cimlog.or.daimlog.or.quaimlog) then
        read (35,'(A50)') dummy_string 
        do i=1,num
          read(35,*)delta(i)
        enddo   
      else
        if( iam .eq. 0 ) write(6,*) 'WARNING - cim info but not a cim run'
        read (35,'(A50)') dummy_string 
        do i=1,num
          read (35,'(A50)') dummy_string 
        enddo    
      endif
    endif

! daim input
    if(daim_input_log .or. quaim_input_log) then 
      if (daimlog.or.quaimlog) then
        read (35,'(A50)') dummy_string 
        do i=1,num
          read(35,*)epsilonx(i),epsilony(i),epsilonz(i)
        enddo    
      else
        if( iam .eq. 0 ) write(6,*) 'WARNING - daim info but not a daim run'
        read (35,'(A50)') dummy_string 
        do i=1,num
          read (35,'(A50)') dummy_string 
        enddo    
      endif 
    endif

! quaim  input
    if(quaim_input_log) then
      if (quaimlog) then
        read (35,'(A50)') dummy_string 
        do i=1,num
          read(35,*)quaimxx(i),quaimyy(i),quaimzz(i)
          read(35,*)quaimxy(i),quaimxz(i),quaimyz(i)
        enddo   
      else 
        if( iam .eq. 0 ) write(6,*) 'WARNING - quaim info but not a quaim run'
        read (35,'(A50)') dummy_string 
        do i=1,num
          read (35,'(A50)') dummy_string 
          read (35,'(A50)') dummy_string 
        enddo    
      endif
    endif

!

    if(.not.dip_input_log.and. .not. conjgradlog .and. .not. rimlog) then
       if( iam .eq. 0 ) write(6,*)'WARNING No dipole data read in but a dipole run ' 
    endif

    if(fullrun_input_log) then
      ntotstp=0
      read (35,'(A50)') dummy_string 
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

    close(35)

    if (iam.eq.0) then 
      write(6,*)
      write(6,*)'**** Run restarted ****'
      write(6,*)
    endif 

  return

END SUBROUTINE
