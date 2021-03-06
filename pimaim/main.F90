PROGRAM main

USE commondata

USE boxdata, ONLY: vg2,vg3

USE recipdata, ONLY: norm_ds,sk_ds,elcall,elsall,emcall,emsall,encall,ensall, &
                     kmaxx,kmaxy,kmaxz,nbin,nspairs, &
                     sk_ds_w

USE lightdata, ONLY: lsint,ncorrtime,ncorrcall,ncorr,elecxu,elecyu,eleczu

use mpipara

IMPLICIT NONE

LOGICAL :: endrun, shutdown

INTEGER :: nnn,i

double precision :: time0, time1, time2

call mpi_init(ierr)
call mpi_comm_size(mpi_comm_world,nprocs,ierr)
call mpi_comm_rank(mpi_comm_world,iam,ierr)

time0 = mpi_wtime()

if (iam.eq.0) then

  write(*,*) '****************************************************************'
  write(*,*) '*                                                              *'
  write(*,*) '*                                                              *'
  write(*,*) '*         ######   ###  #     #     #     ###  #     #         *'
  write(*,*) '*         #     #   #   ##   ##    # #     #   ##   ##         *'
  write(*,*) '*         #     #   #   # # # #   #   #    #   # # # #         *'
  write(*,*) '*         ######    #   #  #  #  #     #   #   #  #  #         *'
  write(*,*) '*         #         #   #     #  #######   #   #     #         *'
  write(*,*) '*         #         #   #     #  #     #   #   #     #         *'
  write(*,*) '*         #        ###  #     #  #     #  ###  #     #         *'
  write(*,*) '*                                                              *'
  write(*,*) '*                           Version 2.1                        *'
  write(*,*) '*                                                              *'
  write(*,*) '*      maintained by Prof G. Watson, TCD (watsong@tcd.ie)      *'
  write(*,*) '*                                                              *'
  write(*,*) '****************************************************************'
  write(*,*) ' '
  write(*,*) '     A molecular dynamics code orgininally developed in the '
  write(*,*) '         group of Prof. Paul Madden at Oxford University'
  write(*,*) ' '
  write(*,*) ' '
  write(*,*) ' Recently updates at TCD'
  write(*,*) '    1) Rewrite of the Ewald cutoff in the correct units and to '
  write(*,*) '       and to ensure equal accuracy in real and reciprical space'
  write(*,*) ' '
  write(*,*) '    2) Rewrite of van der Waals interaction to include correct '
  write(*,*) '       energetics when using damping and the Ewald summation'
  write(*,*) ' '
  write(*,*) '    3) Combination of all energy routines charge, dipole and'
  write(*,*) '       quadrupole into a single routine with conditional '
  write(*,*) '       compilation'
  write(*,*) ' '
  write(*,*) '    4) New compilation systems allowing condition electrostatics,'
  write(*,*) '       van der walls, debugging etc.'
  write(*,*) ' '
  write(*,*) ' Version 2.1 - maintained by Prof G. Watson, TCD (watsong@tcd.ie)'
  write(*,*) ' '
  write(*,*) ' '

endif 

!<--- Parallelization_E
endrun=.false.
shutdown=.false.
nstep=0
nrdfcall=0

vg2=0.d0
vg3=0.d0
pzeta2=0.d0
pzeta3=0.d0
bzeta2=0.d0
bzeta3=0.d0
vpzeta1=0.d0
vpzeta2=0.d0
vpzeta3=0.d0
vbzeta1=0.d0
vbzeta2=0.d0
vbzeta3=0.d0
eps2=0.d0
eps3=0.d0
veps1=0.d0
veps2=0.d0
veps3=0.d0
W=0.d0
Wrec=0.d0
Wgo=0.d0
Wgorec=0.d0
dom=0.d0

num=0
nanion=0
ncation=0

norm_ds=0

call readin                !Read input data

nsp(0)=0



! only need for debye_schere ?
nspairs=0
do i=1,nspec
   nspairs=nspairs+i
enddo

ALLOCATE ( sk_ds(nspairs,0:nbin) )
ALLOCATE ( sk_ds_w(nspairs,0:nbin) )
sk_ds=0.d0


call setup                 !Initialise some constants and settings

! GWW - very shite way of distributing pair interaction - all of the them 
! every i-j irrespective of the cutoff and size of the cell - really stupid 
call numpara

engeff=0.d0

! only need thes eif dipole ! 
xmu=0.d0
ymu=0.d0
zmu=0.d0

! Only need these if quadrupole 
quadxx=0.d0
quadyy=0.d0
quadzz=0.d0
quadxy=0.d0
quadxz=0.d0
quadyz=0.d0
delta=0.d0
epsilonx=0.d0
epsilony=0.d0
epsilonz=0.d0
quaimxx=0.d0
quaimyy=0.d0
quaimzz=0.d0
quaimxy=0.d0
quaimxz=0.d0
quaimyz=0.d0
if(environmentalaimlog)then
   selfeps=0.d0
   selfquaim=0.d0
endif

engft1=0.d0
engft2=0.d0
engft3=0.d0
engft1dotx=0.d0
engft1doty=0.d0
engft1dotz=0.d0
engft2dotx=0.d0
engft2doty=0.d0
engft2dotz=0.d0
engft3dotx=0.d0
engft3doty=0.d0
engft3dotz=0.d0
engft1dotxx=0.d0
engft1dotyy=0.d0
engft1dotzz=0.d0
engft1dotxy=0.d0
engft1dotxz=0.d0
engft1dotyz=0.d0
engft2dotxx=0.d0
engft2dotyy=0.d0
engft2dotzz=0.d0
engft2dotxy=0.d0
engft2dotxz=0.d0
engft2dotyz=0.d0
engft3dotxx=0.d0
engft3dotyy=0.d0
engft3dotzz=0.d0
engft3dotxy=0.d0
engft3dotxz=0.d0
engft3dotyz=0.d0

xk1=0.d0
xk2=0.d0
xk3=0.d0
xk4=0.d0

nstep=1

if(restart) call rstrun      !Restart from a previous run...


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GWW calculation of dynamical matrix and phonons 
  if(dynam) then        !Calculates the dynamical matrix in 
    call dynmat        !order to calculate the phonons...
    call mpi_finalize(ierr)
    stop
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(velinit) call velset  !Sets up initial velocities...

ALLOCATE ( elcall(num,0:kmaxx+1), elsall(num,0:kmaxx+1) )
ALLOCATE ( emcall(num,0:kmaxy+1), emsall(num,0:kmaxy+1) )
ALLOCATE ( encall(num,0:kmaxz+1), ensall(num,0:kmaxz+1) )


! GWW - anothner shite parallelisation over first kvectopr index.
! final process has hardly  any work
! only thing that made it look good was the incorrect reciprical cut off
! which  menat it worked as a very large cube 

   allocate( kmaxy_s(0:kmaxx), kmaxy_e(0:kmaxx) )
   allocate( kmaxz_s(-kmaxy:kmaxy,0:kmaxx), kmaxz_e(-kmaxy:kmaxy,0:kmaxx) )

   call kmaxpara


! GWW - calculate the dipoles, quadrupoles, etc... for the inital
! configuration 
if(environmentalpimlog) then
   firstiter=.true.
   call conjgradpimaim
   firstiter=.false.
else
   if(conjgradaimlog) call conjgradaim
   if(conjgradlog) call conjgrad
endif

!Now to rescale the atom velocities, if appropriate...
  if(rescalelog) then
    if(.not. velinit .and. .not. restart) then
      if( iam .eq. 0 ) then
        write (*,*) 'Cannot rescale from zero temperature!'
      endif
      call mpi_finalize(ierr)
      stop
    endif
    if( iam .eq. 0 ) then
      write(*,*) 'rescale called.' 
    endif
    call rescale
  endif

! Mai loop. - NOT ! 
forfl=.false.

  if( iam .eq. 0 ) then
    write(*,*) 'main loop begun'
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GWW - structural relaxation using conjugate gradients 
  if(relaxconfig) then
    call conjstruct

    shutdown=.true.
    pereng=.true.
    percell=.true.
    pervel=.true.
    if( iam .eq. 0 ) call output

    veldumplog=.false.
    crddumplog=.true.
    chgdumplog=.false.
    fulldumplog=.true.
    fileout='testout.rst'
    if( iam .eq. 0 ) call dump

    call mpi_finalize(ierr)
    stop
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  if( iam .eq. 0 ) then
!    write(*,*) ' start sub ener '
!  endif

  call ener

! GWW - energy of initial configuration only 
!Static energy calculation...
if(nrun.eq.0) then
   shutdown=.true.
   pereng=.true.
   percell=.true.
   pervel=.true.
!---> Parallelization_S
!  call output
   if( iam .eq. 0 ) call output
!<--- Parallelization_E
   veldumplog=.false.
   crddumplog=.true.
   chgdumplog=.false.
   fulldumplog=.true.
   fileout='testout.rst'
!---> Parallelization_S
!  call dump
   if( iam .eq. 0 ) call dump
   call mpi_finalize(ierr)
!<--- Parallelization_E
   stop
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GWW - Start of main MD loop 

  do while(nstep.le.nrun)

    if( iam .eq. 0 ) then
      if(mod(nstep,20).eq.0)then
        write (*,*) nstep
      endif
    endif

! Rescale velocities?
  if(nrscalelog)then
    if(mod(float(nstep),float(nrscale)).eq.0) then
      call rescale
#ifdef debug
        if( iam .eq. 0 ) then
          write(*,*)'Rescale called'
        endif
#endif 
      endif
   endif

!Call routine to move the ions according to NVE,NVT or NPT ensembles...
   call trans_vv

#ifdef debug 
   if( iam .eq. 0 ) then
      do i=1,num
         write(48,*)frrx(i),frry(i),frrz(i)
      enddo
   endif
#endif 

! ****** Periodic parts of main loop ******
!
! Periodic dump file output.
! 
   if(mod(float(nstep),float(npereng)).eq.0) then
!   if(mod(float(nstep),float(npervel)).eq.0) then
      veldumplog=.true.
      crddumplog=.true.
      chgdumplog=.true.
      fulldumplog=.true.
      fileout='testout.rst'
!      if( iam .eq. 0 ) call dump
!   endif
!
! Periodic energies output.
!
!   if(mod(float(nstep),float(npereng)).eq.0) then
      pereng=.true.
!      if( iam .eq. 0 ) call output
!   endif
!
! Periodic velocities output.
!
!   if(mod(float(nstep),float(npervel)).eq.0) then
      pervel=.true.
!      if( iam .eq. 0 ) call output
!   endif
!
! Periodic temperature output.
!
!   if(mod(float(nstep),float(nperfri)).eq.0) then
      perfric=.true.
!      if( iam .eq. 0 ) call output
!   endif

!   if(mod(float(nstep),float(npercell)).eq.0) then
      percell=.true.
!     call output
      if( iam .eq. 0 ) call output
   endif
!===============================================================
! periodic T ramping if required.
   if ((tramplog).and.(mod(float(nstep),float(nsteptramp)).eq.0)) &
      call tramp
   if ((pramplog).and.(mod(float(nstep),float(nsteppramp)).eq.0)) &
      call pramp
!===============================================================
! Periodic radial distribution function calculation.
!
   if(mod(float(nstep),float(nperrdf)).eq.0) then
      call rdf
!!    if( iam .eq. 0 ) call rdf
   endif
!=====================================================================


#ifdef light_scattering 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GWW - light scattering routines
   if (lscalc) then
      if (mod(float(nstep),float(lsint)).eq.0) then

!********** electric fields used in light scattering include SR damping
         elecxu=elecx
         elecyu=elecy
         eleczu=elecz

         call lightfields
         call shortrange
         call lightarray
         call lightcfcalc
!
! update correlation fn array pointers

         ncorrtime=mod(ncorrtime,ncorr)
         ncorrcall=ncorrcall+1
         ncorrtime=ncorrtime+1
     endif
   endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GWW - rather stupid - just break out of loop why set endrun then do a load 
! of stuff it is set and then leave the loop ? 
   if(nstep.eq.nrun) endrun=.true.

   if(endrun) then
      shutdown=.true.
      if( iam .eq. 0 ) call output

      veldumplog=.true.
      crddumplog=.true.
      chgdumplog=.true.
      fulldumplog=.true. 
      fileout='testout.rst'
      if( iam .eq. 0 ) then
        call dump
        call rdfout
#ifdef debye_scherer
        call debye_scherer_out
#endif 
      endif

!  Dump the light scattering correlation functions and the store arrays 
!required to restart the calculation.

#ifdef light_scattering
      if (lscalc) then
         call lightcfdump
      endif
#endif 

!---> Parallelization_S
      if( iam .eq. 0 ) then

       close(21)
       close(22)
       close(23)
       close(25)
       close(26)
       close(30)
       close(58)
       close(59)

       do i=1,nummon
         nnn=0
         close(81+nnn)
         close(82+nnn)
         close(83+nnn)
       enddo    

       endif
     endif
! end of strange sectiopn in main loop which dump stuff if finishing 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! GWW - write HISTORY file 
    if( iam .eq. 0) call output_MABC    

    nstep=nstep+1
  enddo
! GWW - end of main MD loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! GWW - timing for run 
  time1 = mpi_wtime()
  time2 = time1 - time0
  if( iam .eq. 0 ) then
    write(*,'(a,F20.8,a)') ' Elapsed Time : ',time2,' (Sec.)'
  endif

! GWW - close down MIP and stop 
  call mpi_finalize(ierr)
  stop

END PROGRAM
