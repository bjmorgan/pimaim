!************************************************************
!   Reads in the data from the .inpt file.
!************************************************************

SUBROUTINE readin

USE commondata
USE boxdata
USE recipdata, ONLY:nbin,convfac,conv1,bdisp, slens
USE lightdata
!---> Parallelization_S
use mpipara
!---> Parallelization_E

IMPLICIT NONE

INTEGER :: i,j,k,ii,n
CHARACTER(len=40) :: inptfile,FTfile,cffile,lcffile,cellfile

  call date_and_time(values = time_array)
    iseed1 = 1
    iseed2 = 1
    do i = 7, 8
    iseed1 = iseed1*time_array(6)
    end do
    do i = 6, 7
    iseed2 = iseed2*time_array(7)
    end do


!  Read the file containing the simulation details

open(10,file='runtime.inpt',status='old')
 read(10,*)jobname !MABC
 read(10,*)nrun
 read(10,*)trantemp
 read(10,*)nionunit
 read(10,*)nspec

!  Allocate all nspec-dimensional arrays
ALLOCATE ( nsp(0:nspec),chg(nspec),spectype(nspec),amass(nspec), &  
	   polarizablelog(nspec),deformablelog(nspec) ) 
ALLOCATE ( ftalp(nspec,nspec),ftb(nspec,nspec),ftc(nspec,nspec), &
           ftd(nspec,nspec),dddamp(nspec,nspec),dqdamp(nspec,nspec)  )
ALLOCATE ( rvph(nspec,nspec),rvpn(nspec,nspec),rvpr4(nspec,nspec) )
ALLOCATE ( ftbeta(nspec,nspec),ftgamma(nspec,nspec),ftb2(nspec,nspec), &
           ftb3(nspec,nspec)  )
ALLOCATE ( selfB(nspec),selfgam(nspec) )
ALLOCATE ( ftalpx(nspec,nspec), ftbx(nspec,nspec), nrpower(nspec,nspec) )
ALLOCATE ( selfeps1(nspec), selfeps2(nspec), selfC(nspec) )
ALLOCATE ( selfquaim1(nspec), selfH(nspec) )
ALLOCATE ( dampa(nspec,nspec), dampfac(nspec,nspec), nkdamp(nspec,nspec), &
           fgb(nspec,nspec), fgc(nspec,nspec), nkfg(nspec,nspec)  )
ALLOCATE ( alppolar1(nspec), alppolardel(nspec), alppolareff(nspec), &
           Bpolar1(nspec), Cpolar1(nspec), gammapolar1(nspec)  )
ALLOCATE ( nunitcellpos(nspec), cellcoordfile(nspec) )
ALLOCATE ( vartrans(nspec), recamass(nspec), hmass(nspec) )
ALLOCATE ( dddamp2(nspec,nspec),dddamp3(nspec,nspec),dddamp4(nspec,nspec), & 
           dddamp5(nspec,nspec),dddamp6(nspec,nspec),dqdamp2(nspec,nspec), &
           dqdamp3(nspec,nspec),dqdamp4(nspec,nspec),dqdamp5(nspec,nspec), &
           dqdamp6(nspec,nspec),dqdamp7(nspec,nspec),dqdamp8(nspec,nspec) )
ALLOCATE ( bdisp(nspec,nspec) )
ALLOCATE ( rdfpart(0:300,nspec,nspec) )
!---> Parallelization_S
ALLOCATE ( rdftot_w(0:300) )
ALLOCATE ( rdfpart_w(0:300,nspec,nspec) )
rdftot=0.0d0
!---> Parallelization_E
!sgi
alppolar1=0.0d0
alppolardel=0.0d0
alppolareff=0.0d0
rdfpart=0.0d0
Bpolar1=0.0d0
Cpolar1=0.0d0
gammapolar1=0.0d0
!sgi
           
 read(10,*)(spectype(i),i=1,nspec) !MABC
 read(10,*)(nsp(i),i=1,nspec)
 num=0
 do i=1,nspec
    num=num +nsp(i)
 enddo
ALLOCATE (atmnam(num),weight(num),chge(num))
 read(10,*)(chg(i),i=1,nspec)
 read(10,*)(amass(i),i=1,nspec)
 n=0                               !MABC
 do i=1,nspec
        j=nsp(i)
        do k=1,j
          atmnam(n+k)=spectype(i)
          weight(n+k)=amass(i)
          chge(n+k)=chg(i)
        enddo
        n = n + j
 enddo
 read(10,*)(polarizablelog(i),i=1,nspec)
 read(10,*)(deformablelog(i),i=1,nspec)

 read(10,*)dtime
 read(10,'(a)')runtype       !AGUADO Will this work??
    if(runtype(1:3).eq.'rim') rimlog=.true.
    if(runtype(1:6).eq.'dippim') dippimlog=.true.
    if(runtype(1:7).eq.'quadpim') then
       quadpimlog=.true.
       read(10,'(a)')ewtype
          if(ewtype(1:10).eq.'chargequad') chargequadlog=.true.
          if(ewtype(1:7).eq.'dipquad') dipquadlog=.true.
          if(ewtype(1:9).eq.'fullewald') fullewaldlog=.true.
    endif
 read(10,'(a)')runtype2
    if(runtype2(1:3).eq.'epp') epplog=.true.
    if(runtype2(1:3).eq.'cim') cimlog=.true.
    if(runtype2(1:4).eq.'daim') daimlog=.true.
    if(runtype2(1:5).eq.'quaim') quaimlog=.true.
 read(10,*) ooaimlog
 read(10,*) oodamplog
 read(10,*) cluslog
 read(10,*) environmentalpimlog
    if(environmentalpimlog)read (10,*) engtol
 read(10,*) environmentalaimlog
 read(10,*) conjgradlog
    if(conjgradlog) then
       read(10,*)tol
       read(10,*)ftol
    endif
 read(10,*) conjgradaimlog
    if(conjgradaimlog) then
       read(10,*)tolaim
       read(10,*)ftolaim
    endif
 read(10,*) restart
 read(10,*) velinit
 read(10,*) rescalelog
 read(10,*) displace
 read(10,*) moveions
 read(10,*) dynam
 read(10,*) relaxconfig
    if(relaxconfig) then
       read(10,*)tolstruc
       read(10,*)ftolstruc
       read(10,*)relaxcell
       if(relaxcell) then
          read(10,*)pextstruct
          read(10,*)(cellconstraints(i),i=1,6)
       else
          cellconstraints=0
       end if
    endif
 read(10,*)istrj,jstrj,keytrj        !MABC
#ifndef ppfit
 if(iam .eq. 0) then
 open (666,file='HISTORY',status='unknown')
 write(666,'(a150)')jobname
 write(666,'(4i10)') keytrj,2,num, (nrun/jstrj)
 endif
#endif
 read(10,*)npereng
 read(10,*)npervel
 read(10,*)nperfri
 read(10,*)npercell
 read(10,*)nperrdf
 read(10,*)nummon
ALLOCATE ( nmon(nummon) )
 read(10,*)(nmon(i),i=1,nummon)
 read(10,*)etainpt
    if(cluslog)etainpt=1.0d-20
 read(10,*)rcut
 read(10,*)conv1
 read(10,*)convfac
 read(10,*)rcutsr
 read(10,*)nth
    if(nth) then
       read(10,*)relax
    endif
 read(10,*)nrscalelog
    if(nrscalelog) then
       read(10,*)nrscale
    endif
 read(10,*)nib
    if(nib) then
       read(10,*)relaxb
       read(10,*)pext
    endif
 read(10,*)nab
    if(nab) read(10,*)relaxb2
 read(10,*)ortho
! T/p ramping logical
 read(10,*)tramplog
    if (tramplog) then
       read(10,*)nsteptramp
       read(10,*)deltatramp
    endif
 read(10,*)pramplog
    if (pramplog) then
       read(10,*)nsteppramp
       read(10,*)deltapramp
    endif
 read(10,*)ewlog
 read(10,*)verbose
  read(10,*)lscalc
  read(10,*)lseed  !MABC: Seeds for ransvu
  if(lseed) then
    read(10,*) iseed1
    read(10,*) iseed2
  endif
close(10)
    call ranin(iseed1,iseed2)

!    Read the file with the potential parameters
open(10,file='potential.inpt',status='old')

read(10,'(a)')pottype

   if(pottype(1:2).eq.'FT') then
      do i=1,nspec
         do j=i,nspec
            read(10,*)ftalp(i,j)
            read(10,*)ftb(i,j)
            read(10,*)ftc(i,j)
            read(10,*)ftd(i,j)
            read(10,*)dddamp(i,j)
            read(10,*)dqdamp(i,j)

            ftalp(j,i)=ftalp(i,j)
            ftb(j,i)=ftb(i,j)
            ftc(j,i)=ftc(i,j)
            ftd(j,i)=ftd(i,j)
            dddamp(j,i)=dddamp(i,j)
            dqdamp(j,i)=dqdamp(i,j)
         enddo   
      enddo   
   endif

   rvplog=.false.    
   if(pottype(1:3).eq.'RVP') then
      rvplog=.true.    
      do i=1,nspec
         do j=i,nspec
            read(10,*)rvph(i,j)
            read(10,*)rvpn(i,j)
            read(10,*)rvpr4(i,j)
            read(10,*)ftc(i,j)

            ftc(j,i)=ftc(i,j)
            rvph(j,i)=rvph(i,j)
            rvpn(j,i)=rvpn(i,j)
            rvpr4(j,i)=rvpr4(i,j)
         enddo   
      enddo   
   endif

   if((pottype(1:4).eq.'CIM').or.          &
      (pottype(1:4).eq.'DAIM').or.          &
      (pottype(1:5).eq.'QUAIM')) then
      do j=2,nspec
         read(10,*)ftalp(1,j)
         read(10,*)ftbeta(1,j)
         read(10,*)ftgamma(1,j)
         read(10,*)ftb(1,j)
         read(10,*)ftb2(1,j)
         read(10,*)ftb3(1,j)
         read(10,*)ftc(1,j)
         read(10,*)ftd(1,j)
         read(10,*)dddamp(1,j)
         read(10,*)dqdamp(1,j)
         read(10,*)ftalpx(1,j)
         read(10,*)ftbx(1,j)
         !read(10,*)nrpower2

         ftalp(j,1)=ftalp(1,j)
         ftbeta(j,1)=ftbeta(1,j)
         ftgamma(j,1)=ftgamma(1,j)
         ftb(j,1)=ftb(1,j)
         ftb2(j,1)=ftb2(1,j)
         ftb3(j,1)=ftb3(1,j)
         ftc(j,1)=ftc(1,j)
         ftd(j,1)=ftd(1,j)
         dddamp(j,1)=dddamp(1,j)
         dqdamp(j,1)=dqdamp(1,j)
      enddo

      do j=1,nspec
         if(deformablelog(j)) then
            read(10,*)selfB(j)
            read(10,*)selfgam(j)
         endif
      enddo
      read(10,*)extraalpha
      read(10,*)extrab

      do i=1,nspec
         read(10,*)ftalp(i,i)
         read(10,*)ftbeta(i,i)
         read(10,*)ftgamma(i,i)
         read(10,*)ftb(i,i)
         read(10,*)ftb2(i,i)
         read(10,*)ftb3(i,i)
         read(10,*)ftc(i,i)
         read(10,*)ftd(i,i)
         read(10,*)dddamp(i,i)
         read(10,*)dqdamp(i,i)
         read(10,*)ftalpx(i,i)
         read(10,*)ftbx(i,i)
      enddo   

      do i=2,nspec-1
         do j=i+1,nspec
            read(10,*)ftalp(i,j)
            read(10,*)ftbeta(i,j)
            read(10,*)ftgamma(i,j)
            read(10,*)ftb(i,j)
            read(10,*)ftb2(i,j)
            read(10,*)ftb3(i,j)
            read(10,*)ftc(i,j)
            read(10,*)ftd(i,j)
            read(10,*)dddamp(i,j)
            read(10,*)dqdamp(i,j)
            read(10,*)ftalpx(i,j)
            read(10,*)ftbx(i,j)

            ftalp(j,i)=ftalp(i,j)
            ftbeta(j,i)=ftbeta(i,j)
            ftgamma(j,i)=ftgamma(i,j)
            ftb(j,i)=ftb(i,j)
            ftb2(j,i)=ftb2(i,j)
            ftb3(j,i)=ftb3(i,j)
            ftc(j,i)=ftc(i,j)
            ftd(j,i)=ftd(i,j)
            dddamp(j,i)=dddamp(i,j)
            dqdamp(j,i)=dqdamp(i,j)
         enddo   
      enddo   
   endif

   xftlog=.false.    
   if(pottype(1:3).eq.'XFT') then
      xftlog=.true.
      do i=1,nspec
         do j=i,nspec
            read(10,*)ftalp(i,j)
            read(10,*)ftb(i,j)
            read(10,*)ftalpx(i,j)
            read(10,*)ftbx(i,j)
            read(10,*)nrpower(i,j)
            read(10,*)ftc(i,j)
            read(10,*)ftd(i,j)
            read(10,*)dddamp(i,j)
            read(10,*)dqdamp(i,j)

            ftalp(j,i)=ftalp(i,j)
            ftb(j,i)=ftb(i,j)
            ftalpx(j,i)=ftalpx(i,j)
            ftbx(j,i)=ftbx(i,j)
            ftc(j,i)=ftc(i,j)
            ftd(j,i)=ftd(i,j)
            dddamp(j,i)=dddamp(i,j)
            dqdamp(j,i)=dqdamp(i,j)
            nrpower(j,i)=nrpower(i,j)
         enddo   
      enddo   
   endif

   if((pottype(1:4).eq.'DAIM').or.(pottype(1:5).eq.'QUAIM')) then
      do j=1,nspec
         if(deformablelog(j)) then
            read (10,*) selfeps1(j)
            if(environmentalaimlog)read (10,*) selfeps2(j)
            read (10,*) selfC(j)
         endif
      enddo
   endif

   if(pottype(1:5).eq.'QUAIM') then
      do j=1,nspec
         if(deformablelog(j)) then
            read (10,*) selfquaim1(j)
            read (10,*) selfH(j)
         endif
      enddo
   endif

   dampa=0.d0
   dampfac=0.d0
   nkdamp=0
   fgb=0.d0
   fgc=0.d0
   nkfg=0
   if(dippimlog) then
      do i=1,nspec
         if(polarizablelog(i)) then
            read (10,*) alppolar1(i)
            if(environmentalpimlog) read (10,*) alppolardel(i)
            if(environmentalpimlog) read (10,*) alppolareff(i)
            do j=1,nspec
               if(oodamplog) then
                  read (10,*) dampa(j,i)
                  read (10,*) nkdamp(j,i)
                  read (10,*) dampfac(j,i)
               else
                  if(j.ne.i) then
                     read (10,*) dampa(j,i)
                     read (10,*) nkdamp(j,i)
                     read (10,*) dampfac(j,i)
                  endif
               endif
            enddo   
         endif
      enddo   
   endif

   if(quadpimlog) then
      do i=1,nspec
         if(polarizablelog(i)) then
            read (10,*) alppolar1(i)
            if(environmentalpimlog) read (10,*) alppolardel(i)
            if(environmentalpimlog) read (10,*) alppolareff(i)
            read (10,*) Bpolar1(i)
            read (10,*) Cpolar1(i)
            read (10,*) gammapolar1(i)
            do j=1,nspec
               if(oodamplog) then
                  read (10,*) dampa(j,i)
                  read (10,*) nkdamp(j,i)
                  read (10,*) dampfac(j,i)
                  read(10,*)fgb(j,i)
                  read(10,*)nkfg(j,i)
                  read(10,*)fgc(j,i)
               else
                  if(j.ne.i) then
                     read (10,*) dampa(j,i)
                     read (10,*) nkdamp(j,i)
                     read (10,*) dampfac(j,i)
                     read(10,*)fgb(j,i)
                     read(10,*)nkfg(j,i)
                     read(10,*)fgc(j,i)
                  endif
               endif
            enddo   
         endif
      enddo   
   endif

close(10)

open(10,file='cf.inpt',status='old')
   read(10,*)nbin
   read(10,*)msdcalllog
   if(msdcalllog)read(10,*)nmsdcalltime
   do i=1, nspec         !Modified by D. Marrocchelli 11/03/2008
      read(10,*)slens(i)
   enddo
close(10)


if(lscalc) then
!---> Parallelization_S
if( iam .eq. 0 ) then

write(6,*) ' Light scattering implemented - ONLY FOR CUBIC CELL!!!'

endif
!---> Parallelization_E
open(10,file='lightscatt.inpt',status='old')
read(10,*) lsint         ! ls - routines called every ?? steps
read(10,*) ncorr     !length of ls correlation functions
read(10,*)ncfmat2
read(10,*)nkmod2

!  ALLOCATE LIGHT SCAATERING ARRAYS=================================
allocate (txxli(num),tyyli(num),tzzli(num), &
             txyli(num),txzli(num),tyzli(num))
!.........contains DID polarizabilities (num)
                                                                               
allocate (elecxu(num),elecyu(num),eleczu(num))
!.........contains electric fields without the sr terms (num)
                                                                               
allocate (exxls(num),eyyls(num),ezzls(num),   &
           exyls(num),exzls(num),eyzls(num))
allocate (exxsrls(num),eyysrls(num),ezzsrls(num),  &
                    exysrls(num),exzsrls(num),eyzsrls(num))
!.........contains electric field gradients (num)
                                                                               
allocate (srxx(num),sryy(num),srzz(num),   &
                       srxy(num),srxz(num),sryz(num))
!.........contains short-range polarizabilities (num)
                                                                               
allocate ( cohacc(6,4),cohav(6,4),term(6,4))
!              accumulators (nsymmax,ncontribmax)

allocate (lnorm(6,10,0:ncorr))
! normaliser         (nsymmax,nindexmax,0:ncorr)

allocate (scf(6,10,0:ncorr))
!  correlations functions (nsymmax,nindexmax,0:ncorrtimemax)

allocate (storecohacc(6,4,ncorr))
!          storage (nsymmax,ncontribmax,ncorrtimemax)

allocate ( cisokr(6,ncfmat2),cisoki(6,ncfmat2))
!.........raman amplitudes (6,ncfmat2)

allocate ( storeisokr(6,ncfmat2,0:ncorr), &
                     storeisoki(6,ncfmat2,0:ncorr))
!........storage (6,ncfmatmax2,0:ncorrtimemax)

allocate ( scfiso(21,ncfmat2,0:ncorr))
!...............isotropic amplitudes ((21,ncfmatmax2,0:ncorrtimemax)

allocate (xkvec2(3,ncfmat2))
!==============================================================================



do i=1,3
    read(10,*)(xkvec2(i,j),j=1,ncfmat2)
enddo

do i=1,3
   do j=1,ncfmat2
      xkvec2(i,j)=dble(xkvec2(i,j))
   enddo
enddo
!=====================================================================

read(10,*)hyperB
read(10,*)hypergamma
read(10,*)polarundist
do i=1,nspec
   read(10,*) xlipol(i)
enddo
do  i=1,nspec
                                                                                                                 
   do j=i,nspec
      read(10,*)asr(i,j)
      read(10,*)bsr(i,j)
      read(10,*)csr(i,j)
      read(10,*)dsr(i,j)
                                                                                                                 
      asr(j,i)=asr(i,j)
      bsr(j,i)=bsr(i,j)
      csr(j,i)=csr(i,j)
      dsr(j,i)=dsr(i,j)
   enddo                                                                                                          
enddo                                                                                                      
do i=1,nspec
   read(10,*)sig(i)
enddo
close(10)
endif

! end of light scattering section====================================



! The cell vectors are read into the columns of h:
!     ( a_x  b_x  c_x )
!     ( a_y  b_y  c_y )
!     ( a_z  b_z  c_z )
open (10,file='crystal_cell.inpt',status='old')
   read(10,*) h(1,1),h(1,2),h(1,3)
   read(10,*) h(2,1),h(2,2),h(2,3)
   read(10,*) h(3,1),h(3,2),h(3,3)
   read(10,*) boxlenx
   read(10,*) boxleny
   read(10,*) boxlenz
   read(10,*) nunitcellx
   read(10,*) nunitcelly
   read(10,*) nunitcellz
   do i=1,nspec
      read(10,*) nunitcellpos(i)
      read(10,'(a)') cellcoordfile(i)
   enddo           
   read(10,*) a0
   read(10,*) b0
   read(10,*) c0
close(10)

call dcell(h,b)
vol3=boxlenx*boxleny*boxlenz*b(10)
hlab3=h
hlab2=h 
h3(:,1)=h(:,1)*boxlenx/(vol3**(1.0d0/3.0d0))
h3(:,2)=h(:,2)*boxleny/(vol3**(1.0d0/3.0d0))
h3(:,3)=h(:,3)*boxlenz/(vol3**(1.0d0/3.0d0))

call invert(h3,hi3)

if(dynam) then
   open(10,file='dynmat.inpt',status='old')
   read(10,*)amovefac
   read(10,*)tol
   read(10,*)ftol
   read(10,*)nmat1
   read(10,*)nmat2
   do i=1,3
      read(10,*)(amove(i,j),j=1,3)
   enddo    
   close(10)
endif

!---> Parallelization_S
if( iam .eq. 0 ) then

write (6,*)
write (6,*) '**** Run-time parameters read in ****'
write (6,*)

endif
!---> Parallelization_E

! Write out the potential in eV/Angstrom units
#ifndef ppfit
open(999, file = 'ev_ang_potential', status='replace')
if(xftlog .eqv. .true.) then
	do i = 1, nspec
		do j = i, nspec
			write(999,*) 'Short ranged potential for interaction',i,'-',j
			write(999,*) ftalp(i,j)*1.8897261,	'     Yukawa Decay (numerator!)'
			write(999,*) ftb(i,j)*14.399644,	'     Yukawa Prefactor'
			write(999,*) ftalpx(i,j)*3.5710649,	'     Gaussian Decay'
			write(999,*) ftbx(i,j)*27.211383,	'     Gaussian Prefactor'
			write(999,*) nrpower(i,j),			'     Yukawa reciprocal power'
			write(999,*) ftc(i,j)*0.59752683,	'     C6 coefficient'
			write(999,*) ftd(i,j)*0.16732455,	'     C8 coefficient'
			write(999,*) dddamp(i,j)*1.8897261,	'     C6 Tang-Toennies decay'
			write(999,*) dqdamp(i,j)*1.8897261,	'     C8 Tang-Toennies decay'
			write(999,*) ''
		enddo
	enddo
endif
#endif


return
END SUBROUTINE
