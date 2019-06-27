SUBROUTINE setup

use error_function 
USE commondata
USE boxdata, ONLY: boxlenx,boxleny,boxlenz,halfboxminsq,bh,fullhi,fullhit,vg2, &
                   h2,hi2,a0,b0,c0,fullh,h
USE recipdata, ONLY: bdisp, conv1 
USE cgdata
USE lightdata
!---> Parallelization_S
use mpipara
!---> Parallelization_E
IMPLICIT NONE

integer :: ierror 

! local variables
INTEGER :: i,j,k,l,m,ii,jj,kk,nx,ny,nz,idum,nnn,ipoint,jpoint
DOUBLE PRECISION :: ran1,gauss,chgsum,dispsum1,dispsum2
DOUBLE PRECISION :: relaxsq,relaxsqrec,gtran
DOUBLE PRECISION, DIMENSION(nspec,600) :: bx, by, bz   !Modified by D. Marrocchelli 11/03/2008
integer :: iseed                                        ! MABC: These two lines contain elements for control


! local characters
CHARACTER(len=9) indi
CHARACTER(len=8) namex,namey,namez
CHARACTER(len=10) namesrx,namesry,namesrz
CHARACTER(len=10) namexx,nameyy,namezz,namexy,namexz,nameyz
CHARACTER(len=12) namesrxx,namesryy,namesrzz,namesrxy,namesrxz,namesryz
CHARACTER(len=30) filenamex,filenamey,filenamez
CHARACTER(len=30) filenamesrx,filenamesry,filenamesrz
CHARACTER(len=30) filenamexx,filenameyy,filenamezz,filenamexy,filenamexz,filenameyz
CHARACTER(len=30) filenamesrxx,filenamesryy,filenamesrzz        &
            ,filenamesrxy,filenamesrxz,filenamesryz
 
double precision :: dl_tol
indi = '123456789'

if((dippimlog).or.(quadpimlog)) then
!---> Parallelization_S
   if( iam .eq. 0 ) then

   namex = 'dipx.out'
   namey = 'dipy.out'
   namez = 'dipz.out'
   namesrx = 'dipsrx.out'
   namesry = 'dipsry.out'
   namesrz = 'dipsrz.out'
   do i=1,nummon
      nnn=0
      filenamex=namex//indi(i:i)
      filenamey=namey//indi(i:i)
      filenamez=namez//indi(i:i)
      filenamesrx=namesrx//indi(i:i)
      filenamesry=namesry//indi(i:i)
      filenamesrz=namesrz//indi(i:i)
      open(81+nnn,file=filenamex,status='new')
      open(82+nnn,file=filenamey,status='new')
      open(83+nnn,file=filenamez,status='new')
      open(84+nnn,file=filenamesrx,status='new')
      open(85+nnn,file=filenamesry,status='new')
      open(86+nnn,file=filenamesrz,status='new')
   enddo    

   endif
!---> Parallelization_E
endif

if(quadpimlog) then
!---> Parallelization_S
   if( iam .eq. 0 ) then

   namexx = 'quadxx.out'
   nameyy = 'quadyy.out'
   namezz = 'quadzz.out'
   namexy = 'quadxy.out'
   namexz = 'quadxz.out'
   nameyz = 'quadyz.out'
   namesrxx = 'quadsrxx.out'
   namesryy = 'quadsryy.out'
   namesrzz = 'quadsrzz.out'
   namesrxy = 'quadsrxy.out'
   namesrxz = 'quadsrxz.out'
   namesryz = 'quadsryz.out'
   do i=1,nummon
      nnn=0
      filenamexx=namexx//indi(i:i)
      filenameyy=nameyy//indi(i:i)
      filenamezz=namezz//indi(i:i)
      filenamexy=namexy//indi(i:i)
      filenamexz=namexz//indi(i:i)
      filenameyz=nameyz//indi(i:i)
      filenamesrxx=namesrxx//indi(i:i)
      filenamesryy=namesryy//indi(i:i)
      filenamesrzz=namesrzz//indi(i:i)
      filenamesrxy=namesrxy//indi(i:i)
      filenamesrxz=namesrxz//indi(i:i)
      filenamesryz=namesryz//indi(i:i)
      open(71+nnn,file=filenamexx,status='new')
      open(72+nnn,file=filenameyy,status='new')
      open(73+nnn,file=filenamezz,status='new')
      open(74+nnn,file=filenamexy,status='new')
      open(75+nnn,file=filenamexz,status='new')
      open(76+nnn,file=filenameyz,status='new')
      open(91+nnn,file=filenamesrxx,status='new')
      open(92+nnn,file=filenamesryy,status='new')
      open(93+nnn,file=filenamesrzz,status='new')
      open(94+nnn,file=filenamesrxy,status='new')
      open(95+nnn,file=filenamesrxz,status='new')
      open(96+nnn,file=filenamesryz,status='new')
   enddo   

   endif
!---> Parallelization_E
endif

!---> Parallelization_S
if( iam .eq. 0 ) then

open (21,file='eng1.out',status='new')
open (22,file='velocities.out',status='new')
open (23,file='positions.out',status='new')
open (25,file='eng2.out',status='new')
open (30,file='kintemp.out',status='new')
open (32,file='celleng.out',status='new')
open (33,file='engtot.out',status='new')
open (34,file='pzeta.out',status='new')
open (36,file='commomentum.out',status='new')
open (41,file='disp.out',status='new')    
open (44,file='bzeta.out',status='new')

! GWW
! forcefield fitting output
open (46,file='dipoles.out',status='new')
open (47,file='quads.out',status='new')
open (48,file='forces.out',status='new')
open (49,file='polar.out',status='new')


open (51,file='pres_diag.out',status='new')
open (52,file='pressure.out',status='new')
open (53,file='celllens.out',status='new')
open (54,file='cellvol.out',status='new')
open (55,file='cellangles.out',status='new')
open (56,file='diagstress.out',status='new')
open (57,file='xxyyzzstress.out',status='new')
open (58,file='xyxzyzstress.out',status='new')
open (60,file='poscart.out',status='new')
open (61,file='cellbox.out',status='new')
open (62,file='polstress.out',status='new')
open (63,file='coulsrstress.out',status='new')

endif
!---> Parallelization_E

call date_and_time(values=time_array)           !MABC: New seed for random number generator
    !iseed = 1
    !do i = 7, 8
    ! write(*,*) iseed, time_array(i)
     iseed = time_array(7)+time_array(6)
     write(*,*) iseed
    !end do

dummy=ran1(-iseed)

trantkb=boltz*trantemp
amass=amass/(avo*emass)
recamass=1.d0/amass
hmass=amass/2.d0
vartrans=dsqrt(trantkb/amass)
l=0
chgsum=0.d0
dispsum1=0.d0
dispsum2=0.d0
!Allocate all num-dimensional allocatable arrays - even those not used - stupid
!(again) 
ALLOCATE ( x(num),y(num),z(num) )
ALLOCATE ( vx(num),vy(num),vz(num) )
ALLOCATE ( frrx(num),frry(num),frrz(num) )
ALLOCATE ( frrx3(num),frry3(num),frrz3(num) )
ALLOCATE ( xdisp(num),ydisp(num),zdisp(num) )
ALLOCATE ( q(num), ntype(num) )
ALLOCATE ( elecx(num),elecy(num),elecz(num) )
ALLOCATE ( elecxq(num),elecyq(num),eleczq(num) )
ALLOCATE ( elecxsr(num),elecysr(num),eleczsr(num) )
ALLOCATE ( exx(num),eyy(num),ezz(num),exy(num),exz(num),eyz(num) )
ALLOCATE ( exxq(num),eyyq(num),ezzq(num),exyq(num),exzq(num),eyzq(num) )
ALLOCATE ( exxsr(num),eyysr(num),ezzsr(num),exysr(num),exzsr(num),eyzsr(num) )
ALLOCATE ( xk1(num),xk2(num),xk3(num),xk4(num) )
ALLOCATE ( alppolar(num),Bpolar(num),Cpolar(num),gammapolar(num) )
ALLOCATE ( engeff(num),xmu(num),ymu(num),zmu(num),quadxx(num),quadyy(num), &
           quadzz(num),quadxy(num),quadxz(num),quadyz(num),delta(num), &
           epsilonx(num),epsilony(num),epsilonz(num),quaimxx(num), &
           quaimyy(num),quaimzz(num),quaimxy(num),quaimxz(num),quaimyz(num) )
ALLOCATE ( selfeps(num), selfquaim(num) )
ALLOCATE ( reseng(num), dipsqeng(num), dipquadeng(num), quadeng(num) )
ALLOCATE ( srdipx(num),srdipy(num),srdipz(num) )
ALLOCATE ( asdipx(num),asdipy(num),asdipz(num) )
ALLOCATE ( srquadxx(num),srquadyy(num),srquadzz(num), &
           srquadxy(num),srquadxz(num),srquadyz(num) )
ALLOCATE ( asquadxx(num),asquadyy(num),asquadzz(num), &
           asquadxy(num),asquadxz(num),asquadyz(num) )
ALLOCATE ( PCOMAIM(num*10),XICOMAIM(num*10) ) !Maybe this is not good? AGUADO
ALLOCATE ( PCOMSTRUCT(3*num+6),XICOMSTRUCT(3*num+6) )
! The same for two-dimensional arrays...
!---> Memory Reduction_S
!MPIAGUADO-->Transform the erfc and dxsav 2D arrays to 1D arrays of
!            length equal to the number of ion pairs, num2.
!            The 2D array numadr contains the correspondence rule for
!            the transformation between the old (i,j) indexes and the new
!            1D arrays.
num2 = ( num - 1 )* num / 2

!ALLOCATE ( erfc(num,num) )
!ALLOCATE ( dxsav(num,num), dysav(num,num), dzsav(num,num) )
ALLOCATE ( erfc(num2) )
ALLOCATE ( dxsav(num2), dysav(num2), dzsav(num2), STAT=ierror )
if (ierror.ne.0) then
  print*, 'allocation error'
  STOP 999
endif 

!print*,'allocated dxsav ',num, num2 , (num-1)*num/2

ALLOCATE ( numadr(num,num) )
numx = 0
do j = 2, num
   do i = 1, j-1
      numx = numx + 1
      numadr(i,j) = numx
   enddo
enddo
!---> Memmory Reduction_E
!---> Parallelization_S
!MPIAGUADO--> eltmp, sctmp, dimtmp are temporal variables employed to store initially all
!             the contributions to the electric fields, field gradients, etc.
!             This is later transferred to the right variables...
!             Each processor evaluates a portion of eltmp and then sends its
!             result to the other processes. 
ALLOCATE ( eltmp(19*num), elltmp(19*num) )
ALLOCATE ( sctmp(59), scctmp(59) )
ALLOCATE ( dimtmp(nspec,2) )
!<--MPIAGUADO
!<--- Parallelization_E
ALLOCATE ( engft1(num,nspec),engft2(num,nspec),engft3(num,nspec), &
           engft1dotx(num,nspec),engft1doty(num,nspec),engft1dotz(num,nspec), &
           engft2dotx(num,nspec),engft2doty(num,nspec),engft2dotz(num,nspec), &
           engft3dotx(num,nspec),engft3doty(num,nspec),engft3dotz(num,nspec), &
        engft1dotxx(num,nspec),engft1dotyy(num,nspec),engft1dotzz(num,nspec), &
        engft1dotxy(num,nspec),engft1dotxz(num,nspec),engft1dotyz(num,nspec), &
        engft2dotxx(num,nspec),engft2dotyy(num,nspec),engft2dotzz(num,nspec), &
        engft2dotxy(num,nspec),engft2dotxz(num,nspec),engft2dotyz(num,nspec), &
        engft3dotxx(num,nspec),engft3dotyy(num,nspec),engft3dotzz(num,nspec), &
        engft3dotxy(num,nspec),engft3dotxz(num,nspec),engft3dotyz(num,nspec) )
!sgi
reseng=0.0d0
dipsqeng=0.0d0
dipquadeng=0.0d0
quadeng=0.0d0
!sgi
!---> Parallelization_S
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
!<--- Parallelization_E

do i=1,nspec
   if(chg(i).lt.0.0d0)nanion=nanion+nsp(i)
   if(chg(i).gt.0.0d0)ncation=ncation+nsp(i)
   do j=1,nsp(i)
      l=l+1
      ntype(l)=i
      q(l)=chg(i)
      chgsum=chgsum+(chg(i)*chg(i))
      dispsum2=dispsum2+ftc(i,i)
   enddo   
enddo    
do i=1,num
   ipoint=ntype(i)
   do j=1,num
      jpoint=ntype(j)
      dispsum1=dispsum1+ftc(ipoint,jpoint)
   enddo    
enddo    

if(.not.restart) then
   do i=1,nspec
      open (7,file=cellcoordfile(i),status='old')
      do j=1,nunitcellpos(i)
         read (7,*) bx(i,j),by(i,j),bz(i,j)
      enddo    
      close(7)
   enddo   
!
! the .mat files that contain the bx,y,z arrays are in the
! cell-frame and to be visualised, should be transformed into
! the lab-frame.
!
   m=1
   do i=1,nspec
      do j=1,nunitcellpos(i)
         do nx=1,nunitcellx
            do ny=1,nunitcelly
               do nz=1,nunitcellz
                  x(m)=a0*(nx-1)+bx(i,j)*a0
                  y(m)=b0*(ny-1)+by(i,j)*b0
                  z(m)=c0*(nz-1)+bz(i,j)*c0
                  m=m+1
               enddo   
            enddo   
         enddo   
      enddo   
   enddo   
!
! Setup displaced sublattice for phonon investigations:
!
   idum=-iseed
   if(displace) then
!---> Parallelization_S
      if( iam .eq. 0 ) then

      write (6,*) 'Randomly displacing ions off lattice.'

      endif
!---> Parallelization_E
      do i=1,num
         x(i)=x(i)+0.01d0*(ran1(idum+i)-0.5d0)  !MABC: Changed position assignment from using ran1 with
         y(i)=y(i)+0.01d0*(ran1(idum+i)-0.5d0)  ! the same seed to using a different seed for each atom
         z(i)=z(i)+0.01d0*(ran1(idum+i)-0.5d0)  ! and I decreased the magnitude of the displacement
      enddo   
   endif
!
! Enforce periodic boundary conditions.
!
   do i=1,num
      if(x(i).lt.0.0d0) x(i)=x(i)+boxlenx
      if(y(i).lt.0.0d0) y(i)=y(i)+boxleny
      if(z(i).lt.0.0d0) z(i)=z(i)+boxlenz
      if(x(i).gt.boxlenx) x(i)=x(i)-boxlenx
      if(y(i).gt.boxleny) y(i)=y(i)-boxleny
      if(z(i).gt.boxlenz) z(i)=z(i)-boxlenz
   enddo
!
! Output cell-frame coordinates.
!
!---> Parallelization_S
   if( iam .eq. 0 ) then

   open(50,file='cellframe.out',status='new')
   do i=1,num
      write(50,*) x(i),y(i),z(i)
   enddo   
   close(50)
!
! Output coordinates, in lab. frame:
!
   open(50,file='setuppos.out',status='new')
   do i=1,num
      write(50,*) h(1,1)*x(i)+h(1,2)*y(i)+h(1,3)*z(i), &
                  h(2,1)*x(i)+h(2,2)*y(i)+h(2,3)*z(i), &
                  h(3,1)*x(i)+h(3,2)*y(i)+h(3,3)*z(i)
   enddo   
   close(50)

   endif
!---> Parallelization_E

endif

dddamp2=dddamp*dddamp/2.d0
dddamp3=dddamp2*dddamp/3.d0
dddamp4=dddamp3*dddamp/4.d0
dddamp5=dddamp4*dddamp/5.d0
dddamp6=dddamp5*dddamp/6.d0
dqdamp2=dqdamp*dqdamp/2.d0
dqdamp3=dqdamp2*dqdamp/3.d0
dqdamp4=dqdamp3*dqdamp/4.d0
dqdamp5=dqdamp4*dqdamp/5.d0
dqdamp6=dqdamp5*dqdamp/6.d0
dqdamp7=dqdamp6*dqdamp/7.d0
dqdamp8=dqdamp7*dqdamp/8.d0

call boxreset
call dcell(fullh,bh)
halfboxminsq=(0.5d0*dmin1(bh(7),bh(8),bh(9)))**2
rsqmax=rcut**2.0d0
rsqmaxsr=rcutsr**2.0d0

! print*, '1/bh', bh(7), bh(8), bh(9)

fullhit=TRANSPOSE(fullhi)
call dcell(fullhit,bh)

! print*, '1/bhr', 1/bh(7), 1/bh(8), 1/bh(9)

!eta=etainpt/(2.0d0*rcut)
!etasq=eta*eta
!etapi=2.0d0*eta/sqrpi
!etaconst=-1.0d0/(4.0d0*eta*eta)


!eta = sqrt(-log(conv1)) / ( 0.5d0 * rcut )
!etasq=eta*eta
!etapi=2.0d0*eta/sqrpi
!etaconst=-1.0d0/(4.0d0*eta*eta)

eta = sqrt(-log(conv1)) / ( rcut )

! GWW try dlpoly parameters 
!dl_tol = sqrt (-log (conv1 * rcut)) 
!eta = sqrt(-log(conv1 *rcut * tol )) / ( rcut )


! GWW try gulp code - actually what I already had ! 
eta = sqrt(-log(conv1)) / ( rcut )
! rcut = rcut + 5.0 

etasq=eta*eta
etapi=2.0d0*eta/sqrpi
etaconst=-1.0d0/(4.0d0*eta*eta)

! BUT WHAT IS BEST ETA !!!!! GWW

!rcut = rcut + 1.0

! construct erfc table

  call construct_erfc_table 

! GWW - oh dear - not what os done in real space
! whoevery put this in is a real jackass 
!
! Factors needed for the Ewald summation of dispersion.
! Needs to be revised. At the moment it assumes that species
! 1 is always anionic and 2 cationic. It also assumes that
! C++=C+-*C+-/C--
! A. AGUADO. January 2001.
!
bdisp(1,1)=dsqrt(ftc(1,1))
do j=2,nspec
   bdisp(j,j)=ftc(1,j)/dsqrt(ftc(1,1))
enddo
dispalpsq=etasq
dispalp=dsqrt(dispalpsq)
dispalp3=dispalp*dispalpsq
dispalp4=dispalpsq*dispalpsq
dispalp6=dispalpsq*dispalp4

if(.not.restart) call kset

chgcorrec=chgsum*etapi*0.5d0
!print*,'chgcorrec =',chgcorrec 
dispcorrec1=-(pithreehalf/6.0d0)*dispalp3*dispsum1
dispcorrec2=dispsum2*dispalp6/12.0d0
!
! Thermostat parameters.
!
relaxsq=relax*relax
if(nth .or. nib)relaxsqrec=1.0d0/relaxsq
gtran=3.0d0*dble(num-1)*trantkb
gtrankin=2.0d0/(3.0d0*dble(num-1)*boltz)
!
! Initialize mean-squared displacement output file.
!
!---> Parallelization_S
if( iam .eq. 0 ) then

!write(41,*) nmsdcalltime,dtime,nrun  !!MABC: This file gets too big!!
!do i=1,num
!   write(41,*) ntype(i)
!enddo   

endif
!---> Parallelization_E

do i=1,num
   ipoint=ntype(i)
   selfeps(i)=selfeps1(ipoint)
   selfquaim(i)=selfquaim1(ipoint)
enddo

!---> Parallelization_S
if( iam .eq. 0 ) then

write (6,*)
write (6,*) '**** Internal variables set up ****'
write (6,*)

endif
!---> Parallelization_E

free=3.0d0*float(num-1)
W=(free+3.0d0)*trantkb*(relaxb**2.0d0)
Wgo=(free+3.0d0)*trantkb*(relaxb2**2.0d0)/3.0d0
if(nib)then
   Wrec=1.0d0/W
   Wgorec=1.0d0/Wgo
   dom=1.d0
endif
if(.not.nab) then
   Wgo=0.d0
   Wgorec=0.d0
   vg2=0.d0
else
   dom=9.d0
   if(ortho) then
      dom=3.0d0
      do j=1,3
         do k=1,3
            if(j.ne.k) then
               vg2(j,k)=0.d0
               h2(j,k)=0.d0
               hi2(j,k)=0.d0
            endif
         enddo   
      enddo   
   endif
endif
 
if(nth)CUEp=gtran/relaxsqrec
if(nth)CUEp2=CUEp/free
if(nth)CUEprec=1.d0/CUEp
if(nth)CUEp2rec=1.d0/CUEp2
if(.not.nth) then
   CUEprec=0.d0
   CUEp2rec=0.d0
endif
 
if(nib)CUEb=dom*trantkb/relaxsqrec
if(nib)CUEb2=CUEb/dom
if(nib)CUEbrec=1.0d0/CUEb
if(nib)CUEb2rec=1.0d0/CUEb2

if(.not.nth) then
   CUEbrec=0.d0
   CUEb2rec=0.d0
endif
if(.not.nib) then
   CUEbrec=0.d0
   CUEb2rec=0.d0
   CUEb2=0.d0
   CUEb=0.d0
endif

if(lscalc) then
!=======Light scattering=========================
                                                                                                             
bhyperpol3=hyperB/3.0d0
root3=dsqrt(3.0d0)
root12=dsqrt(12.0d0)
root1twodiv3gam=root12*hypergamma/3.0d0
twodiv3gamma=(2.0d0/3.0d0)*hypergamma
root3div3gam=(root3/3.0d0)*hypergamma
gammadiv3=hypergamma/3.0d0
polarundist=1.0d0/polarundist

normav=0
ncorrtime=1
ncorrcall=0

do i=1,6
   do j=1,10
      cohav(i,j)=0.0d0
      do k=0,ncorr
         lnorm(i,j,k)=0
         scf(i,j,k)=0.0d0
      enddo
   enddo
enddo
!=====================================================================
! Multiply the inputed light scattering preexponents by the
! exponential of the sum of the radii to allow scaling

do i=1,nspec
   do j=1,nspec
!---> Parallelization_S
     if( iam .eq. 0 ) then

     write(6,*)'before',asr(i,j),sig(i),sig(j)

     endif
!---> Parallelization_E
     asr(i,j)=asr(i,j)*dexp(csr(i,j)*(sig(i)+sig(j)))
     bsr(i,j)=bsr(i,j)*dexp(dsr(i,j)*(sig(i)+sig(j)))
!---> Parallelization_S
     if( iam .eq. 0 ) then

     write(6,*)'after',asr(i,j)

     endif
!---> Parallelization_E
   enddo
enddo
! Set up array to allow the terms to be calculated.
do ii=1,21
   do jj=1,ncfmat2
      do kk=0,ncorr
         scfiso(ii,jj,kk)=0.0d0
      enddo
   enddo
enddo
endif
!=========================end light scattering=============

! allocate for force check 
ALLOCATE ( gw_force(12*num) )

return
END SUBROUTINE
