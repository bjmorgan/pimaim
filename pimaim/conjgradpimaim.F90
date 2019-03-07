SUBROUTINE conjgradpimaim

use error_function 
USE commondata, ONLY: num,delta,epsilonx,epsilony,epsilonz,quaimxx,quaimyy, &
                      quaimzz,quaimxy,quaimxz,quaimyz,x,y,z,dxsav,dysav, &
                      dzsav,eta,erfc,tolaim,ftolaim,ftol,engpetot,dippimlog, &
                      quadpimlog,conjgradlog,engtol,environmentalaimlog,xmu, &
                      ymu,zmu,quadxx,quadyy,quadzz,quadxy,quadxz,quadyz, &
                      selfeps,selfquaim,alppolar,Bpolar,Cpolar,gammapolar, &
                      verbose,firstiter,ntype,alppolar1,alppolardel, &
                      alppolareff,Bpolar1,Cpolar1,gammapolar1,engeff, &
!---> Memmory Reduction_S
!                     relaxconfig
                      relaxconfig,num2,numx, &
                      cimlog,daimlog,quaimlog,ooaimlog,epplog, &
                      elecxq,elecyq,eleczq,exxq,eyyq,ezzq,exyq,exzq,eyzq, &
                      engeff,sqrpi,etapi,q,ewcorr,elecxsr,elecysr,eleczsr, &
                      exxsr,eyysr,ezzsr,exysr,exzsr,eyzsr
!<--- Memmory Reduction_E
USE recipdata, ONLY: elcall,elsall,emcall,emsall,encall,ensall,kmaxx,kmaxy, &
                     kmaxz
USE boxdata, ONLY: twopiboxx,twopiboxy,twopiboxz,boxlenx,boxleny,boxlenz, &
                   halfboxx,halfboxy,halfboxz,hlab2

!---> Parallelization_S
use mpipara
!<--- Parallelization_E
IMPLICIT NONE

DOUBLE PRECISION :: fac2
DOUBLE PRECISION :: dxcf,dycf,dzcf,ftoltemp,ftolaimtemp,engtemp,drsq,dr
DOUBLE PRECISION, DIMENSION(10*num) :: p
INTEGER :: i,l,j,ipoint

p(1:num)=delta
p(num+1:2*num)=epsilonx
p(2*num+1:3*num)=epsilony
p(3*num+1:4*num)=epsilonz
p(4*num+1:5*num)=quaimxx
p(5*num+1:6*num)=quaimyy
p(6*num+1:7*num)=quaimzz
p(7*num+1:8*num)=quaimxy
p(8*num+1:9*num)=quaimxz
p(9*num+1:10*num)=quaimyz

elcall(:,0)=1.d0
emcall(:,0)=1.d0
encall(:,0)=1.d0
elsall(:,0)=0.d0
emsall(:,0)=0.d0
ensall(:,0)=0.d0
elcall(:,1)=cos(twopiboxx*x)
emcall(:,1)=cos(twopiboxy*y)
encall(:,1)=cos(twopiboxz*z)
elsall(:,1)=sin(twopiboxx*x)
emsall(:,1)=sin(twopiboxy*y)
ensall(:,1)=sin(twopiboxz*z)

do l=2,kmaxx
      elcall(:,l)=elcall(:,l-1)*elcall(:,1)-elsall(:,l-1)*elsall(:,1)
      elsall(:,l)=elsall(:,l-1)*elcall(:,1)+elcall(:,l-1)*elsall(:,1)
enddo   
do l=2,kmaxy
      emcall(:,l)=emcall(:,l-1)*emcall(:,1)-emsall(:,l-1)*emsall(:,1)
      emsall(:,l)=emsall(:,l-1)*emcall(:,1)+emcall(:,l-1)*emsall(:,1)
enddo   
do l=2,kmaxz
      encall(:,l)=encall(:,l-1)*encall(:,1)-ensall(:,l-1)*ensall(:,1)
      ensall(:,l)=ensall(:,l-1)*encall(:,1)+encall(:,l-1)*ensall(:,1)
enddo   

!---> Parallelization_S
!do j=2,num
!   do i=1,j-1
do j=jst,jed
   do i=ist(j),ied(j)
!<--- Parallelization_E
      dxcf=x(i)-x(j)
      dycf=y(i)-y(j)
      dzcf=z(i)-z(j)
      if(abs(dxcf).gt.halfboxx) dxcf=dxcf-sign(boxlenx,dxcf)
      if(abs(dycf).gt.halfboxy) dycf=dycf-sign(boxleny,dycf)
      if(abs(dzcf).gt.halfboxz) dzcf=dzcf-sign(boxlenz,dzcf)
!---> Memmory Reduction_S
      numx = numadr(i,j)
!     dxsav(i,j)=hlab2(1,1)*dxcf+hlab2(1,2)*dycf+hlab2(1,3)*dzcf
!     dysav(i,j)=hlab2(2,1)*dxcf+hlab2(2,2)*dycf+hlab2(2,3)*dzcf
!     dzsav(i,j)=hlab2(3,1)*dxcf+hlab2(3,2)*dycf+hlab2(3,3)*dzcf
!     drsq=dxsav(i,j)*dxsav(i,j)+dysav(i,j)*dysav(i,j)+dzsav(i,j)*dzsav(i,j)
      dxsav(numx)=hlab2(1,1)*dxcf+hlab2(1,2)*dycf+hlab2(1,3)*dzcf
      dysav(numx)=hlab2(2,1)*dxcf+hlab2(2,2)*dycf+hlab2(2,3)*dzcf
      dzsav(numx)=hlab2(3,1)*dxcf+hlab2(3,2)*dycf+hlab2(3,3)*dzcf
      drsq=dxsav(numx)*dxsav(numx)+dysav(numx)*dysav(numx)+dzsav(numx)*dzsav(numx)
!<--- Memmory Reduction_E
      dr=dsqrt(drsq)
!---> Memmory Reduction_S
!     erfc(i,j)=erfunc(eta*dr)
      erfc(numx)=erfunc(eta*dr)
!*<--- Memmory Reduction_E
   write(911,'(i10,4(X,F12.6))') numx,dxsav(numx),dysav(numx),dzsav(numx), dr 
   enddo   
enddo   
!---> Parallelization_S
if( nprocs .ne. 1 ) then
   if( cimlog .or. daimlog .or. quaimlog ) then
      do j=jst2,jed2
         do i=ist2(j),ied2(j)
            numx = numadr(i,j)
            dxcf=x(i)-x(j)
            dycf=y(i)-y(j)
            dzcf=z(i)-z(j)
            if(abs(dxcf).gt.halfboxx) dxcf=dxcf-sign(boxlenx,dxcf)
            if(abs(dycf).gt.halfboxy) dycf=dycf-sign(boxleny,dycf)
            if(abs(dzcf).gt.halfboxz) dzcf=dzcf-sign(boxlenz,dzcf)
            dxsav(numx)=hlab2(1,1)*dxcf+hlab2(1,2)*dycf+hlab2(1,3)*dzcf
            dysav(numx)=hlab2(2,1)*dxcf+hlab2(2,2)*dycf+hlab2(2,3)*dzcf
            dzsav(numx)=hlab2(3,1)*dxcf+hlab2(3,2)*dycf+hlab2(3,3)*dzcf
            drsq=dxsav(numx)*dxsav(numx)+dysav(numx)*dysav(numx)+dzsav(numx)*dzsav(numx)
            dr=dsqrt(drsq)
            erfc(numx)=erfunc(eta*dr)
         enddo
      enddo
   endif
   if( ooaimlog ) then
      if( cimlog .or. daimlog .or. quaimlog ) then
         do j=jst3,jed3
            do i=ist3(j),ied3(j)
               numx = numadr(i,j)
               dxcf=x(i)-x(j)
               dycf=y(i)-y(j)
               dzcf=z(i)-z(j)
               if(abs(dxcf).gt.halfboxx) dxcf=dxcf-sign(boxlenx,dxcf)
               if(abs(dycf).gt.halfboxy) dycf=dycf-sign(boxleny,dycf)
               if(abs(dzcf).gt.halfboxz) dzcf=dzcf-sign(boxlenz,dzcf)
               dxsav(numx)=hlab2(1,1)*dxcf+hlab2(1,2)*dycf+hlab2(1,3)*dzcf
               dysav(numx)=hlab2(2,1)*dxcf+hlab2(2,2)*dycf+hlab2(2,3)*dzcf
               dzsav(numx)=hlab2(3,1)*dxcf+hlab2(3,2)*dycf+hlab2(3,3)*dzcf
               drsq=dxsav(numx)*dxsav(numx)+dysav(numx)*dysav(numx)+dzsav(numx)*dzsav(numx)
               dr=dsqrt(drsq)
               erfc(numx)=erfunc(eta*dr)
            enddo
         enddo
      endif
   elseif( .not. epplog ) then
      do j=jst3,jed3
         do i=ist3(j),ied3(j)
            numx = numadr(i,j)
            dxcf=x(i)-x(j)
            dycf=y(i)-y(j)
            dzcf=z(i)-z(j)
            if(abs(dxcf).gt.halfboxx) dxcf=dxcf-sign(boxlenx,dxcf)
            if(abs(dycf).gt.halfboxy) dycf=dycf-sign(boxleny,dycf)
            if(abs(dzcf).gt.halfboxz) dzcf=dzcf-sign(boxlenz,dzcf)
            dxsav(numx)=hlab2(1,1)*dxcf+hlab2(1,2)*dycf+hlab2(1,3)*dzcf
            dysav(numx)=hlab2(2,1)*dxcf+hlab2(2,2)*dycf+hlab2(2,3)*dzcf
            dzsav(numx)=hlab2(3,1)*dxcf+hlab2(3,2)*dycf+hlab2(3,3)*dzcf
            drsq=dxsav(numx)*dxsav(numx)+dysav(numx)*dysav(numx)+dzsav(numx)*dzsav(numx)
            dr=dsqrt(drsq)
            erfc(numx)=erfunc(eta*dr)
         enddo
      enddo
   endif
   if( .not. epplog ) then
      do j=jst4,jed4
         do i=ist4(j),ied4(j)
            numx = numadr(i,j)
            dxcf=x(i)-x(j)
            dycf=y(i)-y(j)
            dzcf=z(i)-z(j)
            if(abs(dxcf).gt.halfboxx) dxcf=dxcf-sign(boxlenx,dxcf)
            if(abs(dycf).gt.halfboxy) dycf=dycf-sign(boxleny,dycf)
            if(abs(dzcf).gt.halfboxz) dzcf=dzcf-sign(boxlenz,dzcf)
            dxsav(numx)=hlab2(1,1)*dxcf+hlab2(1,2)*dycf+hlab2(1,3)*dzcf
            dysav(numx)=hlab2(2,1)*dxcf+hlab2(2,2)*dycf+hlab2(2,3)*dzcf
            dzsav(numx)=hlab2(3,1)*dxcf+hlab2(3,2)*dycf+hlab2(3,3)*dzcf
            drsq=dxsav(numx)*dxsav(numx)+dysav(numx)*dysav(numx)+dzsav(numx)*dzsav(numx)
            dr=dsqrt(drsq)
            erfc(numx)=erfunc(eta*dr)
         enddo
      enddo
   endif
endif
!<--- Parallelization_E

eltmp=0.d0
call cgfirst1
call cgfirst2

!---> Parallelization_S
!MPIAGUADO--> At this point eltmp has a different value for each node.
!             Here, each processor sends its eltmp information to all
!             other processors, so that in the end all processors have
!             the whole information about electric fields, etc
CALL MPI_ALLREDUCE(eltmp,elltmp,19*num,mpi_double_precision,MPI_SUM, &
                   MPI_COMM_WORLD,ierr)
!<--- Parallelization_E

elecxq=elltmp(1:num)
elecyq=elltmp(num+1:2*num)
eleczq=elltmp(2*num+1:3*num)
exxq=elltmp(3*num+1:4*num)
eyyq=elltmp(4*num+1:5*num)
ezzq=elltmp(5*num+1:6*num)
exyq=elltmp(6*num+1:7*num)
exzq=elltmp(7*num+1:8*num)
eyzq=elltmp(8*num+1:9*num)
engeff=elltmp(9*num+1:10*num)
elecxsr=elltmp(10*num+1:11*num)
elecysr=elltmp(11*num+1:12*num)
eleczsr=elltmp(12*num+1:13*num)
exxsr=elltmp(13*num+1:14*num)
eyysr=elltmp(14*num+1:15*num)
ezzsr=elltmp(15*num+1:16*num)
exysr=elltmp(16*num+1:17*num)
exzsr=elltmp(17*num+1:18*num)
eyzsr=elltmp(18*num+1:19*num)

elecxq=elecxq+elecxsr
elecyq=elecyq+elecysr
eleczq=eleczq+eleczsr

exxq=exxq+exxsr
eyyq=eyyq+eyysr
ezzq=ezzq+ezzsr
exyq=exyq+exysr
exzq=exzq+exzsr
eyzq=eyzq+eyzsr

eleczq=eleczq-ewcorr

!...Self-correction for the electric field gradient (AGUADO)
!...In this way the field gradient (induced by charges)
!...will be traceless
!........fac is 2 times fourpicell (two forces from differentiating double sum)
fac2=4.0d0*(eta**3.0d0)/(3.0d0*sqrpi)
exxq=exxq-fac2*q
eyyq=eyyq-fac2*q
ezzq=ezzq-fac2*q

engeff=engeff-etapi*q

if(firstiter) call cgener
do i=1,num
   ipoint=ntype(i)

   alppolar(i)=alppolar1(ipoint)*(dexp(alppolardel(ipoint)* &
           delta(i))+dexp(-alppolareff(ipoint)*engeff(i)))*0.5d0
   Bpolar(i)=Bpolar1(ipoint)*alppolar(i)
   Cpolar(i)=Cpolar1(ipoint)*alppolar(i)
   gammapolar(i)=gammapolar1(ipoint)*alppolar(i)
enddo

ftolaimtemp=ftolaim
ftoltemp=ftol
tolaim=1.0d-04
ftolaim=1.0d-06
ftol=1.0d-02

1     engtemp=engpetot
if(ftolaim.gt.ftolaimtemp)ftolaim=ftolaim*0.1d0
if(ftol.gt.ftoltemp)ftol=ftol*0.1d0
call frprmnaim(p,num*10)
if(quadpimlog.and.conjgradlog) call dipquadcg
if(dippimlog.and.conjgradlog) call dipcg
!     call cgsr_energy
!---> Parallelization_S
if( iam .eq. 0 ) then

if((.not.relaxconfig).or.verbose) write(6,*)'ENVPIM', abs(engtemp-engpetot),engtol

endif
!<--- Parallelization_E
if(abs(engpetot-engtemp).gt.engtol) goto 1
ftolaim=ftolaimtemp
ftol=ftoltemp

delta=p(1:num)
epsilonx=p(num+1:2*num)
epsilony=p(2*num+1:3*num)
epsilonz=p(3*num+1:4*num)
quaimxx=p(4*num+1:5*num)
quaimyy=p(5*num+1:6*num)
quaimzz=p(6*num+1:7*num)
quaimxy=p(7*num+1:8*num)
quaimxz=p(8*num+1:9*num)
quaimyz=p(9*num+1:10*num)

!---> Parallelization_S
if( iam .eq. 0 ) then

if(verbose)then
   do i=1,num
      if(environmentalaimlog)write(43,*)delta(i),selfeps(i),selfquaim(i)
      write(45,*)xmu(i),ymu(i),zmu(i)
      write(46,*)quadxx(i),quadyy(i),quadzz(i)
      write(47,*)quadxy(i),quadxz(i),quadyz(i)
      write(49,*)real(alppolar(i)),real(Bpolar(i)) &
                ,real(Cpolar(i)),real(gammapolar(i))
   enddo   
endif

if((.not.relaxconfig).or.verbose) then
 write(6,*)
 write(6,*)'**** Anion annealing completed ****'
 write(6,*)
end if

endif
!<--- Parallelization_E

return
END SUBROUTINE
