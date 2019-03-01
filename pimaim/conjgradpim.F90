SUBROUTINE conjgrad

use error_function 
USE commondata, ONLY: num,x,y,z,dxsav,dysav,dzsav,erfc,eta,dippimlog, &
                      quadpimlog,xmu,ymu,zmu,quadxx,quadyy,quadzz,quadxy, &
                      quadxz,quadyz,alppolar,Bpolar,Cpolar,gammapolar,verbose, &
                      engeff,alppolar1,Bpolar1,Cpolar1,gammapolar1, &
!---> Memmory Reduction_S
!                     alppolardel,alppolareff,ntype,delta
                      alppolardel,alppolareff,ntype,delta,num2,numx, &
                      cimlog,daimlog,quaimlog,ooaimlog,epplog, &
                      elecxq,elecyq,eleczq,exxq,eyyq,ezzq,exyq,exzq,eyzq, &
                      engeff,etapi,q,sqrpi,ewcorr,elecxsr,elecysr,eleczsr, &
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
DOUBLE PRECISION :: dxcf,dycf,dzcf,drsq,dr
INTEGER :: i,l,j,ipoint
!
! Setup cos and sin arrays for addition rule.
!
 elcall(:,0)=1.0d0
 emcall(:,0)=1.0d0
 encall(:,0)=1.0d0
 elsall(:,0)=0.0d0
 emsall(:,0)=0.0d0
 ensall(:,0)=0.0d0

 elcall(:,1)=cos(twopiboxx*x)
 emcall(:,1)=cos(twopiboxy*y)
 encall(:,1)=cos(twopiboxz*z)
 elsall(:,1)=sin(twopiboxx*x)
 emsall(:,1)=sin(twopiboxy*y)
 ensall(:,1)=sin(twopiboxz*z)
!
! Calculate all cosines and sines
!
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

!      write(911,'(3(F16.8,2X))')hlab2(1,1),hlab2(1,2),hlab2(1,3)
!      write(911,'(3(F16.8,2X))')hlab2(2,1),hlab2(2,2),hlab2(2,3)
!      write(911,'(3(F16.8,2X))')hlab2(3,1),hlab2(3,2),hlab2(3,3)
do j=jst,jed
   do i=ist(j),ied(j)
      numx = numadr(i,j)
!<--- Parallelization_E
      dxcf=x(i)-x(j)
      dycf=y(i)-y(j)
      dzcf=z(i)-z(j)
!   write(911,'("atom 1",i10,4(X,F12.6))') i,x(i),y(i),z(i)
!   write(911,'("atom 2",i10,4(X,F12.6))') j,x(j),y(j),z(j)
!   write(911,'("dx bef",i10,4(X,F12.6))') numx,dxcf/boxlenx,dycf/boxleny,dzcf/boxlenz 
      if(abs(dxcf).gt.halfboxx) dxcf=dxcf-sign(boxlenx,dxcf)
      if(abs(dycf).gt.halfboxy) dycf=dycf-sign(boxleny,dycf)
      if(abs(dzcf).gt.halfboxz) dzcf=dzcf-sign(boxlenz,dzcf)
!   write(911,'("dx aft",i10,4(X,F12.6))') numx,dxcf/boxlenx,dycf/boxleny,dzcf/boxlenz 
!---> Memmory Reduction_S
!     dxsav(i,j)=hlab2(1,1)*dxcf+hlab2(1,2)*dycf+hlab2(1,3)*dzcf
!     dysav(i,j)=hlab2(2,1)*dxcf+hlab2(2,2)*dycf+hlab2(2,3)*dzcf
!     dzsav(i,j)=hlab2(3,1)*dxcf+hlab2(3,2)*dycf+hlab2(3,3)*dzcf
      dxsav(numx)=hlab2(1,1)*dxcf+hlab2(1,2)*dycf+hlab2(1,3)*dzcf
      dysav(numx)=hlab2(2,1)*dxcf+hlab2(2,2)*dycf+hlab2(2,3)*dzcf
      dzsav(numx)=hlab2(3,1)*dxcf+hlab2(3,2)*dycf+hlab2(3,3)*dzcf

!     drsq=dxsav(i,j)*dxsav(i,j)+dysav(i,j)*dysav(i,j)+dzsav(i,j)*dzsav(i,j)
      drsq=dxsav(numx)*dxsav(numx)+dysav(numx)*dysav(numx)+dzsav(numx)*dzsav(numx)
!<--- Memmory Reduction_E
      dr=dsqrt(drsq)
!---> Memmory Reduction_S
!     erfc(i,j)=erfunc(eta*dr)
      erfc(numx)=erfunc(eta*dr)
!<--- Memmory Reduction_E
!   write(911,'(i10,4(X,F12.6))') numx,dxsav(numx),dysav(numx),dzsav(numx), dr 
   enddo   
enddo   

!---> Parallelization_S
!MPIAGUADO--> This is done separately for the sr_energy part, which is
!             evaluated differently for cation-anion, anion-anion and
!             cation-cation interactions.
!             Nevertheless, this routine is only called in pim runs, so if
!             we would ever want to do an AIM calculation with non polarizable ions
!             (not very probable, I know) this might not work!
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
CALL MPI_ALLREDUCE(eltmp,elltmp,19*num,mpi_double_precision,MPI_SUM,MPI_COMM_WORLD,ierr)
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
alppolar=0.0d0
do i=1,num
   ipoint=ntype(i)

   alppolar(i)=alppolar1(ipoint)*(dexp(alppolardel(ipoint)* &
           delta(i))+dexp(-alppolareff(ipoint)*engeff(i)))*0.5d0
   Bpolar(i)=Bpolar1(ipoint)*alppolar(i)
   Cpolar(i)=Cpolar1(ipoint)*alppolar(i)
   gammapolar(i)=gammapolar1(ipoint)*alppolar(i)
enddo

if(quadpimlog) call dipquadcg
if(dippimlog) call dipcg

!---> Parallelization_S
if( iam .eq. 0 ) then

if(verbose)then
   do i=1,num
      write(45,*)xmu(i),ymu(i),zmu(i)
      write(46,*)quadxx(i),quadyy(i),quadzz(i)
      write(47,*)quadxy(i),quadxz(i),quadyz(i)
      write(49,*)real(alppolar(i)),real(Bpolar(i)) &
                ,real(Cpolar(i)),real(gammapolar(i))
   enddo   
endif

!write(*,*)
!write(*,*)'**** Anion annealing completed ****'
!write(*,*)

endif
!<--- Parallelization_E
return
END SUBROUTINE
