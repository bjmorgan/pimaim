SUBROUTINE dipquadcg

USE commondata, ONLY: num,xmu,ymu,zmu,quadxx,quadyy,quadzz,quadxy,quadxz, &
                      quadyz,alppolar, &
                      Bpolar,Cpolar,gammapolar, &
                      elecx,elecy,elecz,exx,eyy,ezz,exy,exz,eyz, &
                      ftol,xk1,xk2,xk3,xk4,polarizablelog,nspec,nsp,reseng, &
                      dipsqeng,dipquadeng,quadeng,delta,engpetot,verbose, &
                      environmentalaimlog,elecxq,elecyq,eleczq,exxq,eyyq,ezzq, &
                      exyq,exzq,eyzq,exit_code

!---> Parallelization_S
use mpipara
use mpi_spawn
!<--- Parallelization_E

IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(num) :: xmusv,ymusv,zmusv,  &
                    quadxxsv,quadyysv,quadzzsv,quadxysv,quadxzsv,quadyzsv, &
                    brhsx,brhsy,brhsz, &
                    brhsxx,brhsyy,brhszz,brhsxy,brhsxz,brhsyz
                    
DOUBLE PRECISION, DIMENSION(num) :: resx,resy,resz, &
                    resxx,resyy,reszz,resxy,resxz,resyz, &
                    pdirx,pdiry,pdirz, &
                    pdirxx,pdiryy,pdirzz,pdirxy,pdirxz,pdiryz
DOUBLE PRECISION, DIMENSION(num) :: Apx,Apy,Apz,Apxx,Apyy,Apzz,Apxy,Apxz,Apyz, &
                                    elecsq,thetatrace
DOUBLE PRECISION :: res2denom,res2num1,res2num2,cgalpha,cgbeta,pAp
DOUBLE PRECISION :: dipsq,quadsq
INTEGER :: i,j,m,n
INTEGER, PARAMETER :: itmax=2000
      
! the rhs of the equation to be solved are the 'charge induced' dipoles
! and quadrupoles.

elecsq=elecxq*elecxq+elecyq*elecyq+eleczq*eleczq
thetatrace=(exxq+eyyq+ezzq)/3.0d0

brhsx=elecxq*alppolar+Bpolar/3.0d0* &
        (elecxq*exxq+1.5d0*(elecyq*exyq+eleczq*exzq)- &
         0.5d0*elecxq*(eyyq+ezzq))+gammapolar/6.0d0*elecxq*elecsq
brhsy=elecyq*alppolar+Bpolar/3.0d0* &
        (elecyq*eyyq+1.5d0*(elecxq*exyq+eleczq*eyzq)- &
         0.5d0*elecyq*(exxq+ezzq))+gammapolar/6.0d0*elecyq*elecsq
brhsz=eleczq*alppolar+Bpolar/3.0d0* &
        (eleczq*ezzq+1.5d0*(elecxq*exzq+elecyq*eyzq)- &
         0.5d0*eleczq*(exxq+eyyq))+gammapolar/6.0d0*eleczq*elecsq
brhsxx=Cpolar*(exxq-thetatrace) &
         +Bpolar*(0.75d0*elecxq*elecxq-0.25d0*elecsq)
brhsyy=Cpolar*(eyyq-thetatrace) &
         +Bpolar*(0.75d0*elecyq*elecyq-0.25d0*elecsq)
brhszz=Cpolar*(ezzq-thetatrace) &
         +Bpolar*(0.75d0*eleczq*eleczq-0.25d0*elecsq)
brhsxy=Cpolar*exyq+0.75d0*Bpolar*elecxq*elecyq
brhsxz=Cpolar*exzq+0.75d0*Bpolar*elecxq*eleczq
brhsyz=Cpolar*eyzq+0.75d0*Bpolar*elecyq*eleczq

! first calculate the initial residuals
if(.not.environmentalaimlog) call cgener
elecsq=elecx*elecx+elecy*elecy+elecz*elecz
thetatrace=(exx+eyy+ezz)/3.0d0

resx=elecx*alppolar+Bpolar/3.0d0* &
       (elecx*exx+1.5d0*(elecy*exy+elecz*exz)- &
        0.5d0*elecx*(eyy+ezz))+gammapolar/6.0d0*elecx*elecsq-xmu
resy=elecy*alppolar+Bpolar/3.0d0* &
       (elecy*eyy+1.5d0*(elecx*exy+elecz*eyz)- &
        0.5d0*elecy*(exx+ezz))+gammapolar/6.0d0*elecy*elecsq-ymu
resz=elecz*alppolar+Bpolar/3.0d0* &
       (elecz*ezz+1.5d0*(elecx*exz+elecy*eyz)- &
        0.5d0*elecz*(exx+eyy))+gammapolar/6.0d0*elecz*elecsq-zmu
resxx=Cpolar*(exx-thetatrace)-quadxx &
        +Bpolar*(0.75d0*elecx*elecx-0.25d0*elecsq)
resyy=Cpolar*(eyy-thetatrace)-quadyy &
        +Bpolar*(0.75d0*elecy*elecy-0.25d0*elecsq)
reszz=Cpolar*(ezz-thetatrace)-quadzz &
        +Bpolar*(0.75d0*elecz*elecz-0.25d0*elecsq)
resxy=Cpolar*exy-quadxy+0.75d0*Bpolar*elecx*elecy
resxz=Cpolar*exz-quadxz+0.75d0*Bpolar*elecx*elecz
resyz=Cpolar*eyz-quadyz+0.75d0*Bpolar*elecy*elecz

do j=1,itmax
! find the new search direction
   res2denom=SUM(resx*resx+resy*resy+resz*resz)
   res2denom=res2denom+SUM(resxx*resxx+resyy*resyy+reszz*reszz+ &
                       resxy*resxy+resxz*resxz+resyz*resyz)
   if(j.eq.1) then
      pdirx=resx
      pdiry=resy
      pdirz=resz
      pdirxx=resxx
      pdiryy=resyy
      pdirzz=reszz
      pdirxy=resxy
      pdirxz=resxz
      pdiryz=resyz
   else
      pdirx=resx+cgbeta*pdirx
      pdiry=resy+cgbeta*pdiry
      pdirz=resz+cgbeta*pdirz
      pdirxx=resxx+cgbeta*pdirxx
      pdiryy=resyy+cgbeta*pdiryy
      pdirzz=reszz+cgbeta*pdirzz
      pdirxy=resxy+cgbeta*pdirxy
      pdirxz=resxz+cgbeta*pdirxz
      pdiryz=resyz+cgbeta*pdiryz
      if(cgbeta.gt.1.d0)then
!---> Parallelization_S
         if(iam.eq.0) write(6,*)'pim cg reset'
!<--- Parallelization_E
         pdirx=resx+0.5*cgbeta*pdirx
         pdiry=resy+0.5*cgbeta*pdiry
         pdirz=resz+0.5*cgbeta*pdirz
         pdirxx=resxx+0.5*cgbeta*pdirxx
         pdiryy=resyy+0.5*cgbeta*pdiryy
         pdirzz=reszz+0.5*cgbeta*pdirzz
         pdirxy=resxy+0.5*cgbeta*pdirxy
         pdirxz=resxz+0.5*cgbeta*pdirxz
         pdiryz=resyz+0.5*cgbeta*pdiryz
      endif
   endif

! find A.p      
   xmusv=xmu
   ymusv=ymu
   zmusv=zmu
   quadxxsv=quadxx
   quadyysv=quadyy
   quadzzsv=quadzz
   quadxysv=quadxy
   quadxzsv=quadxz
   quadyzsv=quadyz
   xmu=pdirx
   ymu=pdiry
   zmu=pdirz
   quadxx=pdirxx
   quadyy=pdiryy
   quadzz=pdirzz
   quadxy=pdirxy
   quadxz=pdirxz
   quadyz=pdiryz

   call cgener
      
   elecsq=elecx*elecx+elecy*elecy+elecz*elecz
   thetatrace=(exx+eyy+ezz)/3.0d0

   Apx=brhsx-(elecx*alppolar+Bpolar/3.0d0* &
       (elecx*exx+1.5d0*(elecy*exy+elecz*exz)- &
        0.5d0*elecx*(eyy+ezz))+gammapolar/6.0d0*elecx*elecsq-xmu)
   Apy=brhsy-(elecy*alppolar+Bpolar/3.0d0* &
       (elecy*eyy+1.5d0*(elecx*exy+elecz*eyz)- &
        0.5d0*elecy*(exx+ezz))+gammapolar/6.0d0*elecy*elecsq-ymu)
   Apz=brhsz-(elecz*alppolar+Bpolar/3.0d0* &
       (elecz*ezz+1.5d0*(elecx*exz+elecy*eyz)- &
        0.5d0*elecz*(exx+eyy))+gammapolar/6.0d0*elecz*elecsq-zmu)
   Apxx=brhsxx-(Cpolar*(exx-thetatrace)-quadxx &
        +Bpolar*(0.75d0*elecx*elecx-0.25d0*elecsq))
   Apyy=brhsyy-(Cpolar*(eyy-thetatrace)-quadyy &
        +Bpolar*(0.75d0*elecy*elecy-0.25d0*elecsq))
   Apzz=brhszz-(Cpolar*(ezz-thetatrace)-quadzz &
        +Bpolar*(0.75d0*elecz*elecz-0.25d0*elecsq))
   Apxy=brhsxy-(Cpolar*exy-quadxy+0.75d0*Bpolar*elecx*elecy)
   Apxz=brhsxz-(Cpolar*exz-quadxz+0.75d0*Bpolar*elecx*elecz)
   Apyz=brhsyz-(Cpolar*eyz-quadyz+0.75d0*Bpolar*elecy*elecz)

   xmu=xmusv
   ymu=ymusv
   zmu=zmusv
   quadxx=quadxxsv
   quadyy=quadyysv
   quadzz=quadzzsv
   quadxy=quadxysv
   quadxz=quadxzsv
   quadyz=quadyzsv
! now find alpha       
   pAp=SUM(pdirx*Apx+pdiry*Apy+pdirz*Apz)
   pAp=pAp+SUM(pdirxx*Apxx+pdiryy*Apyy+pdirzz*Apzz+ &
               pdirxy*Apxy+pdirxz*Apxz+pdiryz*Apyz)
   cgalpha=res2denom/pAp
! update dipoles and residuals
   xmu=xmu+cgalpha*pdirx
   ymu=ymu+cgalpha*pdiry
   zmu=zmu+cgalpha*pdirz
   quadxx=quadxx+cgalpha*pdirxx
   quadyy=quadyy+cgalpha*pdiryy
   quadzz=quadzz+cgalpha*pdirzz
   quadxy=quadxy+cgalpha*pdirxy
   quadxz=quadxz+cgalpha*pdirxz
   quadyz=quadyz+cgalpha*pdiryz
   call cgener

   res2num2=SUM((resx-cgalpha*Apx)*resx+(resy-cgalpha*Apy)*resy+&
                (resz-cgalpha*Apz)*resz)
   res2num2=res2num2+SUM((resxx-cgalpha*Apxx)*resxx &
                 +(resyy-cgalpha*Apyy)*resyy+(reszz-cgalpha*Apzz)*reszz &
                 +(resxy-cgalpha*Apxy)*resxy+(resxz-cgalpha*Apxz)*resxz &
                 +(resyz-cgalpha*Apyz)*resyz)

   resx=resx-cgalpha*Apx
   resy=resy-cgalpha*Apy
   resz=resz-cgalpha*Apz
   resxx=resxx-cgalpha*Apxx
   resyy=resyy-cgalpha*Apyy
   reszz=reszz-cgalpha*Apzz
   resxy=resxy-cgalpha*Apxy
   resxz=resxz-cgalpha*Apxz
   resyz=resyz-cgalpha*Apyz

   res2num1=SUM(resx*resx+resy*resy+resz*resz)
   res2num1=res2num1+SUM(resxx*resxx+resyy*resyy+reszz*reszz+ &
                         resxy*resxy+resxz*resxz+resyz*resyz)
   cgbeta=max((res2num1-res2num2)/res2denom,0.d0)
!---> Parallelization_S
   if(verbose .and. iam.eq.0) write(6,*)j,dsqrt(res2num1/num),ftol
!<--- Parallelization_E
! can we stop now?      
   if(dsqrt(res2num1/num).lt.ftol) then
      call cgener
! update of contributions to polengtot
      m=1
      do n=1,nspec
         if(polarizablelog(n)) then
            xk1(m:nsp(n)+m-1)=1.0d0/(2.0d0*alppolar(m:nsp(n)+m-1))
            xk3(m:nsp(n)+m-1)=1.0d0/(6.0d0*Cpolar(m:nsp(n)+m-1))
            xk2(m:nsp(n)+m-1)=-Bpolar(m:nsp(n)+m-1)/(4.0d0*alppolar(m:nsp(n)+m-1) &
                              *alppolar(m:nsp(n)+m-1)*Cpolar(m:nsp(n)+m-1))
            xk4(m:nsp(n)+m-1)=(2.0d0*gammapolar(m:nsp(n)+m-1)*Cpolar(m:nsp(n)+m-1)-3.0d0*Bpolar(m:nsp(n)+m-1)* &
                               Bpolar(m:nsp(n)+m-1))/(48.0d0*Cpolar(m:nsp(n)+m-1)* &
                               alppolar(m:nsp(n)+m-1)*alppolar(m:nsp(n)+m-1)* &
                               alppolar(m:nsp(n)+m-1)*alppolar(m:nsp(n)+m-1))
    
            do i=m,nsp(n)+m-1
               dipsq=(xmu(i)*xmu(i))+(ymu(i)*ymu(i))+(zmu(i)*zmu(i))
               reseng(i)=dipsq*xk1(i)
               dipsqeng(i)=dipsq*dipsq*xk4(i)
               dipquadeng(i)=(xmu(i)*quadxx(i)*xmu(i) &
                             +ymu(i)*quadyy(i)*ymu(i)+zmu(i)*quadzz(i)*zmu(i) &
                 +2.0d0*(xmu(i)*quadxy(i)*ymu(i)+xmu(i)*quadxz(i)*zmu(i) &
                        +ymu(i)*quadyz(i)*zmu(i)))*xk2(i)
               quadsq=quadxx(i)*quadxx(i)+quadyy(i)*quadyy(i) &
                     +quadzz(i)*quadzz(i)+2.0d0*(quadxy(i)*quadxy(i) &
                     +quadxz(i)*quadxz(i)+quadyz(i)*quadyz(i))
               quadeng(i)=quadsq*xk3(i)
            enddo      
         endif
         m=m+nsp(n)
      enddo
      return
   endif
enddo
!---> Parallelization_S

exit_code  = 1
#ifndef ppfit
if(iam.eq.0) write(6,*)'cg failed to converge'
call close_down_mpi()
!<--- Parallelization_E
STOP
#endif
RETURN
END SUBROUTINE
