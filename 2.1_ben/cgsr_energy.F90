SUBROUTINE cgsr_energy

USE commondata, ONLY: num,nspec,engft1,engft2,engft3,engft1dotx,engft1doty, &
                      engft1dotz,engft2dotx,engft2doty,engft2dotz, &
                      engft3dotx,engft3doty,engft3dotz,engft1dotxx, &
                      engft1dotyy,engft1dotzz,engft1dotxy,engft1dotxz, &
                      engft1dotyz,engft2dotxx,engft2dotyy,engft2dotzz, &
                      engft2dotxy,engft2dotxz,engft2dotyz,engft3dotxx, &
                      engft3dotyy,engft3dotzz,engft3dotxy,engft3dotxz, &
                      engft3dotyz,cimlog,daimlog,quaimlog,nsp,ntype,dxsav, &
                      dysav,dzsav,rsqmaxsr,ftalp,ftbeta,ftgamma,ftb,ftb2,ftb3, &
                      delta,epsilonx,epsilony,epsilonz,quaimxx,quaimyy, &
                      quaimzz,quaimxy,quaimxz,quaimyz,ooaimlog,nanion,engpetot, &
!---> Memmory Reduction_S
!---> Parallelization_S
!                     deformablelog,extraalpha,extrab
                      deformablelog,extraalpha,extrab,num2,numx
!<--- Parallelization_E
!<--- Memmory Reduction_E
!---> Optimize_S
use mpipara
!<--- Optimize_E

IMPLICIT NONE

!---> Parallelization_S
!INTEGER :: i,j,ipoint,jpoint
INTEGER :: i,j,ipoint,jpoint,kk
!<--- Parallelization_E
DOUBLE PRECISION :: dr,drrec,drsq,drsqrec,engsr,expft, &
                 expftx,expft2,expft3,expft4,dx,dy,dz, &
                 dxx,dyy,dzz,dxy,dxz,dyz,txx,tyy,tzz,txy,txz,tyz, &
                 dxunit,dyunit,dzunit,repsilon,rrtheta, &
                 repsiloni,rrthetai,repsilonj,rrthetaj,repsiloncat,rrthetacat
DOUBLE PRECISION :: engtmp(30*num,nspec),enggtmp(30*num,nspec)

!---> Parallelization_S
sctmp = 0.0d0
engtmp=0.d0
!<--- Parallelization_E

if(cimlog) then
!---> Parallelization_S
!  do j=nsp(1)+1,num
!MPIAGUADO-->Note that now we use jst2,jed2 because here only the
!            cation-anion interactions are evaluated...
   do j = jst2, jed2
!<--- Parallelization_E
      jpoint=ntype(j)
!---> Parallelization_S
!     do i=1,nsp(1)
      do i = ist2(j), ied2(j)
!<--- Parallelization_E
         ipoint=ntype(i)

!---> Memmory Reduction_S
         numx = numadr(i,j)
!        drsq=dxsav(i,j)*dxsav(i,j)+dysav(i,j)*dysav(i,j)+dzsav(i,j)*dzsav(i,j) 
         drsq=dxsav(numx)*dxsav(numx)+dysav(numx)*dysav(numx)+dzsav(numx)*dzsav(numx) 
!<--- Memmory Reduction_E

         if(drsq.ge.rsqmaxsr) CYCLE

         dr=dsqrt(drsq)

         expft=(ftb(ipoint,jpoint))*exp(-ftalp(ipoint,jpoint)* &
                                       (dr-delta(i)-delta(j)))
         expft2=(ftb2(ipoint,jpoint))*exp(-ftbeta(ipoint,jpoint)* &
                                       (dr-delta(i)-delta(j)))
         expft3=(ftb3(ipoint,jpoint))*exp(-ftgamma(ipoint,jpoint)* &
                                       (dr-delta(i)-delta(j)))
         expft4=extrab*exp(-extraalpha*dr)

!---> Parallelization_S
!        engft1(i,jpoint)=engft1(i,jpoint)+expft
!        engft2(i,jpoint)=engft2(i,jpoint)+expft2
!        engft3(i,jpoint)=engft3(i,jpoint)+expft3
         engtmp(i,jpoint)=engtmp(i,jpoint)+expft
         engtmp(num+i,jpoint)=engtmp(num+i,jpoint)+expft2
         engtmp(2*num+i,jpoint)=engtmp(2*num+i,jpoint)+expft3
         if(deformablelog(jpoint))then
!           engft1(j,ipoint)=engft1(j,ipoint)+expft
!           engft2(j,ipoint)=engft2(j,ipoint)+expft2
!           engft3(j,ipoint)=engft3(j,ipoint)+expft3
            engtmp(j,ipoint)=engtmp(j,ipoint)+expft
            engtmp(num+j,ipoint)=engtmp(num+j,ipoint)+expft2
            engtmp(2*num+j,ipoint)=engtmp(2*num+j,ipoint)+expft3
         endif
!        engsr=engsr+expft+expft2+expft3+expft4
         sctmp(1)=sctmp(1)+expft+expft2+expft3+expft4
!<--- Parallelization_E
      enddo   
   enddo   
endif

if(daimlog) then
!---> Parallelization_S
!  do j=nsp(1)+1,num
   do j = jst2, jed2
!<--- Parallelization_E
      jpoint=ntype(j)
!---> Parallelization_S
!     do i=1,nsp(1)
      do i = ist2(j), ied2(j)
!<--- Parallelization_E
         ipoint=ntype(i)

!---> Memmory Reduction_S
         numx = numadr(i,j)
!        drsq=dxsav(i,j)*dxsav(i,j)+dysav(i,j)*dysav(i,j)+dzsav(i,j)*dzsav(i,j)
         drsq=dxsav(numx)*dxsav(numx)+dysav(numx)*dysav(numx)+dzsav(numx)*dzsav(numx)
!<--- Memmory Reduction_E

         if(drsq.ge.rsqmaxsr) CYCLE

         dr=dsqrt(drsq)
         drrec=1.0d0/dr
!---> Memmory Reduction_S
!        dx=dxsav(i,j)
!        dy=dysav(i,j)
!        dz=dzsav(i,j)
         dx=dxsav(numx)
         dy=dysav(numx)
         dz=dzsav(numx)
!<--- Memmory Reduction_E

         repsilon=epsilonx(i)*dx+epsilony(i)*dy+epsilonz(i)*dz
         repsiloncat=epsilonx(j)*dx+epsilony(j)*dy+epsilonz(j)*dz

         repsilon=repsilon*drrec
         repsiloncat=repsiloncat*drrec

         expft=ftb(ipoint,jpoint)*exp(-ftalp(ipoint,jpoint)* &
                    (dr-delta(i)-delta(j)-repsilon+repsiloncat))
         expft2=ftb2(ipoint,jpoint)*exp(-ftbeta(ipoint,jpoint)* &
                    (dr-delta(i)-delta(j)-repsilon+repsiloncat))
         expft3=ftb3(ipoint,jpoint)*exp(-ftgamma(ipoint,jpoint)* &
                    (dr-delta(i)-delta(j)-repsilon+repsiloncat))
         expft4=extrab*exp(-extraalpha*dr)

!---> Parallelization_S
!        engft1(i,jpoint)=engft1(i,jpoint)+expft
!        engft2(i,jpoint)=engft2(i,jpoint)+expft2
!        engft3(i,jpoint)=engft3(i,jpoint)+expft3
         engtmp(i,jpoint)=engtmp(i,jpoint)+expft
         engtmp(num+i,jpoint)=engtmp(num+i,jpoint)+expft2
         engtmp(2*num+i,jpoint)=engtmp(2*num+i,jpoint)+expft3

         dxunit=dx*drrec
         dyunit=dy*drrec
         dzunit=dz*drrec

!        engft1dotx(i,jpoint)=engft1dotx(i,jpoint)+expft*dxunit
!        engft1doty(i,jpoint)=engft1doty(i,jpoint)+expft*dyunit
!        engft1dotz(i,jpoint)=engft1dotz(i,jpoint)+expft*dzunit
!        engft2dotx(i,jpoint)=engft2dotx(i,jpoint)+expft2*dxunit
!        engft2doty(i,jpoint)=engft2doty(i,jpoint)+expft2*dyunit
!        engft2dotz(i,jpoint)=engft2dotz(i,jpoint)+expft2*dzunit
!        engft3dotx(i,jpoint)=engft3dotx(i,jpoint)+expft3*dxunit
!        engft3doty(i,jpoint)=engft3doty(i,jpoint)+expft3*dyunit
!        engft3dotz(i,jpoint)=engft3dotz(i,jpoint)+expft3*dzunit
         engtmp(3*num+i,jpoint)=engtmp(3*num+i,jpoint)+expft*dxunit
         engtmp(4*num+i,jpoint)=engtmp(4*num+i,jpoint)+expft*dyunit
         engtmp(5*num+i,jpoint)=engtmp(5*num+i,jpoint)+expft*dzunit
         engtmp(6*num+i,jpoint)=engtmp(6*num+i,jpoint)+expft2*dxunit
         engtmp(7*num+i,jpoint)=engtmp(7*num+i,jpoint)+expft2*dyunit
         engtmp(8*num+i,jpoint)=engtmp(8*num+i,jpoint)+expft2*dzunit
         engtmp(9*num+i,jpoint)=engtmp(9*num+i,jpoint)+expft3*dxunit
         engtmp(10*num+i,jpoint)=engtmp(10*num+i,jpoint)+expft3*dyunit
         engtmp(11*num+i,jpoint)=engtmp(11*num+i,jpoint)+expft3*dzunit

         if(deformablelog(jpoint))then
!           engft1(j,ipoint)=engft1(j,ipoint)+expft
!           engft2(j,ipoint)=engft2(j,ipoint)+expft2
!           engft3(j,ipoint)=engft3(j,ipoint)+expft3
!           engft1dotx(j,ipoint)=engft1dotx(j,ipoint)-expft*dxunit
!           engft1doty(j,ipoint)=engft1doty(j,ipoint)-expft*dyunit
!           engft1dotz(j,ipoint)=engft1dotz(j,ipoint)-expft*dzunit
!           engft2dotx(j,ipoint)=engft2dotx(j,ipoint)-expft2*dxunit
!           engft2doty(j,ipoint)=engft2doty(j,ipoint)-expft2*dyunit
!           engft2dotz(j,ipoint)=engft2dotz(j,ipoint)-expft2*dzunit
!           engft3dotx(j,ipoint)=engft3dotx(j,ipoint)-expft3*dxunit
!           engft3doty(j,ipoint)=engft3doty(j,ipoint)-expft3*dyunit
!           engft3dotz(j,ipoint)=engft3dotz(j,ipoint)-expft3*dzunit
            engtmp(j,ipoint)=engtmp(j,ipoint)+expft
            engtmp(num+j,ipoint)=engtmp(num+j,ipoint)+expft2
            engtmp(2*num+j,ipoint)=engtmp(2*num+j,ipoint)+expft3
            engtmp(3*num+j,ipoint)=engtmp(3*num+j,ipoint)-expft*dxunit
            engtmp(4*num+j,ipoint)=engtmp(4*num+j,ipoint)-expft*dyunit
            engtmp(5*num+j,ipoint)=engtmp(5*num+j,ipoint)-expft*dzunit
            engtmp(6*num+j,ipoint)=engtmp(6*num+j,ipoint)-expft2*dxunit
            engtmp(7*num+j,ipoint)=engtmp(7*num+j,ipoint)-expft2*dyunit
            engtmp(8*num+j,ipoint)=engtmp(8*num+j,ipoint)-expft2*dzunit
            engtmp(9*num+j,ipoint)=engtmp(9*num+j,ipoint)-expft3*dxunit
            engtmp(10*num+j,ipoint)=engtmp(10*num+j,ipoint)-expft3*dyunit
            engtmp(11*num+j,ipoint)=engtmp(11*num+j,ipoint)-expft3*dzunit
         endif
!        engsr=engsr+expft+expft2+expft3+expft4
         sctmp(1)=sctmp(1)+expft+expft2+expft3+expft4
!<--- Parallelization_E
      enddo   
   enddo   
endif

if(quaimlog) then
!---> Parallelization_S
!  do j=nsp(1)+1,num
   do j = jst2, jed2
!<--- Parallelization_E
      jpoint=ntype(j)
!---> Parallelization_S
!     do i=1,nsp(1)
      do i = ist2(j), ied2(j)
!<--- Parallelization_E
         ipoint=ntype(i)

!---> Memmory Reduction_S
         numx = numadr(i,j)
!        drsq=dxsav(i,j)*dxsav(i,j)+dysav(i,j)*dysav(i,j)+dzsav(i,j)*dzsav(i,j)
         drsq=dxsav(numx)*dxsav(numx)+dysav(numx)*dysav(numx)+dzsav(numx)*dzsav(numx)
!<--- Memmory Reduction_E

         if(drsq.ge.rsqmaxsr) CYCLE

         drsqrec=1.0d0/drsq
         dr=dsqrt(drsq)
         drrec=1.0d0/dr
!---> Memmory Reduction_S
!        dx=dxsav(i,j)
!        dy=dysav(i,j)
!        dz=dzsav(i,j)
         dx=dxsav(numx)
         dy=dysav(numx)
         dz=dzsav(numx)
!<--- Memmory Reduction_E

         repsilon=epsilonx(i)*dx+epsilony(i)*dy+epsilonz(i)*dz
         repsiloncat=epsilonx(j)*dx+epsilony(j)*dy+epsilonz(j)*dz
         repsilon=repsilon*drrec
         repsiloncat=repsiloncat*drrec

         dxx=dx*dx
         dyy=dy*dy
         dzz=dz*dz
         dxy=dx*dy
         dxz=dx*dz
         dyz=dy*dz
         txx=3.0d0*dxx*drsqrec-1.0d0
         tyy=3.0d0*dyy*drsqrec-1.0d0
         tzz=3.0d0*dzz*drsqrec-1.0d0
         txy=3.0d0*dxy*drsqrec
         txz=3.0d0*dxz*drsqrec
         tyz=3.0d0*dyz*drsqrec

         rrtheta=quaimxx(i)*txx+quaimyy(i)*tyy+quaimzz(i)*tzz+ &
              2.0d0*(quaimxy(i)*txy+quaimxz(i)*txz+quaimyz(i)*tyz)
         rrthetacat=quaimxx(j)*txx+quaimyy(j)*tyy+ &
                    quaimzz(j)*tzz+2.0d0*(quaimxy(j)*txy+ &
                    quaimxz(j)*txz+quaimyz(j)*tyz)

         expft=ftb(ipoint,jpoint)*exp(-ftalp(ipoint,jpoint)* &
                (dr-delta(i)-delta(j)-repsilon+repsiloncat-rrtheta-rrthetacat))
         expft2=ftb2(ipoint,jpoint)*exp(-ftbeta(ipoint,jpoint)* &
                (dr-delta(i)-delta(j)-repsilon+repsiloncat-rrtheta-rrthetacat))
         expft3=ftb3(ipoint,jpoint)*exp(-ftgamma(ipoint,jpoint)* &
                (dr-delta(i)-delta(j)-repsilon+repsiloncat-rrtheta-rrthetacat))
         expft4=extrab*exp(-extraalpha*dr)

!---> Parallelization_S
!        engft1(i,jpoint)=engft1(i,jpoint)+expft
!        engft2(i,jpoint)=engft2(i,jpoint)+expft2
!        engft3(i,jpoint)=engft3(i,jpoint)+expft3
         engtmp(i,jpoint)=engtmp(i,jpoint)+expft
         engtmp(num+i,jpoint)=engtmp(num+i,jpoint)+expft2
         engtmp(2*num+i,jpoint)=engtmp(2*num+i,jpoint)+expft3

         dxunit=dx*drrec
         dyunit=dy*drrec
         dzunit=dz*drrec

!        engft1dotx(i,jpoint)=engft1dotx(i,jpoint)+expft*dxunit
!        engft1doty(i,jpoint)=engft1doty(i,jpoint)+expft*dyunit
!        engft1dotz(i,jpoint)=engft1dotz(i,jpoint)+expft*dzunit
!        engft2dotx(i,jpoint)=engft2dotx(i,jpoint)+expft2*dxunit
!        engft2doty(i,jpoint)=engft2doty(i,jpoint)+expft2*dyunit
!        engft2dotz(i,jpoint)=engft2dotz(i,jpoint)+expft2*dzunit
!        engft3dotx(i,jpoint)=engft3dotx(i,jpoint)+expft3*dxunit
!        engft3doty(i,jpoint)=engft3doty(i,jpoint)+expft3*dyunit
!        engft3dotz(i,jpoint)=engft3dotz(i,jpoint)+expft3*dzunit
         engtmp(3*num+i,jpoint)=engtmp(3*num+i,jpoint)+expft*dxunit
         engtmp(4*num+i,jpoint)=engtmp(4*num+i,jpoint)+expft*dyunit
         engtmp(5*num+i,jpoint)=engtmp(5*num+i,jpoint)+expft*dzunit
         engtmp(6*num+i,jpoint)=engtmp(6*num+i,jpoint)+expft2*dxunit
         engtmp(7*num+i,jpoint)=engtmp(7*num+i,jpoint)+expft2*dyunit
         engtmp(8*num+i,jpoint)=engtmp(8*num+i,jpoint)+expft2*dzunit
         engtmp(9*num+i,jpoint)=engtmp(9*num+i,jpoint)+expft3*dxunit
         engtmp(10*num+i,jpoint)=engtmp(10*num+i,jpoint)+expft3*dyunit
         engtmp(11*num+i,jpoint)=engtmp(11*num+i,jpoint)+expft3*dzunit

!        engft1dotxx(i,jpoint)=engft1dotxx(i,jpoint)+expft*txx
!        engft1dotyy(i,jpoint)=engft1dotyy(i,jpoint)+expft*tyy
!        engft1dotzz(i,jpoint)=engft1dotzz(i,jpoint)+expft*tzz
!        engft1dotxy(i,jpoint)=engft1dotxy(i,jpoint)+2.0d0*expft*txy
!        engft1dotxz(i,jpoint)=engft1dotxz(i,jpoint)+2.0d0*expft*txz
!        engft1dotyz(i,jpoint)=engft1dotyz(i,jpoint)+2.0d0*expft*tyz
!        engft2dotxx(i,jpoint)=engft2dotxx(i,jpoint)+expft2*txx
!        engft2dotyy(i,jpoint)=engft2dotyy(i,jpoint)+expft2*tyy
!        engft2dotzz(i,jpoint)=engft2dotzz(i,jpoint)+expft2*tzz
!        engft2dotxy(i,jpoint)=engft2dotxy(i,jpoint)+2.0d0*expft2*txy
!        engft2dotxz(i,jpoint)=engft2dotxz(i,jpoint)+2.0d0*expft2*txz
!        engft2dotyz(i,jpoint)=engft2dotyz(i,jpoint)+2.0d0*expft2*tyz
!        engft3dotxx(i,jpoint)=engft3dotxx(i,jpoint)+expft3*txx
!        engft3dotyy(i,jpoint)=engft3dotyy(i,jpoint)+expft3*tyy
!        engft3dotzz(i,jpoint)=engft3dotzz(i,jpoint)+expft3*tzz
!        engft3dotxy(i,jpoint)=engft3dotxy(i,jpoint)+2.0d0*expft3*txy
!        engft3dotxz(i,jpoint)=engft3dotxz(i,jpoint)+2.0d0*expft3*txz
!        engft3dotyz(i,jpoint)=engft3dotyz(i,jpoint)+2.0d0*expft3*tyz
         engtmp(12*num+i,jpoint)=engtmp(12*num+i,jpoint)+expft*txx
         engtmp(13*num+i,jpoint)=engtmp(13*num+i,jpoint)+expft*tyy
         engtmp(14*num+i,jpoint)=engtmp(14*num+i,jpoint)+expft*tzz
         engtmp(15*num+i,jpoint)=engtmp(15*num+i,jpoint)+2.0d0*expft*txy
         engtmp(16*num+i,jpoint)=engtmp(16*num+i,jpoint)+2.0d0*expft*txz
         engtmp(17*num+i,jpoint)=engtmp(17*num+i,jpoint)+2.0d0*expft*tyz
         engtmp(18*num+i,jpoint)=engtmp(18*num+i,jpoint)+expft2*txx
         engtmp(19*num+i,jpoint)=engtmp(19*num+i,jpoint)+expft2*tyy
         engtmp(20*num+i,jpoint)=engtmp(20*num+i,jpoint)+expft2*tzz
         engtmp(21*num+i,jpoint)=engtmp(21*num+i,jpoint)+2.0d0*expft2*txy
         engtmp(22*num+i,jpoint)=engtmp(22*num+i,jpoint)+2.0d0*expft2*txz
         engtmp(23*num+i,jpoint)=engtmp(23*num+i,jpoint)+2.0d0*expft2*tyz
         engtmp(24*num+i,jpoint)=engtmp(24*num+i,jpoint)+expft3*txx
         engtmp(25*num+i,jpoint)=engtmp(25*num+i,jpoint)+expft3*tyy
         engtmp(26*num+i,jpoint)=engtmp(26*num+i,jpoint)+expft3*tzz
         engtmp(27*num+i,jpoint)=engtmp(27*num+i,jpoint)+2.0d0*expft3*txy
         engtmp(28*num+i,jpoint)=engtmp(28*num+i,jpoint)+2.0d0*expft3*txz
         engtmp(29*num+i,jpoint)=engtmp(29*num+i,jpoint)+2.0d0*expft3*tyz

         if(deformablelog(jpoint))then
!           engft1(j,ipoint)=engft1(j,ipoint)+expft
!           engft2(j,ipoint)=engft2(j,ipoint)+expft2
!           engft3(j,ipoint)=engft3(j,ipoint)+expft3
!           engft1dotx(j,ipoint)=engft1dotx(j,ipoint)-expft*dxunit
!           engft1doty(j,ipoint)=engft1doty(j,ipoint)-expft*dyunit
!           engft1dotz(j,ipoint)=engft1dotz(j,ipoint)-expft*dzunit
!           engft2dotx(j,ipoint)=engft2dotx(j,ipoint)-expft2*dxunit
!           engft2doty(j,ipoint)=engft2doty(j,ipoint)-expft2*dyunit
!           engft2dotz(j,ipoint)=engft2dotz(j,ipoint)-expft2*dzunit
!           engft3dotx(j,ipoint)=engft3dotx(j,ipoint)-expft3*dxunit
!           engft3doty(j,ipoint)=engft3doty(j,ipoint)-expft3*dyunit
!           engft3dotz(j,ipoint)=engft3dotz(j,ipoint)-expft3*dzunit
!           engft1dotxx(j,ipoint)=engft1dotxx(j,ipoint)+expft*txx
!           engft1dotyy(j,ipoint)=engft1dotyy(j,ipoint)+expft*tyy
!           engft1dotzz(j,ipoint)=engft1dotzz(j,ipoint)+expft*tzz
!           engft1dotxy(j,ipoint)=engft1dotxy(j,ipoint)+2.0d0*expft*txy
!           engft1dotxz(j,ipoint)=engft1dotxz(j,ipoint)+2.0d0*expft*txz
!           engft1dotyz(j,ipoint)=engft1dotyz(j,ipoint)+2.0d0*expft*tyz
!           engft2dotxx(j,ipoint)=engft2dotxx(j,ipoint)+expft2*txx
!           engft2dotyy(j,ipoint)=engft2dotyy(j,ipoint)+expft2*tyy
!           engft2dotzz(j,ipoint)=engft2dotzz(j,ipoint)+expft2*tzz
!           engft2dotxy(j,ipoint)=engft2dotxy(j,ipoint)+2.0d0*expft2*txy
!           engft2dotxz(j,ipoint)=engft2dotxz(j,ipoint)+2.0d0*expft2*txz
!           engft2dotyz(j,ipoint)=engft2dotyz(j,ipoint)+2.0d0*expft2*tyz
!           engft3dotxx(j,ipoint)=engft3dotxx(j,ipoint)+expft3*txx
!           engft3dotyy(j,ipoint)=engft3dotyy(j,ipoint)+expft3*tyy
!           engft3dotzz(j,ipoint)=engft3dotzz(j,ipoint)+expft3*tzz
!           engft3dotxy(j,ipoint)=engft3dotxy(j,ipoint)+2.0d0*expft3*txy
!           engft3dotxz(j,ipoint)=engft3dotxz(j,ipoint)+2.0d0*expft3*txz
!           engft3dotyz(j,ipoint)=engft3dotyz(j,ipoint)+2.0d0*expft3*tyz
            engtmp(j,ipoint) =engtmp(j,ipoint) +expft
            engtmp(num+j,ipoint) =engtmp(num+j,ipoint) +expft2
            engtmp(2*num+j,ipoint) =engtmp(2*num+j,ipoint) +expft3
            engtmp(3*num+j,ipoint)=engtmp(3*num+j,ipoint)-expft*dxunit
            engtmp(4*num+j,ipoint)=engtmp(4*num+j,ipoint)-expft*dyunit
            engtmp(5*num+j,ipoint)=engtmp(5*num+j,ipoint)-expft*dzunit
            engtmp(6*num+j,ipoint)=engtmp(6*num+j,ipoint)-expft2*dxunit
            engtmp(7*num+j,ipoint)=engtmp(7*num+j,ipoint)-expft2*dyunit
            engtmp(8*num+j,ipoint)=engtmp(8*num+j,ipoint)-expft2*dzunit
            engtmp(9*num+j,ipoint)=engtmp(9*num+j,ipoint)-expft3*dxunit
            engtmp(10*num+j,ipoint)=engtmp(10*num+j,ipoint)-expft3*dyunit
            engtmp(11*num+j,ipoint)=engtmp(11*num+j,ipoint)-expft3*dzunit
            engtmp(12*num+j,ipoint)=engtmp(12*num+j,ipoint)+expft*txx
            engtmp(13*num+j,ipoint)=engtmp(13*num+j,ipoint)+expft*tyy
            engtmp(14*num+j,ipoint)=engtmp(14*num+j,ipoint)+expft*tzz
            engtmp(15*num+j,ipoint)=engtmp(15*num+j,ipoint)+2.0d0*expft*txy
            engtmp(16*num+j,ipoint)=engtmp(16*num+j,ipoint)+2.0d0*expft*txz
            engtmp(17*num+j,ipoint)=engtmp(17*num+j,ipoint)+2.0d0*expft*tyz
            engtmp(18*num+j,ipoint)=engtmp(18*num+j,ipoint)+expft2*txx
            engtmp(19*num+j,ipoint)=engtmp(19*num+j,ipoint)+expft2*tyy
            engtmp(20*num+j,ipoint)=engtmp(20*num+j,ipoint)+expft2*tzz
            engtmp(21*num+j,ipoint)=engtmp(21*num+j,ipoint)+2.0d0*expft2*txy
            engtmp(22*num+j,ipoint)=engtmp(22*num+j,ipoint)+2.0d0*expft2*txz
            engtmp(23*num+j,ipoint)=engtmp(23*num+j,ipoint)+2.0d0*expft2*tyz
            engtmp(24*num+j,ipoint)=engtmp(24*num+j,ipoint)+expft3*txx
            engtmp(25*num+j,ipoint)=engtmp(25*num+j,ipoint)+expft3*tyy
            engtmp(26*num+j,ipoint)=engtmp(26*num+j,ipoint)+expft3*tzz
            engtmp(27*num+j,ipoint)=engtmp(27*num+j,ipoint)+2.0d0*expft3*txy
            engtmp(28*num+j,ipoint)=engtmp(28*num+j,ipoint)+2.0d0*expft3*txz
            engtmp(29*num+j,ipoint)=engtmp(29*num+j,ipoint)+2.0d0*expft3*tyz
         endif

!        engsr=engsr+expft+expft2+expft3+expft4
         sctmp(1)=sctmp(1)+expft+expft2+expft3+expft4
!<--- Parallelization_E
      enddo   
   enddo   
endif

if(ooaimlog) then
   if(cimlog) then
!---> Parallelization_S
!     do j=2,nanion
      do j = jst3, jed3
!<--- Parallelization_E
         jpoint=ntype(j)
!---> Parallelization_S
!        do i=1,j-1
         do i = ist3(j), ied3(j)
!<--- Parallelization_E
            ipoint=ntype(i)

!---> Memmory Reduction_S
            numx = numadr(i,j)
!           drsq=dxsav(i,j)*dxsav(i,j)+dysav(i,j)*dysav(i,j) &
!               +dzsav(i,j)*dzsav(i,j)
            drsq=dxsav(numx)*dxsav(numx)+dysav(numx)*dysav(numx) &
                +dzsav(numx)*dzsav(numx)
!<--- Memmory Reduction_E

            if(drsq.ge.rsqmaxsr) CYCLE

            dr=dsqrt(drsq)

            expft=ftb(ipoint,jpoint)*exp(-ftalp(ipoint,jpoint) &
                                       *(dr-delta(i)-delta(j)))
            expft2=ftb2(ipoint,jpoint)*exp(-ftbeta(ipoint,jpoint) &
                                       *(dr-delta(i)-delta(j)))
            expft3=ftb3(ipoint,jpoint)*exp(-ftgamma(ipoint,jpoint) &
                                       *(dr-delta(i)-delta(j)))

!---> Parallelization_S
!           engft1(i,jpoint)=engft1(i,jpoint)+expft
!           engft2(i,jpoint)=engft2(i,jpoint)+expft2
!           engft3(i,jpoint)=engft3(i,jpoint)+expft3
!           engft1(j,ipoint)=engft1(j,ipoint)+expft
!           engft2(j,ipoint)=engft2(j,ipoint)+expft2
!           engft3(j,ipoint)=engft3(j,ipoint)+expft3
            engtmp(i,jpoint)=engtmp(i,jpoint)+expft
            engtmp(num+i,jpoint)=engtmp(num+i,jpoint)+expft2
            engtmp(2*num+i,jpoint)=engtmp(2*num+i,jpoint)+expft3
            engtmp(j,ipoint)=engtmp(j,ipoint)+expft
            engtmp(num+j,ipoint)=engtmp(num+j,ipoint)+expft2
            engtmp(2*num+j,ipoint)=engtmp(2*num+j,ipoint)+expft3

!           engsr=engsr+expft+expft2+expft3
            sctmp(1)=sctmp(1)+expft+expft2+expft3
!<--- Parallelization_E
         enddo   
      enddo   
   else if(daimlog) then
!---> Parallelization_S
!     do j=2,nanion
      do j = jst3, jed3
!<--- Parallelization_E
         jpoint=ntype(j)
!---> Parallelization_S
!        do i=1,j-1
         do i = ist3(j), ied3(j)
!<--- Parallelization_E
            ipoint=ntype(i)

!---> Memmory Reduction_S
            numx = numadr(i,j)
!           drsq=dxsav(i,j)*dxsav(i,j)+dysav(i,j)*dysav(i,j) &
!               +dzsav(i,j)*dzsav(i,j)
            drsq=dxsav(numx)*dxsav(numx)+dysav(numx)*dysav(numx) &
                +dzsav(numx)*dzsav(numx)
!<--- Memmory Reduction_E

            if(drsq.ge.rsqmaxsr) CYCLE

            dr=dsqrt(drsq)
            drrec=1.0d0/dr
!---> Memmory Reduction_S
!           dx=dxsav(i,j)
!           dy=dysav(i,j)
!           dz=dzsav(i,j)
            dx=dxsav(numx)
            dy=dysav(numx)
            dz=dzsav(numx)
!<--- Memmory Reduction_E

            repsiloni=epsilonx(i)*dx+epsilony(i)*dy+epsilonz(i)*dz
            repsilonj=epsilonx(j)*dx+epsilony(j)*dy+epsilonz(j)*dz
            repsiloni=repsiloni*drrec
            repsilonj=repsilonj*drrec

            expft=ftb(ipoint,jpoint)*exp(-ftalp(ipoint,jpoint) &
                 *(dr-delta(i)-delta(j)-repsiloni+repsilonj))
            expft2=ftb2(ipoint,jpoint)*exp(-ftbeta(ipoint,jpoint) &
                 *(dr-delta(i)-delta(j)-repsiloni+repsilonj))
            expft3=ftb3(ipoint,jpoint)*exp(-ftgamma(ipoint,jpoint) &
                 *(dr-delta(i)-delta(j)-repsiloni+repsilonj))

!---> Parallelization_S
!           engft1(i,jpoint)=engft1(i,jpoint)+expft
!           engft2(i,jpoint)=engft2(i,jpoint)+expft2
!           engft3(i,jpoint)=engft3(i,jpoint)+expft3
!           engft1(j,ipoint)=engft1(j,ipoint)+expft
!           engft2(j,ipoint)=engft2(j,ipoint)+expft2
!           engft3(j,ipoint)=engft3(j,ipoint)+expft3
            engtmp(i,jpoint)=engtmp(i,jpoint)+expft
            engtmp(num+i,jpoint)=engtmp(num+i,jpoint)+expft2
            engtmp(2*num+i,jpoint)=engtmp(2*num+i,jpoint)+expft3
            engtmp(j,ipoint)=engtmp(j,ipoint)+expft
            engtmp(num+j,ipoint)=engtmp(num+j,ipoint)+expft2
            engtmp(2*num+j,ipoint)=engtmp(2*num+j,ipoint)+expft3

            dxunit=dx*drrec
            dyunit=dy*drrec
            dzunit=dz*drrec

!           engft1dotx(i,jpoint)=engft1dotx(i,jpoint)+expft*dxunit
!           engft1doty(i,jpoint)=engft1doty(i,jpoint)+expft*dyunit
!           engft1dotz(i,jpoint)=engft1dotz(i,jpoint)+expft*dzunit
!           engft2dotx(i,jpoint)=engft2dotx(i,jpoint)+expft2*dxunit
!           engft2doty(i,jpoint)=engft2doty(i,jpoint)+expft2*dyunit
!           engft2dotz(i,jpoint)=engft2dotz(i,jpoint)+expft2*dzunit
!           engft3dotx(i,jpoint)=engft3dotx(i,jpoint)+expft3*dxunit
!           engft3doty(i,jpoint)=engft3doty(i,jpoint)+expft3*dyunit
!           engft3dotz(i,jpoint)=engft3dotz(i,jpoint)+expft3*dzunit
            engtmp(3*num+i,jpoint)=engtmp(3*num+i,jpoint)+expft*dxunit
            engtmp(4*num+i,jpoint)=engtmp(4*num+i,jpoint)+expft*dyunit
            engtmp(5*num+i,jpoint)=engtmp(5*num+i,jpoint)+expft*dzunit
            engtmp(6*num+i,jpoint)=engtmp(6*num+i,jpoint)+expft2*dxunit
            engtmp(7*num+i,jpoint)=engtmp(7*num+i,jpoint)+expft2*dyunit
            engtmp(8*num+i,jpoint)=engtmp(8*num+i,jpoint)+expft2*dzunit
            engtmp(9*num+i,jpoint)=engtmp(9*num+i,jpoint)+expft3*dxunit
            engtmp(10*num+i,jpoint)=engtmp(10*num+i,jpoint)+expft3*dyunit
            engtmp(11*num+i,jpoint)=engtmp(11*num+i,jpoint)+expft3*dzunit

!           engft1dotx(j,ipoint)=engft1dotx(j,ipoint)-expft*dxunit
!           engft1doty(j,ipoint)=engft1doty(j,ipoint)-expft*dyunit
!           engft1dotz(j,ipoint)=engft1dotz(j,ipoint)-expft*dzunit
!           engft2dotx(j,ipoint)=engft2dotx(j,ipoint)-expft2*dxunit
!           engft2doty(j,ipoint)=engft2doty(j,ipoint)-expft2*dyunit
!           engft2dotz(j,ipoint)=engft2dotz(j,ipoint)-expft2*dzunit
!           engft3dotx(j,ipoint)=engft3dotx(j,ipoint)-expft3*dxunit
!           engft3doty(j,ipoint)=engft3doty(j,ipoint)-expft3*dyunit
!           engft3dotz(j,ipoint)=engft3dotz(j,ipoint)-expft3*dzunit
            engtmp(3*num+j,ipoint)=engtmp(3*num+j,ipoint)-expft*dxunit
            engtmp(4*num+j,ipoint)=engtmp(4*num+j,ipoint)-expft*dyunit
            engtmp(5*num+j,ipoint)=engtmp(5*num+j,ipoint)-expft*dzunit
            engtmp(6*num+j,ipoint)=engtmp(6*num+j,ipoint)-expft2*dxunit
            engtmp(7*num+j,ipoint)=engtmp(7*num+j,ipoint)-expft2*dyunit
            engtmp(8*num+j,ipoint)=engtmp(8*num+j,ipoint)-expft2*dzunit
            engtmp(9*num+j,ipoint)=engtmp(9*num+j,ipoint)-expft3*dxunit
            engtmp(10*num+j,ipoint)=engtmp(10*num+j,ipoint)-expft3*dyunit
            engtmp(11*num+j,ipoint)=engtmp(11*num+j,ipoint)-expft3*dzunit

!           engsr=engsr+expft+expft2+expft3
            sctmp(1)=sctmp(1)+expft+expft2+expft3
!<--- Parallelization_E
         enddo    
      enddo    
   else if(quaimlog) then
!---> Parallelization_S
!     do j=2,nanion
      do j = jst3, jed3
!<--- Parallelization_E
         jpoint=ntype(j)
!---> Parallelization_S
!        do i=1,j-1
         do i = ist3(j), ied3(j)
!<--- Parallelization_E
            ipoint=ntype(i)

!---> Memmory Reduction_S
            numx = numadr(i,j)
!           drsq=dxsav(i,j)*dxsav(i,j)+dysav(i,j)*dysav(i,j) &
!               +dzsav(i,j)*dzsav(i,j)
            drsq=dxsav(numx)*dxsav(numx)+dysav(numx)*dysav(numx) &
                +dzsav(numx)*dzsav(numx)
!<--- Memmory Reduction_E

            if(drsq.ge.rsqmaxsr) CYCLE

            drsqrec=1.0d0/drsq
            dr=dsqrt(drsq)
            drrec=1.0d0/dr
!---> Memmory Reduction_S
!           dx=dxsav(i,j)
!           dy=dysav(i,j)
!           dz=dzsav(i,j)
            dx=dxsav(numx)
            dy=dysav(numx)
            dz=dzsav(numx)
!<--- Memmory Reduction_E

            repsiloni=epsilonx(i)*dx+epsilony(i)*dy+epsilonz(i)*dz
            repsilonj=epsilonx(j)*dx+epsilony(j)*dy+epsilonz(j)*dz
            repsiloni=repsiloni*drrec
            repsilonj=repsilonj*drrec

            dxx=dx*dx
            dyy=dy*dy
            dzz=dz*dz
            dxy=dx*dy
            dxz=dx*dz
            dyz=dy*dz
            txx=3.0d0*dxx*drsqrec-1.0d0
            tyy=3.0d0*dyy*drsqrec-1.0d0
            tzz=3.0d0*dzz*drsqrec-1.0d0
            txy=3.0d0*dxy*drsqrec
            txz=3.0d0*dxz*drsqrec
            tyz=3.0d0*dyz*drsqrec

            rrthetai=quaimxx(i)*txx+quaimyy(i)*tyy+quaimzz(i)*tzz+ &
                 2.0d0*(quaimxy(i)*txy+quaimxz(i)*txz+quaimyz(i)*tyz)
            rrthetaj=quaimxx(j)*txx+quaimyy(j)*tyy+quaimzz(j)*tzz+ &
                 2.0d0*(quaimxy(j)*txy+quaimxz(j)*txz+quaimyz(j)*tyz)

            expft=ftb(ipoint,jpoint)*exp(-ftalp(ipoint,jpoint) &
                  *(dr-delta(i)-delta(j)-repsiloni+repsilonj &
                                        -rrthetaj-rrthetai))
            expft2=ftb2(ipoint,jpoint)*exp(-ftbeta(ipoint,jpoint) &
                  *(dr-delta(i)-delta(j)-repsiloni+repsilonj &
                                        -rrthetaj-rrthetai))
            expft3=ftb3(ipoint,jpoint)*exp(-ftgamma(ipoint,jpoint) &
                  *(dr-delta(i)-delta(j)-repsiloni+repsilonj &
                                        -rrthetaj-rrthetai))

!---> Parallelization_S
!           engft1(i,jpoint)=engft1(i,jpoint)+expft
!           engft2(i,jpoint)=engft2(i,jpoint)+expft2
!           engft3(i,jpoint)=engft3(i,jpoint)+expft3
!           engft1(j,ipoint)=engft1(j,ipoint)+expft
!           engft2(j,ipoint)=engft2(j,ipoint)+expft2
!           engft3(j,ipoint)=engft3(j,ipoint)+expft3
            engtmp(i,jpoint)=engtmp(i,jpoint)+expft
            engtmp(num+i,jpoint)=engtmp(num+i,jpoint)+expft2
            engtmp(2*num+i,jpoint)=engtmp(2*num+i,jpoint)+expft3
            engtmp(j,ipoint)=engtmp(j,ipoint)+expft
            engtmp(num+j,ipoint)=engtmp(num+j,ipoint)+expft2
            engtmp(2*num+j,ipoint)=engtmp(2*num+j,ipoint)+expft3

            dxunit=dx*drrec
            dyunit=dy*drrec
            dzunit=dz*drrec

!           engft1dotx(i,jpoint)=engft1dotx(i,jpoint)+expft*dxunit
!           engft1doty(i,jpoint)=engft1doty(i,jpoint)+expft*dyunit
!           engft1dotz(i,jpoint)=engft1dotz(i,jpoint)+expft*dzunit
!           engft2dotx(i,jpoint)=engft2dotx(i,jpoint)+expft2*dxunit
!           engft2doty(i,jpoint)=engft2doty(i,jpoint)+expft2*dyunit
!           engft2dotz(i,jpoint)=engft2dotz(i,jpoint)+expft2*dzunit
!           engft3dotx(i,jpoint)=engft3dotx(i,jpoint)+expft3*dxunit
!           engft3doty(i,jpoint)=engft3doty(i,jpoint)+expft3*dyunit
!           engft3dotz(i,jpoint)=engft3dotz(i,jpoint)+expft3*dzunit
            engtmp(3*num+i,jpoint)=engtmp(3*num+i,jpoint)+expft*dxunit
            engtmp(4*num+i,jpoint)=engtmp(4*num+i,jpoint)+expft*dyunit
            engtmp(5*num+i,jpoint)=engtmp(5*num+i,jpoint)+expft*dzunit
            engtmp(6*num+i,jpoint)=engtmp(6*num+i,jpoint)+expft2*dxunit
            engtmp(7*num+i,jpoint)=engtmp(7*num+i,jpoint)+expft2*dyunit
            engtmp(8*num+i,jpoint)=engtmp(8*num+i,jpoint)+expft2*dzunit
            engtmp(9*num+i,jpoint)=engtmp(9*num+i,jpoint)+expft3*dxunit
            engtmp(10*num+i,jpoint)=engtmp(10*num+i,jpoint)+expft3*dyunit
            engtmp(11*num+i,jpoint)=engtmp(11*num+i,jpoint)+expft3*dzunit

!           engft1dotx(j,ipoint)=engft1dotx(j,ipoint)-expft*dxunit
!           engft1doty(j,ipoint)=engft1doty(j,ipoint)-expft*dyunit
!           engft1dotz(j,ipoint)=engft1dotz(j,ipoint)-expft*dzunit
!           engft2dotx(j,ipoint)=engft2dotx(j,ipoint)-expft2*dxunit
!           engft2doty(j,ipoint)=engft2doty(j,ipoint)-expft2*dyunit
!           engft2dotz(j,ipoint)=engft2dotz(j,ipoint)-expft2*dzunit
!           engft3dotx(j,ipoint)=engft3dotx(j,ipoint)-expft3*dxunit
!           engft3doty(j,ipoint)=engft3doty(j,ipoint)-expft3*dyunit
!           engft3dotz(j,ipoint)=engft3dotz(j,ipoint)-expft3*dzunit
            engtmp(3*num+j,ipoint)=engtmp(3*num+j,ipoint)-expft*dxunit
            engtmp(4*num+j,ipoint)=engtmp(4*num+j,ipoint)-expft*dyunit
            engtmp(5*num+j,ipoint)=engtmp(5*num+j,ipoint)-expft*dzunit
            engtmp(6*num+j,ipoint)=engtmp(6*num+j,ipoint)-expft2*dxunit
            engtmp(7*num+j,ipoint)=engtmp(7*num+j,ipoint)-expft2*dyunit
            engtmp(8*num+j,ipoint)=engtmp(8*num+j,ipoint)-expft2*dzunit
            engtmp(9*num+j,ipoint)=engtmp(9*num+j,ipoint)-expft3*dxunit
            engtmp(10*num+j,ipoint)=engtmp(10*num+j,ipoint)-expft3*dyunit
            engtmp(11*num+j,ipoint)=engtmp(11*num+j,ipoint)-expft3*dzunit

!           engft1dotxx(i,jpoint)=engft1dotxx(i,jpoint)+expft*txx
!           engft1dotyy(i,jpoint)=engft1dotyy(i,jpoint)+expft*tyy
!           engft1dotzz(i,jpoint)=engft1dotzz(i,jpoint)+expft*tzz
!           engft1dotxy(i,jpoint)=engft1dotxy(i,jpoint)+2.0d0*expft*txy
!           engft1dotxz(i,jpoint)=engft1dotxz(i,jpoint)+2.0d0*expft*txz
!           engft1dotyz(i,jpoint)=engft1dotyz(i,jpoint)+2.0d0*expft*tyz
!           engft2dotxx(i,jpoint)=engft2dotxx(i,jpoint)+expft2*txx
!           engft2dotyy(i,jpoint)=engft2dotyy(i,jpoint)+expft2*tyy
!           engft2dotzz(i,jpoint)=engft2dotzz(i,jpoint)+expft2*tzz
!           engft2dotxy(i,jpoint)=engft2dotxy(i,jpoint)+2.0d0*expft2*txy
!           engft2dotxz(i,jpoint)=engft2dotxz(i,jpoint)+2.0d0*expft2*txz
!           engft2dotyz(i,jpoint)=engft2dotyz(i,jpoint)+2.0d0*expft2*tyz
!           engft3dotxx(i,jpoint)=engft3dotxx(i,jpoint)+expft3*txx
!           engft3dotyy(i,jpoint)=engft3dotyy(i,jpoint)+expft3*tyy
!           engft3dotzz(i,jpoint)=engft3dotzz(i,jpoint)+expft3*tzz
!           engft3dotxy(i,jpoint)=engft3dotxy(i,jpoint)+2.0d0*expft3*txy
!           engft3dotxz(i,jpoint)=engft3dotxz(i,jpoint)+2.0d0*expft3*txz
!           engft3dotyz(i,jpoint)=engft3dotyz(i,jpoint)+2.0d0*expft3*tyz
            engtmp(12*num+i,jpoint)=engtmp(12*num+i,jpoint)+expft*txx
            engtmp(13*num+i,jpoint)=engtmp(13*num+i,jpoint)+expft*tyy
            engtmp(14*num+i,jpoint)=engtmp(14*num+i,jpoint)+expft*tzz
            engtmp(15*num+i,jpoint)=engtmp(15*num+i,jpoint)+2.0d0*expft*txy
            engtmp(16*num+i,jpoint)=engtmp(16*num+i,jpoint)+2.0d0*expft*txz
            engtmp(17*num+i,jpoint)=engtmp(17*num+i,jpoint)+2.0d0*expft*tyz
            engtmp(18*num+i,jpoint)=engtmp(18*num+i,jpoint)+expft2*txx
            engtmp(19*num+i,jpoint)=engtmp(19*num+i,jpoint)+expft2*tyy
            engtmp(20*num+i,jpoint)=engtmp(20*num+i,jpoint)+expft2*tzz
            engtmp(21*num+i,jpoint)=engtmp(21*num+i,jpoint)+2.0d0*expft2*txy
            engtmp(22*num+i,jpoint)=engtmp(22*num+i,jpoint)+2.0d0*expft2*txz
            engtmp(23*num+i,jpoint)=engtmp(23*num+i,jpoint)+2.0d0*expft2*tyz
            engtmp(24*num+i,jpoint)=engtmp(24*num+i,jpoint)+expft3*txx
            engtmp(25*num+i,jpoint)=engtmp(25*num+i,jpoint)+expft3*tyy
            engtmp(26*num+i,jpoint)=engtmp(26*num+i,jpoint)+expft3*tzz
            engtmp(27*num+i,jpoint)=engtmp(27*num+i,jpoint)+2.0d0*expft3*txy
            engtmp(28*num+i,jpoint)=engtmp(28*num+i,jpoint)+2.0d0*expft3*txz
            engtmp(29*num+i,jpoint)=engtmp(29*num+i,jpoint)+2.0d0*expft3*tyz

!           engft1dotxx(j,ipoint)=engft1dotxx(j,ipoint)+expft*txx
!           engft1dotyy(j,ipoint)=engft1dotyy(j,ipoint)+expft*tyy
!           engft1dotzz(j,ipoint)=engft1dotzz(j,ipoint)+expft*tzz
!           engft1dotxy(j,ipoint)=engft1dotxy(j,ipoint)+2.0d0*expft*txy
!           engft1dotxz(j,ipoint)=engft1dotxz(j,ipoint)+2.0d0*expft*txz
!           engft1dotyz(j,ipoint)=engft1dotyz(j,ipoint)+2.0d0*expft*tyz
!           engft2dotxx(j,ipoint)=engft2dotxx(j,ipoint)+expft2*txx
!           engft2dotyy(j,ipoint)=engft2dotyy(j,ipoint)+expft2*tyy
!           engft2dotzz(j,ipoint)=engft2dotzz(j,ipoint)+expft2*tzz
!           engft2dotxy(j,ipoint)=engft2dotxy(j,ipoint)+2.0d0*expft2*txy
!           engft2dotxz(j,ipoint)=engft2dotxz(j,ipoint)+2.0d0*expft2*txz
!           engft2dotyz(j,ipoint)=engft2dotyz(j,ipoint)+2.0d0*expft2*tyz
!           engft3dotxx(j,ipoint)=engft3dotxx(j,ipoint)+expft3*txx
!           engft3dotyy(j,ipoint)=engft3dotyy(j,ipoint)+expft3*tyy
!           engft3dotzz(j,ipoint)=engft3dotzz(j,ipoint)+expft3*tzz
!           engft3dotxy(j,ipoint)=engft3dotxy(j,ipoint)+2.0d0*expft3*txy
!           engft3dotxz(j,ipoint)=engft3dotxz(j,ipoint)+2.0d0*expft3*txz
!           engft3dotyz(j,ipoint)=engft3dotyz(j,ipoint)+2.0d0*expft3*tyz
            engtmp(12*num+j,ipoint)=engtmp(12*num+j,ipoint)+expft*txx
            engtmp(13*num+j,ipoint)=engtmp(13*num+j,ipoint)+expft*tyy
            engtmp(14*num+j,ipoint)=engtmp(14*num+j,ipoint)+expft*tzz
            engtmp(15*num+j,ipoint)=engtmp(15*num+j,ipoint)+2.0d0*expft*txy
            engtmp(16*num+j,ipoint)=engtmp(16*num+j,ipoint)+2.0d0*expft*txz
            engtmp(17*num+j,ipoint)=engtmp(17*num+j,ipoint)+2.0d0*expft*tyz
            engtmp(18*num+j,ipoint)=engtmp(18*num+j,ipoint)+expft2*txx
            engtmp(19*num+j,ipoint)=engtmp(19*num+j,ipoint)+expft2*tyy
            engtmp(20*num+j,ipoint)=engtmp(20*num+j,ipoint)+expft2*tzz
            engtmp(21*num+j,ipoint)=engtmp(21*num+j,ipoint)+2.0d0*expft2*txy
            engtmp(22*num+j,ipoint)=engtmp(22*num+j,ipoint)+2.0d0*expft2*txz
            engtmp(23*num+j,ipoint)=engtmp(23*num+j,ipoint)+2.0d0*expft2*tyz
            engtmp(24*num+j,ipoint)=engtmp(24*num+j,ipoint)+expft3*txx
            engtmp(25*num+j,ipoint)=engtmp(25*num+j,ipoint)+expft3*tyy
            engtmp(26*num+j,ipoint)=engtmp(26*num+j,ipoint)+expft3*tzz
            engtmp(27*num+j,ipoint)=engtmp(27*num+j,ipoint)+2.0d0*expft3*txy
            engtmp(28*num+j,ipoint)=engtmp(28*num+j,ipoint)+2.0d0*expft3*txz
            engtmp(29*num+j,ipoint)=engtmp(29*num+j,ipoint)+2.0d0*expft3*tyz

!           engsr=engsr+expft+expft2+expft3
            sctmp(1)=sctmp(1)+expft+expft2+expft3
!<--- Parallelization_E
         enddo   
      enddo   
   endif
endif

!---> Parallelization_S
if(cimlog) then

   CALL MPI_ALLREDUCE(sctmp,scctmp,1,mpi_double_precision,MPI_SUM, &
                      MPI_COMM_WORLD,ierr)
   engsr=scctmp(1)

   CALL MPI_ALLREDUCE(engtmp,enggtmp,30*num*nspec,mpi_double_precision,MPI_SUM, &
                      MPI_COMM_WORLD,ierr)
   engft1(:,:) = enggtmp(1:num,:)
   engft2(:,:) = enggtmp(num+1:2*num,:)
   engft3(:,:) = enggtmp(2*num+1:3*num,:)

else if(daimlog) then

   CALL MPI_ALLREDUCE(sctmp,scctmp,1,mpi_double_precision,MPI_SUM, &
                      MPI_COMM_WORLD,ierr)
   engsr=scctmp(1)

   CALL MPI_ALLREDUCE(engtmp,enggtmp,30*num*nspec,mpi_double_precision,MPI_SUM, &
                      MPI_COMM_WORLD,ierr)
   engft1(:,:) = enggtmp(1:num,:)
   engft2(:,:) = enggtmp(num+1:2*num,:)
   engft3(:,:) = enggtmp(2*num+1:3*num,:)
   engft1dotx(:,:) = enggtmp(3*num+1:4*num,:)
   engft1doty(:,:) = enggtmp(4*num+1:5*num,:)
   engft1dotz(:,:) = enggtmp(5*num+1:6*num,:)
   engft2dotx(:,:) = enggtmp(6*num+1:7*num,:)
   engft2doty(:,:) = enggtmp(7*num+1:8*num,:)
   engft2dotz(:,:) = enggtmp(8*num+1:9*num,:)
   engft3dotx(:,:) = enggtmp(9*num+1:10*num,:)
   engft3doty(:,:) = enggtmp(10*num+1:11*num,:)
   engft3dotz(:,:) = enggtmp(11*num+1:12*num,:)

else if(quaimlog) then

   CALL MPI_ALLREDUCE(sctmp,scctmp,1,mpi_double_precision,MPI_SUM, &
                      MPI_COMM_WORLD,ierr)
   engsr=scctmp(1)

   CALL MPI_ALLREDUCE(engtmp,enggtmp,30*num*nspec,mpi_double_precision,MPI_SUM, &
                      MPI_COMM_WORLD,ierr)
   engft1(:,:) = enggtmp(1:num,:)
   engft2(:,:) = enggtmp(num+1:2*num,:)
   engft3(:,:) = enggtmp(2*num+1:3*num,:)
   engft1dotx(:,:) = enggtmp(3*num+1:4*num,:)
   engft1doty(:,:) = enggtmp(4*num+1:5*num,:)
   engft1dotz(:,:) = enggtmp(5*num+1:6*num,:)
   engft2dotx(:,:) = enggtmp(6*num+1:7*num,:)
   engft2doty(:,:) = enggtmp(7*num+1:8*num,:)
   engft2dotz(:,:) = enggtmp(8*num+1:9*num,:)
   engft3dotx(:,:) = enggtmp(9*num+1:10*num,:)
   engft3doty(:,:) = enggtmp(10*num+1:11*num,:)
   engft3dotz(:,:) = enggtmp(11*num+1:12*num,:)
   engft1dotxx(:,:) = enggtmp(12*num+1:13*num,:)
   engft1dotyy(:,:) = enggtmp(13*num+1:14*num,:)
   engft1dotzz(:,:) = enggtmp(14*num+1:15*num,:)
   engft1dotxy(:,:) = enggtmp(15*num+1:16*num,:)
   engft1dotxz(:,:) = enggtmp(16*num+1:17*num,:)
   engft1dotyz(:,:) = enggtmp(17*num+1:18*num,:)
   engft2dotxx(:,:) = enggtmp(18*num+1:19*num,:)
   engft2dotyy(:,:) = enggtmp(19*num+1:20*num,:)
   engft2dotzz(:,:) = enggtmp(20*num+1:21*num,:)
   engft2dotxy(:,:) = enggtmp(21*num+1:22*num,:)
   engft2dotxz(:,:) = enggtmp(22*num+1:23*num,:)
   engft2dotyz(:,:) = enggtmp(23*num+1:24*num,:)
   engft3dotxx(:,:) = enggtmp(24*num+1:25*num,:)
   engft3dotyy(:,:) = enggtmp(25*num+1:26*num,:)
   engft3dotzz(:,:) = enggtmp(26*num+1:27*num,:)
   engft3dotxy(:,:) = enggtmp(27*num+1:28*num,:)
   engft3dotxz(:,:) = enggtmp(28*num+1:29*num,:)
   engft3dotyz(:,:) = enggtmp(29*num+1:30*num,:)
endif
!<--- Parallelization_E
engpetot=engpetot+engsr

return
END SUBROUTINE
