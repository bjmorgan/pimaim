SUBROUTINE DFUNCAIM(p,xi)

USE commondata, ONLY: nspec,nsp,deformablelog,chg,ftalp,ftbeta,ftgamma, &
                      engft1,engft2,engft3,selfC,selfeps,selfH,selfquaim, &
                      environmentalaimlog,selfeps2,engft1dotx,engft1doty, &
                      engft1dotz,engft2dotx,engft2doty,engft2dotz, &
                      engft3dotx,engft3doty,engft3dotz,engft1dotxx, &
                      engft1dotyy,engft1dotzz,engft1dotxy,engft1dotxz, &
                      engft1dotyz,engft2dotxx,engft2dotyy,engft2dotzz, &
                      engft2dotxy,engft2dotxz,engft2dotyz,engft3dotxx, &
                      engft3dotyy,engft3dotzz,engft3dotxy,engft3dotxz, &
                      engft3dotyz,environmentalpimlog,xmu,ymu,zmu,quadxx, &
                      quadyy,quadzz,quadxy,quadxz,quadyz,xk1,alppolar1, &
                      Bpolar1,Cpolar1,gammapolar1,alppolardel,num,selfB, &
                      selfgam,dippimlog

IMPLICIT NONE

INTEGER :: i,j,k,l
DOUBLE PRECISION :: thetasq,epsilonsq,selfexp,selfengderiv,deltasq, &
                    dipsq,thacetheta3,quadsq,dipquadsq
DOUBLE PRECISION, DIMENSION(num*10), INTENT(IN) :: p
DOUBLE PRECISION, DIMENSION(num*10), INTENT(OUT) :: xi
      
l=1
do j=1,nspec

   if(deformablelog(j))then
      do i=l,nsp(j)+l-1
! CIM component
         if(chg(j).lt.0.0d0)then
            xi(i)=SUM(ftalp(j,:)*engft1(i,:) &
                  +ftbeta(j,:)*engft2(i,:)+ftgamma(j,:)*engft3(i,:))
            xi(i)=xi(i)-(selfgam(j)*selfB(j)*(dexp(-selfgam(j)*p(i)) &
                 -dexp(selfgam(j)*p(i))))
         else
            deltasq=p(i)*p(i)
            selfexp=selfB(j)*dexp(selfgam(j)*selfgam(j)*deltasq)
            selfengderiv=2.0d0*p(i)*selfB(j)*selfgam(j)*selfgam(j)*selfexp
            xi(i)=SUM(ftalp(j,:)*engft1(i,:) &
                  +ftbeta(j,:)*engft2(i,:)+ftgamma(j,:)*engft3(i,:))
            xi(i)=xi(i)+selfengderiv
         endif
! DAIM components
         epsilonsq=p(i+num)*p(i+num)+ &
                   p(i+2*num)*p(i+2*num)+p(i+3*num)*p(i+3*num)
         selfexp=selfC(j)*dexp(selfeps(i)*selfeps(i)*epsilonsq)
         selfengderiv=2.0d0*selfeps(i)*selfeps(i)*selfexp

         if(environmentalaimlog)then
            xi(i)=xi(i)-selfengderiv*epsilonsq*selfeps2(j)
         endif

         xi(i+num)=SUM(ftalp(j,:)*engft1dotx(i,:) &
                  +ftbeta(j,:)*engft2dotx(i,:)+ftgamma(j,:)*engft3dotx(i,:))
         xi(i+2*num)=SUM(ftalp(j,:)*engft1doty(i,:) &
                    +ftbeta(j,:)*engft2doty(i,:)+ftgamma(j,:)*engft3doty(i,:))
         xi(i+3*num)=SUM(ftalp(j,:)*engft1dotz(i,:) &
                    +ftbeta(j,:)*engft2dotz(i,:)+ftgamma(j,:)*engft3dotz(i,:))
         xi(i+num)=xi(i+num)+p(i+num)*selfengderiv
         xi(i+2*num)=xi(i+2*num)+p(i+2*num)*selfengderiv
         xi(i+3*num)=xi(i+3*num)+p(i+3*num)*selfengderiv
! QUAIM components
         thetasq=p(i+4*num)*p(i+4*num)+p(i+5*num)*p(i+5*num)+ &
                 p(i+6*num)*p(i+6*num)+2.0d0*(p(i+7*num)*p(i+7*num)+ &
                        p(i+8*num)*p(i+8*num)+p(i+9*num)*p(i+9*num))
         selfexp=selfH(j)*exp(selfquaim(i)*selfquaim(i)*thetasq)
         selfengderiv=2.0d0*selfquaim(i)*selfquaim(i)*selfexp

         if(environmentalaimlog)then
            xi(i)=xi(i)-selfengderiv*thetasq*selfeps2(j)
         endif

         xi(i+4*num)=SUM(ftalp(j,:)*engft1dotxx(i,:) &
                        +ftbeta(j,:)*engft2dotxx(i,:) &
                        +ftgamma(j,:)*engft3dotxx(i,:))
         xi(i+5*num)=SUM(ftalp(j,:)*engft1dotyy(i,:) &
                        +ftbeta(j,:)*engft2dotyy(i,:) &
                        +ftgamma(j,:)*engft3dotyy(i,:))
         xi(i+6*num)=SUM(ftalp(j,:)*engft1dotzz(i,:) &
                        +ftbeta(j,:)*engft2dotzz(i,:) &
                        +ftgamma(j,:)*engft3dotzz(i,:))
         xi(i+7*num)=SUM(ftalp(j,:)*engft1dotxy(i,:) &
                        +ftbeta(j,:)*engft2dotxy(i,:) &
                        +ftgamma(j,:)*engft3dotxy(i,:))
         xi(i+8*num)=SUM(ftalp(j,:)*engft1dotxz(i,:) &
                        +ftbeta(j,:)*engft2dotxz(i,:) &
                        +ftgamma(j,:)*engft3dotxz(i,:))
         xi(i+9*num)=SUM(ftalp(j,:)*engft1dotyz(i,:) &
                        +ftbeta(j,:)*engft2dotyz(i,:) &
                        +ftgamma(j,:)*engft3dotyz(i,:))
         xi(i+4*num)=xi(i+4*num)+p(i+4*num)*selfengderiv
         xi(i+5*num)=xi(i+5*num)+p(i+5*num)*selfengderiv
         xi(i+6*num)=xi(i+6*num)+p(i+6*num)*selfengderiv
         xi(i+7*num)=xi(i+7*num)+2.0d0*p(i+7*num)*selfengderiv
         xi(i+8*num)=xi(i+8*num)+2.0d0*p(i+8*num)*selfengderiv
         xi(i+9*num)=xi(i+9*num)+2.0d0*p(i+9*num)*selfengderiv

         if(environmentalpimlog)then
! PIM component-influences CIM
            dipsq=xmu(i)*xmu(i)+ymu(i)*ymu(i)+zmu(i)*zmu(i)
            quadsq=quadxx(i)*quadxx(i)+quadyy(i)*quadyy(i)+ &
                   quadzz(i)*quadzz(i)+2.0d0*quadxy(i)*quadxy(i)+ &
                   2.0d0*quadxz(i)*quadxz(i)+2.0d0*quadyz(i)*quadyz(i)
            dipquadsq=xmu(i)*quadxx(i)*xmu(i)+ymu(i)*quadyy(i)*ymu(i)+  &
                      zmu(i)*quadzz(i)*zmu(i)+2.0d0*xmu(i)*quadxy(i)*ymu(i)+ &
                      2.0d0*xmu(i)*quadxz(i)*zmu(i)+ &
                      2.0d0*ymu(i)*quadyz(i)*zmu(i)

            if(dippimlog) then
               selfexp=dipsq
            else
               selfexp=dipsq+quadsq/(3.0d0*Cpolar1(j))- &
                       dipquadsq*xk1(i)*2.0d0*Bpolar1(j)/Cpolar1(j)+ &
                       dipsq*dipsq*xk1(i)*xk1(i)* &
                          (2.0d0*gammapolar1(j)*Cpolar1(j) &
                         -3.0d0*Bpolar1(j)*Bpolar1(j))/(2.0d0*Cpolar1(j))
            endif
            selfengderiv=-alppolar1(j)*alppolardel(j)* &
                          dexp(alppolardel(j)*p(i))*xk1(i)*xk1(i)

            xi(i)=xi(i)+selfexp*selfengderiv
         endif
      enddo    
   endif
   l=l+nsp(j)
enddo   

return
END SUBROUTINE
