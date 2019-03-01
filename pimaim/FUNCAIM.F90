FUNCTION FUNCAIM(p)

USE commondata, ONLY: selfengtot,epsselfengtot,quaimselfengtot,nspec,nsp, &
                      deformablelog,environmentalaimlog,selfeps,selfquaim, &
                      selfeps1,selfeps2,selfquaim1, &
                      selfB,selfC,selfH,selfgam,chg, &
                      environmentalpimlog,ntype,xk1,xk2,xk3,xk4,alppolar, &
                      Bpolar,Cpolar,gammapolar,alppolar1,Bpolar1,Cpolar1, &
                      gammapolar1,alppolardel,alppolareff,engeff,xmu,ymu,zmu, &
                      quadxx,quadyy,quadzz,quadxy,quadxz,quadyz,engpetot,num, &
                      quadpimlog,polarizablelog

IMPLICIT NONE

INTEGER :: i,j,l,ipoint
DOUBLE PRECISION, DIMENSION(num*10), INTENT(IN) :: p
DOUBLE PRECISION :: FUNCAIM
DOUBLE PRECISION :: thetasq,epsilonsq,deltasq,selfexp, &
                    epsselfeng,quaimselfeng,selfeng, &
                    dipsq,quadsq,dipquadsq,pimselfeng

FUNCAIM=0.d0
selfengtot=0.d0
epsselfengtot=0.d0
quaimselfengtot=0.d0

l=1
do j=1,nspec
   if(deformablelog(j)) then
      do i=l,nsp(j)+l-1
         if(environmentalaimlog)then
            selfeps(i)=selfeps1(j)*dexp(-selfeps2(j)*p(i))
            selfquaim(i)=selfquaim1(j)*selfeps(i)
         endif
! CIM component
         if(chg(j).lt.0.0d0)then
            selfeng=selfB(j)*(dexp(-selfgam(j)*p(i))+ &
                              dexp(selfgam(j)*p(i)))
            selfengtot=selfengtot+selfeng
         else
            deltasq=p(i)*p(i)
            selfexp=selfB(j)*dexp(selfgam(j)*selfgam(j)*deltasq)
            selfeng=selfexp-selfB(j)
            selfengtot=selfengtot+selfeng
         endif
! DAIM component
         epsilonsq=p(i+num)*p(i+num)+p(i+2*num)*p(i+2*num)+ &
                   p(i+3*num)*p(i+3*num)
         selfexp=selfC(j)*dexp(selfeps(i)*selfeps(i)*epsilonsq)
         epsselfeng=selfexp-selfC(j)
         epsselfengtot=epsselfengtot+epsselfeng
! QUAIM component
         thetasq=p(i+4*num)*p(i+4*num)+p(i+5*num)*p(i+5*num)+ &
                 p(i+6*num)*p(i+6*num)+2.0d0*(p(i+7*num)*p(i+7*num)+ &
                 p(i+8*num)*p(i+8*num)+p(i+9*num)*p(i+9*num))
         selfexp=selfH(j)*exp(selfquaim(i)*selfquaim(i)*thetasq)
         quaimselfeng=selfexp-selfH(j)
         quaimselfengtot=quaimselfengtot+quaimselfeng

         if(environmentalpimlog)then
            ipoint=ntype(i)
! PIM component
            alppolar(i)=alppolar1(ipoint)*(dexp(alppolardel(ipoint)* &
                       p(i))+dexp(-alppolareff(ipoint)*engeff(i)))*0.5d0
            Bpolar(i)=Bpolar1(ipoint)*alppolar(i)
            Cpolar(i)=Cpolar1(ipoint)*alppolar(i)
            gammapolar(i)=gammapolar1(ipoint)*alppolar(i)
            if(polarizablelog(j)) xk1(i)=1.0d0/(2.0d0*alppolar(i))
            if(quadpimlog)then
               xk3(i)=1.0d0/(6.0d0*Cpolar(i))
               xk2(i)=-Bpolar(i)/(4.0d0*alppolar(i)*alppolar(i)*Cpolar(i))
               xk4(i)=(2.0d0*gammapolar(i)*Cpolar(i)-3.0d0*Bpolar(i)* &
                            Bpolar(i))/(48.0d0*Cpolar(i)* &
                             alppolar(i)*alppolar(i)*alppolar(i)*alppolar(i))
            endif

            dipsq=xmu(i)*xmu(i)+ymu(i)*ymu(i)+zmu(i)*zmu(i)
            quadsq=quadxx(i)*quadxx(i)+quadyy(i)*quadyy(i)+ &
                   quadzz(i)*quadzz(i)+2.0d0*quadxy(i)*quadxy(i)+ &
                   2.0d0*quadxz(i)*quadxz(i)+2.0d0*quadyz(i)*quadyz(i)
            dipquadsq=xmu(i)*quadxx(i)*xmu(i)+ymu(i)*quadyy(i)*ymu(i)+ &
                      zmu(i)*quadzz(i)*zmu(i)+2.0d0*xmu(i)*quadxy(i)*ymu(i)+ &
                      2.0d0*xmu(i)*quadxz(i)*zmu(i)+ &
                      2.0d0*ymu(i)*quadyz(i)*zmu(i)

            pimselfeng=xk1(i)*dipsq+xk4(i)*dipsq*dipsq+ &
                       xk2(i)*dipquadsq+xk3(i)*quadsq

            FUNCAIM=FUNCAIM+pimselfeng
         endif
      enddo   
   endif
   l=l+nsp(j)
enddo    

FUNCAIM=FUNCAIM+engpetot+selfengtot+epsselfengtot+quaimselfengtot

return
END FUNCTION
