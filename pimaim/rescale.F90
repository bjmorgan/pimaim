SUBROUTINE rescale

USE commondata, ONLY: nspec,nsp,vx,vy,vz,hmass,gtrankin,trantemp,num
USE boxdata, ONLY: hlab2

!---> Parallelization_S
use mpipara
!<--- Parallelization_E
IMPLICIT NONE

DOUBLE PRECISION :: xke2,vxtot,vytot,vztot,sprec,vxmean,vymean, &
                    vzmean,vttx,vtty,vttz,tempfac,tkin 
INTEGER :: n,n1,n2,i,j                                

n1=1
n2=1
xke2=0.0d0

do i=1,nspec
   vxtot=0.0d0
   vytot=0.0d0
   vztot=0.0d0
   if(nsp(i).ne.1) then
      do j=1,nsp(i)
         vxtot=vxtot+vx(n1)
         vytot=vytot+vy(n1)
         vztot=vztot+vz(n1)
         n1=n1+1
      enddo   
      sprec=1.0d0/dble(nsp(i))
      vxmean=vxtot*sprec
      vymean=vytot*sprec
      vzmean=vztot*sprec
      do j=1,nsp(i)
         vx(n2)=vx(n2)-vxmean
         vy(n2)=vy(n2)-vymean
         vz(n2)=vz(n2)-vzmean
!
! Noncubic box code:
!
         vttx=hlab2(1,1)*vx(n2)+hlab2(1,2)*vy(n2)+hlab2(1,3)*vz(n2)
         vtty=hlab2(2,1)*vx(n2)+hlab2(2,2)*vy(n2)+hlab2(2,3)*vz(n2)
         vttz=hlab2(3,1)*vx(n2)+hlab2(3,2)*vy(n2)+hlab2(3,3)*vz(n2)

         xke2=xke2+(vttx*vttx+vtty*vtty+vttz*vttz)*hmass(i)
         n2=n2+1
      enddo   
   else
!
! Cope with the case that we only have one
! ion of a given species. (The mean velocity is
! then the actual velocity!)
!
         vttx=hlab2(1,1)*vx(n1)+hlab2(1,2)*vy(n1)+hlab2(1,3)*vz(n1)
         vtty=hlab2(2,1)*vx(n1)+hlab2(2,2)*vy(n1)+hlab2(2,3)*vz(n1)
         vttz=hlab2(3,1)*vx(n1)+hlab2(3,2)*vy(n1)+hlab2(3,3)*vz(n1)

         xke2=xke2+(vttx*vttx+vtty*vtty+vttz*vttz)*hmass(i)
         n1=n1+1
         n2=n2+1
   endif
enddo   
!
! Calculate the kinetic temperature and the rescaling factor.
!
tkin=xke2*gtrankin
tempfac=dsqrt(trantemp/tkin)
n=1
do i=1,nspec
   do j=1,nsp(i)
      vx(n)=vx(n)*tempfac
      vy(n)=vy(n)*tempfac
      vz(n)=vz(n)*tempfac
      n=n+1
   enddo   
enddo   

!---> Parallelization_S
if( iam .eq. 0 ) then

write(6,*)
write(6,*)'**** Velocity rescaling complete ****'
write(6,*)

endif
!<--- Parallelization_E

return
END SUBROUTINE
