SUBROUTINE dynmat

USE commondata, ONLY: num,x,y,z,amove,amovefac,frrx,frry,frrz, &
                      nmat1,nmat2,environmentalpimlog,conjgradaimlog, &
                      conjgradlog
!---> Parallelization_S
use mpipara
!<--- Parallelization_E

IMPLICIT NONE

INTEGER :: i,j,ii,k
DOUBLE PRECISION :: fdiffx,fdiffy,fdiffz

DOUBLE PRECISION, DIMENSION(num) :: xstore,ystore,zstore
DOUBLE PRECISION, DIMENSION(num) :: fstorex,fstorey,fstorez
DOUBLE PRECISION, DIMENSION(3,10,3,10) :: D

xstore=x
ystore=y
zstore=z

! Calculate the dynamical matrix.
do i=nmat1,nmat2
   do j=1,3
      x(i)=xstore(i)+amove(1,j)*amovefac
      y(i)=ystore(i)+amove(2,j)*amovefac
      z(i)=zstore(i)+amove(3,j)*amovefac

      if(environmentalpimlog) then
         call conjgradpimaim
      else
         if(conjgradaimlog) call conjgradaim
         if(conjgradlog) call conjgrad
      endif
! Store the forces for the first move.
      fstorex=frrx
      fstorey=frry
      fstorez=frrz

      x(i)=xstore(i)-amove(1,j)*amovefac
      y(i)=ystore(i)-amove(2,j)*amovefac
      z(i)=zstore(i)-amove(3,j)*amovefac

      if(environmentalpimlog) then
         call conjgradpimaim
      else
         if(conjgradaimlog) call conjgradaim
         if(conjgradlog) call conjgrad
      endif
      do k=1,num
         fdiffx=frrx(k)-fstorex(k)
         fdiffy=frry(k)-fstorey(k)
         fdiffz=frrz(k)-fstorez(k)
         D(j,i,1,k)=0.5d0*fdiffx/amovefac
         D(j,i,2,k)=0.5d0*fdiffy/amovefac
         D(j,i,3,k)=0.5d0*fdiffz/amovefac
      enddo   
   x(i)=xstore(i)
   y(i)=ystore(i)
   z(i)=zstore(i)
   enddo   
enddo   

!---> Parallelization_S
if( iam .eq. 0 ) then

open(50,file='dynmat.dat',status='new',form='unformatted')
write(50)D
close(50)

endif
!<--- Parallelization_E

return
END SUBROUTINE
