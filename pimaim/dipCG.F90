SUBROUTINE dipcg

USE commondata, ONLY: xmu,ymu,zmu,num,nspec,ntype,engeff,alppolar,alppolar1, &
                      alppolardel,alppolareff,elecx,elecy,elecz,ftol,xk1,nsp, &
                      reseng,delta,polarizablelog,verbose,elecxq,elecyq,eleczq
!---> Parallelization_S
use mpipara
!---> Parallelization_E

IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(num) :: xmusv,ymusv,zmusv,brhsx,brhsy,brhsz
DOUBLE PRECISION, DIMENSION(num) :: resx,resy,resz,pdirx,pdiry,pdirz
DOUBLE PRECISION, DIMENSION(num) :: Apx,Apy,Apz
DOUBLE PRECISION :: res2,res2old,cgalpha,cgbeta,pAp,dipsq
INTEGER :: i,j,ipoint,m,n
INTEGER, PARAMETER :: itmax=100
      
! the rhs of the equation to be solved are the 'charge induced' dipoles
brhsx=elecxq*alppolar
brhsy=elecyq*alppolar
brhsz=eleczq*alppolar
! first calculate the initial residuals
call cgener
resx=elecx*alppolar-xmu
resy=elecy*alppolar-ymu
resz=elecz*alppolar-zmu

if(verbose) then
  if( iam .eq. 0 ) then
    write (10,*) ' Initial residues' 
    do i = 1,num
      write (10,'(i6,1x,5(f15.8,2x))') i, resx(i), elecxq(i), elecx(i), alppolar(i), xmu(i) 
    enddo
  endif
endif 

do j=1,itmax
   res2=SUM(resx*resx+resy*resy+resz*resz)
!---> Parallelization_S
   if( iam .eq. 0 ) then

   if(verbose) write(6,'(a,i4,1x,a,f20.14,1x,a,1x,f20.14,1x,a,i4)')'j=',j,'convergence=',dsqrt(res2/real(num,8)),'ftl',ftol,'num=',num

   endif
!---> Parallelization_E
! can we stop now?      
   if(dsqrt(res2/num).lt.ftol) then
      call cgener
! update of contributions to polengtot

      m=1
      do n=1,nspec
         if(polarizablelog(n)) then
            xk1(m:nsp(n)+m-1)=1.0d0/(2.0d0*alppolar(m:nsp(n)+m-1))
            do i=m,nsp(n)+m-1
               dipsq=(xmu(i)*xmu(i))+(ymu(i)*ymu(i))+(zmu(i)*zmu(i))
               reseng(i)=dipsq*xk1(i)
            end do
         endif
         m=m+nsp(n)
      enddo
      return
   endif
! find the new search direction
   if(j.eq.1) then
      pdirx=resx
      pdiry=resy
      pdirz=resz
   else
      cgbeta=res2/res2old
      pdirx=resx+cgbeta*pdirx
      pdiry=resy+cgbeta*pdiry
      pdirz=resz+cgbeta*pdirz
   endif
! finc A.p      
   xmusv=xmu
   ymusv=ymu
   zmusv=zmu
   xmu=pdirx
   ymu=pdiry
   zmu=pdirz

   call cgener
      
   Apx=brhsx-(elecx*alppolar-xmu)
   Apy=brhsy-(elecy*alppolar-ymu)
   Apz=brhsz-(elecz*alppolar-zmu)
   xmu=xmusv
   ymu=ymusv
   zmu=zmusv
! now find alpha       
   pAp=SUM(pdirx*Apx+pdiry*Apy+pdirz*Apz)
   cgalpha=res2/pAp
! update dipoles and residuals
   xmu=xmu+cgalpha*pdirx
   ymu=ymu+cgalpha*pdiry
   zmu=zmu+cgalpha*pdirz
   resx=resx-cgalpha*Apx
   resy=resy-cgalpha*Apy
   resz=resz-cgalpha*Apz

   res2old=res2
enddo
!---> Parallelization_S
if( iam .eq. 0 ) then

    do i = 1,num
      write (10,'(i6,1x,3(f15.8,2x))') i, xmu(i), ymu(i), zmu(i)
    enddo
write(6,*)'cg failed to converge - stopping ' 

endif
!---> Parallelization_E
!---> Parallelization_S
call mpi_finalize(ierr)
!---> Parallelization_E
STOP
RETURN
END SUBROUTINE
