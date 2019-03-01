SUBROUTINE rdf

USE commondata, ONLY: nrdfcall,num,ntype,x,y,z,rdftot,rdfpart, &
!---> Memmory Reduction_S
!                     dxsav,dysav,dzsav
                      dxsav,dysav,dzsav,num2,numx,nspec
!<--- Memmory Reduction_E
USE boxdata, ONLY: cellvol
!---> Parallelization_S
use mpipara
!<--- Parallelization_E

IMPLICIT NONE
 
DOUBLE PRECISION :: dr,dx,dy,dz,dist
INTEGER :: i,ipoint,j,jpoint,nparttype,nbinrdf
DOUBLE PRECISION, DIMENSION(0:300,nspec,nspec) :: tmp
DOUBLE PRECISION, DIMENSION(0:300) :: ttmp

nrdfcall=nrdfcall+1

dr=(cellvol**(1.0d0/3.0d0))/600.0d0

!---> Parallelization_S
rdftot_w=0.0d0
rdfpart_w=0.0d0
!<--- Parallelization_E
!---> Parallelization_S
!do j=2,num
do j=jst,jed
!<--- Parallelization_E
   jpoint=ntype(j)
 
!---> Parallelization_S
!  do i=1,j-1
   do i=ist(j),ied(j)
!<--- Parallelization_E
      ipoint=ntype(i)
      nparttype=ipoint+jpoint-1

!---> Memmory Reduction_S
      numx = numadr(i,j)
!     dx=dxsav(i,j)
!     dy=dysav(i,j)
!     dz=dzsav(i,j)
      dx=dxsav(numx)
      dy=dysav(numx)
      dz=dzsav(numx)
!<--- Memmory Reduction_E

      dist=dx*dx+dy*dy+dz*dz
      dist=dsqrt(dist)
 
      nbinrdf=int((dist)/dr)

      if(nbinrdf.ge.0.and.nbinrdf.le.300) then
!---> Parallelization_S
!        rdftot(nbinrdf)=rdftot(nbinrdf)+2.0d0
         rdftot_w(nbinrdf)=rdftot_w(nbinrdf)+2.0d0

!        rdfpart(nbinrdf,ipoint,jpoint) &
!              =rdfpart(nbinrdf,ipoint,jpoint)+2.0d0
         rdfpart_w(nbinrdf,ipoint,jpoint) &
               =rdfpart_w(nbinrdf,ipoint,jpoint)+2.0d0
!<--- Parallelization_E
      endif
   enddo   
enddo   

!---> Parallelization_S
CALL MPI_ALLREDUCE(rdftot_w,ttmp,301,mpi_double_precision,MPI_SUM, &
                   MPI_COMM_WORLD,ierr)
rdftot = rdftot + ttmp

CALL MPI_ALLREDUCE(rdfpart_w,tmp,301*nspec*nspec,mpi_double_precision,MPI_SUM, &
                   MPI_COMM_WORLD,ierr)
rdfpart = rdfpart + tmp
!<--- Parallelization_E
return
END SUBROUTINE
