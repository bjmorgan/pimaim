!
! Calls correct energy routine for CG minimisations.
!
SUBROUTINE cgener

USE commondata, ONLY: quadpimlog,cluslog,chargequadlog,dipquadlog, &
                      fullewaldlog,dippimlog,eta,sqrpi,exx,eyy,ezz, &
                      exy,exz,eyz,quadxx,quadyy,quadzz,quadxy,quadxz, &
                      quadyz,elecx,elecy,elecz,xmu,ymu,zmu,elecxq, &
                      elecyq,eleczq,exxq,eyyq,ezzq,exyq,exzq,eyzq,num
USE mpipara

IMPLICIT NONE
double precision :: fac2,fac3

if(quadpimlog) then
   if(cluslog) then
      call cgrealE_clus_quadpim
   else
      if(chargequadlog) then
         call cgrecipE_chargequad
         call cgrealE_chargequad
      endif
      if(dipquadlog) then
         call cgrecipE_dipquad
         call cgrealE_dipquad
      endif
      if(fullewaldlog) then
         call cgfullrecip
         call cgfullreal
      endif

      CALL MPI_ALLREDUCE(eltmp,elltmp,9*num,mpi_double_precision,MPI_SUM, &
                         MPI_COMM_WORLD,ierr)

      elecx=elecxq+elltmp(1:num)
      elecy=elecyq+elltmp(num+1:2*num)
      elecz=eleczq+elltmp(2*num+1:3*num)
      exx=exxq+elltmp(3*num+1:4*num)
      eyy=eyyq+elltmp(4*num+1:5*num)
      ezz=ezzq+elltmp(5*num+1:6*num)
      exy=exyq+elltmp(6*num+1:7*num)
      exz=exzq+elltmp(7*num+1:8*num)
      eyz=eyzq+elltmp(8*num+1:9*num)

!............correct the electric field for the mu-mu self-interaction
      fac2=4.0d0*(eta**3.0d0)/(3.0d0*sqrpi)
      elecx=elecx+fac2*xmu
      elecy=elecy+fac2*ymu
      elecz=elecz+fac2*zmu

!...correct the electric field gradient for the quad-quad
!... self-interaction (AGUADO) CHECK...
!... This term does not influence the trace of the field gradient.
      fac3=16.0d0*(eta**5.0d0)/(15.0d0*sqrpi)
      exx=exx+fac3*quadxx
      eyy=eyy+fac3*quadyy
      ezz=ezz+fac3*quadzz
      exy=exy+fac3*quadxy
      exz=exz+fac3*quadxz
      eyz=eyz+fac3*quadyz
   endif
endif

if(dippimlog) then
   if(cluslog) then
      call cgrealE_clus_dippim
   else
      call cgrecipE_dippim
      call cgrealE_dippim

      CALL MPI_ALLREDUCE(eltmp,elltmp,3*num,mpi_double_precision,MPI_SUM, &
                         MPI_COMM_WORLD,ierr)

      elecx=elecxq+elltmp(1:num)
      elecy=elecyq+elltmp(num+1:2*num)
      elecz=eleczq+elltmp(2*num+1:3*num)

!............correct the electric field for the mu-mu self-interaction
      fac2=4.0d0*(eta**3.0d0)/(3.0d0*sqrpi)
      elecx=elecx+fac2*xmu
      elecy=elecy+fac2*ymu
      elecz=elecz+fac2*zmu

   endif
endif

return
END SUBROUTINE
