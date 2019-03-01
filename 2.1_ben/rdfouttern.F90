SUBROUTINE rdfout
 
USE commondata, ONLY: num,fourpi,nrdfcall,nspec,nsp,nionunit,rdftot,rdfpart
USE boxdata, ONLY: cellvol3rec,cellvol

IMPLICIT NONE

DOUBLE PRECISION :: const,drrdf,xtot,rlower,rupper,xnorm,xnormrec
DOUBLE PRECISION, DIMENSION(6) :: xion
INTEGER :: i,j,ii,jj,k

CHARACTER(len=9) indi
CHARACTER(len=7) rdfname
CHARACTER(len=60) filename

indi = '123456789'
rdfname = 'rdf.out'

const=dble(num)*fourpi*cellvol3rec
const=const*dble(nrdfcall*num)

drrdf=(cellvol**(1.0d0/3.0d0))/600.0d0

xtot=0.0d0
do i=1, nspec
xtot=float(nsp(i))/float(nionunit) + xtot
xion(i)=float(nsp(i))/float(nionunit)
enddo    

open(65,file='rdftot.out',status='new')
 
   do i=0,300
      rlower=dble(i)*drrdf
      rupper=rlower+drrdf
      xnorm=const*((rupper**3.0d0)-(rlower**3.0d0))
      xnormrec=1.0d0/xnorm
 
      rdftot(i)=rdftot(i)*xnormrec

      do ii=1,nspec      
         do jj=ii,nspec      
            rdfpart(i,ii,jj)=0.5d0*xtot*xtot/((xion(ii))*(xion(jj))) &
                   *rdfpart(i,ii,jj)/xnorm
            if(ii.eq.jj) then
               rdfpart(i,ii,jj)=rdfpart(i,ii,jj)*2.0d0
            endif
         enddo   
      enddo        
      write(65,*) rlower,rdftot(i)
   enddo   
close(65)

do i=1,nspec
   do j=i,nspec
      filename=rdfname//indi(i:i)//indi(j:j)
      open(65,file=filename,status='new')
         do k=0,300
            rlower=dble(k)*drrdf
            rlower=rlower+(0.5d0*drrdf)
            write(65,*) rlower,rdfpart(k,i,j)
         enddo   
      close(65)
   enddo   
enddo   

write(*,*)
write(*,*) '**** Radial distribution function written out ****'
write(*,*)

return
END SUBROUTINE
