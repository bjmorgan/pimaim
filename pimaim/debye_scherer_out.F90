SUBROUTINE debye_scherer_out

USE commondata, ONLY: nspec,nsp,nrun, ntotstp
USE recipdata, ONLY: rksqmax,norm_ds,nbin,sk_ds,slens

IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(nspec,nspec) :: xijrec
DOUBLE PRECISION :: runrec,fnbin,rk,binwidth,count,skk,totsk(2000),avslens
INTEGER :: i,j,lk,ii,jj,ibin,icount
CHARACTER(len=9) indi
CHARACTER(len=6) skname
CHARACTER(len=30) filename

indi = '123456789'
skname = 'sk.out'

do i=1, nspec
   avslens=avslens+slens(i) !Modified by D. Marrocchelli 11/03/2008
enddo

do i=1,nspec
   do j=i,nspec
      xijrec(i,j)=1.0d0/dsqrt(dble((nsp(i))*(nsp(j))))
   enddo   
enddo   

runrec=1.0d0/dble(nrun-ntotstp)   !Modified by D. Marrocchelli 11/03/2008
fnbin=float(nbin)

lk=1
do ii=1,nspec
   do jj=ii,nspec
#ifndef ppfit
      filename=skname//indi(ii:ii)//indi(jj:jj)
      open(65,file=filename,status='new')
#endif
         do ibin=1,nbin
            rk=dsqrt(rksqmax*(float(ibin)+0.5d0)/fnbin)
            binwidth=dsqrt(float(ibin+1)*rksqmax)- &
                     dsqrt(float(ibin)*rksqmax)
            binwidth=rk**2*binwidth
            count=float(norm_ds(ibin))
!.......NB the 5 here is arbitrary  -- to let in only those
!       bins which are appreciably sampled
!           icount=int(5.0d0*count*runrec)
            icount=int(count*runrec)

            if(icount.ne.0) then
!              skk=sk_ds(lk,ibin)*xijrec(ii,jj)*runrec/binwidth
               skk=sk_ds(lk,ibin)*xijrec(ii,jj)/count
               totsk(ibin)=totsk(ibin)+skk

!Modified by D.Marrocchelli 07/03/2008
           
#ifndef ppfit
               if((ii.eq.nspec).and.(jj.eq.nspec)) then  
                  write(75,*)rk, totsk(ibin)
               endif
               write(65,*)rk,skk 
#endif
	     
            endif
         enddo    
#ifndef ppfit
      close(65)
#endif
      lk=lk+1
   enddo   
enddo 
#ifndef ppfit  
close(75)
#endif
return
END SUBROUTINE
