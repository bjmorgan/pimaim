subroutine lightcfcalc



!............modified for improved isotropic spectrum calculation

!************************************************************
!
!  Auto- and cross- correlate polarizability fluctuations
! due to the 4 mechanisms (SR,DID,asym(B,gamma)).
!
!************************************************************

use lightdata, only: ncfmat2,normav,storeisokr,storeisoki, &
             cisokr,cisoki,scfiso,ncorrcall,ncorr, &
             storecohacc,cohacc,cohav,ncorrtime,lnorm,scf 
implicit none
! *** define our local variables
integer :: i,j,k,L,INDEX,M,NT,III,JJJ,KKK
! *** end of local definitions

!=====================================================================
!  Primary loop.
!  Format of index pointer:
!   __________________________________________
!   | index=1   |   auto    |    sr  |    sr  |   
!   | index=2   |  cross    |    sr  |   did  |   
!   | index=3   |  cross    |    sr  |     B  |   
!   | index=4   |  cross    |    sr  | gamma  |   
!   | index=5   |   auto    |   did  |   did  |   
!   | index=6   |  cross    |   did  |     B  |   
!   | index=7   |  cross    |   did  | gamma  |   
!   | index=8   |   auto    |     B  |     B  |   
!   | index=9   |  cross    |     B  | gamma  |   
!   | index=10  |   auto    | gamma  | gamma  |   
!   ===========================================

!=====================================================================
! Store coherent array in correlation function array.


do l=1,6
   do j=1,ncfmat2

      storeisokr(l,j,ncorrtime)=0.0d0
      storeisoki(l,j,ncorrtime)=0.0d0
   enddo
enddo

do i=1,6

   do j=1,4

      storecohacc(i,j,ncorrtime)=cohacc(i,j)
      cohav(i,j)=cohav(i,j)+cohacc(i,j)
   enddo
enddo

normav=normav+1


do l=1,6
   do j=1,ncfmat2

      storeisokr(l,j,ncorrtime)=cisokr(l,j)
      storeisoki(l,j,ncorrtime)=cisoki(l,j)
    enddo
enddo

!=====================================================================
!****************** Primary loop.
do i=1,6

   index=1

   do j=1,4

      do k=j,4

         do m=1,ncorrtime

            nt=ncorrtime-m
            lnorm(i,index,nt)=lnorm(i,index,nt)+1

            scf(i,index,nt)=scf(i,index,nt)  &
                           +storecohacc(i,j,m)*cohacc(i,k)  &
                           +storecohacc(i,k,m)*cohacc(i,j)

         enddo

         index=index+1

      enddo
   enddo
enddo

do k=1,ncfmat2
   index=1

   do l=1,6

      do j=l,6

         do m=1,ncorrtime

            nt=ncorrtime-m
            scfiso(index,k,nt)=scfiso(index,k,nt) &
                              +storeisokr(l,k,m)*cisokr(j,k) &
                              +storeisokr(j,k,m)*cisokr(l,k) &
                              +storeisoki(l,k,m)*cisoki(j,k) &
                              +storeisoki(j,k,m)*cisoki(l,k)

         enddo
         index=index+1
      enddo
   enddo
enddo
!=====================================================================
!*************** Secondary loop.
 
if (ncorrcall.gt.ncorr) then

   do i=1,6

      index=1

      do j=1,4

         do k=j,4

            do m=ncorrtime+1,ncorr

               nt=ncorrtime-m+ncorr

               lnorm(i,index,nt)=lnorm(i,index,nt)+1

               scf(i,index,nt)=scf(i,index,nt)               &
                            +storecohacc(i,j,m)*cohacc(i,k)  &
                            +storecohacc(i,k,m)*cohacc(i,j)

            enddo
 
            index=index+1

         enddo
       enddo
    enddo

    do k=1,ncfmat2
       
       index=1

       do l=1,6

          do j=l,6

             do m=ncorrtime+1,ncorr

                nt=ncorrtime-m+ncorr

                scfiso(index,k,nt)=scfiso(index,k,nt)       &
                            +storeisokr(l,k,m)*cisokr(j,k)  &
                            +storeisokr(j,k,m)*cisokr(l,k)  &
                            +storeisoki(l,k,m)*cisoki(j,k)  &
                            +storeisoki(j,k,m)*cisoki(l,k)

              enddo
 
              index=index+1

           enddo
        enddo
     enddo
!=====================================================================

endif
 
return

end subroutine

