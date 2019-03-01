subroutine lightcfdump

!*********************************************************************
!
!  Write final correlation function data to file.
!
!************************************************************

use commondata, only: nanion,nspec
use lightdata, only: ncorr,normav,ncfmat2,lnorm, sc,scf,cohav,scfiso

!---> Parallelization_S
use mpipara
!<--- Parallelization_E
implicit none




! *** define our local variables
double precision :: SUM,DUM,DENOM,check
double precision :: THIRD,SIXTH,FOURTH
integer :: i,j,k,l,INDEX,M,I2,J2,NI,NJ,NK,NL
! *** end of local definitions
double precision, dimension(4,4) :: rinpt,rout,rtest
double precision, dimension(6,6) :: pr, acf, ainv, acc
double precision, dimension(3,3) :: threein,threeout
allocate (sc(7,3,0:ncorr))
  
!---> Parallelization_S
if (iam .eq. 0 ) then

open(9,file='lightcf-average.out',status='new')
do i=1,6

!==========write out the anisotropic (k=0) cfs
   do j=1,10

     do k=0,ncorr-1
        if (lnorm(i,j,k).gt.0) then 
           write((600+(10*(j-1))+i),*)k,scf(i,j,k)/float(lnorm(i,j,k))

        endif

      enddo
   enddo

       
    sum=0.0d0
    do j=1,4
       write(9,*)i,j,cohav(i,j)/float(normav)   &
                ,cohav(i,j)/float(normav*nanion)
       sum=sum+cohav(i,j)/float(normav)
    enddo
      
enddo

endif
!<--- Parallelization_E

do k=1,ncfmat2
   index=1
   do i=1,6
      do j=i,6

         do m=0,ncorr-1

            if (lnorm(1,1,m).gt.0) then 

                scfiso(index,k,m)=scfiso(index,k,m)/float(lnorm(1,1,m))
 
             endif

         enddo
         pr(i,j)=scfiso(index,k,0)
         index=index+1
       enddo
    enddo

    do i=1,5
       do j=i+1,6
          pr(j,i)=pr(i,j)
       enddo
    enddo

!---> Parallelization_S
    if (iam .eq. 0 ) then

    do i=1,6
       do j=1,6
          write(9,*) k,i,j,pr(i,j)
       enddo
    enddo

    endif
!<--- Parallelization_E

    do i=1,6
       do j=1,6
          ainv(i,j)=0.0d0
       enddo
    enddo

    if(nspec.eq.4) then
!++++++++++++++++ FOUR SPECIES SECTION
      do i=1,4
         i2=i+2
         do j=1,4
            j2=j+2
            rinpt(i,j)=pr(i2,j2)      
         enddo
!---> Parallelization_S
         if (iam .eq. 0 ) then

         write(9,*)rinpt(i,1),rinpt(i,2),rinpt(i,3),rinpt(i,4)

         endif
!<--- Parallelization_E
      enddo

      call fourinvert(rinpt,rout,4)
!      rtest=matmul(rinpt,rout)
      do i=1,4
         i2=i+2
         do j=1,4
            j2=j+2
            ainv(i2,j2)=rout(i,j)
         enddo
!---> Parallelization_S
         if (iam .eq. 0 ) then

         write(9,*)rout(i,1),rout(i,2),rout(i,3),rout(i,4)
!         write(9,*)rtest(i,1),rtest(i,2),rtest(i,3),rtest(i,4)

         endif
!<--- Parallelization_E
      enddo

   else if(nspec.eq.3) then

!++++++++++++++++ THREE SPECIES SECTION
       do i=1,3
          i2=i+2
          do j=1,3
             j2=j+2
             threein(i,j)=pr(i2,j2)
          enddo
!---> Parallelization_S
          if (iam .eq. 0 ) then

          write(9,*)threein(i,1),threein(i,2),threein(i,3)

          endif
!<--- Parallelization_E
       enddo

!+++++ "dum" contains the determinant, is it OK that we don't keep it?
       call invert(threein,threeout)

       do i=1,3
          i2=i+2
          do j=1,3
             j2=j+2
             ainv(i2,j2)=threeout(i,j)
             check = 0.0d0
             do l=1,3
                check=check+threein(i,l)*threeout(l,j)
             enddo
!---> Parallelization_S
             if (iam .eq. 0 ) then

             write(9,*) check      

             endif
!<--- Parallelization_E
          enddo
!---> Parallelization_S
          if (iam .eq. 0 ) then

          write(9,*)threeout(i,1),threeout(i,2),threeout(i,3)

          endif
!<--- Parallelization_E
       enddo

     else
!+++++++++++TWO SPECIES SECTION
          denom=pr(3,3)*pr(4,4)-pr(3,4)*pr(4,3)
          ainv(3,3)=pr(4,4)/denom
          ainv(4,4)=pr(3,3)/denom
          ainv(3,4)=-pr(4,3)/denom
          ainv(4,3)=ainv(3,4)
     endif

!     Write(6,*) 'before 150 loop'
   
     do m=0,ncorr-1
        index=1
        do i=1,6
           do j=i,6
              acf(i,j)=scfiso(index,k,m)
              index=index+1
           enddo
        enddo

        do i=1,5
           do j=i+1,6
              acf(j,i)=acf(i,j)
           enddo
        enddo
       
        do i=1,2
           do j=i,2
              acc(i,j)=acf(i,j)

              do ni=1,6
                 do nj=1,6
                    acc(i,j)=acc(i,j) -                 &
                         acf(i,ni)*ainv(ni,nj)*pr(nj,j) &
                        -acf(j,ni)*ainv(ni,nj)*pr(nj,i)
                    do nk=1,6
                       do nl=1,6
                          acc(i,j)=acc(i,j) +           &
         acf(ni,nk)*ainv(nk,nl)*pr(nl,j)*ainv(ni,nj)*pr(nj,i)
                       enddo
                    enddo
                 enddo
               enddo
            enddo
         enddo

!............using the scfiso's to store the projected cfs is just  
!..................to save space
         scfiso(3,k,m)=acc(1,1)
         scfiso(4,k,m)=acc(1,2)
         scfiso(5,k,m)=acc(2,2)
         scfiso(8,k,m)=scfiso(12,k,m)
         scfiso(10,k,m)=scfiso(16,k,m)+2.0d0*scfiso(17,k,m)    &
                        +2.0d0*scfiso(18,k,m)+scfiso(19,k,m)   &
                        +2.0d0*scfiso(20,k,m)+scfiso(21,k,m)
!---> Parallelization_S
         if (iam .eq. 0 ) then

         write(8,*)m,scfiso(16,k,m),scfiso(19,k,m),scfiso(21,k,m)
         write(8,*)m,scfiso(17,k,m),scfiso(18,k,m),scfiso(20,k,m)

         endif
!<--- Parallelization_E
     enddo
enddo

 
!     write(6,*) 'before 630 loop'

!.............now collect together the correlation functions for each k 
!  magnitude

do m=0,ncorr-1
   do k=1,3
      do j=1,7

         sc(j,k,m)=0.0d0
      enddo
   enddo

!................the (1,0,0)
   third=1.0d0/3.0d0

!       write(6,*) 'before 660 loop'
   do k=1,3

      sc(1,1,m)=sc(1,1,m)+scfiso(3,k,m)*third

      sc(2,1,m)=sc(2,1,m)+2.0d0*scfiso(4,k,m)*third

      sc(3,1,m)=sc(3,1,m)+scfiso(5,k,m)*third

      sc(4,1,m)=sc(4,1,m)+scfiso(8,k,m)*third

      sc(5,1,m)=sc(5,1,m)+scfiso(10,k,m)*third

      sc(6,1,m)=sc(6,1,m)+scfiso(1,k,m)*third

      sc(7,1,m)=sc(7,1,m)+scfiso(7,k,m)*third

   enddo

!...............the (1,1,0)
        sixth=1.0d0/6.0d0

!       write(6,*) 'before 670 loop'

   do k=4,9

      sc(1,2,m)=sc(1,2,m)+scfiso(3,k,m)*sixth

      sc(2,2,m)=sc(2,2,m)+2.0d0*scfiso(4,k,m)*sixth

      sc(3,2,m)=sc(3,2,m)+scfiso(5,k,m)*sixth

      sc(4,2,m)=sc(4,2,m)+scfiso(8,k,m)*sixth

      sc(5,2,m)=sc(5,2,m)+scfiso(10,k,m)*sixth

      sc(6,2,m)=sc(6,2,m)+scfiso(1,k,m)*sixth

      sc(7,2,m)=sc(7,2,m)+scfiso(7,k,m)*sixth

   enddo

!...........the (1,1,1)
fourth=1.0d0/4.0d0


   do k=10,13

      sc(1,3,m)=sc(1,3,m)+scfiso(3,k,m)*fourth

      sc(2,3,m)=sc(2,3,m)+2.0d0*scfiso(4,k,m)*fourth

      sc(3,3,m)=sc(3,3,m)+scfiso(5,k,m)*fourth

      sc(4,3,m)=sc(4,3,m)+scfiso(8,k,m)*fourth

      sc(5,3,m)=sc(5,3,m)+scfiso(10,k,m)*fourth

      sc(6,3,m)=sc(6,3,m)+scfiso(1,k,m)*fourth

      sc(7,3,m)=sc(7,3,m)+scfiso(7,k,m)*fourth

   enddo
enddo

!---> Parallelization_S
if (iam .eq. 0 ) then

write(6,*) 'about to write out light scattering'

!     Open our light scattering files
open(11,file='lightcf-11.out',status='new')
open(12,file='lightcf-12.out',status='new')
open(13,file='lightcf-13.out',status='new')
open(14,file='lightcf-14.out',status='new')
open(15,file='lightcf-15.out',status='new')
open(16,file='lightcf-16.out',status='new')
open(17,file='lightcf-17.out',status='new')

! SR projected to 11
! (sr x gamma projected) to 12
! gamma-gamma projected  to 13
!   S--(k,w) to 14
!    S++ to 15
!SR unprojected to 16
!gamma-gamma unprojected to 17

 
do m=0,ncorr-1
   write(11,710) m, sc(1,1,m),sc(1,2,m),sc(1,3,m)
   write(12,710) m, sc(2,1,m),sc(2,2,m),sc(2,3,m)
   write(13,710) m, sc(3,1,m),sc(3,2,m),sc(3,3,m)
   write(14,710) m, sc(4,1,m),sc(4,2,m),sc(4,3,m)
   write(15,710) m, sc(5,1,m),sc(5,2,m),sc(5,3,m)
   write(16,710) m, sc(6,1,m),sc(6,2,m),sc(6,3,m)
   write(17,710) m, sc(7,1,m),sc(7,2,m),sc(7,3,m)
enddo

710 format(I6,3E15.5)

!     Close light scattering files
close(9)
close(11)
close(12)
close(13)
close(14)
close(15)
close(16)
close(17)

endif
!<--- Parallelization_E
 
return

end subroutine



