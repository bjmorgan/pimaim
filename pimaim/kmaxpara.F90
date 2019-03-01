!---> Parallelization_S
SUBROUTINE kmaxpara

USE recipdata, ONLY:kmaxx,kmaxy,kmaxz
use mpipara

IMPLICIT NONE
integer :: i,j,k

INTEGER :: ll,mm,nn,mmin,nmin,kmaxt,klmn,klmn_s,klmn_e,kmod,kcnt, &
              mm_s,mm_e,nn_s,nn_e

!MPIAGUADO--> This is similar to the numparam routine but now we
!             distribute k-points instead of ion pairs. Thus, kmaxt is
!             the total number of k-points, klmn the average number of
!             k-points per node, and kmod is used just in case the total
!             number of k-points is not an exact multiple of nprocs, in
!             which case one more k-point is assigned to each processor,
!             starting from the lower end rank=0, until all k-points are
!             distributed.
!             klmn_s and klmn_e denote the starting and end k-values for
!             each node.
kmaxt = kmaxx*( 2*kmaxy + 1 )*( 2*kmaxz + 1 ) + kmaxy*( 2*kmaxz + 1 ) + kmaxz
klmn  = kmaxt / nprocs

klmn_s = klmn*iam + 1
klmn_e = klmn*( iam + 1 )

kmod = mod( kmaxt, nprocs )
if( kmod .ne. 0 ) then
   if( iam .eq.  0 ) then
      klmn_e = klmn_e + 1
   else
      if( iam + 1 .gt. kmod ) then
         klmn_s = klmn_s + kmod
         klmn_e = klmn_e + kmod
      else
         klmn_s = klmn_s + iam
         klmn_e = klmn_e + iam + 1
      endif
   endif
endif
!<--MPIAGUADO

mmin=0
nmin=1
kcnt=0

!MPIAGUADO--> Here it just translates the information contained in klmn_s
!             and klmn_e into the starting and end values, for each processor,
!             of kmaxx, kmaxy, and kmaxz (the "max" here is confusing, but
!             anyway...).
!             Notice that, similarly to the real space part, when the value
!             of kx for a given processor does not coincide with kmaxx_s, then
!             kmaxy_s(kx)=-kmaxy, the minimum value. Also, when ks is not
!             equal to kmaxx_e, then kmaxy_e(kx)=kmaxy, the maximum value.
!             The same happens for the 2D arrays kmaxz_s and kmaxz_e if
!             either kx or ky do not coincide with an edge value.
do ll=0,kmaxx
   do mm=mmin,kmaxy
      do nn=nmin,kmaxz
         kcnt = kcnt + 1
         if( kcnt .eq. klmn_s ) then
            kmaxx_s = ll
            mm_s    = mm
            nn_s    = nn
         endif
         if( kcnt .eq. klmn_e ) then
            kmaxx_e = ll
            mm_e    = mm
            nn_e    = nn
         endif
      enddo
      nmin=-kmaxz
   enddo
   mmin=-kmaxy
enddo

do ll = kmaxx_s , kmaxx_e
   kmaxy_s(ll) = -kmaxy
   kmaxy_e(ll) =  kmaxy
enddo
kmaxy_s(kmaxx_s) = mm_s
kmaxy_e(kmaxx_e) = mm_e

do ll = kmaxx_s, kmaxx_e
   do mm = kmaxy_s(ll), kmaxy_e(ll)
      kmaxz_s(mm,ll) = -kmaxz
      kmaxz_e(mm,ll) =  kmaxz
   enddo
enddo
kmaxz_s(kmaxy_s(kmaxx_s),kmaxx_s) = nn_s
kmaxz_e(kmaxy_e(kmaxx_e),kmaxx_e) = nn_e
!<--MPIAGUADO

#ifdef debug 
   print*,"kmaxt ", iam, kmaxt,klmn, klmn_s, klmn_e 
   print*,"kmaxx ", iam, kmaxx,kmaxy,kmaxz
   print*,"kmaxx ", iam, kmaxx_s,kmaxx_e

!print*, 'check of vectors used' 
!
   do i=kmaxx_s, kmaxx_e 
     do j= kmaxy_s(i), kmaxy_e(i)
       do k = kmaxz_s(j,i), kmaxz_e(j,i)
         print*,'vector', i,j,k
       enddo 
     enddo
  enddo 

#endif 


return
END SUBROUTINE
!<--- Parallelization_E
