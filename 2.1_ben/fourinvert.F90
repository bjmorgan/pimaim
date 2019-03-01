!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!    routine to invert a matrix by Gaussian elimination
!     Ainv=inverse(A)    A(n,n)   Ainv(n,n)
!     Note: sizes maxed at 4 x 4 - re dimension for larger matrices;
!           D matrix is extended.
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine fourinvert( A, Ainv, n )
!---> Parallelization_S
use mpipara
!<--- Parallelization_E
IMPLICIT NONE

!.........local variables
double precision, dimension(4,4):: A,Ainv
double precision, dimension(4,8) ::  D
double precision :: alpha,beta
integer :: n,n2,i,j,k

!  initialize the reduction matrix
n2 = 2*n
do i = 1,n
   do j = 1,n
      D(i,j) = A(i,j)
      d(i,n+j) = 0.0d0
   enddo
   D(i,n+i) = 1.0d0
enddo

!  do the reduction 
do i = 1,n
   alpha = D(i,i)
   if(alpha .ne. 0.) then
      do j = 1,n2
         D(i,j) = D(i,j)/alpha
      enddo

      do k = 1,n
         if((k-i).ne.0) then
            beta = D(k,i)
            do j = 1,n2
               D(k,j) = D(k,j) - beta*D(i,j)
            enddo
         endif
      enddo
   else 
!---> Parallelization_S
      if( iam .eq. 0 ) then

      write(6,*) 'Singular Matrix in Fourinvert'

      endif
!<--- Parallelization_E
      return
   endif
enddo
!  copy result into output matrix
do i = 1,n
   do j = 1,n
      Ainv(i,j) = D(i,j+n)
   enddo
enddo

return
end subroutine
