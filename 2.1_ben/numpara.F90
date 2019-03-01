!---> Parallelization_S
SUBROUTINE numpara

USE commondata, ONLY:num,num_s,num_s2,num_s3,num_s4,nsp,nanion
use mpipara

IMPLICIT NONE

INTEGER :: i,j,numcnt,numave,nmod,stnum,ednum,istw,iedw

allocate( ist(num),  ied(num) )
allocate( ist2(num), ied2(num) )
allocate( ist3(num), ied3(num) )
allocate( ist4(num), ied4(num) )

!MPIAGUADO--> This looks like distributing the total number of ion pairs
!          to the different nodes. If you have, say, 100 ion pairs and
!          10 nodes, stnum and ednum will define intervals [1-10],
!          [11,20],...,[91,100], depending on the iam value (the rank)
!          So "st" and "ed" should stand for "start" and "end", respectively.
numcnt = 0
do j = 2, num
   do i = 1, j-1
      numcnt = numcnt + 1
   enddo
enddo

numave = numcnt / nprocs
stnum = numave*iam + 1
ednum = numave*( iam + 1 )
!<---MPIAGUADO

!MPIAGUADO--> If the number of ion pairs is not an exact multiple of the
!             number of nodes, the prior construction does not cover all
!             ion pairs. The following part fixes this problem. 
!             For example, if you have 100 ion pairs and 8 processors,
!             nmod=100-12*8=4. Then the first 4 processors are each assigned
!             one more ion pair, while the 4 last processors continue to
!             have associated 12 ion pairs each. 
!             We will have some load inbalance... 

nmod = mod( numcnt, nprocs )

if( nmod .ne. 0 ) then
   if( iam .eq. 0 ) then
      ednum = ednum + 1
   else
      if( iam + 1 .gt. nmod ) then
         stnum = stnum + nmod
         ednum = ednum + nmod
      else
         stnum = stnum + iam
         ednum = ednum + iam + 1
      endif
   endif
endif
!<---MPIAGUADO

num_s = stnum - 1

!-->MPIAGUADO--> Determines (jst,ist(jst)) and (jed,ied(jed)). These are the
!             initial and final (j,i) couple of indexes associated to 
!             each node for the distribution of work in the real space
!             routines...
!             In the example above (100 ion pairs and 8 processors), the
!             node with rank equal to zero would take care of those ion
!             pairs from (2,1) to (6,3), the rank=1 node from (6,4) to
!             (8,5), etc...
!             For j values not corresponding to any of the jst(ed) values,
!             the array elements ist(j) are simply set to one and the ied(j)
!             to (j-1). Thus, in our example, ist(3)=1 and ied(4)=3.
!             These elements are needed for the next part of code. Simply,
!             when going from j=2 to j=6 for the rank=0 node, the index i
!             may freely vary when j=3,4,5. This is the reason for the minimum
!             i value being one and the maximum being (j-1)...
numcnt = 0
do j = 2, num
   istw = 0
   iedw = 0

   do i = 1, j-1
      numcnt = numcnt + 1
      if( numcnt .eq. stnum ) then
         jst = j
         istw= i
      endif
      if( numcnt .eq. ednum ) then
         jed = j
         iedw= i
      endif
   enddo

   if( istw .eq. 0 ) then
      ist(j) = 1
   else
      ist(j) = istw
   endif

   if( iedw .eq. 0 ) then
      ied(j) = j-1
   else
      ied(j) = iedw
   endif
enddo
!<--MPIAGUADO

!MPIAGUADO--> Repeates the whole story only for the anion-cation
!             interactions. Defines in this way corresponding couples
!             (jst2,ist2) and (jed2,ied2)...
numcnt = 0
do j = nsp(1)+1, num
   do i = 1, nsp(1)
      numcnt = numcnt + 1
   enddo
enddo

numave = numcnt / nprocs
stnum = numave*iam + 1
ednum = numave*( iam + 1 )

nmod = mod( numcnt, nprocs )

if( nmod .ne. 0 ) then
   if( iam .eq. 0 ) then
      ednum = ednum + 1
   else
      if( iam + 1 .gt. nmod ) then
         stnum = stnum + nmod
         ednum = ednum + nmod
      else
         stnum = stnum + iam
         ednum = ednum + iam + 1
      endif
   endif
endif

num_s2 = stnum - 1

numcnt = 0
do j = nsp(1)+1, num
   istw = 0
   iedw = 0

   do i = 1, nsp(1)
      numcnt = numcnt + 1
      if( numcnt .eq. stnum ) then
         jst2 = j
         istw = i
      endif
      if( numcnt .eq. ednum ) then
         jed2 = j
         iedw = i
      endif
   enddo

   if( istw .eq. 0 ) then
      ist2(j) = 1
   else
      ist2(j) = istw
   endif

   if( iedw .eq. 0 ) then
      ied2(j) = nsp(1)
   else
      ied2(j) = iedw
   endif

enddo
!<--MPIAGUADO

!MPIAGUADO--> And now the same for the anion-anion pairs...
!             leads to the calculation of jst3,jed3,ist3,ied3
numcnt = 0
do j = 2, nanion
   do i = 1, j-1
      numcnt = numcnt + 1
   enddo
enddo

numave = numcnt / nprocs
stnum = numave*iam + 1
ednum = numave*( iam + 1 )

nmod = mod( numcnt, nprocs )

if( nmod .ne. 0 ) then
   if( iam .eq. 0 ) then
      ednum = ednum + 1
   else
      if( iam + 1 .gt. nmod ) then
         stnum = stnum + nmod
         ednum = ednum + nmod
      else
         stnum = stnum + iam
         ednum = ednum + iam + 1
      endif
   endif
endif

num_s3 = stnum - 1

numcnt = 0
do j = 2, nanion
   istw = 0
   iedw = 0

   do i = 1, j-1
      numcnt = numcnt + 1
      if( numcnt .eq. stnum ) then
         jst3 = j
         istw = i
      endif
      if( numcnt .eq. ednum ) then
         jed3 = j
         iedw = i
      endif
   enddo

   if( istw .eq. 0 ) then
      ist3(j) = 1
   else
      ist3(j) = istw
   endif

   if( iedw .eq. 0 ) then
      ied3(j) = j-1
   else
      ied3(j) = iedw
   endif

enddo
!<--MPIAGUADO

!MPIAGUADO--> And now the same for the cation-cation pairs...
!             leads to the calculation of jst4,jed4,ist4,ied4
numcnt = 0
do j = nanion + 2, num
   do i = 1 + nanion, j-1
      numcnt = numcnt + 1
   enddo
enddo

numave = numcnt / nprocs
stnum = numave*iam + 1
ednum = numave*( iam + 1 )

nmod = mod( numcnt, nprocs )

if( nmod .ne. 0 ) then
   if( iam .eq. 0 ) then
      ednum = ednum + 1
   else
      if( iam + 1 .gt. nmod ) then
         stnum = stnum + nmod
         ednum = ednum + nmod
      else
         stnum = stnum + iam
         ednum = ednum + iam + 1
      endif
   endif
endif

num_s4 = stnum - 1

numcnt = 0
do j = nanion + 2, num
   istw = 0
   iedw = 0

   do i = 1 + nanion, j-1
      numcnt = numcnt + 1
      if( numcnt .eq. stnum ) then
         jst4 = j
         istw = i
      endif
      if( numcnt .eq. ednum ) then
         jed4 = j
         iedw = i
      endif
   enddo

   if( istw .eq. 0 ) then
      ist4(j) = 1 + nanion
   else
      ist4(j) = istw
   endif

   if( iedw .eq. 0 ) then
      ied4(j) = j-1
   else
      ied4(j) = iedw
   endif
!<--MPIAGUADO

enddo

return
END SUBROUTINE
!<--- Parallelization_E
