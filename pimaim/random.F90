!
!     It is recommended that the COMMON /RANSET/ .... statement should
!     be included in the user's main program unit.
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                      c
!     *** RANIN ***                                                    c
!                                                                      c
!     Initialises the random number generator using two integer seeds  c
!     IJ and KL,  where 0 <= IJ <= 31328  and  0 <= KL <= 30081.       c
!     It is only called once in a normal run, but is not called if     c
!     a run is restarted using RANRST.                                 c
!                                                                      c
!     Permutations of these seeds give 900 million sequences of        c
!     pseudo-random values, each of length 10**30, and each of which   c
!     is said to be independent of the others.                         c
!                                                                      c
!     This random number generator is due to Marsaglia and Zaman(1987) c
!     with slight modifications by James(1988).                        c
!     It is portable across most 32-bit machines.                      c
!                                                                      c     
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ranin(ij,kl)
!     .. Scalar Arguments ..
      integer ij,kl
!     ..
!     .. Scalars in Common ..
      real cc,cd,cm
      integer i97,j97
!     ..
!     .. Arrays in Common ..
      real cu(97)
!     ..
!     .. Local Scalars ..
      real s,t
      integer i,ii,j,jj,k,l,m
!     ..
!     .. Intrinsic Functions ..
      intrinsic mod
!     ..
!     .. Common blocks ..
      common /ranset/cu,cc,cd,cm,i97,j97
!     ..
      i = mod(ij/177,177) + 2
      j = mod(ij,177) + 2
      k = mod(kl/169,178) + 1
      l = mod(kl,169)
      do 2 ii = 1,97
          s = 0.0
          t = 0.5
          do 3 jj = 1,24
              m = mod(mod(i*j,179)*k,179)
              i = j
              j = k
              k = m
              l = mod(53*l+1,169)
              if (mod(l*m,64).ge.32) s = s + t
    3     t = 0.5*t
    2 cu(ii) = s
      cc = 362436./16777216.
      cd = 7654321./16777216.
      cm = 16777213./16777216.
      i97 = 97
      j97 = 33
      return

      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                      c
!     *** RANSVU ***                                                   c
!                                                                      c
!     Creates a real single value SVU from a uniform distribution,     c
!     where 0.0 <= SVU < 1.0.                                          c
!                                                                      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ransvu(svu)
!     .. Scalar Arguments ..
      real svu
!     ..
!     .. Scalars in Common ..
      real cc,cd,cm
      integer i97,j97
!     ..
!     .. Arrays in Common ..
      real cu(97)
!     ..
!     .. Common blocks ..
      common /ranset/cu,cc,cd,cm,i97,j97
!     ..
      svu = cu(i97) - cu(j97)
      if (svu.lt.0.0) svu = svu + 1.0
      cu(i97) = svu
      i97 = i97 - 1
      if (i97.eq.0) i97 = 97
      j97 = j97 - 1
      if (j97.eq.0) j97 = 97
      cc = cc - cd
      if (cc.lt.0.0) cc = cc + cm
      svu = svu - cc
      if (svu.lt.0.0) svu = svu + 1.0
      return

      end
