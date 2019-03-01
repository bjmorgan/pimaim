!-------------------------------------------------------------------------------------
!
! GWW - module for error function routines implimenting look up table 
!   17/05/2017    first version     
!   30/06/2017    second version  
!                 implimentation of more accurate erfc functions and intrinsic 
!
Module error_function

    implicit none
    double precision, allocatable :: erfc_table(:) 
    double precision :: del_r = 0.00025
    double precision :: recip_del_r 

contains 

!-------------------------------------------------------------------------------------
!
  subroutine construct_erfc_table 
! subroutine to contruct complemnetart error function table for real space ewald 

    implicit none
    double precision :: r 
    integer  ::  maxtable, i 

! GWW code is absolutely fucking stupid - calculating error function for all
! pairs irrespecive of the distance - so calculted for pairs well outside the 
! requirements ! 
! REALLY STUPID 
! For now just make table really large ! 
    r = 0.0

! used to get data back from table - best to do it only once 
    recip_del_r = 1 / del_r 

! should use eta and cutoff ( eta * rcut) - 8.0 should be sufficient 
    maxtable =  int (132.0 / del_r ) 
!    maxtable =  100000 
    allocate (erfc_table(maxtable)) 
    do i = 1, maxtable 
      r = r + del_r
      erfc_table(i) = gerfc(r)

! intrinsic function - looks good with intel - but will be compiler dependent
!      erfc_table(i) = erfc(r)

! origianl erfc function - not very accurate 
!      erfc_table(i) = old_erffunc(r)

    enddo 

  end subroutine 

!-------------------------------------------------------------------------------------
!
! function to provide erfc from table 

  function erfunc(x) 
    implicit none
    double precision :: x, r1, r2, erfunc 
    integer :: ir

    r1 = x * recip_del_r  
    ir = int (r1)     
    r2 = r1 - float(ir) 

    erfunc = erfc_table(ir) - r2*(erfc_table(ir) - erfc_table (ir+1)) 

! alternative use direct calls to erfc function - intrinsic, gerfc, original
!    erfunc = erfc(x) 
!    erfunc = gerfc(x) 
!    erfunc = old_erfunc(x) 

  end function erfunc

!-------------------------------------------------------------------------------------
!
! More accurate function to calculate the error function Calculates complementary error function
! named changed to gerfc to avoid intrisinc fortran finction -which is good ! (intel) 

  function gerfc(x)
  implicit none
  double precision :: gerfc, x 
  integer :: j
  double precision  :: z,sd,sn,xsq,expn,factor
  double precision,save :: p(8),q(8),a(10)
 
  data factor/1.128379167095512d+0/, &
   p/883.47894260850d+0,1549.6793124037d+0,1347.1941340976d+0, &
   723.04000277753d+0,255.50049469496d+0,59.240010112914d+0, &
   8.3765310814197d+0,0.56418955944261d+0/,q/883.47894260850d+0, &
   2546.5785458098d+0,3337.2213699893d+0,2606.7120152651d+0, &
   1333.5699756800d+0,460.28512369160d+0,105.50025439769d+0, &
   14.847012237523d+0/
  data a/ &
   0.10000000000000d+01,-0.33333333333333d+00, 0.10000000000000d+00, &
  -0.23809523809524d-01, 0.46296296296296d-02,-0.75757575757576d-03, &
   0.10683760683761d-03,-0.13227513227513d-04, 0.14589169000934d-05, &
  -0.14503852223151d-06/
!
  z = abs(x)
  xsq = z*z
  if (z.gt.8.0)then
    gerfc = 1.00d0 - sign(1.00d0,x)
  else if (z.gt.0.47) then
    expn = exp(-xsq)
    sn = p(8)
    sd = z + q(8)
    do j = 2,8
      sn = sn*z + p(9-j)
      sd = sd*z + q(9-j)
    enddo
    gerfc = 1.00d0 - sign((1.00d0-sn*expn/sd),x)
  elseif (z.gt.1.0d-15) then
    gerfc = a(10)
    do j = 1,9
      gerfc = gerfc*xsq + a(10-j)
    enddo
    gerfc = 1.00d0 - sign(factor*z*gerfc,x)
  else
    gerfc = 1.00d0
  endif
  return
  end function gerfc 

!-------------------------------------------------------------------------------------
!
! original erfc function - not very accurate 
! 
  FUNCTION old_erfunc(x)

    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: x
    DOUBLE PRECISION :: old_erfunc

      old_erfunc=1.0d0/((1.0d0+x*(0.0705230745d0+x*(0.0422820123d0        &
     +x*(0.0092705272d0+x*(0.0001520143d0+x*(0.0002765672d0+x*      &
     0.0000430638d0))))))**16.0d0)

    return
  END FUNCTION

END MODULE 
