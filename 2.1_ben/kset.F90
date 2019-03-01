SUBROUTINE kset

USE commondata, ONLY: twopi,etasq, eta, rcut 
USE boxdata, ONLY: fullhi,fullhit,bh
USE recipdata, ONLY: rksqmax,conv1,kmaxx,kmaxy,kmaxz,convfac

IMPLICIT NONE

INTEGER :: i,j,kmax
DOUBLE PRECISION :: rkmin
DOUBLE PRECISION :: rkmax 

double precision :: dl_tol, dl_tol1

fullhit=TRANSPOSE(fullhi)

call dcell(fullhit,bh)

rkmin=(twopi*dmin1(dabs(bh(7)),dabs(bh(8)),dabs(bh(9))))

! unit stupid - compared in Ewald to distance - but it is not ! 
!rksqmax=-log(conv1)*4.0d0*etasq/(rkmin*rkmin)
!rksqmax=-log(conv1)*4.0d0*etasq
rksqmax= -log(conv1)*4.0d0*etasq

! gulp
!rksqmax= -log(conv1)*4.0d0*eta


!rksqmax= rksqmax + 1.0 

! seems more sensible here
rksqmax=rksqmax*convfac
!rksqmax=rksqmax*1.2 

! why int ? 
!kmax=int(dsqrt(rksqmax))
rkmax=(dsqrt(rksqmax))

kmaxx=int(rkmax/(twopi*bh(7)))+1
kmaxy=int(rkmax/(twopi*bh(8)))+1
kmaxz=int(rkmax/(twopi*bh(9)))+1

! GWW check to see if this does not work for no-orthogonal
kmaxx = kmaxx + 10
kmaxy = kmaxy + 10
kmaxz = kmaxz + 10

!print *, 'old kmaxx', kmaxx, kmaxy, kmaxz 


! GWW try DLPOLY way
! note DLPOLY uses spme - so doubles the number of vectors 

!dl_tol  = sqrt ((-log (conv1 * rcut)))
!dl_tol1 = sqrt (-log (conv1 * rcut * (2 * dl_tol * eta)**2 ) ) 
!
!kmaxx = 2 * nint (0.25 + ( eta * dl_tol1 / (0.5 * twopi * bh(7))))
!kmaxy = 2 * nint (0.25 + ( eta * dl_tol1 / (0.5 * twopi * bh(8))))
!kmaxz = 2 * nint (0.25 + ( eta * dl_tol1 / (0.5 * twopi * bh(9))))


!kmaxx=int((rkmin/(twopi*bh(7)))*float(kmax))+1
!kmaxy=int((rkmin/(twopi*bh(8)))*float(kmax))+1
!kmaxz=int((rkmin/(twopi*bh(9)))*float(kmax))+1

! GW Stupid place - cut it back after x,y and z vectors calculated !
! rksqmax=rksqmax*convfac

#ifdef debug 
   print*, 'kmax x,y,z', kmaxx, kmaxy, kmaxz 
   print*, 'perpendiculars',bh(7), bh(8), bh(9) 
   print*, 'eta, etasq, conv1', eta, etasq, conv1 
   print*, 'rksqmax ', rksqmax                
#endif 

return
END SUBROUTINE
