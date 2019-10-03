FUNCTION DBRENT(AX,BX,CX,F,DF,TOL,XMINLOC)

use mpipara 
IMPLICIT NONE

INTEGER :: ITER
DOUBLE PRECISION, INTENT(OUT) :: XMINLOC
DOUBLE PRECISION :: AX,BX,CX,F,DF,TOL,DBRENT
DOUBLE PRECISION :: TOL1,TOL2,OLDE,A,B,D,E,U,V,W,X, &
                    D1,D2,DU,DV,DW,DX,FU,FV,FW,FX,U1,U2,XM
INTEGER, PARAMETER :: ITMAX=100
DOUBLE PRECISION, PARAMETER :: ZEPS=1.0d-10
LOGICAL OK1,OK2


A=MIN(AX,CX)
B=MAX(AX,CX)
V=BX
W=V
X=V
D=0.0d0
E=0.0d0
FX=F(X)
DX=DF(X)
FV=FX
FW=FX
DV=DX
DW=DX
DO ITER=1,ITMAX
   XM=0.5d0*(A+B)
! GW the closer to the min the smaller the tolerance ?????
!   TOL1=TOL*ABS(X)+ZEPS
   TOL1 = TOL + ZEPS 
   TOL2=2.0d0*TOL1
   IF(ABS(X-XM).LE.(TOL2-.5d0*(B-A))) EXIT
   IF(ABS(E).GT.TOL1) THEN
     D1=2.*(B-A)
     D2=D1
     IF(DW.NE.DX) D1=(W-X)*DX/(DX-DW)
     IF(DV.NE.DX) D2=(V-X)*DX/(DX-DV)
     U1=X+D1
     U2=X+D2
     OK1=((A-U1)*(U1-B).GT.0.0d0).AND.(DX*D1.LE.0.0d0)
     OK2=((A-U2)*(U2-B).GT.0.0d0).AND.(DX*D2.LE.0.0d0)
     OLDE=E
     E=D
     IF(.NOT.(OK1.OR.OK2))THEN
       GO TO 1
     ELSE IF (OK1.AND.OK2)THEN
       IF(ABS(D1).LT.ABS(D2))THEN
         D=D1
       ELSE
         D=D2
       ENDIF
     ELSE IF (OK1)THEN
       D=D1
     ELSE
       D=D2
     ENDIF
     IF(ABS(D).GT.ABS(0.5d0*OLDE))GO TO 1
     U=X+D
     IF(U-A.LT.TOL2 .OR. B-U.LT.TOL2) D=SIGN(TOL1,XM-X)
     GOTO 2
   ENDIF
1       IF(DX.GE.0.0d0) THEN
     E=A-X
   ELSE
     E=B-X
   ENDIF
   D=0.5*E
2       IF(ABS(D).GE.TOL1) THEN
     U=X+D
     FU=F(U)
   ELSE
     U=X+SIGN(TOL1,D)
     FU=F(U)
     IF(FU.GT.FX) EXIT
   ENDIF
    if (iam.eq.0)  WRITE(6,*) 'In Brent', U, FU  
   DU=DF(U)
   IF(FU.LE.FX) THEN
     IF(U.GE.X) THEN
       A=X
     ELSE
       B=X
     ENDIF
     V=W
     FV=FW
     DV=DW
     W=X
     FW=FX
     DW=DX
     X=U
     FX=FU
     DX=DU
   ELSE
     IF(U.LT.X) THEN
       A=U
     ELSE
       B=U
     ENDIF
     IF(FU.LE.FW .OR. W.EQ.X) THEN
       V=W
       FV=FW
       DV=DW
       W=U
       FW=FU
       DW=DU
     ELSE IF(FU.LE.FV .OR. V.EQ.X .OR. V.EQ.W) THEN
       V=U
       FV=FU
       DV=DU
     ENDIF
   ENDIF
ENDDO    
XMINLOC=X
DBRENT=FX
RETURN
END FUNCTION
