subroutine lightfields

!************************************************************
!************************************************************


use lightdata, only: exxls,eyyls,ezzls,exyls,exzls,eyzls,    &
                     txxli,tyyli,tzzli, &
                     txyli,txzli,tyzli,xlipol,sig,asr,bsr,csr,dsr, &
                     srxx,sryy,srzz,srxy,srxz,sryz
use commondata, only: num,ntype,dxsav,dysav,dzsav,rsqmax,q, &
!---> Memmory Reduction_S
!                     xmu,ymu,zmu
                      xmu,ymu,zmu,num2,numx
!<--- Memmory Reduction_E
!---> Optimize_S
use mpipara
!<--- Optimize_E
! *** define our local variables
implicit none
double precision :: DRSQ,DR,DRSQREC,DRREC,DR3REC,DR7REC,DXUNITSQ
double precision :: DYUNITSQ,DZUNITSQ,SXX,SYY,SZZ,SXY,SXZ,SYZ,TXXTT
double precision :: TYYTT,TZZTT,TXYTT,TXZTT,TYZTT,ASRTERM,BSRTERM
double precision :: XXCORR,YYCORR,ZZCORR,XYCORR,XZCORR,YZCORR
double precision :: THREEDR7,TXXX,TYYY,TZZZ,TXXY,TXYY,TXXZ,TXZZ
double precision :: TYYZ,TYZZ,TXYZ,EXXDIPI,EYYDIPI,EZZDIPI,EXYDIPI
double precision :: EXZDIPI,EYZDIPI,EXXDIPJ,EYYDIPJ,EZZDIPJ,EXYDIPJ
double precision :: EXZDIPJ,EYZDIPJ
integer :: i,j,IPOINT,JPOINT
! *** end of local definitions
 
!================================================
!

do  i=1,num

   exxls(i)=0.0d0
   eyyls(i)=0.0d0
   ezzls(i)=0.0d0
   exyls(i)=0.0d0
   exzls(i)=0.0d0
   eyzls(i)=0.0d0

   txxli(i)=0.0d0
   tyyli(i)=0.0d0
   tzzli(i)=0.0d0
   txyli(i)=0.0d0
   txzli(i)=0.0d0
   tyzli(i)=0.0d0

   srxx(i)=0.0d0
   sryy(i)=0.0d0
   srzz(i)=0.0d0
   srxy(i)=0.0d0
   srxz(i)=0.0d0
   sryz(i)=0.0d0

enddo



do j=2,num

   jpoint=ntype(j)

   do i=1,j-1
      ipoint=ntype(i)
!---> Memmory Reduction_S
      numx = numadr(i,j)
!     drsq=dxsav(i,j)*dxsav(i,j) &
!         +dysav(i,j)*dysav(i,j) &
!         +dzsav(i,j)*dzsav(i,j)
      drsq=dxsav(numx)*dxsav(numx) &
          +dysav(numx)*dysav(numx) &
          +dzsav(numx)*dzsav(numx)
!<--- Memmory Reduction_E

       
!
! implement cut-off:
          if (drsq.ge.rsqmax) CYCLE

             dr=dsqrt(drsq)
             drsqrec=1.0d0/drsq
             drrec=1.0d0/dr
             dr3rec=drrec*drsqrec
             dr7rec=drrec*dr3rec*dr3rec
!
! Dipole-dipole interaction tensor:
!
!---> Memmory Reduction_S
!            dxunitsq=dxsav(i,j)*dxsav(i,j)*drsqrec
!            dyunitsq=dysav(i,j)*dysav(i,j)*drsqrec
!            dzunitsq=dzsav(i,j)*dzsav(i,j)*drsqrec
             dxunitsq=dxsav(numx)*dxsav(numx)*drsqrec
             dyunitsq=dysav(numx)*dysav(numx)*drsqrec
             dzunitsq=dzsav(numx)*dzsav(numx)*drsqrec
!<--- Memmory Reduction_E

             sxx=1.0d0-3.0d0*dxunitsq
             syy=1.0d0-3.0d0*dyunitsq
             szz=1.0d0-3.0d0*dzunitsq
!---> Memmory Reduction_S
!            sxy=-3.0d0*dxsav(i,j)*dysav(i,j)*drsqrec
!            sxz=-3.0d0*dxsav(i,j)*dzsav(i,j)*drsqrec
!            syz=-3.0d0*dzsav(i,j)*dysav(i,j)*drsqrec
             sxy=-3.0d0*dxsav(numx)*dysav(numx)*drsqrec
             sxz=-3.0d0*dxsav(numx)*dzsav(numx)*drsqrec
             syz=-3.0d0*dzsav(numx)*dysav(numx)*drsqrec
!<--- Memmory Reduction_E

!..............DID (NB -not ewalded)
             txxtt=sxx*dr3rec
             tyytt=syy*dr3rec
             tzztt=szz*dr3rec
             txytt=sxy*dr3rec
             txztt=sxz*dr3rec
             tyztt=syz*dr3rec
             txxli(i)=txxli(i)+txxtt*xlipol(ipoint) &
                                *xlipol(jpoint)
             tyyli(i)=tyyli(i)+tyytt*xlipol(ipoint) &
                                 *xlipol(jpoint)
             tzzli(i)=tzzli(i)+tzztt*xlipol(ipoint) &
                                 *xlipol(jpoint)
             txyli(i)=txyli(i)+txytt*xlipol(ipoint) &
                                 *xlipol(jpoint)
             txzli(i)=txzli(i)+txztt*xlipol(ipoint) &
                                 *xlipol(jpoint)
             tyzli(i)=tyzli(i)+tyztt*xlipol(ipoint) &
                                 *xlipol(jpoint)
             txxli(j)=txxli(j)+txxtt*xlipol(ipoint) &
                                 *xlipol(jpoint)
             tyyli(j)=tyyli(j)+tyytt*xlipol(ipoint) &
                                 *xlipol(jpoint)
             tzzli(j)=tzzli(j)+tzztt*xlipol(ipoint) &
                                 *xlipol(jpoint)
             txyli(j)=txyli(j)+txytt*xlipol(ipoint) &
                                 *xlipol(jpoint)
             txzli(j)=txzli(j)+txztt*xlipol(ipoint) &
                                 *xlipol(jpoint)
             tyzli(j)=tyzli(j)+tyztt*xlipol(ipoint) &
                                 *xlipol(jpoint)

!................Field-gradient (NB -- not ewalded)

             exxls(j)=exxls(j)+txxtt*q(i)
             eyyls(j)=eyyls(j)+tyytt*q(i)
             ezzls(j)=ezzls(j)+tzztt*q(i)
             exyls(j)=exyls(j)+txytt*q(i)
             exzls(j)=exzls(j)+txztt*q(i)
             eyzls(j)=eyzls(j)+tyztt*q(i)

             exxls(i)=exxls(i)+txxtt*q(j)
             eyyls(i)=eyyls(i)+tyytt*q(j)
             ezzls(i)=ezzls(i)+tzztt*q(j)
             exyls(i)=exyls(i)+txytt*q(j)
             exzls(i)=exzls(i)+txztt*q(j)
             eyzls(i)=eyzls(i)+tyztt*q(j)

!........... SR light scattering terms.
             asrterm=asr(ipoint,jpoint)*dexp(-csr(ipoint,jpoint)*dr)
             bsrterm=bsr(ipoint,jpoint)*dexp(-dsr(ipoint,jpoint)*dr)

             xxcorr=asrterm-bsrterm*sxx
             yycorr=asrterm-bsrterm*syy
             zzcorr=asrterm-bsrterm*szz
             xycorr=-bsrterm*sxy
             xzcorr=-bsrterm*sxz
             yzcorr=-bsrterm*syz

             srxx(i)=srxx(i)+xxcorr
             sryy(i)=sryy(i)+yycorr
             srzz(i)=srzz(i)+zzcorr
             srxy(i)=srxy(i)+xycorr
             srxz(i)=srxz(i)+xzcorr
             sryz(i)=sryz(i)+yzcorr
             srxx(j)=srxx(j)+xxcorr
             sryy(j)=sryy(j)+yycorr
             srzz(j)=srzz(j)+zzcorr
             srxy(j)=srxy(j)+xycorr
             srxz(j)=srxz(j)+xzcorr
             sryz(j)=sryz(j)+yzcorr

!========================================================
! Dipole-quadrupole interaction tensor (see Buckingham &
! theory book 3, p 135 -).
! There are 10 of these.
!
! this part is not 'ewalded':
!
            threedr7=3.0d0*dr7rec
!---> Memmory Reduction_S
!           txxx=(5.0d0*dxsav(i,j)*dxsav(i,j)*dxsav(i,j) &
!              -3.0d0*dxsav(i,j)*drsq)*threedr7
!           tyyy=(5.0d0*dysav(i,j)*dysav(i,j)*dysav(i,j) &
!               -3.0d0*dysav(i,j)*drsq)*threedr7
!           tzzz=(5.0d0*dzsav(i,j)*dzsav(i,j)*dzsav(i,j) &
!              -3.0d0*dzsav(i,j)*drsq)*threedr7
            txxx=(5.0d0*dxsav(numx)*dxsav(numx)*dxsav(numx) &
               -3.0d0*dxsav(numx)*drsq)*threedr7
            tyyy=(5.0d0*dysav(numx)*dysav(numx)*dysav(numx) &
                -3.0d0*dysav(numx)*drsq)*threedr7
            tzzz=(5.0d0*dzsav(numx)*dzsav(numx)*dzsav(numx) &
               -3.0d0*dzsav(numx)*drsq)*threedr7
!<--- Memmory Reduction_E

!---> Memmory Reduction_S
!           txxy=(5.0d0*dxsav(i,j)*dxsav(i,j)*dysav(i,j) &
!              -drsq*dysav(i,j))*threedr7
!           txyy=(5.0d0*dxsav(i,j)*dysav(i,j)*dysav(i,j) &
!              -drsq*dxsav(i,j))*threedr7
!           txxz=(5.0d0*dxsav(i,j)*dxsav(i,j)*dzsav(i,j) &
!              -drsq*dzsav(i,j))*threedr7
!           txzz=(5.0d0*dxsav(i,j)*dzsav(i,j)*dzsav(i,j) &
!              -drsq*dxsav(i,j))*threedr7
!           tyyz=(5.0d0*dysav(i,j)*dysav(i,j)*dzsav(i,j) &
!              -drsq*dzsav(i,j))*threedr7
!           tyzz=(5.0d0*dysav(i,j)*dzsav(i,j)*dzsav(i,j) &
!              -drsq*dysav(i,j))*threedr7
!
!           txyz=15.0d0*dr7rec*dxsav(i,j)*dysav(i,j)*dzsav(i,j)
            txxy=(5.0d0*dxsav(numx)*dxsav(numx)*dysav(numx) &
               -drsq*dysav(numx))*threedr7
            txyy=(5.0d0*dxsav(numx)*dysav(numx)*dysav(numx) &
               -drsq*dxsav(numx))*threedr7
            txxz=(5.0d0*dxsav(numx)*dxsav(numx)*dzsav(numx) &
               -drsq*dzsav(numx))*threedr7
            txzz=(5.0d0*dxsav(numx)*dzsav(numx)*dzsav(numx) &
               -drsq*dxsav(numx))*threedr7
            tyyz=(5.0d0*dysav(numx)*dysav(numx)*dzsav(numx) &
               -drsq*dzsav(numx))*threedr7
            tyzz=(5.0d0*dysav(numx)*dzsav(numx)*dzsav(numx) &
               -drsq*dysav(numx))*threedr7

            txyz=15.0d0*dr7rec*dxsav(numx)*dysav(numx)*dzsav(numx)
!<--- Memmory Reduction_E

!==================================
!.......Light scattering -- dipole contribution to field-gradients

! 2) Dipole contribution - T_3 mu.
            exxdipi=txxx*xmu(i)+txxy*ymu(i)+txxz*zmu(i)
            eyydipi=tyyy*ymu(i)+txyy*xmu(i)+tyyz*zmu(i)
            ezzdipi=tzzz*zmu(i)+tyzz*ymu(i)+txzz*xmu(i)
            exydipi=txxy*xmu(i)+txyy*ymu(i)+txyz*zmu(i)
            exzdipi=txxz*xmu(i)+txyz*ymu(i)+txzz*zmu(i)
            eyzdipi=txyz*xmu(i)+tyyz*ymu(i)+tyzz*zmu(i)

!            exxls(j)=exxls(j)+exxdipi
!            eyyls(j)=eyyls(j)+eyydipi
!            ezzls(j)=ezzls(j)+ezzdipi
!            exyls(j)=exyls(j)+exydipi
!            exzls(j)=exzls(j)+exzdipi
!            eyzls(j)=eyzls(j)+eyzdipi

            exxdipj=-txxx*xmu(j)-txxy*ymu(j)-txxz*zmu(j)
            eyydipj=-tyyy*ymu(j)-txyy*xmu(j)-tyyz*zmu(j)
            ezzdipj=-tzzz*zmu(j)-tyzz*ymu(j)-txzz*xmu(j)
            exydipj=-txxy*xmu(j)-txyy*ymu(j)-txyz*zmu(j)
            exzdipj=-txxz*xmu(j)-txyz*ymu(j)-txzz*zmu(j)
            eyzdipj=-txyz*xmu(j)-tyyz*ymu(j)-tyzz*zmu(j)

!            exxls(i)=exxls(i)+exxdipj
!            eyyls(i)=eyyls(i)+eyydipj
!            ezzls(i)=ezzls(i)+ezzdipj
!            exyls(i)=exyls(i)+exydipj
!            exzls(i)=exzls(i)+exzdipj
!            eyzls(i)=eyzls(i)+eyzdipj

!        endif
   enddo
enddo

return

end subroutine
