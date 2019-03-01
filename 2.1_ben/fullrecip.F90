SUBROUTINE fullrecip

! GWW - module for erfc lookup table 
  use error_function 

  USE commondata, ONLY: num,q,x,y,z, &
                      twopi,sqrpi,pithreehalf,etasq,chgcorrec, &
                      eta,etaconst,onethird,nspec,nsp, &
                      xmu,ymu, zmu, &
                      quadxx,quadyy,quadzz,quadxy,quadxz,quadyz, &
                      dispcorrec1,dispcorrec2,dispalp, & 
! force check
                      gw_force, & 
                      engpetot, ftc 

  USE recipdata, ONLY: kmaxx,kmaxy,kmaxz,elcall,elsall,emcall,emsall,encall, &
                     ensall,rksqmax,bdisp, &
                     norm_ds,sk_ds,nbin,nspairs, &
                     norm_ds_w,sk_ds_w

  USE boxdata, ONLY: fullhi,cellvol,fourpicell,fac, twopiboxx, twopiboxy, twopiboxz &
                      & ,boxlenx, boxleny, boxlenz

  use mpipara

  IMPLICIT NONE

  INTEGER :: i,j,l,m,n,ll,mm,nn,nmin,mmin,kbin,lk,ilower,iupper
  DOUBLE PRECISION :: rl,rkx1,rky1,rkz1,rm,rkx2,rky2, &
                rkz2,rn,rkx3,rky3,rkz3,xkk,ckcs,ckss, &
                akk,arl,arm,arn,stfac,efac

#ifdef dipole      
  DOUBLE PRECISION :: cx,cy,cz,sx,sy,sz, &
                clx,cly,clz,slx,sly,slz,rdotmus,rdotmuc,&
                dipsum, xmumuengfac ,ddfac,dqfac
#endif 

#ifdef quadrupole 
  DOUBLE PRECISION :: cxx,cyy,czz,sxx,syy,szz,cxy,cxz,cyz,sxy,sxz,syz, &
                arll,armm,arnn,arlm,arln,armn,arijdotsij,arijdotcij, &
                quadsum, quadquadself
  DOUBLE PRECISION :: rkdotcxx,rkdotcyy,rkdotczz,rkdotsxx,rkdotsyy, &
                rkdotszz
#endif 

#ifdef vdw_ewald 
  double precision :: bparam,bparamrec,bparamsq,bparam3rec,erfcc,strfac, &
                drl,drm,drn,forcfac,dispfac,dispckcs,dispckss,sqrxkk, &
                sum_disp_struc
#endif 

! GWW - hope this OR works !
#ifdef vdw_ewald 
  DOUBLE PRECISION, DIMENSION(nspec) :: arsk,aisk
#elif debye_scherer
  DOUBLE PRECISION, DIMENSION(nspec) :: arsk,aisk
#endif 

! GWW - reduce stack size by using allocatable 
  DOUBLE PRECISION, allocatable :: clmall(:),slmall(:),ckcnoqall(:),cksnoqall(:), &
                                    temp(:),temp2(:),temp3(:), temp4(:) 


#ifdef debye_scherer
! GWW - Debye-Scherer - not used 
 DOUBLE PRECISION, DIMENSION(nspairs,0:nbin) :: tmp
 INTEGER, DIMENSION(0:1000) :: ttmp

  norm_ds_w=0
  sk_ds_w=0.0d0
#endif 



 double precision ::  gw_s1, gw_s2, gw_s3 
! GWW - reduce stack size by using allocatable 
! not sure how much time allocatre and deallocate takes ! 
! may be better to have these in one of pimaims so called modules ! (not really)
! could reduce for RIM since does not use temp2 / temp3 
! leave for now - will probably change how these are done  by using real modules 
! to avoid allocation every call 
   allocate (clmall(num), slmall(num))
   allocate (ckcnoqall(num),cksnoqall(num))
   allocate (temp(num),temp2(num),temp3(num),temp4(num)) 

  eltmp=0.d0
  sctmp=0.0d0
  gw_force=0.0d0 

#ifndef dipole 
#ifndef quadrupole
! GWW - need to do this since set in conjgradpim etc for polarisation models
! must change so done at appropriate time fro all systems and not repeating
! code everywehere - which this code seem to do as a virtue 

! Let's add in the initialising stuff for the sine and cosine arrays since it's magically absent
 elcall(:,0)=1.0d0
 emcall(:,0)=1.0d0
 encall(:,0)=1.0d0
 elsall(:,0)=0.0d0
 emsall(:,0)=0.0d0
 ensall(:,0)=0.0d0

 elcall(:,1)=cos(twopiboxx*x)
 emcall(:,1)=cos(twopiboxy*y)
 encall(:,1)=cos(twopiboxz*z)
 elsall(:,1)=sin(twopiboxx*x)
 emsall(:,1)=sin(twopiboxy*y)
 ensall(:,1)=sin(twopiboxz*z)
!
! Calculate all cosines and sines
!
do l=2,kmaxx
   elcall(:,l)=elcall(:,l-1)*elcall(:,1)-elsall(:,l-1)*elsall(:,1)
   elsall(:,l)=elsall(:,l-1)*elcall(:,1)+elcall(:,l-1)*elsall(:,1)
enddo   
do l=2,kmaxy
   emcall(:,l)=emcall(:,l-1)*emcall(:,1)-emsall(:,l-1)*emsall(:,1)
   emsall(:,l)=emsall(:,l-1)*emcall(:,1)+emcall(:,l-1)*emsall(:,1)
enddo   
do l=2,kmaxz
   encall(:,l)=encall(:,l-1)*encall(:,1)-ensall(:,l-1)*ensall(:,1)
   ensall(:,l)=ensall(:,l-1)*encall(:,1)+encall(:,l-1)*ensall(:,1)
enddo   

#endif
#endif 

!  do l = 1,num
!   write(6,('("fract ",i6,3(2X,F12.6))')) l,x(l)/boxlenx, y(l)/boxleny, z(l)/boxlenz 
!  enddo 


! GWW at some point try changing parallel implimentation so that
! it does one loop over stored vector combinations within the cutoff 
! one MPI process per node and X openMP threads
! loop is step number of threads * number of processes
! fisrt X vectors are done by threads on 1st node, next X threads on 2nd node
! etc. 

!  do  ll=0,kmaxx
  do ll=kmaxx_s,kmaxx_e
    l=iabs(ll)
    rl=dble(ll)*twopi
    rkx1=rl*fullhi(1,1)
    rky1=rl*fullhi(1,2)
    rkz1=rl*fullhi(1,3)

!    do mm=mmin,kmaxy
    do mm=kmaxy_s(ll),kmaxy_e(ll)
      m=iabs(mm)
      rm=dble(mm)*twopi
      rkx2=rkx1+rm*fullhi(2,1)
      rky2=rky1+rm*fullhi(2,2)
      rkz2=rkz1+rm*fullhi(2,3)

      if(mm.ge.0) then
        clmall(:)=elcall(:,l)*emcall(:,m)-elsall(:,l)*emsall(:,m)
        slmall(:)=elsall(:,l)*emcall(:,m)+emsall(:,m)*elcall(:,l)
      else
        clmall(:)=elcall(:,l)*emcall(:,m)+elsall(:,l)*emsall(:,m)
        slmall(:)=elsall(:,l)*emcall(:,m)-emsall(:,m)*elcall(:,l)
      endif

!     do nn=nmin,kmaxz
      do nn=kmaxz_s(mm,ll),kmaxz_e(mm,ll)
         n=iabs(nn)
         rn=dble(nn)*twopi
         rkx3=rkx2+rn*fullhi(3,1)
         rky3=rky2+rn*fullhi(3,2)
         rkz3=rkz2+rn*fullhi(3,3)
         xkk=rkx3*rkx3+rky3*rky3+rkz3*rkz3

#ifdef vdw_ewald 
! GWW-  length of recip vector (h) 
         sqrxkk=dsqrt(xkk)

!  GWW - dispalp = eta - not sure why the seperate variable ? 
!  b = h / 2 eta 
         bparam=sqrxkk/(2.0d0*dispalp)
! 1/b
         bparamrec=1.0d0/bparam
! b^2
         bparamsq=bparam*bparam
! 1/b^3 
         bparam3rec=bparamrec*bparamrec*bparamrec
         erfcc=erfunc(bparam)
#endif 

         if(xkk.le.rksqmax) then

            if(nn.ge.0) then
               ckcnoqall(:)=clmall(:)*encall(:,n)-slmall(:)*ensall(:,n)
               cksnoqall(:)=slmall(:)*encall(:,n)+clmall(:)*ensall(:,n)
            else
               ckcnoqall(:)=clmall(:)*encall(:,n)+slmall(:)*ensall(:,n)
               cksnoqall(:)=slmall(:)*encall(:,n)-clmall(:)*ensall(:,n)
            endif

! GWW - this looks a bit complicated but is quite simple
! vdw Ewald needs dispckcs etc. so needs all the following sections
! Debye-Scherer only needs arsk etc. so needs part of the loop
#ifdef vdw_ewald
            dispckcs=0.d0
            dispckss=0.d0

            do i=1,nspec
               ilower=SUM(nsp(0:i-1))+1
               iupper=SUM(nsp(0:i))
               arsk(i)=SUM(ckcnoqall(ilower:iupper))
               aisk(i)=SUM(cksnoqall(ilower:iupper))
! GWW not needed for non rule based VDW 
!               dispckcs=dispckcs+bdisp(i,i)*arsk(i)
!               dispckss=dispckss+bdisp(i,i)*aisk(i)
            enddo 

! GWW - summ of structure factors with seperate i,j VDW parameters (but pair
! type 
          sum_disp_struc = 0.d0 
          do i = 1, nspec
            do j = 1, nspec 
                sum_disp_struc = sum_disp_struc + ftc(i,j) * ( arsk(i) * arsk(j) + aisk(i) * aisk(j) ) 
            enddo
          enddo
   

#elif debye_scherer
            do i=1,nspec
               ilower=SUM(nsp(0:i-1))+1
               iupper=SUM(nsp(0:i))
               arsk(i)=SUM(ckcnoqall(ilower:iupper))
               aisk(i)=SUM(cksnoqall(ilower:iupper))
            enddo 
#endif

! GWW - Charges
            ckcs=SUM(ckcnoqall*q)
            ckss=SUM(cksnoqall*q)

! GWW - dipole
#ifdef dipole
            cx=SUM(ckcnoqall*xmu)
            cy=SUM(ckcnoqall*ymu)
            cz=SUM(ckcnoqall*zmu)

            sx=SUM(cksnoqall*xmu)
            sy=SUM(cksnoqall*ymu)
            sz=SUM(cksnoqall*zmu)
#endif 

! GWW - quadropoles 
#ifdef quadrupole 
            cxx=SUM(ckcnoqall*quadxx)
            cyy=SUM(ckcnoqall*quadyy)
            czz=SUM(ckcnoqall*quadzz)
            cxy=SUM(ckcnoqall*quadxy)
            cxz=SUM(ckcnoqall*quadxz)
            cyz=SUM(ckcnoqall*quadyz)

            sxx=SUM(cksnoqall*quadxx)
            syy=SUM(cksnoqall*quadyy)
            szz=SUM(cksnoqall*quadzz)
            sxy=SUM(cksnoqall*quadxy)
            sxz=SUM(cksnoqall*quadxz)
            syz=SUM(cksnoqall*quadyz)
#endif 


! GWW - Debye-Scherer - not usally used to put in an ifdef 
#ifdef debye_scherer 
            kbin=int((xkk/rksqmax)*float(nbin))
!           norm_ds(kbin)=norm_ds(kbin)+1
             norm_ds_w(kbin)=norm_ds_w(kbin)+1

             lk=1
             do i=1,nspec
               do j=i,nspec
!                sk_ds(lk,kbin)=sk_ds(lk,kbin)+(arsk(i)*arsk(j)) &
!                              +(aisk(i)*aisk(j))
                  sk_ds_w(lk,kbin)=sk_ds_w(lk,kbin)+(arsk(i)*arsk(j)) &
                                  +(aisk(i)*aisk(j))
                  lk=lk+1
               enddo
             enddo
#endif 


!
! GWW - charge-charge energy:
!
            akk=exp(etaconst*xkk)/xkk
            efac=akk*(ckss*ckss+ckcs*ckcs)
!           engpetot=engpetot+efac
            sctmp(55)=sctmp(55)+efac
!
! N.B. a self-interaction term needs to be subtracted at the end:
! this is chgcorrec, calculated in setup.f
!


! GWW - dispersion energy (reciprocal space part):
! unforunately uses different C parameters to real space - based 
! on rules and does not include the damping which is done on the
! real space  component. 
! It was only correct if the C terms obey the rules and no damping ! 
! modified to be correct  - sum_dis_struc calculated above used 
!


#ifdef vdw_ewald 
            forcfac=-(pithreehalf/(6.0d0*cellvol))*xkk*sqrxkk* &
                    (sqrpi*erfcc+exp(-bparamsq)* &
                    (0.5d0*bparam3rec-bparamrec))
! old rule based VDW 
!            dispfac=(dispckss*dispckss+dispckcs*dispckcs)
!            dispfac=forcfac*0.5d0*dispfac 
! new VDW 
            dispfac=forcfac*0.5d0*sum_disp_struc 
            
! dispacc
            sctmp(56)=sctmp(56)+dispfac

# endif 


! Use the cartesian components of k:
!
#ifdef vdw_ewald 
            drl=forcfac*rkx3
            drm=forcfac*rky3
            drn=forcfac*rkz3
#endif

            arl=akk*rkx3
            arm=akk*rky3
            arn=akk*rkz3


#ifdef quadrupole 
            arll=arl*rkx3
            armm=arm*rky3
            arnn=arn*rkz3
            arlm=arl*rky3
            arln=arl*rkz3
            armn=arm*rkz3

            arijdotsij=arll*sxx+armm*syy+arnn*szz+ &
                       2.0d0*(arlm*sxy+arln*sxz+armn*syz)
            arijdotcij=arll*cxx+armm*cyy+arnn*czz+ &
                       2.0d0*(arlm*cxy+arln*cxz+armn*cyz)
#endif 

!
! GWW - stress tensor 
!
! GWW - Charge - charge stress factor for stress tensor 
            stfac=(4.0d0*etasq+xkk)/(4.0d0*etasq*xkk)

! stcxx, stcyy, stczz, stcxy, stcxz, stcyz
            sctmp(7) =sctmp(7) +efac*(1.0d0-2.0d0*rkx3*rkx3*stfac)
            sctmp(10)=sctmp(10)+efac*(1.0d0-2.0d0*rky3*rky3*stfac)
            sctmp(12)=sctmp(12)+efac*(1.0d0-2.0d0*rkz3*rkz3*stfac)
            sctmp(8) =sctmp(8) -efac*2.0d0*rkx3*rky3*stfac
            sctmp(9) =sctmp(9) -efac*2.0d0*rkx3*rkz3*stfac
            sctmp(11)=sctmp(11)-efac*2.0d0*rky3*rkz3*stfac

#ifdef vdw_ewald 
! GWW - ewald vdw contribution to the  stress tensor 
            strfac=-(pithreehalf/(12.0d0*cellvol))* &
                   3.0d0*sqrxkk*(sqrpi*erfcc- &
                   bparamrec*exp(-bparamsq))* &
! old rule based VDW
!                   (dispckcs*dispckcs+dispckss*dispckss)
! new VDW 
                   sum_disp_struc 

! stsrxx, stsryy, stsrzz, stsrxy, stsrxz, stsryz
            sctmp(1)=sctmp(1)+dispfac+strfac*rkx3*rkx3
            sctmp(4)=sctmp(4)+dispfac+strfac*rky3*rky3
            sctmp(6)=sctmp(6)+dispfac+strfac*rkz3*rkz3
            sctmp(2)=sctmp(2)+strfac*rkx3*rky3
            sctmp(3)=sctmp(3)+strfac*rkx3*rkz3
            sctmp(5)=sctmp(5)+strfac*rky3*rkz3
#endif 

#ifdef dipole
            clx=cx*rkx3
            cly=cy*rky3
            clz=cz*rkz3
            slx=sx*rkx3
            sly=sy*rky3
            slz=sz*rkz3

            rdotmuc=clx+cly+clz
            rdotmus=slx+sly+slz
            ddfac=rdotmuc*rdotmuc+rdotmus*rdotmus

!
! GWW - dipole-dipole contributions to stress tensor
!
! stp2xx, stp2yy, stp2zz, stp2xy, stp2xz, stp2yz
            sctmp(19)=sctmp(19)+akk*(ddfac*(1.0d0-2.0d0*rkx3*rkx3*stfac) &
                  +2.0d0*(clx*rdotmuc+slx*rdotmus))
            sctmp(22)=sctmp(22)+akk*(ddfac*(1.0d0-2.0d0*rky3*rky3*stfac) &
                  +2.0d0*(cly*rdotmuc+sly*rdotmus))
            sctmp(24)=sctmp(24)+akk*(ddfac*(1.0d0-2.0d0*rkz3*rkz3*stfac) &
                  +2.0d0*(clz*rdotmuc+slz*rdotmus))
            sctmp(20)=sctmp(20)+akk*(-ddfac*2.0d0*rkx3*rky3*stfac &
                  +(rkx3*cy+rky3*cx)*rdotmuc+(rkx3*sy+rky3*sx)*rdotmus)
            sctmp(21)=sctmp(21)+akk*(-ddfac*2.0d0*rkx3*rkz3*stfac &
                  +(rkx3*cz+rkz3*cx)*rdotmuc+(rkx3*sz+rkz3*sx)*rdotmus)
            sctmp(23)=sctmp(23)+akk*(-ddfac*2.0d0*rky3*rkz3*stfac &
                  +(rky3*cz+rkz3*cy)*rdotmuc+(rky3*sz+rkz3*sy)*rdotmus)
!
! GWW - charge-dipole contributions to stress tensor
!
            dqfac=-2.0d0*akk*(ckcs*rdotmus-ckss*rdotmuc)

! stpxx stpyy, stpzz, stpxy, stpxz, stpyz
            sctmp(13)=sctmp(13)+2.0d0*arl*(cx*ckss-sx*ckcs) &
                 +dqfac*(1.0d0-2.0d0*rkx3*rkx3*stfac)
            sctmp(16)=sctmp(16)+2.0d0*arm*(cy*ckss-sy*ckcs) &
                 +dqfac*(1.0d0-2.0d0*rky3*rky3*stfac)
            sctmp(18)=sctmp(18)+2.0d0*arn*(cz*ckss-sz*ckcs) &
                 +dqfac*(1.0d0-2.0d0*rkz3*rkz3*stfac)
            sctmp(14)=sctmp(14)+arl*(cy*ckss-sy*ckcs) &
                 +arm*(cx*ckss-sx*ckcs)-dqfac*2.0d0*rkx3*rky3*stfac
            sctmp(15)=sctmp(15)+arl*(cz*ckss-sz*ckcs) &
                 +arn*(cx*ckss-sx*ckcs)-dqfac*2.0d0*rkx3*rkz3*stfac
            sctmp(17)=sctmp(17)+arm*(cz*ckss-sz*ckcs) &
                 +arn*(cy*ckss-sy*ckcs)-dqfac*2.0d0*rky3*rkz3*stfac
#endif 

! GWW - force and electric field due to charge - charge interaction 
!
            temp=cksnoqall*ckcs-ckcnoqall*ckss

#ifdef dipole
! GWW  charge - dip energy 
!
! -mu.field gives charge-dipole energy:
! extra factor of 2 at end over and above q-q and mu-mu energies
! since there are two of these terms in the energy equation:
!
            temp2=arl*xmu+arm*ymu+arn*zmu
! qmueng
            sctmp(57)=sctmp(57)-2.d0*SUM(temp*temp2)
#endif 

! GWW - and now the charge - charge forces 
            temp=temp*q

! frrx, frry, frrz
            eltmp(1:num)=eltmp(1:num)+arl*temp
            eltmp(num+1:2*num)=eltmp(num+1:2*num)+arm*temp
            eltmp(2*num+1:3*num)=eltmp(2*num+1:3*num)+arn*temp
!

#ifdef dipole
! GWW - force and electric field for charge - dipole interaction
!
! -mu.field gives dipole-dipole energy:
! multiplied by same factor as q-q energy at end.
!
            temp=-ckcnoqall*rdotmuc-cksnoqall*rdotmus
! xmumueng
            sctmp(57)=sctmp(57)-SUM(temp*temp2)

            temp=temp*q

! q-mu frrx, frry, frrz
            eltmp(1:num)=eltmp(1:num)+arl*temp
            eltmp(num+1:2*num)=eltmp(num+1:2*num)+arm*temp
            eltmp(2*num+1:3*num)=eltmp(2*num+1:3*num)+arn*temp

! GWW - forces and electric field  due to dipole - dipole interaction 
!
            temp=cksnoqall*rdotmuc-ckcnoqall*rdotmus
            temp3=rkx3*xmu+rky3*ymu+rkz3*zmu
            temp=temp*temp3

! mu-mu frrx, frry, frrz
            eltmp(1:num)=eltmp(1:num)+arl*temp
            eltmp(num+1:2*num)=eltmp(num+1:2*num)+arm*temp
            eltmp(2*num+1:3*num)=eltmp(2*num+1:3*num)+arn*temp

! GWW - dipole - charge  force and electric field 

            temp=ckcnoqall*ckcs+cksnoqall*ckss
            temp4=temp*temp3

! mu-q  frrx, frry, frrz
            eltmp(1:num)=eltmp(1:num)+arl*temp4
            eltmp(num+1:2*num)=eltmp(num+1:2*num)+arm*temp4
            eltmp(2*num+1:3*num)=eltmp(2*num+1:3*num)+arn*temp4

            gw_force(3*num+1:4*num)=gw_force(3*num+1:4*num)+arl*temp4
            gw_force(4*num+1:5*num)=gw_force(4*num+1:5*num)+arm*temp4
            gw_force(5*num+1:6*num)=gw_force(5*num+1:6*num)+arn*temp4
#endif 



#ifdef quadrupole
!
! -quad*field gradient gives quadrupole-charge energy 
!
            temp2=arll*quadxx+armm*quadyy+arnn*quadzz+2.0d0*( &
                  arlm*quadxy+arln*quadxz+armn*quadyz)

! qquadeng - factor of 2 acounts for quad-q and q-quad 
            sctmp(58)=sctmp(58)-onethird*2.d0*SUM(temp2*temp)

! GWW - forces and electric field due to various quadrupole interaction 

! GWW - force and electronc field for charge - quadrupole interaction
!
! -mu.field gives dipole-quadrupole energy
! This is the same as minus the product of quadi by the electric
! field gradient created by the dipoles (derived above).
!
            temp=onethird*(ckcnoqall*arijdotsij-cksnoqall*arijdotcij)

! dipquadengrec -  mu-quad energy  + quad-mu  (hence factor of 2) 
            sctmp(58)=sctmp(58)-2.d0*SUM(temp*temp3)


! GWW - charge - quadrupole forces 
            temp=temp*q

! q-quad  frrx, frry, frrz
            eltmp(1:num)=eltmp(1:num)+rkx3*temp
            eltmp(num+1:2*num)=eltmp(num+1:2*num)+rky3*temp
            eltmp(2*num+1:3*num)=eltmp(2*num+1:3*num)+rkz3*temp
!
! GWW  - forces and electric field due quad-charge 
! This involves the gradient of the field gradient, which is irrelevant
! for the quadrupole forces (it would enter the charge-octupole energy )
!
            temp=onethird*temp2*(ckcnoqall*ckss-cksnoqall*ckcs)

!           frrx=frrx+rkx3*temp
!           frry=frry+rky3*temp
!           frrz=frrz+rkz3*temp
! quad-q  frrx, frry, frrz
            eltmp(1:num)=eltmp(1:num)+rkx3*temp
            eltmp(num+1:2*num)=eltmp(num+1:2*num)+rky3*temp
            eltmp(2*num+1:3*num)=eltmp(2*num+1:3*num)+rkz3*temp


! GWW dipole - quadrupole forces and electric fields
!
! -onethird*quadi*egrad gives quadrupole-quadrupole energy 
!
            temp=-onethird*(ckcnoqall*arijdotcij+cksnoqall*arijdotsij)
! quadquadengrec
            sctmp(58)=sctmp(58)-onethird*SUM(temp*temp2)/akk

            temp=temp*temp3

! mu-quad  frrx, frry, frrz
            eltmp(1:num)=eltmp(1:num)+rkx3*temp
            eltmp(num+1:2*num)=eltmp(num+1:2*num)+rky3*temp
            eltmp(2*num+1:3*num)=eltmp(2*num+1:3*num)+rkz3*temp

            temp=onethird*temp2*(ckcnoqall*rdotmuc+cksnoqall*rdotmus)

! quad-mu  frrx, frry, frrz
            eltmp(1:num)=eltmp(1:num)+rkx3*temp
            eltmp(num+1:2*num)=eltmp(num+1:2*num)+rky3*temp
            eltmp(2*num+1:3*num)=eltmp(2*num+1:3*num)+rkz3*temp

! GWW - quadrupole - quadrupole forces 
!
            temp=onethird*onethird*temp2*(cksnoqall*arijdotcij-ckcnoqall*arijdotsij)/akk

! quad-quad  frrx, frry, frrz
            eltmp(1:num)=eltmp(1:num)+rkx3*temp
            eltmp(num+1:2*num)=eltmp(num+1:2*num)+rky3*temp
            eltmp(2*num+1:3*num)=eltmp(2*num+1:3*num)+rkz3*temp
#endif 

#ifdef vdw_ewald
!
! GWW - forces for reciprical space dispersion 
!
!            temp=cksnoqall*dispckcs-ckcnoqall*dispckss
!            do i=1,nspec
!               ilower=SUM(nsp(0:i-1))+1
!               iupper=SUM(nsp(0:i))
!               temp(ilower:iupper)=temp(ilower:iupper)*bdisp(i,i)
!            enddo
!            print*,' old VDW force temp ', temp 


! new VDWE forces not using the rules for C paramters 

            temp = 0.d0 
            do i=1 ,nspec 
               ilower=SUM(nsp(0:i-1))+1
               iupper=SUM(nsp(0:i))
               do j = 1, nspec 
               temp(ilower:iupper)= temp(ilower:iupper) + ftc(i,j) * ( cksnoqall(ilower:iupper) * arsk(j) - ckcnoqall(ilower:iupper) * aisk(j) ) 
               enddo
            enddo 


!            print*,' new VDW force temp ', temp 


! Ewald vdw  dispfrrx, dispfrry, dispfrrz
            eltmp(3*num+1:4*num)=eltmp(3*num+1:4*num)+drl*temp
            eltmp(4*num+1:5*num)=eltmp(4*num+1:5*num)+drm*temp
            eltmp(4*num+1:6*num)=eltmp(5*num+1:6*num)+drn*temp
#endif 

#ifdef quadrupole 
! GWW - Stress tensor for the charge - quadrupole interaction 
!
            rkdotcxx=rkx3*cxx+rky3*cxy+rkz3*cxz
            rkdotcyy=rkx3*cxy+rky3*cyy+rkz3*cyz
            rkdotczz=rkx3*cxz+rky3*cyz+rkz3*czz
            rkdotsxx=rkx3*sxx+rky3*sxy+rkz3*sxz
            rkdotsyy=rkx3*sxy+rky3*syy+rkz3*syz
            rkdotszz=rkx3*sxz+rky3*syz+rkz3*szz

!---> Parallelization_S
!           stpqquadxx=stpqquadxx+onethird* &
            sctmp(25) =sctmp(25) +onethird* &
                       (ckcs*arijdotcij+ckss*arijdotsij)* &
                       (1.0d0-2.0d0*rkx3*rkx3*stfac)- &
                       2.0d0*onethird*arl*(ckcs*rkdotcxx+ckss*rkdotsxx)
!           stpqquadyy=stpqquadyy+onethird* &
            sctmp(28) =sctmp(28) +onethird* &
                       (ckcs*arijdotcij+ckss*arijdotsij)* &
                       (1.0d0-2.0d0*rky3*rky3*stfac)- &
                       2.0d0*onethird*arm*(ckcs*rkdotcyy+ckss*rkdotsyy)
!           stpqquadzz=stpqquadzz+onethird* &
            sctmp(30) =sctmp(30) +onethird* &
                       (ckcs*arijdotcij+ckss*arijdotsij)* &
                       (1.0d0-2.0d0*rkz3*rkz3*stfac)- &
                       2.0d0*onethird*arn*(ckcs*rkdotczz+ckss*rkdotszz)
!           stpqquadxy=stpqquadxy-onethird* &
            sctmp(26) =sctmp(26) -onethird* &
                       (ckcs*arijdotcij+ckss*arijdotsij)* &
                       2.0d0*rkx3*rky3*stfac-onethird*(arl* &
                       (ckcs*rkdotcyy+ckss*rkdotsyy)+arm* &
                       (ckcs*rkdotcxx+ckss*rkdotsxx))
!           stpqquadxz=stpqquadxz-onethird* &
            sctmp(27) =sctmp(27) -onethird* &
                       (ckcs*arijdotcij+ckss*arijdotsij)* &
                       2.0d0*rkx3*rkz3*stfac-onethird*(arl* &
                       (ckcs*rkdotczz+ckss*rkdotszz)+arn* &
                       (ckcs*rkdotcxx+ckss*rkdotsxx))
!           stpqquadyz=stpqquadyz-onethird* &
            sctmp(29) =sctmp(29) -onethird* &
                       (ckcs*arijdotcij+ckss*arijdotsij)* &
                       2.0d0*rky3*rkz3*stfac-onethird*(arm* &
                       (ckcs*rkdotczz+ckss*rkdotszz)+arn* &
                       (ckcs*rkdotcyy+ckss*rkdotsyy))
 
! GWW - Stress tensor for the dipole - quadrupole interaction 
 
!           stpdipquadxx=stpdipquadxx+onethird* &
            sctmp(31)   =sctmp(31)   +onethird* &
                         (rdotmus*arijdotcij-rdotmuc*arijdotsij)* &
                         (1.0d0-2.0d0*rkx3*rkx3*stfac)+rkx3*onethird*( &
                          (cx*arijdotsij-sx*arijdotcij)+2.0d0*akk* &
                          (rdotmuc*rkdotsxx-rdotmus*rkdotcxx))
!           stpdipquadyy=stpdipquadyy+onethird* &
            sctmp(34)   =sctmp(34)   +onethird* &
                         (rdotmus*arijdotcij-rdotmuc*arijdotsij)* &
                         (1.0d0-2.0d0*rky3*rky3*stfac)+rky3*onethird*( &
                          (cy*arijdotsij-sy*arijdotcij)+2.0d0*akk* &
                          (rdotmuc*rkdotsyy-rdotmus*rkdotcyy))
!           stpdipquadzz=stpdipquadzz+onethird* &
            sctmp(36)   =sctmp(36)   +onethird* &
                         (rdotmus*arijdotcij-rdotmuc*arijdotsij)* &
                         (1.0d0-2.0d0*rkz3*rkz3*stfac)+rkz3*onethird*( &
                          (cz*arijdotsij-sz*arijdotcij)+2.0d0*akk* &
                          (rdotmuc*rkdotszz-rdotmus*rkdotczz))
!           stpdipquadxy=stpdipquadxy-2.0d0*onethird* &
            sctmp(32)   =sctmp(32)   -2.0d0*onethird* &
                         (rdotmus*arijdotcij-rdotmuc*arijdotsij)* &
                         rkx3*rky3*stfac+rkx3*onethird*0.5d0*( &
                          (cy*arijdotsij-sy*arijdotcij)+akk* &
                          (rdotmuc*rkdotsyy-rdotmus*rkdotcyy))+ &
                         rky3*onethird*0.5d0*( &
                          (cx*arijdotsij-sx*arijdotcij)+akk* &
                          (rdotmuc*rkdotsxx-rdotmus*rkdotcxx))
!           stpdipquadxz=stpdipquadxz-2.0d0*onethird* &
            sctmp(33)   =sctmp(33)   -2.0d0*onethird* &
                         (rdotmus*arijdotcij-rdotmuc*arijdotsij)* &
                         rkx3*rkz3*stfac+rkx3*onethird*0.5d0*( &
                          (cz*arijdotsij-sz*arijdotcij)+akk* &
                          (rdotmuc*rkdotszz-rdotmus*rkdotczz))+ &
                         rkz3*onethird*0.5d0*( &
                          (cx*arijdotsij-sx*arijdotcij)+akk* &
                          (rdotmuc*rkdotsxx-rdotmus*rkdotcxx))
!           stpdipquadyz=stpdipquadyz-2.0d0*onethird* &
            sctmp(35)   =sctmp(35)   -2.0d0*onethird* &
                         (rdotmus*arijdotcij-rdotmuc*arijdotsij)* &
                         rky3*rkz3*stfac+rky3*onethird*0.5d0*( &
                          (cz*arijdotsij-sz*arijdotcij)+akk* &
                          (rdotmuc*rkdotszz-rdotmus*rkdotczz))+ &
                         rkz3*onethird*0.5d0*( &
                          (cy*arijdotsij-sy*arijdotcij)+akk* &
                          (rdotmuc*rkdotsyy-rdotmus*rkdotcyy)) 

! GWW - Stress tensor for the quadrupole - quadrupole interaction 

!           stpquadquadxx=stpquadquadxx+onethird*onethird* &
            sctmp(37)    =sctmp(37)    +onethird*onethird* &
                          (arijdotcij**2.0d0+arijdotsij**2.0d0)* &
                          (1.0d0-2.0d0*rkx3*rkx3*stfac)+ &
                          4.0d0*onethird*onethird*rkx3* &
                          (rkdotcxx*arijdotcij+rkdotsxx*arijdotsij)
!           stpquadquadyy=stpquadquadyy+onethird*onethird* &
            sctmp(40)    =sctmp(40)    +onethird*onethird* &
                          (arijdotcij**2.0d0+arijdotsij**2.0d0)* &
                          (1.0d0-2.0d0*rky3*rky3*stfac)+ &
                          4.0d0*onethird*onethird*rky3* &
                          (rkdotcyy*arijdotcij+rkdotsyy*arijdotsij)
!           stpquadquadzz=stpquadquadzz+onethird*onethird* &
            sctmp(42)    =sctmp(42)    +onethird*onethird* &
                          (arijdotcij**2.0d0+arijdotsij**2.0d0)* &
                          (1.0d0-2.0d0*rkz3*rkz3*stfac)+ &
                          4.0d0*onethird*onethird*rkz3* &
                          (rkdotczz*arijdotcij+rkdotszz*arijdotsij)
!           stpquadquadxy=stpquadquadxy-onethird*onethird* &
            sctmp(38)    =sctmp(38)    -onethird*onethird* &
                          (arijdotcij**2.0d0+arijdotsij**2.0d0)* &
                                (2.0d0*rkx3*rky3*stfac)+ &
                          2.0d0*onethird*onethird*(rkx3* &
                          (rkdotcyy*arijdotcij+rkdotsyy*arijdotsij)+ &
                           rky3*(rkdotcxx*arijdotcij+rkdotsxx*arijdotsij))
!           stpquadquadxz=stpquadquadxz-onethird*onethird* &
            sctmp(39)    =sctmp(39)    -onethird*onethird* &
                          (arijdotcij**2.0d0+arijdotsij**2.0d0)* &
                                (2.0d0*rkx3*rkz3*stfac)+ &
                          2.0d0*onethird*onethird*(rkx3* &
                          (rkdotczz*arijdotcij+rkdotszz*arijdotsij)+ &
                           rkz3*(rkdotcxx*arijdotcij+rkdotsxx*arijdotsij))
!           stpquadquadyz=stpquadquadyz-onethird*onethird* &
            sctmp(41)    =sctmp(41)    -onethird*onethird* &
                          (arijdotcij**2.0d0+arijdotsij**2.0d0)* &
                                (2.0d0*rky3*rkz3*stfac)+ &
                          2.0d0*onethird*onethird*(rky3* &
                          (rkdotczz*arijdotcij+rkdotszz*arijdotsij)+ &
                           rkz3*(rkdotcyy*arijdotcij+rkdotsyy*arijdotsij))
#endif 

         endif

      enddo   
!     nmin=-kmaxz
   enddo   
!  mmin=-kmaxy
enddo   
!
!
eltmp(1:3*num)=eltmp(1:3*num)*fac
gw_force(1:9*num)=gw_force(1:9*num)*fac

sctmp(55)=sctmp(55)*fourpicell
sctmp(57)=sctmp(57)*fourpicell
sctmp(58)=sctmp(58)*fourpicell

! scale sctmp for all components - even those zero when not being done (e.g.
! 37-42 - is for quad-quad terms. 
sctmp(7:42)=sctmp(7:42)*fourpicell
sctmp(25:36)=sctmp(25:36)*2.d0

! MS modif 31/03
!sctmp(1)=sctmp(1)+dispcorrec1/cellvol
!sctmp(4)=sctmp(4)+dispcorrec1/cellvol
!sctmp(6)=sctmp(6)+dispcorrec1/cellvol
! MS endmodif 31/03

!print*,'Recip E = ',sctmp(55)

! GWW - chgcorrec term (charge-charge self-interaction)is calculated in setup.
  if (iam.eq.0) sctmp(55) = sctmp(55) -chgcorrec

!print*,'Recip E after self term = ',sctmp(55)

#ifdef vdw_ewald
! GWW - Ewald dispersion correction (self energy)
  if (iam.eq.0) sctmp(56) = sctmp(56) +(dispcorrec1/cellvol)+dispcorrec2
#endif 
 
#ifdef dipole      
! GWW - dipole-dipole self-interaction
  if (iam.eq.0) then 
    dipsum=SUM(xmu*xmu+ymu*ymu+zmu*zmu)
    xmumuengfac=2.0d0*(eta**3.0d0)/(3.0d0*sqrpi)*dipsum
    sctmp(57)=sctmp(57)-xmumuengfac
  endif
#endif 


#ifdef quadrupole 
! GWW - quadrupole-quadrupole self-interaction 
  if (iam.eq.0) then 
    quadsum=SUM(quadxx*quadxx+quadyy*quadyy+quadzz*quadzz+ &
          2.0d0*(quadxy*quadxy+quadxz*quadxz+quadyz*quadyz))
    quadquadself=8.0d0*(eta**5.0d0)/(45.0d0*sqrpi)*quadsum
    sctmp(58)=sctmp(58)-quadquadself
  endif 
#endif 
   
#ifdef debye_scherer
! GWW - Stuff for Debye-Scherer 
 CALL MPI_ALLREDUCE(norm_ds_w,ttmp,1001,mpi_integer,MPI_SUM, &
                    MPI_COMM_WORLD,ierr)
 norm_ds = norm_ds + ttmp
 
 CALL MPI_ALLREDUCE(sk_ds_w,tmp,nspairs*(nbin+1),mpi_double_precision,MPI_SUM, &
                    MPI_COMM_WORLD,ierr)
 sk_ds = sk_ds + tmp
#endif 

   deallocate (clmall, slmall)
   deallocate (ckcnoqall,cksnoqall)
   deallocate (temp,temp2,temp3,temp4)



return
END SUBROUTINE
