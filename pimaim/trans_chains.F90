SUBROUTINE trans_vv

USE commondata, ONLY: num,forfl,frrx,frry,frrz, &
          vpzeta1,vpzeta2,vpzeta3,nth,  &
          pzeta2,pzeta3,nib,vbzeta1,vbzeta2,vbzeta3,bzeta2,bzeta3, &
          veps1,veps2,veps3,eps2,eps3,vol3,ntype,vx,vy,vz,x,y,z,  &
          amass,pint2,stsrxx,stsryy,stsrzz,stsrxy,stsrxz,stsryz,  &
          stcxx,stcyy,stczz,stcxy,stcxz,stcyz,stpxx,stpyy,stpzz,  &
          stpxy,stpxz,stpyz,stp2xx,stp2yy,stp2zz,stp2xy,stp2xz,  &
          stp2yz,stpsrxx,stpsryy,stpsrzz,stpsrxy,stpsrxz,stpsryz,  &
          stqsrxx,stqsryy,stqsrzz,stqsrxy,stqsrxz,stqsryz,  &
          stpqquadxx,stpqquadyy,stpqquadzz,stpqquadxy,stpqquadxz,  &
          stpqquadyz,stpdipquadxx,stpdipquadyy,stpdipquadzz,  &
          stpdipquadxy,stpdipquadxz,stpdipquadyz,stpquadquadxx,  &
          stpquadquadyy,stpquadquadzz,stpquadquadxy,stpquadquadxz, &
          stpquadquadyz,stewzz,ortho,trantkb,free,CUEp,CUEp2,  &
          CUEprec,CUEp2rec,CUEb,CUEbrec,CUEb2,CUEb2rec,dtime,dom,  &
          W,Wgo,Wrec,Wgorec,pext,rcut,rsqmax,recamass,  &
          environmentalpimlog,conjgradaimlog,conjgradlog,Pcomx,  &
          Pcomy,Pcomz,tranke,nspec,nsp,tcell,tvol,tbzeta,tpzeta,  &
          PeeVee,pzeta,bzeta,tke,gtrankin,nstep,msdcalllog,  &
          nmsdcalltime,moveions,nab,hmass,frrx3,frry3,frrz3 ,engpetot, &
          xdisp,ydisp,zdisp
USE boxdata, ONLY: vg2,vg3,h2,h3,hi2,hi3,hlab2,hlab3 &
         ,cellvol,cellvol3rec,bee,boxlenx,boxleny,boxlenz,hi,h,fullh,fullhi 
!---> Parallelization_S
use mpipara
!<--- Parallelization_E
!---> Debug_S
use recipdata, only:elcall,elsall,emcall,emsall,encall,ensall,kmaxx,kmaxy,kmaxz
!<--- Debug_E

IMPLICIT NONE

INTEGER :: i,j,k,l,m,ir,ionpoint,icnt,ion,ipoint
INTEGER, DIMENSION(3,3) :: iden
!---> Debug_S
INTEGER :: kmaxx_w,kmaxy_w,kmaxz_w
!<--- Debug_E
DOUBLE PRECISION :: xkesq,xkesq2,rmassvol,pinst2,trvg2,Feps2,volfac, &
     Vcomx,Vcomy,Vcomz,totmass,COMke,trace,trvg3,Feps3,converge, &
     pinst3,tranko,vol2
DOUBLE PRECISION, DIMENSION(num) :: vtx,vty,vtz,vxni,vyni,vzni, &
     vxold,vyold,vzold,rxn,ryn,rzn,tx,ty,tz, &
     vxn,vyn,vzn,ax,ay,az,ax3,ay3,az3,vmag,rxnn,rynn, &
     rznn,comx,comy,comz,comx2,comy2,comz2
DOUBLE PRECISION, DIMENSION(5) :: gpzeta2,gpzeta3,gbzeta2,gbzeta3
DOUBLE PRECISION, DIMENSION(3,3) :: riden,hcom,com1,pint3,pint3c, &
              com3,Fgo2,Fgo3,com2,vg3sq,vg2sq,comexp,hhi,hh


! GWW resetting comke does not look correct - varible to check total KE / mv
double precision newCOMke 

if(.not. forfl) then
   xdisp=0.d0
   ydisp=0.d0
   zdisp=0.d0
end if
if(forfl) then
   frrx=frrx3
   frry=frry3
   frrz=frrz3
endif

forfl=.true.


vpzeta1=vpzeta2
vpzeta2=vpzeta3
pzeta2=pzeta3
gpzeta2=0.0d0
gpzeta3=0.0d0
if(.not.nth)then
   vpzeta1=0.d0
   vpzeta2=0.d0
   vpzeta3=0.d0
   pzeta2=0.d0
endif

vbzeta1=vbzeta2
vbzeta2=vbzeta3
bzeta2=bzeta3
gbzeta2=0.0d0
gbzeta3=0.0d0


! GWW - vhy is this zeroed for nab ? 
if((.not.nth).or.(.not.nib))then
   vbzeta1=0.d0
   vbzeta2=0.d0
   vbzeta3=0.d0
   bzeta2=0.d0
endif

veps1=veps2
veps2=veps3
eps2=eps3
vol2=vol3

! GWW - vhy is this zeroed for nab ? 
if(.not.nib)then
   veps1=0.d0
   veps2=0.d0
   veps3=0.d0
endif

vg2=vg3
riden=0.0d0
iden=0
do j=1,3
   riden(j,j)=1.0d0
   iden(j,j)=1
enddo   

h2=h3
hi2=hi3
hlab2=hlab3
  
! vtx,y,z are the velocities in the Cartesian frame.
vtx=hlab2(1,1)*vx+hlab2(1,2)*vy+hlab2(1,3)*vz
vty=hlab2(2,1)*vx+hlab2(2,2)*vy+hlab2(2,3)*vz
vtz=hlab2(3,1)*vx+hlab2(3,2)*vy+hlab2(3,3)*vz
! positions in the Cartesian frame.
tx=hlab2(1,1)*x+hlab2(1,2)*y+hlab2(1,3)*z
ty=hlab2(2,1)*x+hlab2(2,2)*y+hlab2(2,3)*z
tz=hlab2(3,1)*x+hlab2(3,2)*y+hlab2(3,3)*z
! 'xkesq2' is the k.e. at time t.
xkesq2=0.0d0
do i=1,num
   ionpoint=ntype(i)
   xkesq2=xkesq2+(amass(ionpoint)*(vtx(i)*vtx(i)+ &
               vty(i)*vty(i)+vtz(i)*vtz(i)))         
enddo   


!     Obtain pressure at time t=0
pint2=0.0d0
do i=1,num
   ionpoint=ntype(i)
   rmassvol=amass(ionpoint)/cellvol
   pint2(1,1)=pint2(1,1)+ rmassvol*vtx(i)*vtx(i)
   pint2(1,2)=pint2(1,2)+ rmassvol*vtx(i)*vty(i)
   pint2(1,3)=pint2(1,3)+ rmassvol*vtx(i)*vtz(i)
   pint2(2,2)=pint2(2,2)+ rmassvol*vty(i)*vty(i)
   pint2(2,3)=pint2(2,3)+ rmassvol*vty(i)*vtz(i)
   pint2(3,3)=pint2(3,3)+ rmassvol*vtz(i)*vtz(i)
enddo   

pint2(1,1)=pint2(1,1)+(((3.0d0*(stpxx+ &
                 stcxx+stsrxx+stp2xx+stpsrxx+stqsrxx+ &
                 stpqquadxx+stpdipquadxx+stpquadquadxx)) &
                 *cellvol3rec))
pint2(2,2)=pint2(2,2)+(((3.0d0*(stpyy+ &
                 stcyy+stsryy+stp2yy+stpsryy+stqsryy+ &
                 stpqquadyy+stpdipquadyy+stpquadquadyy)) &
                 *cellvol3rec))
pint2(3,3)=pint2(3,3)+(((3.0d0*(stpzz+ &
                 stczz+stsrzz+stp2zz+stpsrzz+stqsrzz+ &
                 stpqquadzz+stpdipquadzz+stpquadquadzz)) &
                 *cellvol3rec))+stewzz
pint2(1,2)=pint2(1,2)+(((3.0d0*(stpxy+ &
                 stcxy+stsrxy+stp2xy+stpsrxy+stqsrxy+ &
                 stpqquadxy+stpdipquadxy+stpquadquadxy)) &
                 *cellvol3rec))
pint2(1,3)=pint2(1,3)+(((3.0d0*(stpxz+ &
                 stcxz+stsrxz+stp2xz+stpsrxz+stqsrxz+ &
                 stpqquadxz+stpdipquadxz+stpquadquadxz)) &
                 *cellvol3rec))
pint2(2,3)=pint2(2,3)+(((3.0d0*(stpyz+ &
                 stcyz+stsryz+stp2yz+stpsryz+stqsryz+ &
                 stpqquadxz+stpdipquadxz+stpquadquadxz)) &
                 *cellvol3rec))
pint2(2,1)=pint2(1,2)
pint2(3,1)=pint2(1,3)
pint2(3,2)=pint2(2,3)

pinst2=(pint2(1,1)+pint2(2,2)+pint2(3,3))/3.0d0

trvg2=0.0d0
do j=1,3
   do k=1,3
      trvg2=trvg2+(vg2(j,k)*vg2(k,j))
   enddo   
enddo   
!.......update particle zetas
gpzeta2(1)=(CUEprec)*((xkesq2)-(free*trantkb))-vpzeta2(1)*vpzeta2(2)
gpzeta2(2)=(CUEp2rec)*(CUEp*vpzeta2(1)**2-trantkb)-vpzeta2(2)*vpzeta2(3)
gpzeta2(3)=(CUEp2rec)*(CUEp2*vpzeta2(2)**2-trantkb)-vpzeta2(3)*vpzeta2(4)
gpzeta2(4)=(CUEp2rec)*(CUEp2*vpzeta2(3)**2-trantkb)-vpzeta2(4)*vpzeta2(5)
gpzeta2(5)=(CUEp2rec)*(CUEp2*vpzeta2(4)**2-trantkb)
     
pzeta3=pzeta2+vpzeta2*dtime+gpzeta2*(dtime**2)/2.0d0
!.......update barostat zetas
gbzeta2(1)=(CUEbrec)*( W*veps2**2+Wgo*trvg2-dom*trantkb)-vbzeta2(1)*vbzeta2(2)
gbzeta2(2)=(CUEb2rec)*(CUEb*vbzeta2(1)**2-trantkb)-vbzeta2(2)*vbzeta2(3)
gbzeta2(3)=(CUEb2rec)*(CUEb2*vbzeta2(2)**2-trantkb)-vbzeta2(3)*vbzeta2(4)
gbzeta2(4)=(CUEb2rec)*(CUEb2*vbzeta2(3)**2-trantkb)-vbzeta2(4)*vbzeta2(5)
gbzeta2(5)=(CUEb2rec)*(CUEb2*vbzeta2(4)**2-trantkb) 

bzeta3=bzeta2+vbzeta2*dtime+gbzeta2*(dtime**2)/2.0d0
!...........update the volume
Feps2=((3.0d0)*vol2*(pinst2-pext))+((3.0d0/free)*xkesq2)
   
eps3=eps2+(dtime*veps2)+(((dtime**2)/2.0d0)* &
       ((Feps2*Wrec)-(veps2*vbzeta2(1))))

volfac=exp(eps3-eps2)
vol3=vol2*volfac**3

!.......invert h0(t=0)
call invert(h2,hi2)
!........update the cell matrix

! isotropic term 
Fgo2=vol2*(pint2-(pext*riden))  &
    -((vol2/3.0d0)*(((pint2(1,1)-(riden(1,1)*pext))+  &
     (pint2(2,2)-(riden(2,2)*pext))+(pint2(3,3)-(riden(3,3)*pext)))*riden))

! anisotropic term 
vg2sq=MATMUL(vg2,vg2)
hcom=riden+(vg2*dtime)+(((dtime**2  &
    /2.0d0)*(Fgo2*Wgorec)+(vg2sq)-(vg2*vbzeta2(1))))

!............update h3, which is Martyna's h_0(Delta_t)
h3=MATMUL(hcom,h2)
if(ortho)then
   h3(1,2)=0.d0
   h3(1,3)=0.d0
   h3(2,3)=0.d0
   h3(2,1)=0.d0
   h3(3,1)=0.d0
   h3(3,2)=0.d0
endif

! Impose constraint that |h(t+dt)| is equal to 1
call shake
!.......invert h(t+dt)
call invert(h3,hi3)
!............hh corresponds to the matrix of cell-side vectors,
!............it is also known as fullh elsewhere in the program
hh=(vol3**(1.0d0/3.0d0))*h3

call dcell(hh,bee)

bee(4:6)=acos(bee(4:6))
boxlenx=bee(1)
boxleny=bee(2)
boxlenz=bee(3)
!..........hlab is the matrix which transforms from  cell to cartesian
!.......i.e. it is the matrix of unit vectors along the cell sides
!........it is the same as the matrix h which appears elsewhere in the program
hlab3(:,1)=hh(:,1)/boxlenx
hlab3(:,2)=hh(:,2)/boxleny
hlab3(:,3)=hh(:,3)/boxlenz

h=hlab3


! Update number density
! boxreset redetermines all the cell-dependent parameters....
!....AND ALSO INVERTS THE CELL MATRIX h -> hi as WELL AS fullh -> fullhi
!...i.e. the above is a highly circuitous route (but does some useful things!)
call boxreset

!if(nib)then
!rcut=dmin1(boxlenx,boxleny,boxlenz)/2.0d0
!rsqmax=rcut**2.0d0
!endif

!---> Debug_S
kmaxx_w=kmaxx
kmaxy_w=kmaxy
kmaxz_w=kmaxz
!<--- Debug_E

call kset

!---> Debug_S
if( kmaxx .gt. kmaxx_w ) then
   deallocate(elcall,elsall)
   ALLOCATE ( elcall(num,0:kmaxx+1), elsall(num,0:kmaxx+1) )
endif
if( kmaxy .gt. kmaxy_w ) then
   deallocate(emcall,emsall)
   ALLOCATE ( emcall(num,0:kmaxy+1), emsall(num,0:kmaxy+1) )
endif
if( kmaxz .gt. kmaxz_w ) then
   deallocate(encall,ensall)
   ALLOCATE ( encall(num,0:kmaxz+1), ensall(num,0:kmaxz+1) )
endif
!<--- Debug_E
!---> Parallelization_S
if( kmaxx .ne. kmaxx_w .or. kmaxy .ne. kmaxy_w .or. kmaxz .ne. kmaxz_w) then
   deallocate (kmaxy_s,kmaxy_e,kmaxz_s,kmaxz_e)
   allocate( kmaxy_s(0:kmaxx), kmaxy_e(0:kmaxx) )
   allocate( kmaxz_s(-kmaxy:kmaxy,0:kmaxx), kmaxz_e(-kmaxy:kmaxy,0:kmaxx) )

   call kmaxpara
endif
!<--- Parallelization_E

comexp=volfac*MATMUL(h3,hi2)
! GWW - why zeroed for nab ? 
if(.not.nib) then
   comexp=0.d0
   comexp(1,1)=1.0d0
   comexp(2,2)=1.0d0
   comexp(3,3)=1.0d0
endif

do i=1,num
   ionpoint=ntype(i)

   ax(i)=frrx(i)*recamass(ionpoint)
   ay(i)=frry(i)*recamass(ionpoint)
   az(i)=frrz(i)*recamass(ionpoint)
enddo
!     calculate positions at t+dt
rxn=tx+(dtime*vtx)+((0.5d0*dtime**2)*(ax-(vpzeta2(1)*vtx) &
    -(2.0d0*((vg2(1,1)*vtx)+(vg2(1,2)*vty)+ &
    (vg2(1,3)*vtz)))-((2.0d0+(3.0d0/free))*vtx*veps2)))
ryn=ty+(dtime*vty)+((0.5d0*dtime**2)*(ay-(vpzeta2(1)*vty) &
    -(2.0d0*((vg2(2,1)*vtx)+(vg2(2,2)*vty)+ &
    (vg2(2,3)*vtz)))-((2.0d0+(3.0d0/free))*vty*veps2)))
rzn=tz+(dtime*vtz)+((0.5d0*dtime**2)*(az-(vpzeta2(1)*vtz) &
    -(2.0d0*((vg2(3,1)*vtx)+(vg2(3,2)*vty)+ &
    (vg2(3,3)*vtz)))-((2.0d0+(3.0d0/free))*vtz*veps2)))

! GWW at this point rxn is the carteisan coordinates 
rxnn=rxn
rynn=ryn
rznn=rzn

rxn=(comexp(1,1)*rxnn)+(comexp(1,2)*rynn)+(comexp(1,3)*rznn)
ryn=(comexp(2,1)*rxnn)+(comexp(2,2)*rynn)+(comexp(2,3)*rznn)
rzn=(comexp(3,1)*rxnn)+(comexp(3,2)*rynn)+(comexp(3,3)*rznn)

rxnn=rxn
rynn=ryn
rznn=rzn

! displacements for 'disp.out'. In Cartesiam frame.
xdisp=xdisp+(rxn-tx)
ydisp=ydisp+(ryn-ty)
zdisp=zdisp+(rzn-tz)

! GWW turn rxn into distance along vectors - not fractional not cartesian
! really stupid using the same variable for different representations 
!.........rotate new ionic positions back to the cell frame ready for
!..... the update step
rxn=rxnn*hi(1,1)+rynn*hi(1,2)+rznn*hi(1,3)
ryn=rxnn*hi(2,1)+rynn*hi(2,2)+rznn*hi(2,3)
rzn=rxnn*hi(3,1)+rynn*hi(3,2)+rznn*hi(3,3)
!     forces at time t+dt are required 
if(moveions) then
   x=rxn
   y=ryn
   z=rzn
   do i=1,num
      if(x(i).lt.0.0d0) x(i)=x(i)+boxlenx
      if(y(i).lt.0.0d0) y(i)=y(i)+boxleny
      if(z(i).lt.0.0d0) z(i)=z(i)+boxlenz
      if(x(i).gt.boxlenx) x(i)=x(i)-boxlenx
      if(y(i).gt.boxleny) y(i)=y(i)-boxleny
      if(z(i).gt.boxlenz) z(i)=z(i)-boxlenz
   enddo   
endif
     
if(environmentalpimlog) then
   call conjgradpimaim
else
   if(conjgradaimlog) call conjgradaim
   if(conjgradlog) then
	 call conjgrad
	endif
endif


call ener

frrx3=frrx
frry3=frry
frrz3=frrz

do i=1,num
   ionpoint=ntype(i)

   ax3(i)=frrx3(i)*recamass(ionpoint)
   ay3(i)=frry3(i)*recamass(ionpoint)
   az3(i)=frrz3(i)*recamass(ionpoint)
enddo   
!       Initial guess for vzeta3 and veps3 (from Verlet)
vpzeta3=vpzeta1+(2.0d0*gpzeta2*dtime)
vbzeta3=vbzeta1+(2.0d0*gbzeta2*dtime)
veps3=veps1+(2.0d0*dtime)*((Feps2*Wrec)-(veps2*vbzeta2(1)))
!....update the configurational part of the stress tensor
pint3c(1,1)=((3.0d0*(stpxx+ &
                 stcxx+stsrxx+stp2xx+stpsrxx+stqsrxx+ &
                 stpqquadxx+stpdipquadxx+stpquadquadxx)) &
                 *cellvol3rec)
pint3c(2,2)=((3.0d0*(stpyy+ &
                 stcyy+stsryy+stp2yy+stpsryy+stqsryy+ &
                 stpqquadyy+stpdipquadyy+stpquadquadyy)) &
                 *cellvol3rec)
pint3c(3,3)=((3.0d0*(stpzz+ &
                 stczz+stsrzz+stp2zz+stpsrzz+stqsrzz+ &
                 stpqquadzz+stpdipquadzz+stpquadquadzz)) &
                 *cellvol3rec)+stewzz
pint3c(1,2)=((3.0d0*(stpxy+ &
                 stcxy+stsrxy+stp2xy+stpsrxy+stqsrxy+ &
                 stpqquadxy+stpdipquadxy+stpquadquadxy)) &
                 *cellvol3rec)
pint3c(1,3)=((3.0d0*(stpxz+ &
                 stcxz+stsrxz+stp2xz+stpsrxz+stqsrxz+ &
                 stpqquadxz+stpdipquadxz+stpquadquadxz)) &
                 *cellvol3rec)
pint3c(2,3)=((3.0d0*(stpyz+ &
                 stcyz+stsryz+stp2yz+stpsryz+stqsryz+ &
                 stpqquadyz+stpdipquadyz+stpquadquadyz)) &
                 *cellvol3rec)
pint3c(2,1)=pint3c(1,2)
pint3c(3,1)=pint3c(1,3)
pint3c(3,2)=pint3c(2,3)
!..........parts of velocity updates at dt=0 -> vxni, vyni etc
comx=vtx+((0.5d0*dtime)*(ax-(vtx*vpzeta2(1)) &
    -(2.0d0*((vg2(1,1)*vtx)+(vg2(1,2)*vty) &
    +(vg2(1,3)*vtz)))-((2.0d0+(3.0d0/free))*veps2*vtx)))
comy=vty+((0.5d0*dtime)*(ay-(vty*vpzeta2(1)) &
    -(2.0d0*((vg2(2,1)*vtx)+(vg2(2,2)*vty) &
    +(vg2(2,3)*vtz)))-((2.0d0+(3.0d0/free))*veps2*vty)))
comz=vtz+((0.5d0*dtime)*(az-(vtz*vpzeta2(1)) &
    -(2.0d0*((vg2(3,1)*vtx)+(vg2(3,2)*vty) &
    +(vg2(3,3)*vtz)))-((2.0d0+(3.0d0/free))*veps2*vtz)))
 
vxni=((comexp(1,1)*comx)+(comexp(1,2)*comy)+(comexp(1,3)*comz))
vyni=((comexp(2,1)*comx)+(comexp(2,2)*comy)+(comexp(2,3)*comz))
vzni=((comexp(3,1)*comx)+(comexp(3,2)*comy)+(comexp(3,3)*comz))
!.........initial guess at the remaining parts of velocities --
!.........uses vg and v from t=0 rather than delta t -> vxn etc
comx2=(0.5d0*dtime)*((ax3)-(vtx*vpzeta3(1)) &
    -(2.0d0*((vg2(1,1)*vtx)+(vg2(1,2)*vty) &
    +(vg2(1,3)*vtz)))-((2.0d0+(3.0d0/free))*vtx*veps3))
comy2=(0.5d0*dtime)*((ay3)-(vtx*vpzeta3(1)) &
    -(2.0d0*((vg2(2,1)*vtx)+(vg2(2,2)*vty) &
    +(vg2(2,3)*vtz)))-((2.0d0+(3.0d0/free))*vty*veps3))
comz2=(0.5d0*dtime)*((az3)-(vtz*vpzeta3(1)) &
    -(2.0d0*((vg2(3,1)*vtx)+(vg2(3,2)*vty) &
    +(vg2(3,3)*vtz)))-((2.0d0+(3.0d0/free))*vtz*veps3))

vxn=vxni+comx2
vyn=vyni+comy2
vzn=vzni+comz2

hhi=MATMUL(h2,hi3)
com1=(Fgo2*Wgorec)+vg2sq-(vbzeta2(1)*vg2)
com2=vg2+(dtime*com1/2.0d0)
com3=MATMUL(com2,hhi)
vg3=com3+(dtime*com1/2.0d0)
!.......make an initial guess at vg3
trace=vg3(1,1)+vg3(2,2)+vg3(3,3)
vg3(1,1)=vg3(1,1)-(1.0d0/3.0d0)*trace
vg3(2,2)=vg3(2,2)-(1.0d0/3.0d0)*trace
vg3(3,3)=vg3(3,3)-(1.0d0/3.0d0)*trace

! GWW why not just not do the calculation if constant volume ????
! instead of doing it with zero parameters ! 
if(.not.nib) then
vg3=0.d0
veps2=0.d0
veps3=0.d0
endif


!     begin iteration loop for velocities
icnt=0
380 continue
icnt=icnt+1
vxold=vxn
vyold=vyn
vzold=vzn
!     obtain pressure at time t+dt
pint3=0.0d0

do i=1,num
   ionpoint=ntype(i)
   rmassvol=amass(ionpoint)/cellvol
   pint3(1,1)=pint3(1,1)+ rmassvol*vxn(i)*vxn(i)
   pint3(1,2)=pint3(1,2)+ rmassvol*vxn(i)*vyn(i)
   pint3(1,3)=pint3(1,3)+ rmassvol*vxn(i)*vzn(i)
   pint3(2,2)=pint3(2,2)+ rmassvol*vyn(i)*vyn(i)
   pint3(2,3)=pint3(2,3)+ rmassvol*vyn(i)*vzn(i)
   pint3(3,3)=pint3(3,3)+ rmassvol*vzn(i)*vzn(i)
enddo   

pint3=pint3+pint3c
pint3(2,1)=pint3(1,2)
pint3(3,1)=pint3(1,3)
pint3(3,2)=pint3(2,3)
 
pinst3=(pint3(1,1)+pint3(2,2)+pint3(3,3))/3.0d0
!......get new force on cell shape
vg3sq=MATMUL(vg3,vg3)

Fgo3=vol3*(pint3-(pext*riden)) &
        -((vol3/3.0d0)*(((pint3(1,1)-(riden(1,1)*pext)) &
                        +(pint3(2,2)-(riden(2,2)*pext)) &
                        +(pint3(3,3)-(riden(3,3)*pext)))*riden))
!.......update cell shape velocity 
vg3=com3+(Fgo3*Wgorec+vg3sq-vg3*vbzeta3(1))*(dtime/2.0d0)
 if(ortho)then
    vg3(1,2)=0.d0
    vg3(1,3)=0.d0
    vg3(2,3)=0.d0
    vg3(2,1)=0.d0
    vg3(3,1)=0.d0
    vg3(3,2)=0.d0
 endif

! GWW - why do this for nab ? 
trace=vg3(1,1)+vg3(2,2)+vg3(3,3)
vg3(1,1)=vg3(1,1)-(1.0d0/3.0d0)*trace
vg3(2,2)=vg3(2,2)-(1.0d0/3.0d0)*trace
vg3(3,3)=vg3(3,3)-(1.0d0/3.0d0)*trace

! GWW - why zeroed for nab ? 
if(.not.nib) then
vg3=0.d0
veps2=0.d0
veps3=0.d0
endif
!...........now update the velocities
xkesq=0.0d0
Pcomx=0.0d0
Pcomy=0.0d0
Pcomz=0.0d0
comx2=(0.5d0*dtime)*((ax3)-(vxn*vpzeta3(1)) &
       -(2.0d0*((vg3(1,1)*vxn)+(vg3(1,2)*vyn) &
       +(vg3(1,3)*vzn)))-((2.0d0+(3.0d0/free))*vxn*veps3))
comy2=(0.5d0*dtime)*((ay3)-(vyn*vpzeta3(1)) &
       -(2.0d0*((vg3(2,1)*vxn)+(vg3(2,2)*vyn) &
       +(vg3(2,3)*vzn)))-((2.0d0+(3.0d0/free))*vyn*veps3))
comz2=(0.5d0*dtime)*((az3)-(vzn*vpzeta3(1)) &
       -(2.0d0*((vg3(3,1)*vxn)+(vg3(3,2)*vyn) &
       +(vg3(3,3)*vzn)))-((2.0d0+(3.0d0/free))*vzn*veps3))
vxn=vxni+comx2
vyn=vyni+comy2
vzn=vzni+comz2
do i=1,num
   ionpoint=ntype(i)
 
   xkesq=xkesq+amass(ionpoint)*(vxn(i)*vxn(i)+ &
               vyn(i)*vyn(i)+vzn(i)*vzn(i))

   Pcomx=Pcomx+vxn(i)*amass(ionpoint)
   Pcomy=Pcomy+vyn(i)*amass(ionpoint)
   Pcomz=Pcomz+vzn(i)*amass(ionpoint)
enddo   

! GWW lots of COM resets later - but uses Ke - with tranlational energy ? 
! try correcting velecoties before using ! 
!
!totmass=SUM(amass*nsp(1:nspec))
!
!COMke=Pcomx**2+Pcomy**2+Pcomz**2
!COMke=COMke/(2.0d0*totmass)
!      
!if(COMke.gt.free*trantkb*1.0d-30) then
!   Vcomx=Pcomx/totmass
!   Vcomy=Pcomy/totmass
!   Vcomz=Pcomz/totmass
!
!   Pcomx=0.0d0
!   Pcomy=0.0d0
!   Pcomz=0.0d0
!   xkesq=0.0d0 
!   vxn=vxn-Vcomx
!   vyn=vyn-Vcomy
!   vzn=vzn-Vcomz
!   do i=1,num
!      ipoint=ntype(i)
!      Pcomx=Pcomx+vxn(i)*amass(ipoint)
!      Pcomy=Pcomy+vyn(i)*amass(ipoint)
!      Pcomz=Pcomz+vzn(i)*amass(ipoint)
!   xkesq=xkesq+amass(ionpoint)*(vxn(i)*vxn(i)+ &
!               vyn(i)*vyn(i)+vzn(i)*vzn(i))
!   enddo   
! 
!newCOMke=Pcomx**2+Pcomy**2+Pcomz**2
!newCOMke=newCOMke/(2.0d0*totmass)
!
!   if( iam .eq. 0 ) then
!
!   write(6,111) COMke,free*trantkb*1.0d-30 &
!                   ,Pcomx,Pcomy,Pcomz, newCOMke 
!  111  format('iCOM reset: COMke=',E12.6, &
!       ' Target=',E12.6,' mv(x,y,z)=',3(E12.6,X),' newCOMke=',E12.6) 
!
!   endif
!endif         
!
! GWW end of COM reset in iterative loop

!     Update the friction coefficient data
vg3sq=MATMUL(vg3,vg3)
trvg3=vg3sq(1,1)+vg3sq(2,2)+vg3sq(3,3)
!............update particle  zeta
gpzeta3(1)=(CUEprec)*((xkesq)-(free*trantkb))-vpzeta3(1)*vpzeta3(2)
gpzeta3(2)=(CUEp2rec)*(CUEp*vpzeta3(1)**2-trantkb)-vpzeta3(2)*vpzeta3(3)
gpzeta3(3)=(CUEp2rec)*(CUEp2*vpzeta3(2)**2-trantkb)-vpzeta3(3)*vpzeta3(4)
gpzeta3(4)=(CUEp2rec)*(CUEp2*vpzeta3(3)**2-trantkb)-vpzeta3(4)*vpzeta3(5)
gpzeta3(5)=(CUEp2rec)*(CUEp2*vpzeta3(4)**2-trantkb)
vpzeta3=vpzeta2+0.5d0*dtime*(gpzeta2+gpzeta3)
!............  update barostat zeta 
gbzeta3(1)=(CUEbrec)*(W*veps3**2+Wgo*trvg3-dom*trantkb)-vbzeta3(1)*vbzeta3(2)
vbzeta3(1)=vbzeta2(1)+0.5d0*dtime*(gbzeta2(1)+gbzeta3(1))
gbzeta3(2)=(CUEb2rec)*(CUEb*vbzeta3(1)**2-trantkb)-vbzeta3(2)*vbzeta3(3)
vbzeta3(2)=vbzeta2(2)+0.5d0*dtime*(gbzeta2(2)+gbzeta3(2))
gbzeta3(3)=(CUEb2rec)*(CUEb2*vbzeta3(2)**2-trantkb)-vbzeta3(3)*vbzeta3(4)
vbzeta3(3)=vbzeta2(3)+0.5d0*dtime*(gbzeta2(3)+gbzeta3(3))
gbzeta3(4)=(CUEb2rec)*(CUEb2*vbzeta3(3)**2-trantkb)-vbzeta3(4)*vbzeta3(5)
vbzeta3(4)=vbzeta2(4)+0.5d0*dtime*(gbzeta2(4)+gbzeta3(4))
gbzeta3(5)=(CUEb2rec)*(CUEb2*vbzeta3(4)**2-trantkb)
vbzeta3(5)=vbzeta2(5)+0.5d0*dtime*(gbzeta2(5)+gbzeta3(5))
!.............update epsilon
Feps3=3.0d0*vol3*(pinst3-pext)+((3.0d0/free)*xkesq)

veps3=veps2+(dtime/2.0d0)*(Feps2*Wrec- &
        veps2*vbzeta2(1))+(dtime/2.0d0)*(Feps3*Wrec-veps3*vbzeta3(1))
!.......this is a pretty crude convergence indicator
!.......might be better to test the individual velocities
converge=0.0d0 
do i=1,num
   ion=ntype(i)
   converge=converge+amass(ion)*((vxn(i)-vxold(i))**2 &
           +(vyn(i)-vyold(i))**2+(vzn(i)-vzold(i))**2)
enddo
if(converge.gt.(free*trantkb*1.0d-17)) then
!if(converge.gt.(free*trantkb*1.0d-30)) then
!---> Parallelization_S
!   if( iam .eq. 0 ) then
!
!      if (mod(nstep,20).eq.0) then
!         write(6,*)"convergence", nstep, icnt,converge,free*trantkb*1.0d-17
!      endif
!
!   endif
!<--- Parallelization_E
   if(icnt.gt.50) then
!---> Parallelization_S
      if( iam .eq. 0 ) then

      write(6,*)'Convergence not in under 50 steps. Abort.'

      endif
      call mpi_finalize(ierr)
!<--- Parallelization_E
      stop
   endif
   goto 380
endif

totmass=SUM(amass*nsp(1:nspec))

COMke=Pcomx**2+Pcomy**2+Pcomz**2
COMke=COMke/(2.0d0*totmass)
      
if(COMke.gt.free*trantkb*1.0d-30) then
   Vcomx=Pcomx/totmass
   Vcomy=Pcomy/totmass
   Vcomz=Pcomz/totmass

   Pcomx=0.0d0
   Pcomy=0.0d0
   Pcomz=0.0d0
   vxn=vxn-Vcomx
   vyn=vyn-Vcomy
   vzn=vzn-Vcomz
   do i=1,num
      ipoint=ntype(i)
      Pcomx=Pcomx+vxn(i)*amass(ipoint)
      Pcomy=Pcomy+vyn(i)*amass(ipoint)
      Pcomz=Pcomz+vzn(i)*amass(ipoint)
   enddo   
 
newCOMke=Pcomx**2+Pcomy**2+Pcomz**2
newCOMke=newCOMke/(2.0d0*totmass)

!---> Parallelization_S
   if( iam .eq. 0 ) then

   write(6,110) COMke,free*trantkb*1.0d-30 &
                   ,Pcomx,Pcomy,Pcomz, newCOMke 
  110  format('COM reset: COMke=',E12.6, &
       ' Target=',E12.6,' mv(x,y,z)=',3(E12.6,X),' newCOMke=',E12.6) 

   endif
!<--- Parallelization_E
endif         
!====================================================================
! Calculate the kinetic energy at time t. (i.e vmag() is the mean
! of the leapfroged v & vn velocities.)
!
!....  Note, these velocities are now in the cartesian frame
!.....and current with the positions and dipoles. i.e. do the energy
!......conservation at this point
vtx=vxn
vty=vyn
vtz=vzn
vmag=vxn*vxn+vyn*vyn+vzn*vzn
vxn=hi(1,1)*vtx+hi(1,2)*vty+hi(1,3)*vtz
vyn=hi(2,1)*vtx+hi(2,2)*vty+hi(2,3)*vtz
vzn=hi(3,1)*vtx+hi(3,2)*vty+hi(3,3)*vtz
if (moveions) then
  vx=vxn
  vy=vyn
  vz=vzn
end if

!!! !              VXN NOW CONTAIN NEW VELOCITIES IN THE CELL FRAME --
l=1
tranko=tranke
tranke=0.0d0
do i=1,nspec
   do j=1,nsp(i)
      tranke=tranke+(hmass(i)*vmag(l))
      l=l+1
   enddo   
enddo 

tcell=(Wgo/2.0d0)*trvg3
tvol=(W/2.0d0)*veps3**2
tpzeta=(CUEp/2.0d0)*vpzeta3(1)**2
tbzeta=(CUEb/2.0d0)*vbzeta3(1)**2
tpzeta=tpzeta+(CUEp2/2.d0)*SUM(vpzeta3(2:5)*vpzeta3(2:5))
tbzeta=tbzeta+(CUEb2/2.d0)*SUM(vbzeta3(2:5)*vbzeta3(2:5))

PeeVee= pext*vol3

pzeta=free*trantkb*pzeta3(1)
bzeta=dom*trantkb*bzeta3(1)
pzeta=pzeta+trantkb*SUM(pzeta3(2:5))
bzeta=bzeta+trantkb*SUM(bzeta3(2:5))
!...just add the PE to these variables -- should give the conserved quantity
!......Pcomx, Pcomy, Pcomz are available to monitor the total momentum
tke=tranke*gtrankin
!Dump msd information if required
if(msdcalllog) then
   if(mod(float(nstep),float(nmsdcalltime)).eq.0) then
!---> Parallelization_S
      if( iam .eq. 0 ) then       

      do i=1,num
         write(41,*)xdisp(i),ydisp(i),zdisp(i)
      enddo

      endif
!<--- Parallelization_E
      xdisp=0.d0
      ydisp=0.d0
      zdisp=0.d0
   endif
endif

return
END SUBROUTINE
