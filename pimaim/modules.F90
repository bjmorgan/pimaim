MODULE COMMONDATA

IMPLICIT NONE
SAVE

! extra array for checking forces

double precision, allocatable :: gw_force(:) 
! Constants

DOUBLE PRECISION, PARAMETER :: PI=3.14159265358979323d0
DOUBLE PRECISION, PARAMETER :: PISQ=PI*PI
DOUBLE PRECISION, PARAMETER :: SQRPI=1.7724538509055160249d0
DOUBLE PRECISION, PARAMETER :: PITHREEHALF=SQRPI*PI
DOUBLE PRECISION, PARAMETER :: TWOPI=2.D0*PI 
DOUBLE PRECISION, PARAMETER :: FOURPI=4.D0*PI
DOUBLE PRECISION, PARAMETER :: EIGHTPI=8.D0*PI
DOUBLE PRECISION, PARAMETER :: FOURPISQ=4.D0*PISQ
DOUBLE PRECISION, PARAMETER :: ONETHIRD=1.D0/3.D0

! Conversion of units 

DOUBLE PRECISION, PARAMETER ::  &
        emass=9.109534d-28,     &      !Electron mass in grams
        boltz=3.166829689d-6,   &      !(Inverse Kelvins)
        avo=6.022045d023               !Avogadro's number

! Functions

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: X,Y,Z     ! particle coordinates
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: VX,VY,VZ  ! particle velocities
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: FRRX,FRRY,FRRZ  ! forces on particles 
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: FRRX3,FRRY3,FRRZ3  ! forces on particles 
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: xdisp,ydisp,zdisp !mean square displacement
! in fact they are the internal coordinates, ranging from 0 to 1
! ie   r=x a1 + y a2 + z a3  NOT r=x i + y j + z k
!---> Memmory Reduction_S
!DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: dxsav,dysav,dzsav !Interatomic distances
!DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: erfc   !Error function (Change name??)
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: dxsav,dysav,dzsav !Interatomic distances
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: erfc   !Error function (Change name??)
!---> Memmory Reduction_E
DOUBLE PRECISION :: stsrxx,stsrxy,stsrxz,stsryy,stsryz,stsrzz,stcxx,stcxy,stcxz,stcyy,stcyz,stczz
DOUBLE PRECISION :: stpxx,stpxy,stpxz,stpyy,stpyz,stpzz
DOUBLE PRECISION :: stp2xx,stp2xy,stp2xz,stp2yy,stp2yz,stp2zz,stpsrxx,stpsrxy,stpsrxz,stpsryy,stpsryz,stpsrzz
DOUBLE PRECISION :: stewxx,stewxy,stewxz,stewyy,stewyz,stewzz,stewyx,stewzx,stewzy
DOUBLE PRECISION :: stpqquadxx,stpqquadyy,stpqquadzz,stpqquadxy,stpqquadxz,stpqquadyz
DOUBLE PRECISION :: stpdipquadxx,stpdipquadyy,stpdipquadzz,stpdipquadxy,stpdipquadxz,stpdipquadyz
DOUBLE PRECISION :: stpquadquadxx,stpquadquadyy,stpquadquadzz,stpquadquadxy,stpquadquadxz,stpquadquadyz
DOUBLE PRECISION :: stqsrxx,stqsryy,stqsrzz,stqsrxy,stqsrxz,stqsryz,pext,pextstruct
                         !Components of the stress tensor (Change in future?)
DOUBLE PRECISION, DIMENSION(3,3) :: pint2   !Internal pressure tensor

DOUBLE PRECISION :: dtime   !Time step

INTEGER, ALLOCATABLE, DIMENSION(:) :: ntype !1 for 1st comp, 2 for 2nd comp ...
INTEGER, ALLOCATABLE, DIMENSION(:) :: nsp   !Number of ions of each species...
INTEGER, ALLOCATABLE, DIMENSION(:) :: nunitcellpos !Number of ions of each species in the unit cell
INTEGER, ALLOCATABLE, DIMENSION(:) :: nmon   !targets ions for which to show detailed information...
!---> Memmory Reduction_S
!INTEGER :: num,nanion,ncation,nspec,nummon,nunitcellx,nunitcelly,nunitcellz,nionunit
INTEGER :: num,nanion,ncation,nspec,nummon,nunitcellx,nunitcelly,nunitcellz,nionunit,num2,numx,num_s,num_s2,num_s3,num_s4
INTEGER :: istrj,jstrj,keytrj                          	!MABC: DL_POLY-like trajectory output control
!<--- Memmory Reduction_E

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: elecx,elecy,elecz !Electric fields
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: elecxq,elecyq,eleczq !Electric fields
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: elecxsr,elecysr,eleczsr !Electric fields (short-range contribution)
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: exx,eyy,ezz,exy,exz,eyz !Electric field gradients
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: exxq,eyyq,ezzq,exyq,exzq,eyzq !Electric field gradients
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: exxsr,eyysr,ezzsr,exysr,exzsr,eyzsr !Electric field gradients (short-range contribution)

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: weight,amass,hmass,recamass  !MABC: weight contains individual species mass
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: chg,q,chge                 !MABC: chge contains individual species charge
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: xmu,ymu,zmu     !dipoles
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: srdipx,srdipy,srdipz    !SR dipoles
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: asdipx,asdipy,asdipz    !LR dipoles
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: quadxx,quadyy,quadzz,quadxy,quadxz,quadyz !quadrupoles
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: srquadxx,srquadyy,srquadzz,srquadxy,srquadxz,srquadyz    !SR quadrupoles
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: asquadxx,asquadyy,asquadzz,asquadxy,asquadxz,asquadyz    !LR quadrupoles

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: xk1,xk2,xk3,xk4,engeff !polarizabilites
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: alppolar,Bpolar,Cpolar,gammapolar
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: alppolar1,alppolardel,alppolareff
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: Bpolar1,Cpolar1,gammapolar1

DOUBLE PRECISION :: engpetot, tranke, tke, trantemp !Potential and kinetic energies, and temperatures
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: reseng,quadeng,dipquadeng,dipsqeng !PIM
DOUBLE PRECISION :: selfengtot,epsselfengtot,quaimselfengtot !CIM/AIM selfenergies

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: delta           !CIM distortion
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: epsilonx,epsilony,epsilonz !DAIM distortion
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: quaimxx,quaimyy,quaimzz,quaimxy,quaimxz,quaimyz     !QUAIM distortion

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: ftalp,ftbeta,ftgamma,ftb,ftb2,ftb3
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: rvph,rvpn,rvpr4,ftalpx,ftbx
INTEGER         , ALLOCATABLE, DIMENSION(:,:) :: nrpower
INTEGER                                       :: nrpower2
                                              ! Short-range repulsion potential
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: engft1,engft2,engft3
                                !CIM CG "Forces"
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: engft1dotx,engft1doty,engft1dotz
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: engft2dotx,engft2doty,engft2dotz
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: engft3dotx,engft3doty,engft3dotz
                                !DAIM CG "forces"
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: engft1dotxx,engft1dotyy,engft1dotzz
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: engft1dotxy,engft1dotxz,engft1dotyz
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: engft2dotxx,engft2dotyy,engft2dotzz
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: engft2dotxy,engft2dotxz,engft2dotyz
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: engft3dotxx,engft3dotyy,engft3dotzz
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: engft3dotxy,engft3dotxz,engft3dotyz
                                !QUAIM CG "Forces"

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: ftc,ftd         !(dispersion potential)
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: dddamp,dddamp2,dddamp3,dddamp4,dddamp5,dddamp6    !Dipole-dipole dispersion damping
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: dqdamp,dqdamp2,dqdamp3,dqdamp4,dqdamp5,dqdamp6,dqdamp7,dqdamp8    !Dipole-quadrupole dispersion damping

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: dampa,dampfac,fgb,fgc !polarization damping
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: nkdamp,nkfg            !polarization damping

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: selfgam,selfB !CIM/AIM selfenergy parameters
DOUBLE PRECISION :: extraalpha,extrab
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: selfeps,selfC,selfeps1,selfeps2
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: selfquaim,selfH,selfquaim1

INTEGER, DIMENSION(100) :: nstpbrk
INTEGER :: nmat1,nmat2
INTEGER :: nmsdcalltime,nperrdf,nrdfcall,nrscale,npereng,npervel,nperfri,npercell
INTEGER :: nsofar,ntotstp,nrun,nstep
INTEGER, DIMENSION(6) :: cellconstraints !for structural relaxation of cell
DOUBLE PRECISION :: trantkb, gtrankin, amovefac
REAL    :: dummy 
INTEGER,DIMENSION(1) :: iseed1,iseed2                                        ! MABC: These three lines contain elements for control
INTEGER, DIMENSION(8) :: time_array                     ! of random number generator. 
LOGICAL :: lseed
DOUBLE PRECISION, DIMENSION(3,3) :: amove
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: vartrans
DOUBLE PRECISION :: chgcorrec,dispcorrec1,dispcorrec2   !Ewald selfenergies
DOUBLE PRECISION :: etainpt,eta,etapi,etaconst,etasq    !Ewald smearing parameters
DOUBLE PRECISION :: rsqmax,rsqmaxsr,rcut,rcutsr       !Several cutoffs
DOUBLE PRECISION :: dispalpsq,dispalp,dispalp3,dispalp4,dispalp6  !For Dispersion Ewald summation

CHARACTER(len=20), ALLOCATABLE, DIMENSION(:) :: cellcoordfile   !X.mat, etc... Will this work??
CHARACTER(len=8), ALLOCATABLE, DIMENSION(:) :: atmnam,spectype           !MABC: arrays with atom/species name
CHARACTER(len=20) :: fileout
CHARACTER(len=150) :: jobname
CHARACTER(len=10) :: runtype,runtype2,pottype,ewtype,rstchar  !MABC

DOUBLE PRECISION, DIMENSION(0:300) :: rdftot      !Total and partial g(r)'s
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: rdfpart

DOUBLE PRECISION :: relax,free,eps1,eps2,eps3,veps1,veps2,veps3,vol3,relaxb,relaxb2 !Thermostat-Barostat
DOUBLE PRECISION, DIMENSION(5) :: vpzeta1,vpzeta2,vpzeta3,pzeta2,pzeta3    !Thermostat
DOUBLE PRECISION, DIMENSION(5) :: vbzeta1,vbzeta2,vbzeta3,bzeta2,bzeta3    !barostat
DOUBLE PRECISION :: tcell,tvol,tbzeta,tpzeta,PeeVee,pzeta,bzeta,Pcomx,Pcomy,Pcomz
DOUBLE PRECISION :: CUEp,CUEp2,CUEprec,CUEp2rec,dom,CUEb,CUEb2,CUEbrec,CUEb2rec,Wgo,W,Wrec,Wgorec

INTEGER :: nsteptramp, nsteppramp ! temperature/pressure ramp
DOUBLE PRECISION :: deltatramp, deltapramp, ewcorr

DOUBLE PRECISION :: tol,ftol,tolaim,ftolaim,tolstruc,ftolstruc,engtol !Tolerances

LOGICAL :: rim,dippimlog,quadpimlog,cimlog,daimlog,quaimlog,rimlog,epplog,xftlog,rvplog !Potential function
LOGICAL :: ooaimlog,oodamplog,cluslog,environmentalaimlog,environmentalpimlog !Levels of theory
LOGICAL, ALLOCATABLE, DIMENSION(:) :: polarizablelog,deformablelog !Level of theory for each species
LOGICAL :: displace,restart,rescalelog,velinit,relaxconfig,dynam,moveions,relaxcell     !Simulation options
LOGICAL :: conjgradlog,msdcalllog,ewlog,conjgradaimlog
LOGICAL :: chargequadlog,dipquadlog,fullewaldlog    !Level of Ewald summation
LOGICAL :: veldumplog,crddumplog,chgdumplog,fulldumplog,rdfcall   !Write to disk
LOGICAL :: pereng,perfric,pervel,percell,nrscalelog             !Write to disk
LOGICAL :: forfl,nth,nib,nab,ortho      !Variable cell
LOGICAL :: verbose  !to control the amount of information written to disk
LOGICAL :: firstiter !To optimize the number of calls to cgener...
LOGICAL :: tramplog,pramplog ! T/p ramp
LOGICAL :: lscalc   !Do the light scattering calculation

LOGICAL :: ppfit_dipoles, ppfit_forces, ppfit_stresses ! does ppfit want forces, dipoles and stresses
INTEGER :: exit_code ! for broadcasting exit info to ppfit - lets it know if cg failed
INTEGER :: ppfit_stop
REAL,DIMENSION(6) :: stress_tensor

END MODULE COMMONDATA

!----------------------------------------!

MODULE LIGHTDATA
IMPLICIT NONE
SAVE

!===============================================
INTEGER :: lsint,ncorrtime,ncfmat2,nkmod2,ncorr,ncorrcall

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: xkvec2
double precision, allocatable,dimension(:) :: txxli,tyyli,tzzli, &
             txyli,txzli,tyzli
!.........contains DID polarizabilities (num)

double precision, allocatable, dimension(:) :: elecxu,elecyu,eleczu
!.........contains electric fields without the sr terms (num)

double precision, allocatable, dimension(:) :: exxls,eyyls,ezzls,  &
                                       exyls,exzls,eyzls
double precision, allocatable, dimension(:) :: exxsrls,eyysrls,ezzsrls,  &
                    exysrls,exzsrls,eyzsrls
!.........contains electric field gradients (num)

double precision, allocatable, dimension(:) :: srxx,sryy,srzz,   &
                       srxy,srxz,sryz
!.........contains short-range polarizabilities (num)
 
double precision, allocatable, dimension(:,:) :: cohacc,cohav,term
!              accumulators (nsymmax,ncontribmax)
integer, allocatable, dimension(:,:,:) :: lnorm
! normaliser         (nsymmax,nindexmax,0:ncorrtimemax)
double precision, allocatable, dimension(:,:,:) :: scf
!  correlations functions (nsymmax,nindexmax,0:ncorrtimemax)
double precision, allocatable, dimension(:,:,:) :: sc
!  correlations functions (7,3,0:ncorrtimemax)
double precision, allocatable, dimension(:,:,:) :: storecohacc
!          storage (nsymmax,ncontribmax,ncorrtimemax)
double precision, allocatable, dimension(:,:) :: cisokr,cisoki
!.........raman amplitudes (6,ncfmat)
double precision, allocatable, dimension(:,:,:) :: storeisokr, storeisoki
!........storage (6,ncfmatmax2,0:ncorrtimemax)
double precision, allocatable, dimension(:,:,:) :: scfiso
!...............isotropic amplitudes ((21,ncfmatmax2,0:ncorrtimemax)
!================================================
! Double precisions in common.
double precision :: hyperB,hypergamma,polarundist 

double precision :: root3,root12,root1twodiv3gam,twodiv3gamma &
                  ,root3div3gam,gammadiv3,bhyperpol3

double precision, dimension (4) :: sig,xlipol
double precision, dimension (4,4) :: asr,bsr,csr,dsr      

integer :: normav

!.................end of light scattering

end module lightdata
!===================================================================================



MODULE BOXDATA

IMPLICIT NONE
SAVE

DOUBLE PRECISION, DIMENSION(3,3) :: h,hi,fullhi,fullhit,fullh !Several cell matrices
DOUBLE PRECISION, DIMENSION(3,3) :: h2, hi2, h3, hi3, vg2, vg3, hlab2, hlab3       
DOUBLE PRECISION, DIMENSION(10) :: bh, b, bee

DOUBLE PRECISION :: boxlenx,boxleny,boxlenz,halfboxx,halfboxy,halfboxz
DOUBLE PRECISION :: halfboxxrec,halfboxyrec,halfboxzrec,halfboxminsq
DOUBLE PRECISION :: twopiboxx,twopiboxy,twopiboxz

DOUBLE PRECISION :: a0,b0,c0,cellvol,fourpicell,cellvol3rec,fac

END MODULE BOXDATA

!----------------------------------------!

MODULE RECIPDATA

IMPLICIT NONE
SAVE

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: elcall,emcall,encall,elsall,emsall,ensall
                                     !sines and cosines for recipE (Change in future??)
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: bdisp !for recipr. dispersion energy

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: sk_ds !Structure factor
INTEGER :: nbin !Structure factor related
INTEGER, DIMENSION(0:1000) :: norm_ds !Structure factor related
!---> Parallelization_S
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: sk_ds_w
INTEGER, DIMENSION(0:1000) :: norm_ds_w
INTEGER :: nspairs
!<--- Parallelization_E

DOUBLE PRECISION :: qmueng,xmumueng,qquadeng,dipquadengrec,quadquadengrec
                                !Reciprocal space energy contributions
INTEGER :: kmaxx,kmaxy,kmaxz,ksqmax    !Reciprocal space maximum vectors...
DOUBLE PRECISION :: conv1,convfac  !Convergence factors for reciprocal space summations...
DOUBLE PRECISION :: rksqmax     !Reciprocal space cutoff
DOUBLE PRECISION :: slens(10) !Scattering lens !Added by D. Marrocchelli 05/03/2008

END MODULE RECIPDATA
!----------------------------------------!

MODULE CGDATA    

IMPLICIT NONE
SAVE

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: PCOMAIM,XICOMAIM
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: PCOMSTRUCT,XICOMSTRUCT

END MODULE CGDATA

!---> Parallelization_S
module mpipara

implicit none
include 'mpif.h'
save

integer :: nprocs, iam, ierr
integer :: kmaxx_s,kmaxx_e
integer :: jst,    jed
integer :: jst2,   jed2
integer :: jst3,   jed3
integer :: jst4,   jed4
integer, allocatable, dimension(:)   :: kmaxy_s,kmaxy_e
integer, allocatable, dimension(:,:) :: kmaxz_s,kmaxz_e
integer, allocatable, dimension(:)   :: ist,ied
integer, allocatable, dimension(:)   :: ist2,ied2
integer, allocatable, dimension(:)   :: ist3,ied3
integer, allocatable, dimension(:)   :: ist4,ied4
integer, allocatable, dimension(:,:) :: numadr
double precision, allocatable, dimension(:)     :: eltmp,elltmp
double precision, allocatable, dimension(:)     :: sctmp,scctmp
double precision, allocatable, dimension(:,:)   :: dimtmp
double precision, allocatable, dimension(:)     :: rdftot_w
double precision, allocatable, dimension(:,:,:) :: rdfpart_w
CHARACTER(len=256)::hostname, parentdir
INTEGER:: hosterror, hostlength

end module mpipara
!<--- Parallelization_E
