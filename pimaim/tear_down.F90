MODULE tear_down

! Contains routines for use if pimaim is
! spawned by the MPI Spawn command, used by
! ppfit

IMPLICIT NONE

CONTAINS


    SUBROUTINE free_arrays

    use error_function
    use mpipara
    USE commondata
    USE boxdata
    USE recipdata, ONLY: bdisp, sk_ds, sk_ds_w,elcall,elsall,&
                        emcall, emsall, encall, ensall,kmaxx,kmaxy,kmaxz 
                       
    USE lightdata
    USE cgdata


    IMPLICIT NONE


    DEALLOCATE ( x,y,z )
    DEALLOCATE ( vx,vy,vz )
    DEALLOCATE ( frrx,frry,frrz )
    DEALLOCATE ( frrx3,frry3,frrz3 )
    DEALLOCATE ( xdisp,ydisp,zdisp )
    DEALLOCATE ( q, ntype )
    DEALLOCATE ( elecx,elecy,elecz )
    DEALLOCATE ( elecxq,elecyq,eleczq )
    DEALLOCATE ( elecxsr,elecysr,eleczsr )
    DEALLOCATE ( exx,eyy,ezz,exy,exz,eyz )
    DEALLOCATE ( exxq,eyyq,ezzq,exyq,exzq,eyzq )
    DEALLOCATE ( exxsr,eyysr,ezzsr,exysr,exzsr,eyzsr )
    DEALLOCATE ( xk1,xk2,xk3,xk4 )
    DEALLOCATE ( alppolar,Bpolar,Cpolar,gammapolar )
    DEALLOCATE ( engeff,xmu,ymu,zmu,quadxx,quadyy, &
                quadzz,quadxy,quadxz,quadyz,delta, &
                epsilonx,epsilony,epsilonz,quaimxx, &
                quaimyy,quaimzz,quaimxy,quaimxz,quaimyz )
    DEALLOCATE ( selfeps, selfquaim )
    DEALLOCATE ( reseng, dipsqeng, dipquadeng, quadeng )
    DEALLOCATE ( srdipx,srdipy,srdipz )
    DEALLOCATE ( asdipx,asdipy,asdipz )
    DEALLOCATE ( srquadxx,srquadyy,srquadzz, &
               srquadxy,srquadxz,srquadyz )
    DEALLOCATE ( asquadxx,asquadyy,asquadzz, &
               asquadxy,asquadxz,asquadyz )
    DEALLOCATE ( PCOMAIM,XICOMAIM ) !Maybe this is not good? AGUADO
    DEALLOCATE ( PCOMSTRUCT,XICOMSTRUCT)

    DEALLOCATE ( erfc )
    DEALLOCATE ( dxsav, dysav, dzsav )

    DEALLOCATE ( engft1,engft2,engft3, &
                engft1dotx,engft1doty,engft1dotz, &
                engft2dotx,engft2doty,engft2dotz, &
                engft3dotx,engft3doty,engft3dotz, &
                engft1dotxx,engft1dotyy,engft1dotzz, &
                engft1dotxy,engft1dotxz,engft1dotyz, &
                engft2dotxx,engft2dotyy,engft2dotzz, &
                engft2dotxy,engft2dotxz,engft2dotyz, &
                engft3dotxx,engft3dotyy,engft3dotzz, &
                engft3dotxy,engft3dotxz,engft3dotyz )

    DEALLOCATE ( nsp,chg,spectype,amass, &
             polarizablelog,deformablelog )
    DEALLOCATE ( ftalp,ftb,ftc, &
             ftd,dddamp,dqdamp  )
    DEALLOCATE ( rvph,rvpn,rvpr4 )
    DEALLOCATE ( ftbeta,ftgamma,ftb2, ftb3  )
    DEALLOCATE ( selfB,selfgam )
    DEALLOCATE ( ftalpx, ftbx, nrpower )
    DEALLOCATE ( selfeps1, selfeps2, selfC )
    DEALLOCATE ( selfquaim1, selfH )
    DEALLOCATE ( dampa, dampfac, nkdamp, &
               fgb, fgc, nkfg  )
    DEALLOCATE ( alppolar1, alppolardel, alppolareff, &
               Bpolar1, Cpolar1, gammapolar1  )
    DEALLOCATE ( nunitcellpos, cellcoordfile )
    DEALLOCATE ( vartrans, recamass, hmass )
    DEALLOCATE ( dddamp2,dddamp3,dddamp4, &
               dddamp5,dddamp6,dqdamp2, &
               dqdamp3,dqdamp4,dqdamp5, &
               dqdamp6,dqdamp7,dqdamp8 )
    DEALLOCATE ( bdisp )
    DEALLOCATE ( rdfpart )
    DEALLOCATE ( rdftot_w )
    DEALLOCATE ( rdfpart_w )
    DEALLOCATE ( nmon, atmnam, weight, chge)
    DEALLOCATE ( numadr)
    DEALLOCATE ( eltmp, elltmp )
    DEALLOCATE ( sctmp, scctmp )
    DEALLOCATE ( dimtmp )
    DEALLOCATE (erfc_table)
    DEALLOCATE (gw_force )

    if(lscalc) then
        DEALLOCATE (txxli,tyyli,tzzli, &
                  txyli,txzli,tyzli)
        DEALLOCATE (elecxu,elecyu,eleczu)
        DEALLOCATE (exxls,eyyls,ezzls,   &
                   exyls,exzls,eyzls)
        DEALLOCATE (exxsrls,eyysrls,ezzsrls,  &
                    exysrls,exzsrls,eyzsrls)
        DEALLOCATE (srxx,sryy,srzz,   &
                       srxy,srxz,sryz)
        DEALLOCATE ( cohacc,cohav,term)
        DEALLOCATE (lnorm)
        DEALLOCATE (scf)
        DEALLOCATE (storecohacc)
        DEALLOCATE ( cisokr,cisoki)
        DEALLOCATE ( storeisokr, storeisoki)
        DEALLOCATE ( scfiso)
        DEALLOCATE (xkvec2)
    end if

    DEALLOCATE( ist,  ied )
    DEALLOCATE( ist2, ied2 )
    DEALLOCATE( ist3, ied3 )
    DEALLOCATE( ist4, ied4 )
    DEALLOCATE( sk_ds, sk_ds_w )
    
    DEALLOCATE ( elcall, elsall )
    DEALLOCATE ( emcall, emsall )
    DEALLOCATE ( encall, ensall )

    DEallocate( kmaxy_s, kmaxy_e )
    DEallocate( kmaxz_s, kmaxz_e )



    !remove this for ppfit version
    !if (iam .eq. 0)  call system("rm *out*")
    !!if (iam .eq. 0) call system("rm ev_ang_potential")
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
   
end subroutine

END MODULE tear_down

