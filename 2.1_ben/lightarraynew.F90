subroutine lightarray

!........This routine contains modifications which will lead to
!........a better isotropic spectrum

!************************************************************
!
!  Following Madden & (Board,O`Sullivan,Fowler,Seddon) we are looking
! at the following 6 symmetries:
!              i)  xy                   [anisotrophic]
!             ii)  xz                   [anisotrophic]
!            iii)  yz                   [anisotrophic]
!             iv)  xx + yy + zz         [isotropic]
!              v)  xx - yy              [anisotropic]
!             iv)  zz - yy - xx         [anisotropic]
!  We consider the following contributions:
!              i)  Short range
!             ii)  Dipole-induced-dipole
!            iii)  B-terms
!             iv)  Gamma-terms.
!
!************************************************************


use commondata, only: num,nanion,nsp,x,y,z
use lightdata, only: hyperB,hypergamma,polarundist,ncfmat2, & 
                 root3,root12,root1twodiv3gam,twodiv3gamma &
                  ,root3div3gam,gammadiv3,bhyperpol3, &
                     exxls,eyyls,ezzls,exyls,exzls,eyzls,srxx,sryy,srzz, &
                     srxy,srxz,sryz,elecxu,elecyu,eleczu,term, &
                  txxli,tyyli,tzzli,txyli,txzli,tyzli,cisokr,cisoki, &
                   cohacc,xkvec2
use boxdata, only: twopiboxx,twopiboxy,twopiboxz

implicit none

! *** define our local variables
double precision :: XKVECRI,CSL,SL
integer :: i,j,k,ISYM,ICONTTYPE,L
! *** end of local definitions



!================================================

do i=1,nanion

!=====================================================================
!  Apply field gradient term hyperpolarizability (B) to the field
! gradient.
! [ bhyperbol3 = B/3 ]
   exxls(i)=exxls(i)*bhyperpol3
   eyyls(i)=eyyls(i)*bhyperpol3
   ezzls(i)=ezzls(i)*bhyperpol3
   exyls(i)=exyls(i)*bhyperpol3
   exzls(i)=exzls(i)*bhyperpol3
   eyzls(i)=eyzls(i)*bhyperpol3
enddo 

!=====================================================================
!  Zero cation contributions to B-term, gamma term and SR term

do i=nanion+1,num
   exxls(i)=0.0d0
   eyyls(i)=0.0d0
   ezzls(i)=0.0d0
   exyls(i)=0.0d0
   exzls(i)=0.0d0
   eyzls(i)=0.0d0

   srxx(i)=0.0d0
   sryy(i)=0.0d0
   srzz(i)=0.0d0
   srxy(i)=0.0d0
   srxz(i)=0.0d0
   sryz(i)=0.0d0

   elecxu(i)=0.0d0
   elecyu(i)=0.0d0
   eleczu(i)=0.0d0
enddo

!=====================================================================
!  Zero coherent term accumulator.
!  isym is the number of symmetries being reported. 
!  iconttype is the number of contributions calculated.
!write(6,*) ' from lightarray'
                                                                                
!write(6,*)  exxls(1),eyyls(1),ezzls(1)
!write(6,*)  exyls(1),exzls(1),eyzls(1)
!write(6,*) ' '
!write(6,*) txxli(1),tyyli(1),tzzli(1)
!write(6,*) txyli(1),txzli(1),tyzli(1)

!write(6,*) 'elecxu(1)',elecxu(1),elecyu(1),eleczu(1)

do isym=1,6
   do iconttype=1,4

      cohacc(isym,iconttype)=0.0d0

   enddo
enddo

do  l=1,6
    do j=1,ncfmat2
            
       cisokr(l,j)=0.0d0
       cisoki(l,j)=0.0d0

    enddo
enddo

!****************************************************************


do i=1,num

!***********SR terms

   term(1,1)=root12*srxy(i)
   term(2,1)=root12*srxz(i)
   term(3,1)=root12*sryz(i)
   term(4,1)=srxx(i)+sryy(i)+srzz(i)
   term(5,1)=root3*(srxx(i)-sryy(i))
   term(6,1)=2.0d0*srzz(i)-srxx(i)-sryy(i)

!*********** DID terms.
   term(1,2)=-root12*txyli(i)
   term(2,2)=-root12*txzli(i)
   term(3,2)=-root12*tyzli(i)
   term(4,2)=-(txxli(i)+tyyli(i)+tzzli(i))
   term(5,2)=-(root3*(txxli(i)-tyyli(i)))
   term(6,2)=-(2.0d0*tzzli(i)-txxli(i)-tyyli(i))

!************ B terms.
   term(1,3)=root12*exyls(i)
   term(2,3)=root12*exzls(i)
   term(3,3)=root12*eyzls(i)
   term(4,3)=exxls(i)+eyyls(i)+ezzls(i)
   term(5,3)=root3*(exxls(i)-eyyls(i))
   term(6,3)=2.0d0*ezzls(i)-exxls(i)-eyyls(i)

!************ Gamma terms.
!  [ root1twodiv3gam = (sqrt(12)/3)*gamma) ]
!  [ twodiv3gam = (2/3)*gamma) ]
!  [ root3div3gam = (sqrt(3)/3)*gamma) ]
!  [ gammadiv3 = (1/3)*gamma) ]
   term(1,4)=root1twodiv3gam*elecxu(i)*elecyu(i)
   term(2,4)=root1twodiv3gam*elecxu(i)*eleczu(i)
   term(3,4)=root1twodiv3gam*elecyu(i)*eleczu(i)
   term(4,4)=twodiv3gamma*(elecxu(i)*elecxu(i) &
           +elecyu(i)*elecyu(i)+eleczu(i)*eleczu(i))
   term(5,4)=root3div3gam*(elecxu(i)*elecxu(i) &
           -elecyu(i)*elecyu(i))
   term(6,4)=gammadiv3*(2.0d0*eleczu(i)*eleczu(i) &
            -elecxu(i)*elecxu(i) -elecyu(i)*elecyu(i))

!=====================================================================
!  Accumulate terms in cohacc.
   do j=1,6
      do k=1,4

         cohacc(j,k)=cohacc(j,k)+term(j,k)

      enddo
   enddo 




!..... begin the new isotropic calculation. 
!......... Note the isotropic polarizability terms 
!........are always on anions --- hence

!............cisokr(1,..) contains the short-range polarizability density
!............cisokr(2,..)  contains gamma term
!.............cisokr(3,...) contains the anion density 
!............cisokr(4,...) contains the first cation density 
!............cisokr(5,...) contains the second cation density 
!............cisokr(6,...) contains the third cation density
   if(i.le.nanion) then
      do j=1,ncfmat2

         xkvecri=twopiboxx*xkvec2(1,j)*x(i) &
                +twopiboxy*xkvec2(2,j)*y(i) &
                +twopiboxz*xkvec2(3,j)*z(i)
         csl=cos(xkvecri)
         sl=sin(xkvecri)
           
         cisokr(1,j)=cisokr(1,j)+csl*term(4,1)
         cisoki(1,j)=cisoki(1,j)+sl*term(4,1)
         cisokr(2,j)=cisokr(2,j)+csl*term(4,4)
         cisoki(2,j)=cisoki(2,j)+sl*term(4,4)
         cisokr(3,j)=cisokr(3,j)+csl
         cisoki(3,j)=cisoki(3,j)+sl
      enddo

   else if(i.le.(nsp(2)+nanion)) then
       do j=1,ncfmat2

          xkvecri=twopiboxx*xkvec2(1,j)*x(i) &
                 +twopiboxy*xkvec2(2,j)*y(i) &
                 +twopiboxz*xkvec2(3,j)*z(i)
  
          csl=cos(xkvecri)
          sl=sin(xkvecri)
         
          cisokr(4,j)=cisokr(4,j)+csl
          cisoki(4,j)=cisoki(4,j)+sl
       enddo

   else if(i.le.(nsp(3)+nsp(2)+nanion)) then
       do j=1,ncfmat2

          xkvecri=twopiboxx*xkvec2(1,j)*x(i) &
                 +twopiboxy*xkvec2(2,j)*y(i) &
                 +twopiboxz*xkvec2(3,j)*z(i)
  
          csl=cos(xkvecri)
          sl=sin(xkvecri)
           
          cisokr(5,j)=cisokr(5,j)+csl
          cisoki(5,j)=cisoki(5,j)+sl
       enddo

   else
       do  j=1,ncfmat2

           xkvecri=twopiboxx*xkvec2(1,j)*x(i) &
                   +twopiboxy*xkvec2(2,j)*y(i) &
                   +twopiboxz*xkvec2(3,j)*z(i)

            csl=cos(xkvecri)
            sl=sin(xkvecri)
           
            cisokr(6,j)=cisokr(6,j)+csl
            cisoki(6,j)=cisoki(6,j)+sl
        enddo
   endif

enddo
!write(6,*) 'cohacc(1,1)',  cohacc(1,1)
!write(6,*)  'cisokr(1,1),cisoki(1,1)',cisokr(1,1),cisoki(1,1)

return

end subroutine



