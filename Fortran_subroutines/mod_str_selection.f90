!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module mod_str_selection
use commondat_mod
use quality_mod
use mod_symmetry_cal_pi
use mod_symmetry_cal_sig
use mod_rumer_structures
use mod_quality_factor
use mod_qult_str_arrange
use mod_write_rumer_xmi
use mod_vector_rep 
use mod_write_symm_xmi_new
use mod_main_bond_cal
use mod_check_str_bond
use mod_eq_dstr_set
use rum_rad_mod
use mod_All_Rumer_set
use final_str_mod
implicit none

contains
subroutine str_selection(astr,nl,m4,perm_nstr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer::i1,i2,i3,i4,i5,i6,i7,ii,nl,str1_cnt,cnt,m4,m8,m9,lp_cnt,&
elporb,perm_nstr,total_str,rep,allowed_nstr
integer::str_cnt(15000), lps(50)
integer::qulsymm(15000),symqq(15000),tqlty1,bqlty1,sqlty1, factorial
integer, pointer::astr(:,:), str1(:,:), str2(:,:)
integer, pointer:: fvec(:,:),q_fac(:)


if(.not. allocated(rumer))then
   allocate(rumer(MaxStrOepo))
endif
if(.not. allocated(rumer_rad))then
   allocate(rumer_rad(MaxStrOepo))
endif

print*,'enter str_selection'
tqlty1=0
bqlty1=0
sqlty1=0

!**********************************************************************************************************
!!!!!!!!!!!!!!!! If number of permisible structures and available structures are same START !!!!!!!!!!!!!!!
!**********************************************************************************************************

if(perm_nstr.eq.m4)then

   call rumer_structures(nl,astr,m4)
   call quality_factor(nl,astr,m4,quality_fac,str_quality_1,str_quality_2,bondq, mbondq, pref_radical)
   call write_structures_file(str1_cnt, str2)!,q_fac)!, mbondq, pref_radical)

   if(input_flg.eq.1)call nnat_bond_cal(nl,astr,m4,bondq)
   if(input_flg.eq.0)call nnat_bond_cal_2(nl,astr,m4,bondq)
   tqlty=0
   bqlty=0
   sqlty=0
   do m8=1,m4
      if(niao.eq.0)then
         write(9,900)str_quality_1(m8),bondq(m8),str_quality_2(m8),'|',(astr(m8,m9),m9=1,nae)
      endif
      if(niao.gt.1)then
         write(9,901)str_quality_1(m8),bondq(m8),str_quality_2(m8),'|',1,':',niao,(astr(m8,m9),m9=1,nae)
      endif
      if(niao.eq.1)then
         write(9,901)str_quality_1(m8),bondq(m8),str_quality_2(m8),'|',1,1,(astr(m8,m9),m9=1,nae)
      endif
      tqlty=tqlty+str_quality_1(m8)
      bqlty=bqlty+bondq(m8)
      sqlty=sqlty+str_quality_2(m8)
   enddo
   write(9,914)'qualities:','intra_bond','sym_break','nn_bond'
   write(9,915)tqlty,sqlty,bqlty

900 format(I3,x,I3,x,I3,x,a,3x,25I4)
901 format(I3,x,I3,x,I3,x,a,3x,I1,a,I2,x,25I4)
914 format(a,x,a,x,a,x,a)
915 format(15x,I3,7x,I3,7x,I3)

return
endif
!**********************************************************************************************************
!!!!!!!!!!!!!!!! If number of permisible structures and available structures are same END !!!!!!!!!!!!!!!
!**********************************************************************************************************

allocate(str1(MaxStrOepo, nae))
allocate(str2(MaxStrOepo, nae))
!allocate(symsc(MaxStrOepo))
cnt=0
do i5=1,10000
str_cnt(i5)=0
enddo

loop1:do i1=1,m4
   str1 = 0
   str2 = 0

   if(repl.ne.0)then
      do i4=1,repl
         rep=0
         do i3=1,nae
            if(repl_struc(i4,i3).eq.astr(i1,i3))then
               rep=rep+1
            endif
         enddo 
         if(rep.eq.nae) cycle loop1
      enddo
   endif

   str1_cnt=0

   do i5=1,cnt
      if(str_cnt(i5).eq.i1) cycle loop1
   enddo

   cnt=cnt+1
   str1_cnt=str1_cnt+1

   do i3=1,nae
      str1(str1_cnt,i3)=astr(i1,i3)
   enddo 

   symqq(str1_cnt)=symq(i1)
   qulsymm(str1_cnt)=qulsym(i1)

   lp_cnt=0

   do i3=1,nl*2,2
      lp_cnt=lp_cnt+1
      lps(lp_cnt)=astr(i1,i3)
   enddo

   str_cnt(cnt)=i1

   loop2:do i2=i1,m4
      if (repl.ne.0)then
         do i4=1,repl
            rep=0
            do i3=1,nae
               if (repl_struc(i4,i3).eq.astr(i2,i3))then
                  rep=rep+1
               endif
            enddo 
            if (rep.eq.nae) cycle loop2
         enddo
      endif

      do i5=1,cnt
         if (str_cnt(i5).eq.i2) cycle loop2
      enddo

      if (nl.ne.0)then
         ii=0
         do i6=1,nl*2,2
            do i7=1,nl*2,2
               if (astr(i1,i6).eq.astr(i2,i7))then
                  ii=ii+1
               endif
            enddo
         enddo
         if (ii.ne.nl) cycle loop2
      endif


      loop3:do i3=nl*2+1,nae
         do i4=nl*2+1,nae
            if(astr(i1,i3).eq.astr(i2,i4)) cycle loop3
         enddo
         cycle loop2
      enddo loop3

      cnt=cnt+1
      str1_cnt=str1_cnt+1

      do i3=1,nae
         str1(str1_cnt,i3)=astr(i2,i3)
      enddo

      symqq(str1_cnt)=symq(i2)
      qulsymm(str1_cnt)=qulsym(i2)
      str_cnt(cnt)=i2
      elporb=nae-nl*2
      total_str=factorial(elporb-1-mod(elporb,2))
      call wigner(elporb, allowed_nstr)
      if (flg_ion.eq.1)then
         if (total_str.eq.cnt) exit
      endif
      if (flg_cov.eq.1)then
         if (total_str*elporb.eq.cnt) exit
      endif
   enddo loop2

!***********************************************************************
!!!!!!!!!!!!!!!! Rumer Structures selection  Start !!!!!!!!!!!!!!!!!!!!!
!***********************************************************************
   if (nl.ne.0)write(7,*)'                 lone pair =',(lps(i3),i3=1,lp_cnt)
   if (flg1.eq.1) then
     if (nfset.eq.0) then
        call quality_factor(nl,str1,str1_cnt,quality_fac,str_quality_1,str_quality_2,bondq, mbondq, pref_radical)
!        call qult_str_arrange(nl,str1,str1_cnt,str2)
        call rumer_structures(nl,str1,str1_cnt)
        call write_structures_file(str1_cnt,str1)
        if(symm.eq.1)then
          if(sig_sym_flg.eq.1)call symmetry_cal_sig(nl, str1, str1_cnt, symq, nssym)
          if(sig_sym_flg.ne.1)call symmetry_cal_pi(nl, str1, str1_cnt, symq, nssym)
        endif
        call write_rumer_xmi(nl,str1,str1_cnt)
     endif

     if (nfset.eq.1.or.nfset.eq.2) then
        call quality_factor(nl,str1,str1_cnt,quality_fac,str_quality_1,str_quality_2,bondq, mbondq, pref_radical)
!        call qult_str_arrange(nl,str1,str1_cnt,str2)
        call rumer_structures(nl,str1,str1_cnt)
        call write_structures_file(str1_cnt,str1)
        if(symm.eq.1)then
          if(sig_sym_flg.eq.1)call symmetry_cal_sig(nl, str1, str1_cnt, symq, nssym)
          if(sig_sym_flg.ne.1)call symmetry_cal_pi(nl, str1, str1_cnt, symq, nssym)
        endif
        !print*,'symq',sig_sym_flg,'|',(symq(i7),i7=1,str1_cnt)
        !stop

        call All_Rumer_set(str1_cnt,str1,nl)
     endif
   endif
!***********************************************************************
!!!!!!!!!!!!!!!! Rumer Structures selection Ends !!!!!!!!!!!!!!!!!!!!!
!***********************************************************************
!***********************************************************************

!***********************************************************************
!!!!!!!!!!!!!!!! Covalent Structures selection Starts !!!!!!!!!!!!!!!!!!
!***********************************************************************

   if (flg1.eq.0.and.flg_cov.eq.1.and.flg_ion.eq.0)then
      if (symm.eq.1) then
         call quality_factor(nl,str1,str1_cnt,quality_fac,str_quality_1,str_quality_2,bondq, mbondq, pref_radical)
         if(sig_sym_flg.eq.1)call symmetry_cal_sig(nl,str1,str1_cnt,symq,nssym)
         if(sig_sym_flg.ne.1)call symmetry_cal_pi(nl,str1,str1_cnt,symq,nssym)
         !print*,'symq',sig_sym_flg,'|',(symq(i7),i7=1,str1_cnt)
         !stop
         call qult_str_arrange(nl,str1,str1_cnt,str2)
         call vector_rep(nl,str2,str1_cnt,fvec)
         call rumer_structures(nl,str2,str1_cnt)
         call write_structures_file(str1_cnt,str2)
         call write_symm_xmi_new(nl,allowed_nstr,str2,str1_cnt)
      endif

      if (symm.eq.0.and.nfset.ne.4)then
         call quality_factor(nl,str1,str1_cnt,quality_fac,str_quality_1,str_quality_2,bondq, mbondq, pref_radical)
         call qult_str_arrange(nl,str1,str1_cnt,str2)
         call rumer_structures(nl,str2,str1_cnt)
         call write_structures_file(str1_cnt, str2)
         call vector_rep(nl,str2,str1_cnt,fvec)
         call main_bond_cal(nl,str2,str1_cnt,mbondq)
      endif

      if (symm.eq.0.and.nfset.eq.4)then
         call quality_factor(nl,str1,str1_cnt,quality_fac,str_quality_1,str_quality_2,bondq, mbondq, pref_radical)
         call qult_str_arrange(nl,str1,str1_cnt,str2)
         call rumer_structures(nl,str2,str1_cnt)
         call write_structures_file(str1_cnt, str2)
         call main_bond_cal(nl,str2,str1_cnt,mbondq)
         call check_str_bond(nl,str2,str1_cnt)
         call eq_dstr_set(str1_cnt,nl,allowed_nstr,str2)
      endif
   endif

!***********************************************************************
!!!!!!!!!!!!!!!! Covalent Structures selection Ends !!!!!!!!!!!!!!!!!!
!***********************************************************************

   tqlty1=tqlty1+tqlty
   bqlty1=bqlty1+bqlty
   sqlty1=sqlty1+sqlty
enddo loop1

if (nfset.eq.1.or.flg1.eq.1) then
   write(9,910)'qualities:',' intra_bond','sym_break','nn_bond'
   write(9,911)tqlty1,sqlty1,bqlty1
   write(9,*)'    '
endif

910 format(a,x,a,x,a,x,a)
911 format(15x,I3,7x,I5,7x,I3)

deallocate(str1)
deallocate(str2)

print*,'exit str_selection'
return
end subroutine str_selection

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_structures_file(str1_cnt, str2)!, mbondq, pref_radical)
integer:: m8, m9, str1_cnt
!integer::mbondq(15000), pref_radical(15000)
character(len=5)::rumstr
integer, pointer:: str2(:,:)!, q_fac(:)!, mbondq(:), pref_radical(:)

print*,'enter write_structures_file'
if (nfset.eq.1.or.nfset.eq.2) then
  strplp = str1_cnt
  if(.not. allocated(f_str))then
    allocate(f_str(strplp, nae))
  f_str = 0
  f_str = str2
  endif
endif


do m8=1,str1_cnt
   if(rumer(m8)*rumer_rad(m8).eq.1)rumstr='Rumer'
   if(rumer(m8)*rumer_rad(m8).eq.0)rumstr='-'
   if(niao.ne.0)write(7,301)'structure',m8,'[',str_quality_1(m8),bondq(m8),str_quality_2(m8),&
     mbondq(m8),pref_radical(m8),']','{',quality_fac(m8),'}',rumstr,'1:',niao,(str2(m8,m9),m9=1,nae)
   if(niao.eq.0)write(7,302)'structure',m8,'[',str_quality_1(m8),bondq(m8),str_quality_2(m8),&
     mbondq(m8),pref_radical(m8),']','{',quality_fac(m8),'}',rumstr,(str2(m8,m9),m9=1,nae)
enddo

write(7,*)'----------------------------  END FILE  ----------------------------'
!close(7)
301 format(a,2x,I0,3x,a,2x,I0,2x,I0,2x,I0,2x,I0,2x,I0,2x,a,3x,a,2x,I0,2x,a,14x,a,14x,a,I0,1x,*(I0, 1x))
302 format(a,2x,I0,3x,I0,2x,I0,2x,I0,2x,I0,2x,I0,3x,a,I0,a,14x,a,14x,*(I0, 1x))
print*,'exit write_structures_file'

return
end subroutine write_structures_file

end module mod_str_selection
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
