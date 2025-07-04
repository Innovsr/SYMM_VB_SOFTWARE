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
use str_module
implicit none

contains
subroutine str_selection(astr,nl,m4,perm_nstr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer::i,i1,i2,i3,i4,i5,i6,i7,ii,nl,str1_cnt,cnt,m4,m8,m9,lp_cnt
integer::elporb,perm_nstr,total_str,rep,allowed_nstr, sf1, sf2, flg2
integer::str_cnt(15000), lps(50)
integer::qulsymm(15000),symqq(15000),tqlty1,bqlty1,sqlty1, factorial
integer, pointer::astr(:,:), str1(:,:), str2(:,:)
integer, pointer:: fvec(:,:)!,q_fac(:)
!integer, allocatable::str12(:,:), col9(:)


tqlty1=0
bqlty1=0
sqlty1=0

!**********************************************************************************************************
!!!!!!!!!!!!!!!! If number of permisible structures and available structures are same then  !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!   they are not send to the independency check. They are directly printed   !!!!!!!!!!!!!!!
!**********************************************************************************************************

print*,'perm_nstr, m4',perm_nstr, m4
flg2=1 ! this flag instruct 'space num will be written in write_rumer_xmi subroutine'
if(perm_nstr.eq.m4)then
   flg2=0 ! this flag instruct 'space num will not be written in write_rumer_xmi subroutine'
   if(flg_cov.eq.1)then
     cov_space=cov_space+perm_nstr
     write(5,913)'Cov_Space',cov_space-perm_nstr+1,'-',cov_space+1
     write(9+u1,914)'============= Cov Space num:',cov_space-perm_nstr+1,'-',cov_space+1,'============='
   endif
   if(flg_ion.eq.1)then
     ion_space=ion_space+perm_nstr
     write(5,913)'Ion_Space',ion_space-perm_nstr+1,'-',ion_space+1
     write(9+u1,914)'============= Ion Space num:',ion_space-perm_nstr+1,'-',ion_space+1,'============='
   endif

   call rumer_structures(nl,astr,m4)
   call quality_factor(nl,astr,m4)
   call write_structures_file(m4, astr)

   if(input_flg.eq.1)call nnat_bond_cal(nl,astr,m4,bondq)
   if(input_flg.eq.0)call nnat_bond_cal_2(nl,astr,m4,bondq)
   tqlty=0
   bqlty=0
   sqlty=0
   do m8 = 1, m4
     do i=1,nae
       str12(m8,i)=astr(m8,i)
     enddo
     col9(m8) = m8
     qq10(m8)=quality_fac(m8)
     qq11(m8)=str_quality_1(m8)
     qq12(m8)=str_quality_2(m8)
     bondq14(m8)=bondq(m8)
   enddo
   set_number = 0
   if (flg1.eq.0) call write_output(sf1, sf2, perm_nstr, str12, col9)
   if (flg1.eq.1) call write_rumer_xmi(nl, astr, perm_nstr,flg2)

deallocate(rumer)
deallocate(rumer_rad)

913 format(a,2x,I0,2x,a,2x,I0)
914 format(a,2x,I0,2x,a,2x,I0,2x,a)

return
endif


!**********************************************************************************************************
!!!!!!!!!!!!!!!! If number of permisible structures and available structures are not same   !!!!!!!!!!!!!!!
!**********************************************************************************************************


allocate(str1(MaxStrOepo, nae))
allocate(str2(MaxStrOepo, nae))
!allocate(symsc(MaxStrOepo))
cnt=0
do i5=1,10000
str_cnt(i5)=0
enddo

loop1:do i1=1,m4
   if(.not. allocated(symq))then
     allocate(symq(m4))
     symq = 0
   endif
   print*,'loop1',i1
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

   print*,'flg1',flg1
   if (flg1.eq.1) then
     if (nfset.eq.0) then
        print*,'Rumer'
        call quality_factor(nl,str1,str1_cnt)
!        call qult_str_arrange(nl,str1,str1_cnt,str2)
        call rumer_structures(nl,str1,str1_cnt)
        call write_structures_file(str1_cnt,str1)
        if(symm.eq.1)then
          if(sig_sym_flg.eq.1)call symmetry_cal_sig(nl, str1, str1_cnt, symq, nssym)
          if(sig_sym_flg.ne.1)call symmetry_cal_pi(nl, str1, str1_cnt, symq, nssym)
        endif
        call write_rumer_xmi(nl,str1,str1_cnt,flg2)
        deallocate(rumer)
        deallocate(rumer_rad)
        deallocate(symq)
     endif

     if (nfset.eq.1.or.nfset.eq.2) then
        print*,'Rumer2'
        call quality_factor(nl,str1,str1_cnt)!,quality_fac,str_quality_1,str_quality_2,bondq, mbondq, pref_radical)
!        call qult_str_arrange(nl,str1,str1_cnt,str2)
        call rumer_structures(nl,str1,str1_cnt)
        call write_structures_file(str1_cnt,str1)
        if(symm.eq.1)then
          print*,'Symm_str_select',symm
          if(sig_sym_flg.eq.1)call symmetry_cal_sig(nl, str1, str1_cnt, symq, nssym)
          if(sig_sym_flg.ne.1)call symmetry_cal_pi(nl, str1, str1_cnt, symq, nssym)
        endif
        !print*,'symq',sig_sym_flg,'|',(symq(i7),i7=1,str1_cnt)
        !stop

        call All_Rumer_set(str1_cnt,str1,nl)
        deallocate(rumer)
        deallocate(rumer_rad)
        deallocate(symq)
     endif
   endif
!***********************************************************************
!!!!!!!!!!!!!!!! Rumer Structures selection Ends !!!!!!!!!!!!!!!!!!!!!
!***********************************************************************
!***********************************************************************

!***********************************************************************
!!!!!!!!!!!!!!!! Chem-Inst Structures selection Starts !!!!!!!!!!!!!!!!!!
!***********************************************************************

   if (flg1.eq.0)then
      if (symm.eq.1.and.asymm.eq.0) then
         print*,'chem-Inst symm'
         call quality_factor(nl,str1,str1_cnt)
         if(sig_sym_flg.eq.1)call symmetry_cal_sig(nl,str1,str1_cnt,symq,nssym)
         if(sig_sym_flg.ne.1)call symmetry_cal_pi(nl,str1,str1_cnt,symq,nssym)
         !print*,'symq',sig_sym_flg,'|',(symq(i7),i7=1,str1_cnt)
         !stop
         print*,'qult_str_arr in'
         call qult_str_arrange(nl,str1,str1_cnt,str2)
         print*,'qult_str_arr out'
         call vector_rep(nl,str2,str1_cnt,fvec)
         print*,'vector_rep_out'
         call rumer_structures(nl,str2,str1_cnt)
         print*,'rumer_structures out'
         call write_structures_file(str1_cnt,str2)
         print*,'write_structures_file out'
         call write_symm_xmi_new(nl,allowed_nstr,str2,str1_cnt)
         print*,'write_symm_xmi_new out'
         deallocate(rumer)
         deallocate(rumer_rad)
         deallocate(symq)
         print*,'i am here'
      endif

      if (symm.eq.0.and.nfset.ne.4.and.asymm.eq.0)then
         print*,'chem-Inst'
         call quality_factor(nl,str1,str1_cnt)!,quality_fac,str_quality_1,str_quality_2,bondq, mbondq, pref_radical)
         call qult_str_arrange(nl,str1,str1_cnt,str2)
         call rumer_structures(nl,str2,str1_cnt)
         call write_structures_file(str1_cnt, str2)
         call vector_rep(nl,str2,str1_cnt,fvec)
         call main_bond_cal(nl,str2,str1_cnt,mbondq)
         deallocate(rumer)
         deallocate(rumer_rad)
         deallocate(symq)
      endif

      if (symm.eq.0.and.nfset.eq.4.and.asymm.eq.0)then
         call quality_factor(nl,str1,str1_cnt)!,quality_fac,str_quality_1,str_quality_2,bondq, mbondq, pref_radical)
         call qult_str_arrange(nl,str1,str1_cnt,str2)
         call rumer_structures(nl,str2,str1_cnt)
         call write_structures_file(str1_cnt, str2)
         call main_bond_cal(nl,str2,str1_cnt,mbondq)
         call check_str_bond(nl,str2,str1_cnt)
         call eq_dstr_set(str1_cnt,nl,allowed_nstr,str2)
         deallocate(rumer)
         deallocate(rumer_rad)
         deallocate(symq)
      endif

      if (symm.eq.0.and.asymm.eq.1) then
         print*,'sourav inside the cheminst assymmetric str generation'
         call quality_factor(nl,str1,str1_cnt)!,quality_fac,str_quality_1,str_quality_2,bondq, mbondq, pref_radical)
         call qult_str_arrange(nl,str1,str1_cnt,str2)
         print*,'qult_str_arr out'
         call vector_rep(nl,str2,str1_cnt,fvec)
         call rumer_structures(nl,str2,str1_cnt)
         call write_structures_file(str1_cnt,str2)
         call write_symm_xmi_new(nl,allowed_nstr,str2,str1_cnt)
         deallocate(rumer)
         deallocate(rumer_rad)
         deallocate(symq)
      endif
   endif

!***********************************************************************
!!!!!!!!!!!!!!!! Covalent Structures selection Ends !!!!!!!!!!!!!!!!!!
!***********************************************************************

print*,'sourav1'
   tqlty1=tqlty1+tqlty
   bqlty1=bqlty1+bqlty
   sqlty1=sqlty1+sqlty
enddo loop1

print*,'sourav2'
if (nfset.eq.1.or.flg1.eq.1) then
   write(9,910)'qualities:',' intra_bond','sym_break','nn_bond'
   write(9,911)tqlty1,sqlty1,bqlty1
   write(9,*)'    '
endif
print*,'sourav3'

910 format(a,x,a,x,a,x,a)
911 format(15x,I3,7x,I5,7x,I3)

deallocate(str1)
deallocate(str2)

print*,'exit str_selection'
return
end subroutine str_selection

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_structures_file(str1_cnt, str2)
integer:: m8, m9, str1_cnt
character(len=5)::rumstr
integer, pointer:: str2(:,:)

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
