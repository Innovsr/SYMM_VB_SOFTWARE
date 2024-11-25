!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine str_selection(astr,nl,m4,perm_nstr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

common/quality/str_quality_1,str_quality_2,bondq,tqlty,bqlty,sqlty,tnqs,nssym,qulsym,symq,&
sigsym,tnqs_sig

integer::i,j,i1,i2,i3,i4,i5,i6,i7,ii,nl,str1_cnt,cnt,m4,m8,m9,m10,lp_cnt,tnqs,&
wig2,elporb,perm_nstr,total_str,mbondq(15000),q_fac(15000),q_fac1(15000),ijk,n,nssym,rep
integer::str_cnt(15000),astr(15000,20),str1(15000,20),str2(15000,20),lps(50),sigsym(15000),tnqs_sig&
,quality_fac(15000),str_quality_1(15000),str_quality_2(15000),strno(1000)&
,rumer(15000),rumer_rad(15000),bondq(15000),fvec(15000,1000),qulsymm(15000),&
symq(15000),symqq(15000),&
tndet,tqlty,tqlty1,bqlty,sqlty,bqlty1,sqlty1,qulsym(15000)
real*8::factorial,symsc(15000)
character(5)::rumstr


print*,'enter str_selection'
print*,'ovopt',ovopt,nfset
tqlty1=0
bqlty1=0
sqlty1=0
!do i1=1,m4
!write(*,305),i1,(astr(i1,i2),i2=1,nae)
!enddo
!print*,perm_nstr,nl,m4

!**********************************************************************************************************
!!!!!!!!!!!!!!!! If number of permisible structures and available structures are same START !!!!!!!!!!!!!!!
!**********************************************************************************************************

if(perm_nstr.eq.m4)then
call rumer_structures(nl,astr,m4,rumer,rumer_rad)
call quality_factor(nl,astr,m4,quality_fac,str_quality_1,str_quality_2,bondq)


write(7,*)'all',m4,' structures are permisible among below '
!n=0
!do ijk=nssym,1,-1
do m8=1,m4
!if(symq(m8).eq.ijk)then
!n=n+1
if(rumer(m8)*rumer_rad(m8).eq.1)rumstr='R'
if(rumer(m8)*rumer_rad(m8).eq.0)rumstr='-'
write(7,301)'str',m8,')','[',str_quality_1(m8),',',str_quality_2(m8),']','{',quality_fac(m8),'}',&
rumstr,(astr(m8,m9),m9=1,nae)
!write(7,305)(fvec(m8,m9),m9=1,ndet)
!write(*,305)(str1(m8,m9),m9=1,nae)
!endif
enddo
!enddo

301 format(a,x,I5,a,x,a,I3,a,I3,x,a,a,I3,x,a,x,a,x,30I5)
305 format(30I5)

!if(nnnatom.ne.0)then
if(input_flg.eq.1)call nnat_bond_cal(nl,astr,m4,bondq)
if(input_flg.eq.0)call nnat_bond_cal_2(nl,astr,m4,bondq)
!endif
tqlty=0
bqlty=0
sqlty=0
do m8=1,m4
if(niao.eq.0)then
write(9,900)str_quality_1(m8),bondq(m8),str_quality_2(m8),'|',(astr(m8,m9),m9=1,nae)
endif
!if(niao.gt.1.and.nnnatom.ne.0)then
if(niao.gt.1)then
write(9,901)str_quality_1(m8),bondq(m8),str_quality_2(m8),'|',1,':',niao,(astr(m8,m9),m9=1,nae)
endif
!if(niao.eq.1.and.nnnatom.ne.0)then
if(niao.eq.1)then
write(9,901)str_quality_1(m8),bondq(m8),str_quality_2(m8),'|',1,1,(astr(m8,m9),m9=1,nae)
endif
!if(niao.ne.0.and.nnnatom.eq.0)then
!write(9,909),quality_fac(m8),1,':',niao,(astr(m8,m9),m9=1,nae)
!endif
tqlty=tqlty+str_quality_1(m8)
bqlty=bqlty+bondq(m8)
sqlty=sqlty+str_quality_2(m8)
enddo
write(9,914)'qualities:','intra_bond','sym_break','nn_bond'
write(9,915)tqlty,sqlty,bqlty

900 format(I3,x,I3,x,I3,x,a,3x,25I4)
901 format(I3,x,I3,x,I3,x,a,3x,I1,a,I2,x,25I4)
909 format(I3,x,I1,a,I1,x,25I4)
914 format(a,x,a,x,a,x,a)
915 format(15x,I3,7x,I3,7x,I3)
goto 374
endif
!**********************************************************************************************************
!!!!!!!!!!!!!!!! If number of permisible structures and available structures are same END !!!!!!!!!!!!!!!
!**********************************************************************************************************

cnt=0
do i5=1,10000
str_cnt(i5)=0
enddo

do i1=1,m4
do m8=1,10000
do m9=1,nae
str1(m8,m9)=0
enddo
enddo

if(repl.ne.0)then
do i4=1,repl
rep=0
do i3=1,nae
if(repl_struc(i4,i3).eq.astr(i1,i3))then
rep=rep+1
endif
enddo 
if(rep.eq.nae)goto 373
enddo
endif


str1_cnt=0

do i5=1,cnt
if(str_cnt(i5).eq.i1)goto 373
enddo

cnt=cnt+1
str1_cnt=str1_cnt+1

do i3=1,nae
str1(str1_cnt,i3)=astr(i1,i3)
enddo 

!print*,'sourav r',str1_cnt,i1
symqq(str1_cnt)=symq(i1)
qulsymm(str1_cnt)=qulsym(i1)
!print*,'symqq(str1_cnt),qulsymm(str1_cnt)',symqq(str1_cnt),qulsymm(str1_cnt)
!print*,symqq(str1_cnt)

lp_cnt=0

do i3=1,nl*2,2
lp_cnt=lp_cnt+1
lps(lp_cnt)=astr(i1,i3)
enddo

!write(*,305),str1_cnt,str1_cnt,(str1(str1_cnt,m9),m9=1,nae)
str_cnt(cnt)=i1

do i2=i1,m4

if(repl.ne.0)then
do i4=1,repl
rep=0
!print*,(repl_struc(i4,i3),i3=1,nae)
do i3=1,nae
if(repl_struc(i4,i3).eq.astr(i2,i3))then
rep=rep+1
endif
enddo 
if(rep.eq.nae)goto 372
enddo
endif

do i5=1,cnt
if(str_cnt(i5).eq.i2)goto 372
enddo
!write(*,305),i2,(astr(i2,m9),m9=1,nae)


if(nl.ne.0)then
ii=0
do i6=1,nl*2,2
do i7=1,nl*2,2
if(astr(i1,i6).eq.astr(i2,i7))then
ii=ii+1
endif
enddo
enddo
if(ii.ne.nl)goto 372
endif


do i3=nl*2+1,nae
do i4=nl*2+1,nae
if(astr(i1,i3).eq.astr(i2,i4))goto 371
enddo
goto 372
371 enddo

cnt=cnt+1
str1_cnt=str1_cnt+1

do i3=1,nae
str1(str1_cnt,i3)=astr(i2,i3)
enddo

symqq(str1_cnt)=symq(i2)
qulsymm(str1_cnt)=qulsym(i2)
!print*,'symqq(str1_cnt),qulsymm(str1_cnt)',symqq(str1_cnt),qulsymm(str1_cnt)

str_cnt(cnt)=i2

elporb=nae-nl*2
total_str=factorial(elporb-1-mod(elporb,2))
!print*,'total_str',total_str,elporb,nl,nae
!print*,'cnt',nl,total_str,elporb
if(flg_ion.eq.1)then
if(total_str.eq.cnt)goto 375
endif
if(flg_cov.eq.1)then
if(total_str*elporb.eq.cnt)goto 375
endif
372 enddo 

!elporb=nae-nl*2
375 call wigner(elporb,wig2)
write(7,*)wig2,' structures are permisible among below ',str1_cnt
if(nl.ne.0)write(7,*)'                 lone pair =',(lps(i3),i3=1,lp_cnt)

call rumer_structures(nl,str1,str1_cnt,rumer,rumer_rad)
call quality_factor(nl,str1,str1_cnt,quality_fac,str_quality_1,str_quality_2,bondq)


!n=0
!do ijk=nssym,1,-1
do m8=1,str1_cnt
!if(symq(m8).eq.ijk)then
!n=n+1
!write(7,*),rumer(m8),rumer_rad(m8),rumer(m8)*rumer_rad(m8)
if(rumer(m8)*rumer_rad(m8).eq.1)rumstr='R'
if(rumer(m8)*rumer_rad(m8).eq.0)rumstr='-'
write(7,301)'str',m8,')','[',str_quality_1(m8),',',str_quality_2(m8),']','{',quality_fac(m8),'}',&
rumstr,(str1(m8,m9),m9=1,nae)
!write(7,305)(fvec(m8,m9),m9=1,ndet)
!write(*,305)(str1(m8,m9),m9=1,nae)
!endif
enddo
!write(7,*)
!enddo


!***********************************************************************
!!!!!!!!!!!!!!!! Rumer Structures selection  Start !!!!!!!!!!!!!!!!!!!!!
!***********************************************************************

if(flg1.eq.1)then

call quality_factor(nl,str1,str1_cnt,quality_fac,str_quality_1,str_quality_2,bondq)

!do i=1,str1_cnt
!write(*,231),(str1(i,j),j=1,nae),str_quality_1(i),str_quality_2(i),bondq(i),quality_fac(i)
!enddo


!do i=1,str1_cnt
!write(*,231),(str2(i,j),j=1,nae),rumer(i),rumer_rad(i)
!enddo

call qult_str_arrange(nl,str1,str1_cnt,quality_fac,str2,q_fac)
call rumer_structures(nl,str2,str1_cnt,rumer,rumer_rad)
call write_rumer_xmi(nl,str2,str1_cnt,rumer,rumer_rad,quality_fac)
endif
print*,'sourav1'
!***********************************************************************
!!!!!!!!!!!!!!!! Rumer Structures selection Ends !!!!!!!!!!!!!!!!!!!!!
!***********************************************************************
!***********************************************************************

!***********************************************************************
!!!!!!!!!!!!!!!! Covalent Structures selection Starts !!!!!!!!!!!!!!!!!!
!***********************************************************************

!print*,' flg_ion,flg_cov',flg_ion,flg_cov
if(flg1.eq.0.and.flg_cov.eq.1.and.flg_ion.eq.0)then
!print*,'sourav1'

if(symm.eq.1) then
call quality_factor(nl,str1,str1_cnt,quality_fac,str_quality_1,str_quality_2,bondq)


if(sig_sym_flg.eq.1)call symmetry_cal_sig(nl,str1,str1_cnt,symsc,symq,nssym)
if(sig_sym_flg.ne.1)call symmetry_cal_pi(nl,str1,str1_cnt,symsc,symq,nssym)
call qult_str_arrange(nl,str1,str1_cnt,quality_fac,str2,q_fac)
do i=1,str1_cnt
write(*,231)(str2(i,j),j=1,nae)
enddo

call vector_rep(nl,str2,str1_cnt,fvec)
call write_symm_xmi_new(nl,wig2,str2,str1_cnt,q_fac)

endif

if(symm.eq.0)then
if(nfset.ne.4)then
call quality_factor(nl,str1,str1_cnt,quality_fac,str_quality_1,str_quality_2,bondq)

call qult_str_arrange(nl,str1,str1_cnt,quality_fac,str2,q_fac)

call vector_rep(nl,str2,str1_cnt,fvec)
print*,'ovopt',ovopt,niao
call main_bond_cal(nl,str2,str1_cnt,mbondq)
print*,'niao:str_sel',niao

endif

if(nfset.eq.4)then
print*,'str1_cnt',str1_cnt
!stop
call quality_factor(nl,str1,str1_cnt,quality_fac,str_quality_1,str_quality_2,bondq)

call qult_str_arrange(nl,str1,str1_cnt,quality_fac,str2,q_fac)
print*,'ovopt',ovopt,niao
call main_bond_cal(nl,str2,str1_cnt,mbondq)
call check_str_bond(nl,wig2,str2,str1_cnt)
print*,'niao:str_sel',niao
call eq_dstr_set(str1_cnt,nl,wig2,str2)
endif
endif

endif

!***********************************************************************
!!!!!!!!!!!!!!!! Covalent Structures selection Ends !!!!!!!!!!!!!!!!!!
!***********************************************************************
tqlty1=tqlty1+tqlty
bqlty1=bqlty1+bqlty
sqlty1=sqlty1+sqlty
373 enddo
231 format(30I3)

if(nfset.eq.1.or.flg1.eq.1) then
write(9,910)'qualities:',' intra_bond','sym_break','nn_bond'
write(9,911)tqlty1,sqlty1,bqlty1
write(9,*)'    '
endif

910 format(a,x,a,x,a,x,a)
911 format(15x,I3,7x,I5,7x,I3)

print*,'exit str_selection'
374 return
end subroutine str_selection

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
