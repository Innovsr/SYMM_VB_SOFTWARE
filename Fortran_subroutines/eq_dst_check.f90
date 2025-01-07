!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine eq_dst_check(nts,strsln,fail,nl,str2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none
integer::totbnd,totseorb,str_bnd(15000,1000),str_rad(500,20),bond_sl(10000,2),num_orb(20)
integer::i,i1,j,m,m1,l,add_edst(500),n1,n2,n3,quality_fac(15000),str_quality_1(15000),&
str_quality_2(15000),nssym,qulsym(15000),symq(15000),tqlty,bqlty,sqlty,tnqs,sigsym(15000)&
,tnqs_sig,fail,strsln(1000),nl,tnbd,nts,nbd,str2(15000,20),bondq(15000),strflg,num_str,numpstr,det_inv
integer::str3(15000,20),fvec(15000,1000),strsln1(1000)
real*8::factorial
common /chek/totbnd,totseorb,str_bnd,str_rad,bond_sl,num_orb
common/quality/str_quality_1,str_quality_2,bondq,tqlty,bqlty,sqlty,tnqs,nssym,qulsym,symq,&
sigsym,tnqs_sig
common/nst/num_str,numpstr


fail=0
do l=1,totbnd
   add_edst(l)=0
   do i=1,nts
      add_edst(l)=add_edst(l)+str_bnd(strsln(i),l)
   enddo
enddo

nbd=((nao-nl-nlast)/2)*numpstr
n1=nao-nl
n2=nao-nl-2
tnbd=(int(factorial(n1)/(factorial(n2)))/2)

if(nbd.gt.tnbd)n3=nbd/tnbd
if(nbd.le.tnbd)n3=tnbd/nbd

do i=1,totbnd
   if(add_edst(i).gt.n3)fail=1
enddo

if(nts.eq.numpstr.and.fail.eq.0)then
   do i=1,numpstr
      strsln1(i)=i
      do i1=1,nae
         str3(i,i1)=str2(strsln(i),i1)
      enddo
   enddo


   call vector_rep(nl,str3,numpstr,fvec)
   
   if(tnbd.ne.nbd)call mat_ind(nl,nts,numpstr,strsln1,fail,det_inv)
   
   if(fail.eq.0)then
      num_str=num_str+1
      do m=1,numpstr
         if(niao.eq.0)then
            write(9,900)str_quality_1(strsln(m)),bondq(strsln(m)),str_quality_2(strsln(m)),'|',(str2(strsln(m),m1),m1=1,nae)
         endif
         if(niao.gt.1)then
            write(9,901)str_quality_1(strsln(m)),bondq(strsln(m)),str_quality_2(strsln(m)),'|',1,':',niao,&
                    (str2(strsln(m),m1),m1=1,nae)
         endif
         if(niao.eq.1)then
            write(9,909)str_quality_1(strsln(m)),bondq(strsln(m)),str_quality_2(strsln(m)),'|',1,1,(str2(strsln(m),m1),m1=1,nae)
         endif
      enddo
   
   fail=1

   write(9,*)''
   write(9,*)'set number',num_str
   
   write(9,*)''
   write(9,*)''

   endif
endif

if(num_str.eq.1000) then
   stop
endif

900 format(I3,2x,I3,2x,I3,2x,a,x,25I4)
901 format(I3,2x,I3,2x,I3,2x,a,x,I2,a,I2,x,25I4)
909 format(I3,2x,I3,2x,I3,2x,a,x,I3,I3,x,25I4)
231 format(a,2x,50I3)
102 format (50I5)

print*,'exit eq_dstr_set'
return
end subroutine eq_dst_check
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
