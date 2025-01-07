!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine qult_str_arrange(nl,str2,n,q_fac,str1,q_fac1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

common/quality/str_quality_1,str_quality_2,bondq,tqlty,bqlty,sqlty,tnqs,nssym,qulsym,symq,&
sigsym,tnqs_sig

integer::i1,i5,i6,i7,l,l1,k,m,bnd,nsymstr,ll,nstsymset,sym_str_num_final(1000),sym_str_qual(1000)
integer::sym_str_sl(1000,100),sym_str_sl1(1000),sym_str_sl_final(1000,100),p,nstrt1,&
sym_str_num(1000),sym_str_num1(1000),val,mq_fac,msymq,mloopsymsc,ind,sl,sll
integer::i,i2,i3,i4,i9,m16,m18,m19,q_fac(15000),q_fac1(15000),k6,k7,k8,k9,st_ct(15000),st_ct1(15000),ii,kk
integer::str1(15000,20),str2(15000,20),n,jj,j,j1,nqul1,nl,iii,x,y,z,iij,i8,q_cnt(1000)&
,str_quality_1(15000),str_quality_2(15000),bondq(15000),sigsym(15000),tnqs_sig&
,tqlty,bqlty,sqlty,tnqs,nssym,qulsym(15000),symq(15000),symq_final(15000),ij,ijk,qt,qtg,q&
,mat(100,100),cnt(15000),qt_max,q_max,lll,sym_str_qual1(1000,100),str_cnt1(1000000),nssym1,mm
integer,dimension(:,:),allocatable::str3
integer,dimension(:),allocatable::qual4,qual3,str_cnt,q_fac3,q_fac4,q_fac5,q_fac6,q_fac7,q_fac2,bdq,bdq1,qul1
integer,dimension(:),allocatable::stq1,stq2,stq3,stq4,sl_num,sl_qfac,sl_qfac1

!set_order desides if the symmetric sets arranged according to quality or size
! set_order=0 means quality
! set_order=1 means  smaller to larger sets
! set_order=2 means larger to smaller sets
print*,'set_order',set_order
!stop
print*,'enter qult_str_arrange'
x=50000
y=20
z=1000

allocate(q_fac2(x))
allocate(q_fac3(x))
allocate(q_fac4(x))
allocate(q_fac5(x))
allocate(q_fac6(x))
allocate(q_fac7(x))
allocate(qul1(x))
allocate(qual3(x))
allocate(qual4(x))
allocate(bdq(x))
allocate(bdq1(x))
allocate(str_cnt(z))
allocate(stq1(x))
allocate(stq2(x))
allocate(stq3(x))
allocate(stq4(x))
allocate(sl_num(x))
allocate(sl_qfac(x))
allocate(sl_qfac1(x))

q_fac3(x)=0
q_fac4(x)=0
q_fac5(x)=0
q_fac6(x)=0
q_fac7(x)=0
q_fac2(x)=0
qul1(x)=0
qual3(x)=0
qual4(x)=0
bdq(x)=0
bdq1(x)=0
str_cnt(z)=0
stq1(x)=0
stq2(x)=0
stq3(x)=0
stq4(x)=0
sl_num(x)=0
sl_qfac(x)=0
sl_qfac1(x)=0

do i=1,15000
st_ct1(i)=0
enddo

print*,'nstrtIIIIIIIIIIIIIIIII',nstrt

jj=1
do m19=1,n
write(*,231)'str2:qarng',m19,(str2(m19,m18),m18=1,nae),q_fac(m19),symq(m19)
!print*,'str2',m19,q_fac(m19)
if(m19.eq.1)qul1(1)=q_fac(1)
j=jj
do i=1,j
if(qul1(i).eq.q_fac(m19))goto 373
enddo
jj=jj+1
qul1(jj)=q_fac(m19)
!print*,'qul1',qul1(jj)
373 enddo
nqul1=jj

print*,'nqul1',nqul1


do k6=1,15000
qual3(k6)=0
enddo

print*,'check1'
do k6=1,nqul1
do k7=1,nqul1
do k8=1,k6
if(qual3(k8).eq.qul1(k7))goto 389
enddo
if(qul1(k7).gt.qual3(k6))then
qual3(k6)=qul1(k7)
print*,'**',qual3(k6)
endif
389 enddo
enddo
print*,'check2'

k7=0
do k6=jj,1,-1
k7=k7+1
qual4(k7)=qual3(k6)
enddo

print*,'check3'

do m18=1,15000
do m19=1,nae
str1(m18,m19)=0
enddo
enddo

print*,n,noqult
!!!!! for the non-symmetric system .. if user dont want to take the symmetry in the calculations !!!!!!!!! 
if(symm.eq.0)then
i4=0
do m18=1,nqul1
do m19=1,n
!print*,q_fac(m19),qult(m18)
if(q_fac(m19).eq.qual4(m18))then 
i4=i4+1
do i3=1,nae
str1(i4,i3)=str2(m19,i3)
enddo
q_fac1(i4)=q_fac(m19)
stq1(i4)=str_quality_1(m19)
stq2(i4)=str_quality_2(m19)
bdq(i4)=bondq(m19)
!lfst1(i4)=lfst(m19)
endif
enddo
enddo


endif

!do m18=1,n
!print*,'sourav1',m18,(str2(m18,m19),m19=1,nae),qulsym(m18),symq(m18),q_fac(m18),tnqs,nssym,nqul1
!enddo

!!!!! for the symmetric system .. if user opted to take the symmetry calculations !!!!!!!!! 
if(symm.ne.0)then
do i=1,n
print*,'i,qfac',i,q_fac(i),symq(i)
enddo

print*,'check3'
ind=0
do i=1,n
if(ind.lt.q_fac(i))ind=q_fac(i)
enddo
mq_fac=ind

ind=0
do i=1,n
if(ind.lt.symq(i))ind=symq(i)
enddo
msymq=ind

ind=0
do i=1,n
if(ind.lt.loopsymsc(i))ind=loopsymsc(i)
enddo
mloopsymsc=ind
print*,'mq_fac,msymq,mloopsymsc',mq_fac,msymq,mloopsymsc

bnd=(nae-nl*2-nlast)/2
m=0
ll=0
mm=0
301 do i=1,mq_fac
if(mloopsymsc.ne.0)then
do k=1,mloopsymsc
do j=1,msymq
lll=1
l=0
ll=ll+1
do i1=1,n
do i2=1,mm
if(str_cnt1(i2).eq.i1)goto 307
enddo
if(q_fac(i1).eq.i)then
if(loopsymsc(i1).eq.k)then
if(symq(i1).eq.j)then
lll=0
m=m+1
mm=mm+1
print*,'mmmmm',m,i1,i,j,ll,loopsymsc(i1)
l=l+1

str_cnt1(mm)=i1
sym_str_sl(ll,l)=i1
sym_str_qual(ll)=q_fac(i1)
sym_str_qual1(ll,l)=q_fac(i1)
endif
endif
endif
307 enddo
if(lll.eq.0)sym_str_num(ll)=l

if(lll.eq.1)ll=ll-1
if(m.eq.n)then

nssym1=ll
goto 101
endif
enddo
enddo
endif
if(mloopsymsc.eq.0)then
do j=1,msymq
lll=1
l=0
ll=ll+1
do i1=1,n
do i2=1,mm
if(str_cnt1(i2).eq.i1)goto 407
enddo
if(q_fac(i1).eq.i)then
if(symq(i1).eq.j)then
lll=0
m=m+1
mm=mm+1
print*,'mmmmm',m,i1,i,j,ll,loopsymsc(i1)
l=l+1

str_cnt1(mm)=i1
sym_str_sl(ll,l)=i1
sym_str_qual(ll)=q_fac(i1)
sym_str_qual1(ll,l)=q_fac(i1)
endif
endif
407 enddo
if(lll.eq.0)sym_str_num(ll)=l

if(lll.eq.1)ll=ll-1
if(m.eq.n)then

nssym1=ll
goto 101
endif
enddo
endif
enddo
print*,'check4',m,n
goto 301

101 do i2=1,nssym1
print*,(sym_str_sl(i2,i3),i3=1,sym_str_num(i2)),'|',(sym_str_qual1(i2,i3),i3=1,sym_str_num(i2)),'|',sym_str_qual(i2)
enddo

!print*,'sourav_is_here'
!do i2=1,nssym1
!print*,'sym_str_num',sym_str_num(i2)
!enddo


if(symtype.eq.'check')then
call sym_check(nl,str2,n,sym_str_sl,sym_str_num,nssym1)
endif

!nstsymset=0
!ll=0
!if(nstrt.ne.0)then
!do i1=1,nstrt
!print*,(strt_struc(i1,i5),i5=1,nae)
!do i3=1,n
!!print*,'i3333',i3
!m=0
!do i4=nl*2+1,nl*2+bnd*2,2
!do i7=nl*2+1,nl*2+bnd*2,2
!l=0
!do i5=i4,i4+1
!do i8=i7,i7+1
!!print*,'str6,str2',str6(i1,i5),str2(i3,i8)
!if(strt_struc(i1,i5).eq.str2(i3,i8))then
!l=l+1
!goto 109
!endif
!enddo
!109 enddo
!!print*,'llll',l
!if(l.eq.2)then
!m=m+1
!goto 111
!endif
!enddo
!goto 112
!111 enddo
!!print*,'mmmm',m,bnd
!if(m.eq.bnd)then
!!print*,'i3i3i3',i3
!do i2=1,nssym1
!do i4=1,sym_str_num(i2)
!if(sym_str_sl(i2,i4).eq.i3)then
!ll=ll+1
!sym_str_sl1(ll)=i2
!print*,'sym_str_sl1(ll)',sym_str_sl1(ll)
!endif
!enddo
!enddo
!endif
!112 enddo
!enddo
!!!!! nstsymset=number of starting symmetry set
!nstsymset=ll

!print*,'ll',ll,sym_str_sl1(1)
!
!
!do i2=1,nssym1
!write(*,*),'sym_str_sl',(sym_str_sl(i2,i3),i3=1,sym_str_num(i2))
!enddo
!
!do i4=1,nstsymset
!do i3=1,sym_str_num(sym_str_sl1(i4))
!sym_str_sl_final(i4,i3)=sym_str_sl(sym_str_sl1(i4),i3)
!enddo
!sym_str_num_final(i4)=sym_str_num(sym_str_sl1(i4))
!enddo 

!i9=0
!do i4=1,nstsymset
!do i3=1,sym_str_num(sym_str_sl1(i4))
!i9=i9+1
!str_cnt(i9)=sym_str_sl(sym_str_sl1(i4),i3)
!symq_final(i9)=symq(sym_str_sl(sym_str_sl1(i4),i3))
!do j=1,nae
!str1(i9,j)=str2(sym_str_sl(sym_str_sl1(i4),i3),j)
!enddo
!q_fac1(i9)=q_fac(sym_str_sl(sym_str_sl1(i4),i3))
!stq1(i9)=str_quality_1(sym_str_sl(sym_str_sl1(i4),i3))
!stq2(i9)=str_quality_2(sym_str_sl(sym_str_sl1(i4),i3))
!bdq(i9)=bondq(sym_str_sl(sym_str_sl1(i4),i3))
!enddo
!enddo



!endif

l=0
do i=1,1000
do i2=1,nssym1
!if(nstrt.ne.0)then
!do i4=1,nstsymset
!if(i2.eq.sym_str_sl1(i4))goto 115
!enddo
!endif
if(sym_str_qual(i2).eq.i)then
l=l+1
sym_str_num_final(l)=sym_str_num(i2)
!print*,'sym_str_num_final(l)',sym_str_num_final(l)
do i3=1,sym_str_num(i2)
sym_str_sl_final(l,i3)=sym_str_sl(i2,i3)
!print*,'sym_str_sl_final(l,i3)',sym_str_sl_final(l,i3)
enddo
endif
115 enddo
enddo
!endif

!do i2=1,nssym1
!print*,'sym_str_num_final(i2)',sym_str_num_final(i2)
!enddo

if(set_order.eq.1.or.set_order.eq.2)then
do l=1,nssym1
sym_str_num(l)=sym_str_num_final(l)
do i3=1,sym_str_num(l)
sym_str_sl(l,i3)=sym_str_sl_final(l,i3)
enddo
enddo

endif


!!!bigger set is going top
if(set_order.eq.2)then
l=ll
do i=100,1,-1
do i2=1,nssym1
!if(nstrt.ne.0)then
!do i4=1,nstsymset
!if(i2.eq.sym_str_sl1(i4))goto 113
!enddo
!endif
if(sym_str_num(i2).eq.i)then
l=l+1
sym_str_num_final(l)=sym_str_num(i2)
do i3=1,sym_str_num(i2)
sym_str_sl_final(l,i3)=sym_str_sl(i2,i3)
enddo
endif
113 enddo
enddo
endif


!!!smaller sets is going top
if(set_order.eq.1)then

l=ll
do i=1,100
do i2=1,nssym1
!if(nstrt.ne.0)then
!do i4=1,nstsymset
!if(i2.eq.sym_str_sl1(i4))goto 114
!enddo
!endif
if(sym_str_num(i2).eq.i)then
l=l+1
sym_str_num_final(l)=sym_str_num(i2)
do i3=1,sym_str_num(i2)
sym_str_sl_final(l,i3)=sym_str_sl(i2,i3)
enddo
endif
114 enddo
enddo

endif


print*,'*********************************'
print*,'nssym1',nssym1
do i2=1,nssym1
print*,'i2',i2
print*,'sym_str_num_final(i2)',sym_str_num_final(i2)
write(*,906)(sym_str_sl_final(i2,i3),i3=1,sym_str_num_final(i2))
enddo

906 format (50I3)

i9=0
do i1=1,nssym1
do i2=1,sym_str_num_final(i1)
i9=i9+1
str_cnt(i9)=sym_str_sl_final(i1,i2)
do j=1,nae
str1(i9,j)=str2(sym_str_sl_final(i1,i2),j)
enddo
q_fac1(i9)=i1
!q_fac1(i9)=q_fac(sym_str_sl_final(i1,i2))
stq1(i9)=str_quality_1(sym_str_sl_final(i1,i2))
stq2(i9)=str_quality_2(sym_str_sl_final(i1,i2))
bdq(i9)=bondq(sym_str_sl_final(i1,i2))
enddo
enddo

do i9=1,n
write(*,231)'str2str2',i9,(str1(i9,m18),m18=1,nae),q_fac1(i9),bdq(i9)
enddo
231 format(a,30I3)
endif
do i=1,n
str_quality_1(i)=stq1(i)
str_quality_2(i)=stq2(i)
bondq(i)=bdq(i)
enddo

do i=1,n
write(*,231)'symmstr',i,(str1(i,m18),m18=1,nae),q_fac1(i),bdq(i),str_quality_1(i),str_quality_2(i)
enddo

!!! if user want some symmetric set must be availeble stats bellow !!!

allocate(str3(x,y))
str3(x,y)=0
if(nstrt.ne.0.and.symm.eq.1)then
do j=1,n
do k=1,nstrt
if(nl.ne.0)then
sl=0
do i=1,nl*2,2
do i1=1,nl*2,2
if(str1(j,i).eq.strt_struc(k,i1))then
sl=sl+1
endif
enddo
enddo
if (sl.eq.nl)goto 372
goto 471
endif

372 sll=0
do i=nl*2+1,nae-nlast,2
do i1=nl*2+1,nae-nlast,2
sl=0
do i3=i,i+1
do i4=i1,i1+1
if(str1(j,i3).eq.strt_struc(k,i4))then
sl=sl+1
if(sl.eq.2)then
sll=sll+1
goto 473
endif
endif
enddo
enddo
enddo
473 enddo
if(sll.eq.(nae-nlast-nl*2)/2)goto 374
goto 471

374 if(nlast.ne.0)then
sl=0
do i=nae-(nlast-1),nae
do i1=nae-(nlast-1),nae
if(str1(j,i).eq.strt_struc(k,i1))then
sl=sl+1
endif
enddo
enddo
if(sl.eq.nlast)goto 375
goto 471
endif

375 sl_qfac(k)=q_fac1(j)
print*,'sl_qfac',sl_qfac(k)
471 enddo
enddo

l=1
sl_qfac1(l)=sl_qfac(l)
do i=1,nstrt
p=sl_qfac(i)
if(sl_qfac(i-1).eq.p.or.i.eq.1) goto 477
l=l+1
sl_qfac1(l)=p
print*,'sl_qfac(i)',sl_qfac(i)
477 enddo
nstrt1=l

print*,'lllll',l
do i=1,nstrt1
print*,'sl_qfac1(l)',sl_qfac1(i)
enddo

j1=0
l=0
do j=1,n
do i=1,nstrt1
if(sl_qfac1(i).eq.q_fac1(j))then
l=l+1
sl_num(l)=j
print*,'sl_num(l)',sl_num(l)
goto 376
endif
enddo
!j1=j1+1
!print*,'jjjjjjjjj',j,j1
!do i1=1,nae
!str3(j1,i1)=str1(j,i1)
!enddo
!q_fac3(j1)=q_fac1(j)
!stq1(j1)=str_quality_1(j)
!stq2(j1)=str_quality_2(j)
!bdq(j1)=bondq(j)
376 enddo

do i=1,j1
print*,'i,str3(i,i1)',i,(str3(i,i1),i1=1,nae)
enddo
do j=1,n
do i=1,l
if(sl_num(i).ne.j)goto 377
j1=j1+1
do i1=1,nae
str3(j1,i1)=str1(j,i1)
enddo
q_fac3(j1)=q_fac1(j)
stq1(j1)=str_quality_1(j)
stq2(j1)=str_quality_2(j)
bdq(j1)=bondq(j)
377 enddo
enddo

do j=1,n
do i=1,l
if(sl_num(i).eq.j)goto 378
enddo
j1=j1+1
do i1=1,nae
str3(j1,i1)=str1(j,i1)
enddo
q_fac3(j1)=q_fac1(j)
stq1(j1)=str_quality_1(j)
stq2(j1)=str_quality_2(j)
bdq(j1)=bondq(j)
378 enddo


do j=1,n
print*,'j,rearranged',j,(str3(j,i1),i1=1,nae),q_fac3(j)
enddo
do i=1,n
do i2=1,nae
str1(i,i2)=str3(i,i2)
enddo
q_fac1(i)=q_fac3(i)
str_quality_1(i)=stq1(i)
str_quality_2(i)=stq2(i)
bondq(i)=bdq(i)
enddo
endif
deallocate(str3)
!stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! below part works when user wish to have some structures always in the top of the list !!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(str3(x,y))
str3(x,y)=0

if(nstrt.ne.0.and.symm.eq.0)then
jj=0
if(nl.ne.0)then
jj=0
do i3=1,nstrt
ii=0
do i=1,nl*2,2 
do i2=1,nl*2,2
if(strt_struc(i3,i).eq.str1(1,i2))then
ii=ii+1
endif
enddo
enddo
if(ii.eq.nl)then
jj=jj+1
st_ct(jj)=i3
endif
enddo
endif

print*,'***********nl,jj***************',nl,jj
if(jj.ne.0)then
kk=0
do i=1,n
!print*,'str1',(str1(i,m19),m19=1,nae)
do k8=1,jj
iii=0
k9=st_ct(k8)
!print*,'strt_struc',(strt_struc(k9,m19),m19=1,nae)
do i2=nl*2+1,nae-nlast,2
do k6=nl*2+1,nae-nlast,2
ii=0
do i4=i2,i2+1
do i3=k6,k6+1
!print*,'k9,i3,i4,strt_struc(k9,i3),str1(n,i4)',k9,i3,i4,strt_struc(k9,i3),str1(n,i4)
if(strt_struc(k9,i3).eq.str1(i,i4))then
ii=ii+1
!print*,'ii',ii
endif
enddo
enddo
if(ii.eq.2)then
iii=iii+ii
!print*,'ii,iii**',ii,iii,nae-nl*2-nlast
if(iii.eq.nae-nl*2-nlast)then
kk=kk+1
st_ct1(kk)=i
!print*,'st_ct1(kk)',st_ct1(kk)
if(kk.eq.jj)goto 340
goto 360
endif
goto 350
endif
enddo
350 enddo
enddo
360 enddo


340 jj=0
do i=1,kk
jj=jj+1
do m18=1,nae
str3(jj,m18)=str1(st_ct1(i),m18)
enddo
q_fac2(jj)=q_fac1(st_ct1(i))
stq3(jj)=stq1(st_ct1(i))
stq4(jj)=stq2(st_ct1(i))
bdq1(jj)=bdq(st_ct1(i))
enddo

do i=1,n
do i2=1,kk
if(i.eq.st_ct1(i2))goto 370
enddo
jj=jj+1
do m18=1,nae
str3(jj,m18)=str1(i,m18)
enddo
q_fac2(jj)=q_fac1(i)
stq3(jj)=stq1(i)
stq4(jj)=stq2(i)
bdq1(jj)=bdq(i)
370 enddo

print*,'souravooooo'
do i=1,n
q_fac1(i)=0
str_quality_1(i)=0
str_quality_2(i)=0
bondq(i)=0
do i2=1,nae
str1(i,i2)=0
enddo
enddo

do i=1,n
q_fac1(i)=q_fac2(i)
str_quality_1(i)=stq3(i)
str_quality_2(i)=stq4(i)
bondq(i)=bdq1(i)
do i2=1,nae
str1(i,i2)=str3(i,i2)
enddo
enddo

endif


if(nl.eq.0)then
kk=0
do i=1,n
!print*,'str1',(str1(i,m19),m19=1,nae)
do k9=1,nstrt
iii=0
print*,'strt_struc',(strt_struc(k9,m19),m19=1,nae)
do i2=1,nae-nlast,2
do k6=1,nae-nlast,2
ii=0
do i4=i2,i2+1
do i3=k6,k6+1
!print*,'k9,i3,i4,strt_struc(k9,i3),str1(n,i4)',k9,i3,i4,strt_struc(k9,i3),str1(n,i4)
if(strt_struc(k9,i3).eq.str1(i,i4))then
ii=ii+1
!print*,'ii',ii
endif
enddo
enddo
if(ii.eq.2)then
iii=iii+ii
!print*,'ii,iii**',ii,iii,nae-nl*2-nlast
if(iii.eq.nae-nlast)then
kk=kk+1
st_ct1(kk)=i
!print*,'st_ct1(kk)',st_ct1(kk)
if(kk.eq.nstrt)goto 341
goto 361
endif
goto 351
endif
enddo
351 enddo
enddo
361 enddo


341 jj=0
do i=1,kk
jj=jj+1
do m18=1,nae
str3(jj,m18)=str1(st_ct1(i),m18)
enddo
q_fac2(jj)=q_fac1(st_ct1(i))
stq3(jj)=stq1(st_ct1(i))
stq4(jj)=stq2(st_ct1(i))
bdq1(jj)=bdq(st_ct1(i))
enddo

do i=1,n
do i2=1,kk
if(i.eq.st_ct1(i2))goto 371
enddo
jj=jj+1
do m18=1,nae
str3(jj,m18)=str1(i,m18)
enddo
q_fac2(jj)=q_fac1(i)
stq3(jj)=stq1(i)
stq4(jj)=stq2(i)
bdq1(jj)=bdq(i)
371 enddo

print*,'souravooooo'
do i=1,n
q_fac1(i)=0
str_quality_1(i)=0
str_quality_2(i)=0
bondq(i)=0
do i2=1,nae
str1(i,i2)=0
enddo
enddo

do i=1,n
q_fac1(i)=q_fac2(i)
str_quality_1(i)=stq3(i)
str_quality_2(i)=stq4(i)
bondq(i)=bdq1(i)
do i2=1,nae
str1(i,i2)=str3(i,i2)
enddo
enddo

endif


endif

do m18=1,n
Print*,'str1*******',m18,(str1(m18,m19),m19=1,nae),q_fac1(m18),bondq(m18)
!,str_quality_1(m18),str_quality_2(m18),bondq(m18),qulsym(m18),symq(m18)
!print*,'q_fac1',q_fac1(m18)
enddo
!do i=1,kk
!print*,'st_ct1(kk)',st_ct1(i),kk
!enddo
!print*,'i4',i4
deallocate(str3)
deallocate(q_fac3)
deallocate(q_fac5)
deallocate(q_fac6)
deallocate(q_fac7)
deallocate(q_fac2)
deallocate(qul1)
deallocate(qual3)
deallocate(qual4)
deallocate(bdq)
deallocate(bdq1)
deallocate(str_cnt)
deallocate(stq1)
deallocate(stq2)
deallocate(stq3)
deallocate(stq4)
deallocate(sl_num)
deallocate(sl_qfac)
deallocate(sl_qfac1)
print*,'exit qult_str_arrange'

!stop
return
end subroutine qult_str_arrange
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
