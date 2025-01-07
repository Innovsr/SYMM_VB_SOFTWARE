!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! All possible Covalent structures are generated here it has three parts 
!!! lone-paires, covalent-bonds, and unpaired or radicals, each part generated 
!!! separately 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cov_struc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use mpi
use commondat
implicit none

common/quality/str_quality_1,str_quality_2,bondq,tqlty,bqlty,sqlty,tnqs,nssym,qulsym,symq,&
sigsym,tnqs_sig

integer::ijk,ij,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,perm_nstr,nssym,totset,&
i14,i15,i16,i17,i18,i19,i20,i21,k1,k2,k3,k4,k5,STDOUT,a,c,b,d,e,f,nnn,elporb,x,y
integer::i,j,k,l,j1,m,j2,ii,jj,iii,n1,j3,j4,j5,wig2,tnqs,lonep,bonds,allp,i_local,i_global
integer::alstr,fullcovstr(15000,20),symq(15000),sigsym(15000),tnqs_sig&
,qulsym(15000),str_quality_1(15000),str_quality_2(15000),bondq(15000),tqlty,bqlty,sqlty,symqq(15000)
real*8::factorial,symsc(15000),T1,T2
integer::ierr, rank, size, i_start, i_end
integer,dimension(:,:),allocatable::strc,fstr,strct,num
integer,dimension(:),allocatable::n,nn


!!! memory allocation !!!!!!!
x=100000
y=100

allocate(strc(x,y))
allocate(strct(x,y))
allocate(n(x))
allocate(nn(x))
allocate(num(x,y))
!!! initialization of arrays !!!
strc(x,y)=0
strct(x,y)=0
n(x)=0
nn(x)=0
num(x,y)=0
!!! ionic and covalent flags
flg_ion=0
flg_cov=1
!!!! new structures generation part starts !!!!!
vacorb=nao-nae
lonep=nae-nao
if(vacorb.gt.0)lonep=0
bonds=(nae-lonep*2-nlast)/2
d=nao
e=2
f=d-e
c=factorial(d)/(factorial(e)*factorial(f))
!!!! production of the set of bonded orbitals strats !!!!!
i4=0
do i1=1,nao-1
  do i2=i1+1,nao
    i4=i4+1
    do i3=1,2
    i5=i2
      if(i3.eq.1)i5=i1
    strc(i4,i3)=i5
    n(i4)=i1
    enddo
  enddo
enddo
totset=i4
!!!! production of the set of bonded orbitals ends !!!!!
!!!! production of the covalent bonding part start !!!!!
!!!! nested loops generates all combinations of bons !!!
!!!! for the covalent past of the structures; 7 nested !
!!!! loops can generate upto seven bond structures !!!!!

  ! Divide the `i1` loop among processes
  i_start = rank * (totset / size) + 1
  i_end = min((rank + 1) * (totset / size), totset)
  !local variable : i_local & global variable: i_global

i_local=0
j=0
do i1=i_start, i_end
j=1
nn(1)=i1
  if(j.eq.bonds)then
  i_local=i_local+1
  k4=0
    do k2=1,bonds
      do k5=1,2
      k4=k4+1
      strct(i_local,k4)=strc(nn(k2),k5)
      enddo
    enddo
!  goto 121
  cycle
  endif

do i2=i1+1,totset
j=2
l=0
  do k=1,2
    if(strc(nn(1),k).eq.n(i2))then
    l=l+1
    endif
  enddo
  if(l.eq.0)goto 142
!    goto 122
  cycle
142 l=0 
  do k=1,2
    do k1=1,2
      if(strc(i2,k).eq.strc(nn(1),k1))then
      l=l+1
      endif
    enddo
  enddo
  if(l.eq.0)goto 162
!  goto 122
cycle
162 nn(2)=i2
  if(j.eq.bonds)then
    i=i+1
    k4=0
    do k2=1,bonds
      do k5=1,2
        k4=k4+1
        strct(i,k4)=strc(nn(k2),k5)
      enddo
    enddo
!goto 122
    cycle
  endif

do i3=i2+1,totset
  j=3
  l=0
  do k1=1,j-1
    do k=1,2
      if(strc(nn(k1),k).eq.n(i3))then
        l=l+1
      endif
    enddo
  enddo
if(l.eq.0)goto 143
!goto 123
  cycle
143 l=0 
  do k2=1,j-1
    do k=1,2
      do k1=1,2
        if(strc(i3,k).eq.strc(nn(k2),k1))then
          l=l+1
        endif
      enddo
    enddo
  enddo
if(l.eq.0)goto 163
!goto 123
cycle
163 nn(3)=i3
  if(j.eq.bonds)then
    i=i+1
    k4=0
    do k2=1,bonds
      do k5=1,2
        k4=k4+1
        strct(i,k4)=strc(nn(k2),k5)
      enddo
    enddo
!   goto 123
   cycle
  endif

do i4=i3+1,totset
  j=4
  l=0
  do k1=1,j-1
    do k=1,2
      if(strc(nn(k1),k).eq.n(i4))then
        l=l+1
      endif
    enddo
  enddo
if(l.eq.0)goto 144
!goto 124
cycle
144 l=0 
  do k2=1,j-1
    do k=1,2
      do k1=1,2
        if(strc(i4,k).eq.strc(nn(k2),k1))then
          l=l+1
        endif
      enddo
    enddo
  enddo
if(l.eq.0)goto 164
!goto 124
cycle
164 nn(4)=i4
  if(j.eq.bonds)then
    i=i+1
    k4=0
    do k2=1,bonds
      do k5=1,2
        k4=k4+1
        strct(i,k4)=strc(nn(k2),k5)
      enddo
    enddo
!goto 124
  cycle
  endif


do i5=i4+1,totset
  j=5
  l=l+1
  do k1=1,j-1
    do k=1,2
      if(strc(nn(k1),k).eq.n(i5))then
        l=l+1
      endif
    enddo
  enddo
if(l.eq.0)goto 145
!goto 125
  cycle
145 l=0 
  do k2=1,j-1
    do k=1,2
      do k1=1,2
        if(strc(i5,k).eq.strc(nn(k2),k1))then
          l=l+1
        endif
      enddo
    enddo
  enddo
if(l.eq.0)goto 165
!goto 125
  cycle
165 nn(5)=i5
  if(j.eq.bonds)then
    i=i+1
    k4=0
    do k2=1,bonds
      do k5=1,2
        k4=k4+1
        strct(i,k4)=strc(nn(k2),k5)
      enddo
    enddo
!goto 125
  cycle
  endif

do i6=i5+1,totset
  j=6
  l=0
  do k1=1,j-1
    do k=1,2
      if(strc(nn(k1),k).eq.n(i6))then
        l=l+1
      endif
    enddo
  enddo
if(l.eq.0)goto 146
!goto 126
  cycle
146 l=0 
  do k2=1,j-1
    do k=1,2
      do k1=1,2
        if(strc(i6,k).eq.strc(nn(k2),k1))then
          l=l+1
        endif
      enddo
    enddo
  enddo
if(l.eq.0)goto 166
!goto 126
  cycle
166 nn(6)=i6
  if(j.eq.bonds)then
    i=i+1
    k4=0
    do k2=1,bonds
      do k5=1,2
        k4=k4+1
        strct(i,k4)=strc(nn(k2),k5)
      enddo
    enddo
!goto 126
  cycle
  endif

do i7=i6+1,totset
  j=7
  l=0
  do k1=1,j-1
    do k=1,2
      if(strc(nn(k1),k).eq.n(i7))then
        l=l+1
      endif
    enddo
  enddo
if(l.eq.0)goto 147
!goto 127
  cycle
  147 l=0 
  do k2=1,j-1
    do k=1,2
      do k1=1,2
        if(strc(i7,k).eq.strc(nn(k2),k1))then
          l=l+1
        endif
      enddo
    enddo
  enddo
if(l.eq.0)goto 167
!goto 127
  cycle
  167 nn(7)=i7
  if(j.eq.bonds)then
    i=i+1
    k4=0
    do k2=1,bonds
      do k5=1,2
        k4=k4+1
        strct(i,k4)=strc(nn(k2),k5)
      enddo
    enddo
!goto 127
  cycle
  endif

127 enddo
126 enddo
125 enddo
124 enddo
123 enddo
122 enddo
121 enddo
!!!! production of the covalent bonding part of the structures ends !!!!!
!!!! bond part of the structures stored in 'strct' matrix

! Collect results from all processes
  call MPI_REDUCE(i_local, i_global, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
   ! Output results from the master process
  if (rank == 0) then
    print *, "Total structures processed:", i_global
  end if


alstr=i_global ! total number of possible structures revealed

strc(x,y)=0
n(x)=0
nn(x)=0


do i1=1,alstr
do k4=1,bonds*2
strc(i1,k4)=strct(i1,k4)
enddo
enddo

!!!!! production of the lone pairs and radical part of the structures starts !!!!!!
if(nlast.ne.0.or.lonep.ne.0)then
i=0
  do ii=1,alstr
  k5=0
  i2=0
   do i1=1,nao
     do k1=1,bonds*2
       if(i1.eq.strct(ii,k1)) goto 530
     enddo
i2=i2+1
num(ii,i2)=i1

530 enddo
allp=i2
!!!!! production of the lone pair part of the structures starts !!!!!!
!!!!! 10 nested loops can creat 10 lone paires in each structures !!!!
!!!!! due to combinations of the lone paires the number of structures!
!!!!! may increase.

if(lonep.ne.0)then


  ! Divide the `i1` loop among processes
  i_start = rank * (allp / size) + 1
  i_end = min((rank + 1) * (allp / size), allp)
  !local variable : i_local & global variable: i_global

j=0
m=0
  do i1=1,i_start, i_end
  j=1
  n(1)=i1
    if (j.eq.lonep) then
    i_local=i_local+1
    j4=0
      do j1=lonep,1,-1
        do j2=1,2
        j4=j4+1
        strc(i_local,j4)=num(ii,n(j1))
        enddo
      enddo
    j2=0
      do j3=1,bonds*2
      j2=j3+j4
      strc(i_local,j2)=strct(ii,j3)
      enddo
    goto 399
    endif

do i2=i1+1,allp
j=2
  if(i2.eq.n(1))goto 401
  n(2)=i2
    if (j.eq.lonep) then
    i=i+1
    j4=0
      do j1=1,lonep
        do j2=1,2
        j4=j4+1
        strc(i,j4)=num(ii,n(j1))
        enddo
      enddo
    j2=0
      do j3=1,bonds*2
      j2=j3+j4
      strc(i,j2)=strct(ii,j3)
      enddo
      goto 401
    endif

do i3=i2+1,allp
j=3
 do k=1,2
  if(i3.eq.n(k))goto 402
 enddo
 n(3)=i3
   if (j.eq.lonep) then
     i=i+1
     j4=0
     do j1=1,lonep
       do j2=1,2
         j4=j4+1
         strc(i,j4)=num(ii,n(j1))
       enddo
     enddo
     j2=0
     do j3=1,bonds*2
       j2=j3+j4
       strc(i,j2)=strct(ii,j3)
     enddo
     goto 402
   endif

do i4=i3+1,allp
j=4
do k=1,3
if(i4.eq.n(k))goto 403
enddo
n(4)=i4
if (j.eq.lonep) then
i=i+1
j4=0
do j1=1,lonep
do j2=1,2
j4=j4+1
strc(i,j4)=num(ii,n(j1))
enddo
enddo
j2=0
do j3=1,bonds*2
j2=j3+j4
strc(i,j2)=strct(ii,j3)
enddo
goto 403
endif

do i5=i4+1,allp
j=5
do k=1,4
if(i5.eq.n(k))goto 404
enddo
n(5)=i5
if (j.eq.lonep) then
i=i+1
j4=0
do j1=1,lonep
do j2=1,2
j4=j4+1
strc(i,j4)=num(ii,n(j1))
enddo
enddo
j2=0
do j3=1,bonds*2
j2=j3+j4
strc(i,j2)=strct(ii,j3)
enddo
goto 404
endif

do i6=i5+1,allp
j=6
do k=1,5
if(i6.eq.n(k))goto 405
enddo
n(6)=i6
if (j.eq.lonep) then
i=i+1
j4=0
do j1=1,lonep
do j2=1,2
j4=j4+1
strc(i,j4)=num(ii,n(j1))
enddo
enddo
j2=0
do j3=1,bonds*2
j2=j3+j4
strc(i,j2)=strct(ii,j3)
enddo
goto 405
endif

do i7=i6+1,allp
j=7
do k=1,6
if(i7.eq.n(k))goto 406
enddo
n(7)=i7
if (j.eq.lonep) then
i=i+1
j4=0
do j1=1,lonep
do j2=1,2
j4=j4+1
strc(i,j4)=num(ii,n(j1))
enddo
enddo
j2=0
do j3=1,bonds*2
j2=j3+j4
strc(i,j2)=strct(ii,j3)
enddo
goto 406
endif

do i8=i7+1,allp
j=8
do k=1,7
if(i8.eq.n(k))goto 407
enddo
n(8)=i8
if (j.eq.lonep) then
i=i+1
j4=0
do j1=1,lonep
do j2=1,2
j4=j4+1
strc(i,j4)=num(ii,n(j1))
enddo
enddo
j2=0
do j3=1,bonds*2
j2=j3+j4
strc(i,j2)=strct(ii,j3)
enddo
goto 407
endif

do i9=i8+1,allp
j=9
do k=1,8
if(i9.eq.n(k))goto 408
enddo
n(9)=i9
if (j.eq.lonep) then
i=i+1
j4=0
do j1=1,lonep
do j2=1,2
j4=j4+1
strc(i,j4)=num(ii,n(j1))
enddo
enddo
j2=0
do j3=1,bonds*2
j2=j3+j4
strc(i,j2)=strct(ii,j3)
enddo
goto 408
endif

do i10=i9+1,allp
j=10
do k=1,9
if(i10.eq.n(k))goto 409
enddo
n(10)=i10
if (j.eq.lonep) then
i=i+1
j4=0
do j1=1,lonep
do j2=1,2
j4=j4+1
strc(i,j4)=num(ii,n(j1))
enddo
enddo
j2=0
do j3=1,bonds*2
j2=j3+j4
strc(i,j2)=strct(ii,j3)
enddo
goto 409
endif

409 enddo
408 enddo
407 enddo
406 enddo
405 enddo
404 enddo
403 enddo
402 enddo
401 enddo
399 enddo
endif
enddo

! Collect results from all processes
  call MPI_REDUCE(i_local, i_global, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
   ! Output results from the master process
  if (rank == 0) then
    print *, "Total structures processed:", i_global
  end if

if (lonep.ne.0) alstr=i_global ! new number of structures revealed
!!!!! production of the radical part of the structures starts !!!!!!
if(nlast.ne.0)then
  do i=1,alstr
  j4=nae-nlast
    do i1=1,nae
      do i2=1,nae-nlast
        if(strc(i,i2).eq.i1)goto 555
      enddo
      j4=j4+1
      strc(i,j4)=i1
555 enddo
enddo
endif
endif
!!!!! production of the lone pairs and radical part of the structures ends !!!!!!
deallocate(num)
990 format(30I3)
656 deallocate(strct)
!!!!!!! new structures generation part ends !!!!!!!!!!!!!!!!!!!!!!!
d=nao
e=nlp
f=nao-nlp
elporb=nao-nlp
c=factorial(d)/(factorial(e)*factorial(f))
call wigner(elporb,wig2)
deallocate(n)
deallocate(nn)

!!! header of the 'structures.dat' file is written here
write(7,*)'  '
write(7,*)'                  covalent structures   '
write(7,*)'                  -------------------   '
write(7,*)'  '

write(7,307)' You have',alstr,' covalent structures among'
write(7,307)'them',wig2*c,' covalent structures are permissible'
write(7,*)' '
write(7,*)'      [] and {} in front of the structures specifies the quality of the structures.'
write(7,*)'1st number in [] specifies intra-atomic bond quality and 2nd number specifies ' 
write(7,*)'symmetry breaking bond quality. Where {} specifies the overall quality of the structures'
write(7,*)'*************************************************************** '
307 format(2x,a,I7,a)


do i=1,alstr
do i1=1,nae
fullcovstr(i,i1)=strc(i,i1)+niao
enddo
enddo

deallocate(strc)

if(symm.eq.1) then
if(sig_sym_flg.eq.1)call symmetry_cal_sig(nlp,fullcovstr,alstr,symsc,symq,nssym)
if(sig_sym_flg.ne.1)call symmetry_cal_pi(nlp,fullcovstr,alstr,symsc,symq,nssym)

nnn=0
do ij=1,tnqs
do ijk=nssym,1,-1
do j2=1,alstr
if(symq(j2).eq.ijk)then
nnn=nnn+1
if(niao.eq.0)then
write(7,900)'cov structure',nnn,')',qulsym(j2),(fullcovstr(j2,j1),j1=1,nao-nlp+nlp*2)
endif
if(niao.ne.0)then
write(7,901)'cov structure',nnn,')',qulsym(j2),1,':',niao,(fullcovstr(j2,j1),j1=1,nao-nlp+nlp*2)
endif
endif
210 enddo
write(7,*)'******************************************************'
enddo
enddo
else

do j2=1,alstr
if(niao.eq.0)then
write(7,902)'cov structure',j2,')',(fullcovstr(j2,j1),j1=1,nao-nlp+nlp*2)
endif
if(niao.ne.0)then
write(7,903)'cov structure',j2,')',1,':',niao,(fullcovstr(j2,j1),j1=1,nao-nlp+nlp*2)
endif

enddo
write(7,*)'******************************************************'
endif
900 format(a,I5,a,x,I3,3x,30(I5))
901 format(a,I5,a,x,I3,3x,I3,x,a,I3,30(I5))
902 format(a,I5,a,x,30(I5))
903 format(a,I5,a,x,I3,x,a,I3,30(I5))


if(nfset.eq.3)write(23,*)'********** RUMER STRUCTURES ******************'
if(nfset.eq.5)write(10,*)'**********ALL POSSIBLE RUMER STRUCTURES ******************'
if(flg1.eq.1)then
write(10,*)'********** RUMER STRUCTURES ******************'
endif
if(flg1.eq.0.and.nfset.ne.5)then
write(10,*)'******* CHEM. QUAL. STRUCTURES **************'
endif

perm_nstr=wig2*c
write(10,*)perm_nstr,' covalent structures' 
write(23,*)perm_nstr,' covalent structures' 
if(nfset.ne.5)write(10,913)'IAB','NNB','SBB','fqual'
write(23,913)'IAB','NNB','SBB','fqual'
913 format(x,a,x,a,x,a,x,a)
call str_selection(fullcovstr,nlp,alstr,perm_nstr)
print*,'exit cove_struc'
return
end subroutine cov_struc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
