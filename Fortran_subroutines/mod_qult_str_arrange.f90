! This subroutine arrange the structures according to their quality values 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module mod_qult_str_arrange
use commondat_mod
use quality_mod
!use symm_check
implicit none

contains
subroutine qult_str_arrange(nl,str2,n,str1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer::i1,l,k,m,bnd,ll,p,nstrt1,mq_fac,msymq,mloopsymsc,ind,sl,sll
integer::i,i2,i3,i4,i9,m18,m19,k6,k7,k8,k9,ii,kk,n,jj,j,j1,nqul1,nl
integer::lll,nssym1,mm,iii,x,y
integer, allocatable::str3(:,:), sym_str_sl(:,:), sym_str_sl_final(:,:)
integer, allocatable::qual4(:),qual3(:),str_cnt(:),q_fac3(:),q_fac4(:),q_fac5(:),q_fac6(:),q_fac1(:)
integer, allocatable::stq1(:), stq2(:), stq3(:), stq4(:), sl_num(:), sl_qfac(:), sl_qfac1(:), str_cnt1(:)
integer, allocatable::st_ct(:), st_ct1(:), sym_str_num(:),q_fac7(:),q_fac2(:),bdq(:),bdq1(:),qul1(:)
integer, allocatable::sym_str_num_final(:), sym_str_qual(:), sym_str_qual1(:,:)
integer, pointer::str1(:,:),str2(:,:)

!set_order desides if the symmetric sets arranged according to quality or size
! set_order=0 means quality
! set_order=1 means  smaller to larger sets
! set_order=2 means larger to smaller sets

print*,'set_order',set_order
print*,'enter qult_str_arrange'
x=MaxStrOepo
y=nae

allocate(sym_str_qual1(x,y))
allocate(sym_str_num_final(x))
allocate(sym_str_qual(x))
allocate(sym_str_num(x))
allocate(sym_str_sl(x,y))
allocate(sym_str_sl_final(x,y))
allocate(st_ct(x))
allocate(st_ct1(x))
allocate(str_cnt1(x))
allocate(q_fac1(x))
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
allocate(str_cnt(x))
allocate(stq1(x))
allocate(stq2(x))
allocate(stq3(x))
allocate(stq4(x))
allocate(sl_num(x))
allocate(sl_qfac(x))
allocate(sl_qfac1(x))

sym_str_qual1=0
sym_str_num_final=0
sym_str_qual=0
sym_str_num=0
sym_str_sl_final=0
sym_str_sl=0
st_ct=0
st_ct1=0
str_cnt1=0
q_fac1=0
q_fac2=0
q_fac3=0
q_fac4=0
q_fac5=0
q_fac6=0
q_fac7=0
qul1=0
qual3=0
qual4=0
bdq=0
bdq1=0
str_cnt=0
stq1=0
stq2=0
stq3=0
stq4=0
sl_num=0
sl_qfac=0
sl_qfac1=0


jj=1
loop1:do m19=1,n
  !write(*,231)'str2:qarng',m19,(str2(m19,m18),m18=1,nae),q_fac(m19),symq(m19)
  if(m19.eq.1)qul1(1)=quality_fac(1)
  j=jj
  do i=1,j
    if(qul1(i).eq.quality_fac(m19)) cycle loop1
  enddo
  jj=jj+1
  qul1(jj)=quality_fac(m19)
enddo loop1
nqul1=jj

!print*,'nqul1',nqul1

  qual3=0

do k6=1,nqul1
  loop2:do k7=1,nqul1
    do k8=1,k6
      if(qual3(k8).eq.qul1(k7)) cycle loop2
    enddo
    if(qul1(k7).gt.qual3(k6))then
      qual3(k6)=qul1(k7)
    endif
  enddo loop2
enddo

k7=0
do k6=jj,1,-1
  k7=k7+1
  qual4(k7)=qual3(k6)
enddo

str1 = 0
!!!!! for the non-symmetric system .. if user dont want to take the symmetry in the calculations !!!!!!!!! 
if(symm.eq.0)then
  i4=0
  do m18=1,nqul1
    do m19=1,n
      if(quality_fac(m19).eq.qual4(m18))then 
        i4=i4+1
        do i3=1,nae
          str1(i4,i3)=str2(m19,i3)
        enddo
        q_fac1(i4)=quality_fac(m19)
        stq1(i4)=str_quality_1(m19)
        stq2(i4)=str_quality_2(m19)
        bdq(i4)=bondq(m19)
      endif
    enddo
  enddo
endif

!!!!! for the symmetric system .. if user opted to take the symmetry calculations !!!!!!!!! 
if(symm.ne.0)then
do i=1,n
  print*,'i,qfac',i,quality_fac(i),symq(i)
enddo

  mq_fac = 0
  msymq = 0
  mloopsymsc = 0
  
  ind=0
  do i=1,n
    if(ind.lt.quality_fac(i))ind=quality_fac(i)
  enddo
  mq_fac=ind  ! mq_fac : maximum quality factor
  
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
  
  !print*,'mq_fac,msymq,mloopsymsc',mq_fac,msymq,mloopsymsc
  !do i = 1,n
  !  print*,'symq',symq(i), q_fac(i),loopsymsc(i)
  !enddo
  
  
  bnd=(nae-nl*2-nlast)/2
  m=0
  ll=0
  mm=0
  do
    do i=1,mq_fac
      if(mloopsymsc.ne.0)then
        do k=1,mloopsymsc
          do j=1,msymq
            lll=1
            l=0
            ll=ll+1
            loop3:do i1=1,n
              do i2=1,mm
                if(str_cnt1(i2).eq.i1) cycle loop3
              enddo
              if (quality_fac(i1).eq.i)then
                if (loopsymsc(i1).eq.k)then
                  if (symq(i1).eq.j)then
                    lll=0
                    m=m+1
                    mm=mm+1
                    !print*,'mmmmm',m,i1,i,j,ll,loopsymsc(i1)
                    l=l+1
                    str_cnt1(mm)=i1
                    sym_str_sl(ll,l)=i1
                    sym_str_qual(ll)=quality_fac(i1)
                    sym_str_qual1(ll,l)=quality_fac(i1)
                  endif
                endif
              endif
            enddo loop3
            if(lll.eq.0)sym_str_num(ll)=l
    
            if(lll.eq.1)ll=ll-1
            if(m.eq.n)then
              nssym1=ll
              goto 101
            endif
          enddo
        enddo
      endif
    
      if (mloopsymsc.eq.0)then
        do j=1,msymq
          lll=1
          l=0
          ll=ll+1
          loop4:do i1=1,n
            do i2=1,mm
              if(str_cnt1(i2).eq.i1) cycle loop4
            enddo
            if (quality_fac(i1).eq.i)then
              if (symq(i1).eq.j)then
                lll=0
                m=m+1
                mm=mm+1
                l=l+1
                str_cnt1(mm)=i1
                sym_str_sl(ll,l)=i1
                sym_str_qual(ll)=quality_fac(i1)
                sym_str_qual1(ll,l)=quality_fac(i1)
              endif
            endif
          enddo loop4
          if (lll.eq.0)sym_str_num(ll)=l
    
          if (lll.eq.1)ll=ll-1
          if (m.eq.n)then
            nssym1=ll
            goto 101
            !exit
          endif
        enddo
      endif
    enddo
  enddo
  
  !if(symtype.eq.'check')then
  !  call sym_check(nl,str2,n,sym_str_sl,sym_str_num,nssym1)
  !endif
  
  !goto 301
  !101 do i2=1,nssym1
  !print*,(sym_str_sl(i2,i3),i3=1,sym_str_num(i2)),'|',(sym_str_qual1(i2,i3),i3=1,sym_str_num(i2)),'|',sym_str_qual(i2)
  !enddo
  
  
101  l=0
  do i=1,1000
    do i2=1,nssym1
      if(sym_str_qual(i2).eq.i)then
        l=l+1
        sym_str_num_final(l)=sym_str_num(i2)
        do i3=1,sym_str_num(i2)
          sym_str_sl_final(l,i3)=sym_str_sl(i2,i3)
        enddo
      endif
    enddo
  enddo
  
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
        if(sym_str_num(i2).eq.i)then
          l=l+1
          sym_str_num_final(l)=sym_str_num(i2)
          do i3=1,sym_str_num(i2)
            sym_str_sl_final(l,i3)=sym_str_sl(i2,i3)
          enddo
        endif
      enddo
    enddo
  endif
  
  
  !!!smaller sets is going top
  if(set_order.eq.1)then
    l=ll
    do i=1,100
      do i2=1,nssym1
        if(sym_str_num(i2).eq.i)then
          l=l+1
          sym_str_num_final(l)=sym_str_num(i2)
          do i3=1,sym_str_num(i2)
            sym_str_sl_final(l,i3)=sym_str_sl(i2,i3)
          enddo
        endif
      enddo
    enddo
  endif
  
  
  !print*,'nssym1',nssym1
  !do i2=1,nssym1
  !print*,'i2',i2
  !print*,'sym_str_num_final(i2)',sym_str_num_final(i2)
  !write(*,906)(sym_str_sl_final(i2,i3),i3=1,sym_str_num_final(i2))
  !enddo
  
 ! 906 format (50I3)
  
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
  
  !do i9=1,n
  !write(*,231)'str2str2',i9,(str1(i9,m18),m18=1,nae),q_fac1(i9),bdq(i9)
  !enddo
  !231 format(a,30I3)
endif
do i=1,n
  quality_fac(i)=q_fac1(i)
  str_quality_1(i)=stq1(i)
  str_quality_2(i)=stq2(i)
  bondq(i)=bdq(i)
enddo

!do i=1,n
!write(*,231)'symmstr',i,(str1(i,m18),m18=1,nae),q_fac1(i),bdq(i),str_quality_1(i),str_quality_2(i)
!enddo

!!! if user want some symmetric set must be availeble stats bellow !!!

allocate(str3(x,y))
str3(x,y)=0
if(nstrt.ne.0.and.symm.eq.1)then
  do j=1,n
    loop6:do k=1,nstrt
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
        !goto 471
        cycle loop6
      endif
372 sll=0
      loop5:do i=nl*2+1,nae-nlast,2
        do i1=nl*2+1,nae-nlast,2
          sl=0
          do i3=i,i+1
            do i4=i1,i1+1
              if(str1(j,i3).eq.strt_struc(k,i4))then
                sl=sl+1
                if(sl.eq.2)then
                  sll=sll+1
                  cycle loop5
                endif
              endif
            enddo
          enddo
        enddo
      enddo loop5
      if(sll.eq.(nae-nlast-nl*2)/2)goto 374
      !goto 471
      cycle loop6
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
        !goto 471
        cycle loop6
      endif

375 sl_qfac(k)=q_fac1(j)
    enddo loop6
  enddo

  l=1
  sl_qfac1(l)=sl_qfac(l)
  loop7:do i=1,nstrt
    p=sl_qfac(i)
    if(sl_qfac(i-1).eq.p.or.i.eq.1) cycle loop7
    l=l+1
    sl_qfac1(l)=p
  enddo loop7
  nstrt1=l

!print*,'lllll',l
!do i=1,nstrt1
!print*,'sl_qfac1(l)',sl_qfac1(i)
!enddo

  j1=0
  l=0
  loop8:do j=1,n
    do i=1,nstrt1
      if(sl_qfac1(i).eq.q_fac1(j))then
        l=l+1
        sl_num(l)=j
        !goto 376
        cycle loop8
      endif
    enddo
  enddo loop8

!do i=1,j1
!print*,'i,str3(i,i1)',i,(str3(i,i1),i1=1,nae)
!enddo
  do j=1,n
    loop9:do i=1,l
      if(sl_num(i).ne.j) cycle loop9
      j1=j1+1
      do i1=1,nae
        str3(j1,i1)=str1(j,i1)
      enddo
      q_fac3(j1)=q_fac1(j)
      stq1(j1)=str_quality_1(j)
      stq2(j1)=str_quality_2(j)
      bdq(j1)=bondq(j)
    enddo loop9
  enddo

  loop10:do j=1,n
    do i=1,l
      if(sl_num(i).eq.j) cycle loop10
    enddo
    j1=j1+1
    do i1=1,nae
      str3(j1,i1)=str1(j,i1)
    enddo
    q_fac3(j1)=q_fac1(j)
    stq1(j1)=str_quality_1(j)
    stq2(j1)=str_quality_2(j)
    bdq(j1)=bondq(j)
  enddo loop10


!do j=1,n
!print*,'j,rearranged',j,(str3(j,i1),i1=1,nae),q_fac3(j)
!enddo
  do i=1,n
    do i2=1,nae
      str1(i,i2)=str3(i,i2)
    enddo
    q_fac1(i)=q_fac3(i)
    quality_fac(i) = q_fac3(i)
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

  if(jj.ne.0)then
    kk=0
    loop11:do i=1,n
      do k8=1,jj
        iii=0
        k9=st_ct(k8)
        loop12:do i2=nl*2+1,nae-nlast,2
          do k6=nl*2+1,nae-nlast,2
            ii=0
            do i4=i2,i2+1
              do i3=k6,k6+1
                if(strt_struc(k9,i3).eq.str1(i,i4))then
                  ii=ii+1
                endif
              enddo
            enddo
            if(ii.eq.2)then
              iii=iii+ii
              if(iii.eq.nae-nl*2-nlast)then
                kk=kk+1
                st_ct1(kk)=i
                if(kk.eq.jj)goto 340
                !goto 360
                cycle loop11
              endif
              !goto 350
              cycle loop12
            endif
          enddo
        enddo loop12
      enddo
    enddo loop11


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

    loop13:do i=1,n
      do i2=1,kk
        if(i.eq.st_ct1(i2)) cycle loop13
      enddo
      jj=jj+1
      do m18=1,nae
        str3(jj,m18)=str1(i,m18)
      enddo
      q_fac2(jj)=q_fac1(i)
      stq3(jj)=stq1(i)
      stq4(jj)=stq2(i)
      bdq1(jj)=bdq(i)
    enddo loop13

    do i=1,n
      q_fac1(i)=0
      quality_fac(i)=0
      str_quality_1(i)=0
      str_quality_2(i)=0
      bondq(i)=0
      do i2=1,nae
        str1(i,i2)=0
      enddo
    enddo

    do i=1,n
      q_fac1(i)=q_fac2(i)
      quality_fac(i)=q_fac2(i)
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
    loop14:do i=1,n
      do k9=1,nstrt
        iii=0
        loop15:do i2=1,nae-nlast,2
          do k6=1,nae-nlast,2
            ii=0
            do i4=i2,i2+1
              do i3=k6,k6+1
                if(strt_struc(k9,i3).eq.str1(i,i4))then
                  ii=ii+1
                endif
              enddo
            enddo
            if(ii.eq.2)then
              iii=iii+ii
              if(iii.eq.nae-nlast)then
                kk=kk+1
                st_ct1(kk)=i
                if(kk.eq.nstrt)goto 341
                  !goto 361
                  cycle loop14
              endif
              !goto 351
              cycle loop15
            endif
          enddo
        enddo loop15
      enddo
    enddo loop14


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

    loop16:do i=1,n
      do i2=1,kk
        if(i.eq.st_ct1(i2)) cycle loop16
      enddo
      jj=jj+1
      do m18=1,nae
        str3(jj,m18)=str1(i,m18)
      enddo
      q_fac2(jj)=q_fac1(i)
      stq3(jj)=stq1(i)
      stq4(jj)=stq2(i)
      bdq1(jj)=bdq(i)
    enddo loop16

!print*,'souravooooo'
    do i=1,n
      q_fac1(i)=0
      quality_fac(i)=0
      str_quality_1(i)=0
      str_quality_2(i)=0
      bondq(i)=0
      do i2=1,nae
        str1(i,i2)=0
      enddo
    enddo

    do i=1,n
      q_fac1(i)=q_fac2(i)
      quality_fac(i) = q_fac2(i)
      str_quality_1(i)=stq3(i)
      str_quality_2(i)=stq4(i)
      bondq(i)=bdq1(i)
      do i2=1,nae
        str1(i,i2)=str3(i,i2)
      enddo
    enddo
  endif
endif

!do m18=1,n
!Print*,'str:quality_arrange',m18,(str1(m18,m19),m19=1,nae),q_fac1(m18),bondq(m18) !&
!!,str_quality_1(m18),str_quality_2(m18),bondq(m18),qulsym(m18),symq(m18)
!!!print*,'q_fac1',q_fac1(m18)
!enddo
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
end module mod_qult_str_arrange
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
