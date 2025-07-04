!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module mod_symmetry_cal_sig
use commondat_mod
use commondat1_mod
use loops_mod
use mod_symm_loops
use mod_nnat_bond_sig  
use mod_coordination_val
implicit none


contains
subroutine symmetry_cal_sig(nl,str1,ncqs,symq,nssym)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer::m19,m18,i,i1,i2,i3,i4,i5,nl,ncqs,k,k1,k2,n1,n2,n3,n4,n5,n6,l1
integer::nd,k3,k4,k5,n,nssym,j,jj,ii2,nsig,nnatom,nn4,nn5,nna
integer, allocatable::sig_orb(:),sig_orb_1(:)
real*8, allocatable::score(:,:),new_score(:,:)
real*8, allocatable::lone_score(:),order(:),ssym(:),stsymsc(:),bndscore(:)
integer::nnscore,ipnum,loop_score,loopsc,ncrow,bonds
integer, allocatable::nn_count(:), symq(:), coordination_mat(:,:)
real*8::tscore,rscore,strscore,scr,coord_score,lpna,losc,n7,n8
integer::numbond,numbond1
integer, pointer:: str1(:,:)
!double precision, pointer::stsymsc(:)
!real*8, pointer::stsymsc(:)


print*,'enter symmetry_cal_sig',ncqs

!! iabd--> inactive bonds associated with the atoms
!! ialp--> inactive lone pairs associated with the atoms
!! score()--> scoring of the atoms depending on the the 'iabd', 'ialp' and atomic number (at_number)

open(unit=15,file='symsig.temp',status='unknown')

if(.not. allocated(stsymsc))then
allocate(stsymsc(MaxStrOepo))
stsymsc = 0
endif
bonds=(nae-nl*2-nlast)/2
if (bonds.eq.0) goto 500
!!! initialise the score array, it stores the scores of each atoms
allocate(score(atom,2))

score = 0.0

!!! atom scoring starts here; scoring criteria: number of inactive bonds,
!!! number of inactive lone paires, number of charge if any associated with atoms
i4=0
do i3=1,atom
  lpna=0.0
  k=i3
  if(input_flg.eq.1)then
    do i1=1,niabd
      if(k.eq.iabd(i1))then
        lpna=lpna+1.0/prime_num(2)
      endif
    enddo
    do i2=1,nialp
      if(k.eq.ialp(i2))then
        lpna=lpna+1/prime_num(1)
      endif
    enddo
    do i2=1,niach
      if(k.eq.iach(i2))then
        lpna=lpna+1/prime_num(3)
      endif
    enddo
  else

  lpna=lpna+biasval(i3)+dist_nnat(i3)
  endif
  score(i3,2)=lpna+at_num(i3)
print*,'score:symme',score(i3,2),at_num(i3),biasval(i3),dist_nnat(i3)
enddo

if (.not. allocated(nnat_bond_new))then
   allocate(nnat_bond_new(nnnatom+nao, 2))
   nnat_bond_new = 0
endif

!! strating generating extended connectivity matrix require for sigma-pi system
n6=0
n7=0
n8=0
do i2=1,nnnatom
  do i1=1,2
    nnat_bond_new(i2,i1)=nnat_bond(i2,i1)+nao+niao
  enddo
enddo

!print*,'nnnatom',nnnatom

allocate(tot_orb((nnnatom+nao)*2))

n5=0
do i2=1,nnnatom
  loop1:do i1=1,2
    do i3=1,n5
      if(nnat_bond_new(i2,i1).eq.tot_orb(i3)) then
        cycle loop1
      endif
    enddo      
    n5=n5+1
    tot_orb(n5)=nnat_bond_new(i2,i1)
  enddo loop1
enddo      

allocate(nn_count(atom))

do i1=1,atom
  n1=0
  do i3=1,nnnatom
    do i4=1,2
      if(active_atoms(i1).eq.nnat_bond(i3,i4))then
        n1=n1+1
      endif
    enddo
  enddo
  nn_count(i1)=n1
enddo

deallocate(nn_count)

allocate(sig_orb(nao))
if (.not. allocated(sig_orb_1))then
  allocate(sig_orb_1(nao+atom))
  sig_orb_1 = 0
endif

n2=nnnatom
n4=0
n3=0
n6=0
do i3=1,atom
  n1=0
  do i4=1,atn(active_atoms(i3))
    do k2=1,syn(nsym)
      if(atoset(active_atoms(i3),i4).eq.atsymset(nsym,k2))then
        n1=n1+1
        n6=n6+1
        sig_orb(n1)=atoset(active_atoms(i3),i4)
        sig_orb_1(n6)=sig_orb(n1)
      endif
    enddo
 enddo
  if(n1.gt.0)then
    n3=n3+1
    k=0
    loop2:do i2=1,nnnatom
      do i1=1,2
        if(i3.eq.nnat_bond(i2,i1))then
          k=k+1
          if(k.gt.n1) then
            cycle loop2
          endif
          n4=n4+1
          nnat_bond_new(i2,i1)=sig_orb(k)
        endif
      enddo      
    enddo loop2
    do i1=1,n1
      n2=n2+1
      n5=n5+1
      nnat_bond_new(n2,1)=niao+nao+i3
      nnat_bond_new(n2,2)=sig_orb(i1)
      tot_orb(n5)=sig_orb(i1)
    enddo
    if(n1.gt.2)then
      do i1=1,n1
        if(i1.ne.n1)ii2=i1+1
        if(i1.eq.n1)ii2=1
          n2=n2+1
          nnat_bond_new(n2,1)=sig_orb(i1)
          nnat_bond_new(n2,2)=sig_orb(ii2)
      enddo
    endif
  endif
enddo
!!!!!!!!!!!!

deallocate(sig_orb)
nnatom=n2
nsig=n6
!!!! nsig=number of sigma orbital
!!!! nnatom= number of connected pair of orbitals in the extended CM

natom=n5
!!!! natom=number of extended orbitals

!do i2=1,nnnatom
!print*,'nnat_bond',(nnat_bond(i2,i3),i3=1,2)
!enddo
!do i2=1,nnatom
!print*,'nnat_bond_new',(nnat_bond_new(i2,i3),i3=1,2),nnatom
!enddo
!!print*,'sig_orb_1(n4)',(sig_orb_1(i2),i2=1,n4)
!do i2=1,n5
!print*,'tot_orb',tot_orb(i2),n5,nnatom
!enddo

call nnat_bond_sig(nnatom,sig_orb_1,nsig,nna)

!do i2=1,nna
!print*,'nnat_bond_new_new',(nnat_bond_new_1(i2,i3),i3=1,2),nna
!enddo

allocate(nn_group(natom,10))
allocate(sl_group(natom,10))
allocate(nelimt(natom))
allocate(full_nn_group(natom*10))

m18=0
m19=0
do i3=1,natom
  k=tot_orb(i3)
  m19=0
  do i2=1,nna
    do i1=1,2
      if(k.eq.nnat_bond_new_1(i2,i1))then
        if(i1.eq.1)i5=2
        if(i1.eq.2)i5=1
          m19=m19+1
          m18=m18+1
          nn_group(i3,m19)=nnat_bond_new_1(i2,i5)
          sl_group(i3,m19)=m18
      endif
    enddo
  enddo
nelimt(i3)=m19
enddo


!!! 'nn_group' specify orbitals in tot_orb array are connected with which other
!!! orbitals

!do i3=1,n5
!print*,'*nn_group*',nelimt(i3),(nn_group(i3,i4),i4=1,nelimt(i3))
!enddo
!do i3=1,n5
!print*,'*sl_group*',(sl_group(i3,i4),i4=1,nelimt(i3))
!enddo


m19=0
do i3=1,natom
  do i4=1,nelimt(i3)
    m19=m19+1
    full_nn_group(m19)=nn_group(i3,i4)
  enddo
enddo
fullgrp=m19

!write(*,*)'full_nn_g',m19,(full_nn_group(i3),i3=1,m19)
!endif
!!!!! Scoring of the structures has been started from here !!!
if (.not. allocated(loop_score_row))then
allocate(loop_score_row(20))
loop_score_row = 0
endif
if (.not. allocated(coordination_mat)) then
  allocate(coordination_mat((nae-nlast)/2, nao+natom))
  coordination_mat = 0
endif

do i=1,ncqs
  loopsc=0

  allocate(new_score(atom,2))
  new_score = 0

  do i3=1,atom
    new_score(i3,2)=score(i3,2)
  enddo

!  write(*,231)(str1(i,i1),i1=1,nae)
!  231 format(*(I0, 1x))

!! if the structures have active lone pairs and or radicals the associated atoms are being scored in 'new_score()'.
!! Scoring of the bonds of the structures started from here
!! the score of the bond is given by the score of the atoms involeved in the


  strscore=0.0
  ipnum=3
  numbond=0
  numbond1=0
  do i1=1+nl*2,nae-nlast,2
    tscore=0.0
    rscore=0.0
    n1=0
    n2=0
    n3=0
    n4=0
    do i2=i1,i1+1
      do i3=1,atom
        do i4=1,atn(active_atoms(i3))
          if(str1(i,i2).eq.atoset(active_atoms(i3),i4))then
            if(i2.eq.i1) then
              n1=active_atoms(i3)
            endif
            if(i2.eq.i1+1) then
              n2=active_atoms(i3)
            endif
          endif
        enddo
      enddo
    enddo
    n4=str1(i,i1)
    n5=str1(i,i1+1)
    do i2=1,nsym
      do i3=1,syn(i2)
        if(n4.eq.atsymset(i2,i3))nn4=i2
        if(n5.eq.atsymset(i2,i3))nn5=i2
      enddo
    enddo

  !! different type of scors given to different type of bonding
  !! three types of them are pi-pi, sigma-sigma, pi-sigma
  !! two types of calculations available 'loose' and 'tight'
  !! loose: consider pi(x) and pi(y) are same symmetry
  !! tight: consider pi(x) and pi(y) are different
  !! so they have different type of scoring

  symtype = 'loose' ! I am setting 'loose' as default as i think 'tight' option is not effective

    if(nsym.eq.1)then
      if(symtype.eq.'tight'.or.symtype.eq.'loose')then
        if(nn4.eq.1.and.nn5.eq.1)scr=1.0
      endif
    endif
    if(nsym.eq.2)then
      if(symtype.eq.'tight'.or.symtype.eq.'loose')then
        if(nn4.eq.2.and.nn5.eq.2)scr=1.0
        if(nn4.eq.1.and.nn5.eq.1)scr=2.0
        if(nn4.eq.1.and.nn5.eq.2)scr=4.0
        if(nn4.eq.2.and.nn5.eq.1)scr=4.0
        print*,'scr',scr
      endif
    endif
    if(nsym.eq.3)then
      if(symtype.eq.'tight')then
        if(nn4.eq.3.and.nn5.eq.3)scr=1.0
        if(nn4.eq.1.and.nn5.eq.1)scr=2.0
        if(nn4.eq.2.and.nn5.eq.2)scr=2.0
        if(nn4.eq.1.and.nn5.eq.2)scr=3.0
        if(nn4.eq.2.and.nn5.eq.1)scr=3.0
        if(nn4.eq.1.and.nn5.eq.3)scr=4.0
        if(nn4.eq.3.and.nn5.eq.1)scr=4.0
        if(nn4.eq.2.and.nn5.eq.3)scr=5.0
        if(nn4.eq.3.and.nn5.eq.2)scr=5.0
      endif
      if(symtype.eq.'loose')then
        if(nn4.eq.3.and.nn5.eq.3)scr=1.0
        if(nn4.eq.1.and.nn5.eq.1)scr=2.0
        if(nn4.eq.2.and.nn5.eq.2)scr=2.0
        if(nn4.eq.1.and.nn5.eq.2)scr=2.0
        if(nn4.eq.2.and.nn5.eq.1)scr=2.0
        if(nn4.eq.1.and.nn5.eq.3)scr=3.0
        if(nn4.eq.3.and.nn5.eq.1)scr=3.0
        if(nn4.eq.2.and.nn5.eq.3)scr=3.0
        if(nn4.eq.3.and.nn5.eq.2)scr=3.0
      endif
    endif
    tscore=new_score(n1,2)+new_score(n2,2)+scr/(new_score(n1,2)+new_score(n2,2))
    !print*,'ttttscore',tscore, new_score(n1,2), new_score(n2,2), scr
    nnscore=0.0
    nd=0
    do i2=1,nactorb
      if(str1(i,i1).eq.active_orbs(i2))then
        k1=atm_nb_orbs(i2)
      endif
      if(str1(i,i1+1).eq.active_orbs(i2))then
        k2=atm_nb_orbs(i2)
      endif
    enddo
    tscore=tscore+dist_mat(k1,k2)

    n1=0
    n2=0
    n3=0
    n4=0
    loop3:do i2=i1,i1+1
      do i3=1,nsig
        if(sig_orb_1(i3).eq.str1(i,i2))then
          n3=n3+1
          if(i2.eq.i1)n1=str1(i,i2)
          if(i2.eq.i1+1)n2=str1(i,i2)
          cycle loop3
        endif
      enddo
      do i3=1,atom
        do i4=1,atn(active_atoms(i3))
          if(str1(i,i2).eq.atoset(active_atoms(i3),i4))then
            if(i2.eq.i1) then
              n1=i3+nao+niao
            endif
            if(i2.eq.i1+1) then
              n2=i3+nao+niao
            endif
          endif
        enddo
      enddo
    enddo loop3


    if(n1.gt.n2)then
      k1=n2
      k2=n1
    else
      k1=n1
      k2=n2
    endif
    do i2=1,natom
      if(k1.eq.tot_orb(i2))l1=i2
    enddo
    do i2=1,nelimt(l1)
      if(k2.eq.nn_group(l1,i2))then
        n=1
        ncrow=1
        loop_score=1
        numbond=numbond+1
        loop_score_row(numbond)=loop_score+1
        coordination_mat(numbond,1)=k1
        coordination_mat(numbond,2)=k2
        goto 555
      endif
    enddo

    !allocate(connectivity_row(ncrow, ncrow))
    !! calculating connectivity score of each bonds with 'symm_loops'
    if(k1.ne.k2)call symm_loops(k1,k2,loop_score,ncrow)


    do i3=1,ncrow
      numbond=numbond+1
      loop_score_row(numbond)=loop_score+1
      do i2=1,loop_score+1
        coordination_mat(numbond,i2)=connectivity_row(i3,i2)
      enddo
    enddo

!    deallocate(connectivity_row)
    !print*,'i am here1'
    if(k1.eq.k2)loop_score=0
555 if(loop_score.ne.0)tscore=tscore+1.0/loop_score
    rscore=1.0/tscore
    !print*,'tscore',tscore

    if (.not. allocated(bndscore))then
      allocate(bndscore(nao))
      bndscore = 0.0
    endif
    do i2=1,ncrow
      numbond1=numbond1+1
      bndscore(numbond1)=rscore
    enddo

!    print*,'loop_score',loop_score

!    if (loop_score.gt.88)then
!    print*,'loop_score',loop_score
!    stop
!    endif
    strscore=strscore+rscore
    if(loop_score.ne.0)loopsc=loopsc+int(prime_num(loop_score))
    if(loop_score.eq.0)loopsc=loopsc+1000

  enddo
  
  deallocate(new_score)

!    print*,'i am here2'
!    print*,'numbond',numbond, bndscore
!do i1=1,numbond
!print*,'loop_score_row',loop_score_row(i1)
!print*,'coordination_mat',(coordination_mat(i1,i2),i2=1,loop_score_row(i1))
!enddo

  !! calculating coordinatiom score with coordination_val
  call coordination_val(coordination_mat,numbond,bndscore,coord_score)
  if(.not. allocated(stsymsc))then
  allocate(stsymsc(MaxStrOepo))
  stsymsc = 0
  endif
  stsymsc(i)=strscore+coord_score
  print*,'stsymsc',stsymsc(i)
  loopsymsc(i)=loopsc

!  write(15,*)stsymsc(i)
enddo

rewind(15)

deallocate(score)
!allocate(stsymsc1(ncqs))

do i=1,ncqs
  write(*,231)'structures','>',(str1(i,i1),i1=1,nae)
  write(*,*)stsymsc(i)
enddo
231 format(a,2x,a,1x,*(I0, 1x))
!909 format (F10.6)

!stop
deallocate(tot_orb)
deallocate(nn_group)
deallocate(sl_group)
deallocate(nelimt)
deallocate(full_nn_group)
500 allocate(ssym(ncqs))
ssym = 0

jj=1
loop4:do m19=1,ncqs
  if(m19.eq.1)ssym(1)=stsymsc(1)
    j=jj
    do i=1,j
      if(ssym(i).eq.stsymsc(m19)) then
        cycle loop4
      endif
    enddo
    jj=jj+1
    ssym(i)=stsymsc(m19)
    !print*,'ssym',ssym(i)
enddo loop4
nssym=jj
!print*,'i am here',ncqs, nssym
!stop

! sort out the structures according their scores
allocate(order(ncqs))
order = 0.0

do k3=1,nssym
  loop5:do k4=1,nssym
    do k5=1,k3
      if(order(k5).eq.ssym(k4)) then
        cycle loop5
      endif
    enddo
    if(ssym(k4).gt.order(k3))then
      order(k3)=ssym(k4)
      !print*,'order',order(k3)
    endif
  enddo loop5
enddo

!print*,'order',(order(k3),k3=1,nssym)
! scores converring into integers
allocate(lone_score(ncqs))
!allocate(symq(ncqs))

jj=0
do i=1,nssym
  jj=jj+1
  do j=1,ncqs
    if(order(i).eq.stsymsc(j))then
      if(nlast.ne.0)then
        losc=0.0
        do i1=nae-nlast,nae
          loop6:do i2=1,syn(1)
            if(str1(j,i1).eq.atsymset(1,i2)) then
              losc=losc+1.0/prime_num(str1(j,i1))
              exit loop6
            endif
          enddo loop6
        enddo
        lone_score(j)=losc
      endif
      symq(j)=jj
      !print*,'symq(j)',symq(j)
    endif
  enddo
enddo

if(.not. allocated(symm_groups))then
   allocate(symm_groups(ncqs, ncqs))
   symm_groups = 0
endif
if(.not. allocated(symm_groups_count)) then
   allocate(symm_groups_count(ncqs))
   symm_groups_count = 0
endif
jj = 0
do i = 1, ncqs
  if (jj.lt.symq(i)) jj = symq(i)
enddo
symm_maxval = jj

do i = 1, symm_maxval
   jj = 0
   do j = 1, ncqs
      if (i.eq.symq(j)) then
         jj = jj + 1
         symm_groups(i,jj)= j
      endif
   enddo
   symm_groups_count(i) = jj
enddo

do i = 1, symm_maxval
   print*,'symm_groups_count_sig',(symm_groups(i, j), j = 1, symm_groups_count(i))
enddo

!stop
!print*,'symq:sig',(symq(j),j=1,ncqs)
CALL SYSTEM ("rm symsig.temp")

deallocate(order)
deallocate(lone_score)
deallocate(ssym)

!deallocate(stsymsc(MaxStrOepo))
return
end subroutine symmetry_cal_sig
end module mod_symmetry_cal_sig
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
