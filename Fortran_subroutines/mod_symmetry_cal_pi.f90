!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! this subroutine designed to identify symmetric structure groups for the systems !
!! where only pi type active orbitals exists, If a system posses both sigma and    !
!! pi orbitals then need symmetry_cal_sig sunroutine.                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module mod_symmetry_cal_pi
use commondat_mod
use commondat1_mod
use loops_mod
use mod_symm_loops
use mod_coordination_val 
implicit none

contains
subroutine symmetry_cal_pi(nl,str1,ncqs,symq,nssym)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer::m19,m18,i,i1,i2,i3,i4,i5,nl,ncqs,k,k1,k2,n1,n2,n3,n4,n5,n6,&
nd,k3,k4,k5,n,nssym,j,jj,kkk,n11,n12,n17,n18,n9,b3,b4
integer, allocatable::nn_count(:)
integer::nnscore,ipnum,b1,b2,l1
double precision::tscore,strscore,lpna
double precision, allocatable:: ssym(:),order(:)
integer::numbond,loop_score,ncrow, bonds
real*8::coord_score
real*8, allocatable::stsymsc(:),bndscore(:),score(:),new_score(:)
integer, allocatable::symq(:)
integer, allocatable::coordination_mat(:,:)
integer, pointer:: str1(:,:)
integer, allocatable::connectivity_row(:,:)
!real*8, pointer::stsymsc(:)

print*,'--- symmetry_cal_pi ---'
do i= 1, ncqs
print*,'str1',(str1(i, j),j=1,nae)
enddo
open(unit=15,file='sympi.temp',status='unknown')

allocate(loopsymsc(ncqs))
loopsymsc=0
if(.not. allocated(stsymsc))then
allocate(stsymsc(MaxStrOepo))
stsymsc = 0
endif
bonds=(nae-nl*2-nlast)/2
if (bonds.eq.0) goto 500

!!! initialise the score array, it stores the scores of each atoms
allocate(score(atom))
score = 0.0

!!! atom scoring starts here; scoring criteria: number of inactive bonds,
!!! number of inactive lone paires, number of charge if any associated with atoms

!! if any inactive biase (one or cluster of atoms) connected to any atom taken as a criteria
!  lpna=biasval(i3)+dist_nnat(i3)
!  endif

lpna = 0.0
do i3=1,atom
  !! if any inactive biase (one or cluster of atoms) connected to any atom it is taken as 'biasval'
  !! the covalent distances of an active atom with all others bonded active atoms taken as sum. 'dist_nnat'
  !! is the matrix of this distances.
  !! atomic number is a unique criteria of atoms
  score(i3)=at_num(i3)+biasval(i3)+dist_nnat(i3)
!  print*,i3,score(i3),lpna,at_num(i3)
enddo

allocate(tot_orb(nnnatom*2))

do i=1,nnnatom
print*,'nnat_bond',(nnat_bond(i,j),j=1,2)
enddo
!! preparation of connetivity matrix 'nnat_bond' starts 
n5=0
loop1:do i2=1,nnnatom
  loop2:do i1=1,2
    loop3:do i3=1,n5
      if(nnat_bond(i2,i1).eq.tot_orb(i3)) then
        cycle loop2
      endif
    enddo loop3
    n5=n5+1
    tot_orb(n5)=nnat_bond(i2,i1)
  enddo loop2
enddo loop1
natom=n5

allocate(nn_count(atom))
!print*,'active_atoms',(active_atoms(i),i=1,atom)
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
allocate(nn_group(natom,10))
allocate(sl_group(natom,10))
allocate(nelimt(natom))
allocate(full_nn_group(natom*10))

m18=0
m19=0
do i3=1,natom
 k=tot_orb(i3)
 m19=0
  do i2=1,nnnatom
   do i1=1,2
    if(k.eq.nnat_bond(i2,i1))then
      if(i1.eq.1)i5=2
      if(i1.eq.2)i5=1
      m19=m19+1
      m18=m18+1
      nn_group(i3,m19)=nnat_bond(i2,i5)
      sl_group(i3,m19)=m18
    endif
   enddo
  enddo
 nelimt(i3)=m19
enddo

!!! 'nn_group' specify orbitals in tot_orb array are connected with which other
!!! orbitals


m19=0
 do i3=1,natom
  do i4=1,nelimt(i3)
   m19=m19+1
   full_nn_group(m19)=nn_group(i3,i4)
  enddo
 enddo
 fullgrp=m19

do i1 =1,5
print*,'atoset',(atoset(i1,j),j=1,5)
enddo
do i1 =1,5
print*,'atsymset',(atsymset(i1,j),j=1,5)
enddo
print*,'full_nn_group',(full_nn_group(i1),i1=1,fullgrp)
do i3=1,n5
print*,'nn_group',nelimt(i3),(nn_group(i3, i4), i4 = 1, nelimt(i3))
enddo
do i3 = 1, n5
print*,'sl_group',(sl_group(i3, i4), i4 = 1, nelimt(i3))
enddo

!print*,'ncqs',ncqs,atom
do i=1,ncqs
  allocate(new_score(atom))
  do i3=1,atom
   new_score(i3)=score(i3)
   !print*,'score',i,i3,score(i3)
  enddo
  write(*,233),'structures',(str1(i,i3),i3=1,nae)
233 format(a,20I3)


  !! Scoring of the bonds of the structures started from here
  !! the score of the bond is given by the score of the atoms involeved in the bond
  
  strscore=0.0
  ipnum=0
  numbond=0
  
  !print*,'nl',nl,nae,nlast
  
  if (.not. allocated(loop_score_row))then
    allocate(loop_score_row(20))
    loop_score_row = 0
  endif
  
  do i1=1+nl*2,nae-nlast,2
    tscore=0.0
    n1=0
    n2=0
    n3=0
    n4=0
    n5=0
    n6=0
    n11=0
    n12=0
    n17=0
    n18=0
    n9=0
    kkk=0
      do i2=i1,i1+1
        b1=str1(i,i1)
        b2=str1(i,i1+1)
        do i3=1,atom
          do i4=1,atn(active_atoms(i3))
            if(str1(i,i2).eq.atoset(active_atoms(i3),i4))then
              if(i2.eq.i1) then
                n1=i3
                n11=active_atoms(i3)
                do k1=1,nsym
                  do k2=1,syn(k1)
                    if(atoset(active_atoms(i3),i4).eq.atsymset(k1,k2))then
                      n3=k1
                      n18=atsymset(k1,k2)
                    endif
                  enddo
                enddo
              endif
              if(i2.eq.i1+1) then
              n2=i3
              n12=active_atoms(i3)
                do k1=1,nsym
                  do k2=1,syn(k1)
                    if(atoset(active_atoms(i3),i4).eq.atsymset(k1,k2))then
                      n4=k1
                      n5=atsymset(k1,k2)
                    endif
                  enddo
                enddo
              endif
            endif
          enddo
        enddo
      enddo
    
    if(n3.eq.n4)then
      do k1=1,nsym
        do k2=1,syn(k1)
          if(n5.eq.atsymset(k1,k2))then
            n6=at_sym(k1)
          endif
        enddo
      enddo
    if(n6.eq.1)kkk=ipnum+2
    if(n6.eq.2)kkk=ipnum+1
      tscore=new_score(n1)+new_score(n2)+kkk/(new_score(n1)+new_score(n2))
      !print*,'tscore1',tscore
    endif
    
    
    if(n3.ne.n4)then
      do k1=1,nsym
        do k2=1,syn(k1)
          if(n5.eq.atsymset(k1,k2))then
            n6=at_sym(k1)
          endif
        enddo
      enddo
    do k1=1,nsym
      do k2=1,syn(k1)
        if(n18.eq.atsymset(k1,k2))then
          n9=at_sym(k1)
        endif
      enddo
    enddo
    
    if(symtype.eq.'loose')then
      if(n6.eq.1.and.n9.eq.1)kkk=ipnum+2
        if(n6.ne.n9)kkk=ipnum+3
      endif
    if(symtype.eq.'tight')then
      if(n6.eq.1.and.n9.eq.1)then
      kkk=ipnum+3
      else
      if(n3.eq.1.or.n4.eq.1)kkk=ipnum+4
      if(n3.eq.2.or.n4.eq.2)kkk=ipnum+5
      endif
    endif
    tscore=new_score(n1)+new_score(n2)+kkk/(new_score(n1)+new_score(n2))
    !print*,'tscore2',tscore
    endif
  
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
    print*,'ttscore3,b1,b2',tscore,b1,b2

    print*,'active_atoms',(active_atoms(i3),i3=1,atom)
    do i3=1,atom
      print*,'atoset',(atoset(active_atoms(i3),i4),i4=1,atn(active_atoms(i3)))
    enddo
    do i3=1,atom
     do i4=1,atn(active_atoms(i3))
      if(b1.eq.atoset(active_atoms(i3),i4))then
        b3=active_atoms(i3)
      endif
     enddo
    enddo
    
    do i3=1,atom
     do i4=1,atn(active_atoms(i3))
      if(b2.eq.atoset(active_atoms(i3),i4))then
       b4=active_atoms(i3)
      endif
     enddo
    enddo
    print*,'b3,b4,natom *****',b3,b4,natom 
    do i2=1,natom
      if(b3.eq.tot_orb(i2))l1=i2
    enddo
    print*,'l1',l1
      if (.not. allocated(bndscore))then
        allocate(bndscore(nao))
        bndscore = 0.0
      endif
    if (.not. allocated(coordination_mat)) then
      allocate(coordination_mat(2**((nae-nlast)/2), nao+natom))
      coordination_mat = 0
    endif
    do i2=1,nelimt(l1)
      if(b4.eq.nn_group(l1,i2))then
        print*,'sourav is here'
        loop_score=1
        numbond=numbond+1
        loop_score_row(numbond)=loop_score+1
        coordination_mat(numbond,1)=b3
        coordination_mat(numbond,2)=b4
        bndscore(numbond)=tscore
        goto 555
      endif
    enddo
    allocate(connectivity_row(nae, nae))
    
    print*,'i am here 0,b3,b4',b3,b4
    call symm_loops(b3,b4,loop_score,ncrow, connectivity_row)
    
    print*,'ncrow',ncrow
    do i3=1,ncrow
    print*,'loop_score 1',loop_score
    print*,'connectivity_row',(connectivity_row(i3,i2),i2=1,loop_score)
    enddo
    print*,'i am here 1'
    
    do i3=1,ncrow
      numbond=numbond+1
      loop_score_row(numbond)=loop_score+1
      do i2=1,loop_score+1
        coordination_mat(numbond,i2)=connectivity_row(i3,i2)
      enddo
      bndscore(numbond)=tscore
    enddo
    
    print*,'i am here 2'
    deallocate(connectivity_row)
555 strscore=strscore+tscore
  enddo
  !print*,'numbond',numbond
  !do i1=1,numbond
  !print*,'loop_score_row',loop_score_row(i1)
  !print*,'coordination_mat',(coordination_mat(i1,i2),i2=1,loop_score_row(i1))
  !enddo
  
  call coordination_val(coordination_mat,numbond,bndscore,coord_score)
  !stop
  
  print*,'sourav_pi_1',strscore
  stsymsc(i)=real(1.0/strscore+coord_score, kind=4)
!  write(15,*)stsymsc(i)
  deallocate(new_score)
enddo
deallocate(score)

    print*,'i am here 3'
!rewind(15)

!allocate(stsymsc1(ncqs))

!do i=1,ncqs
!  read(15,105)stsymsc1(i)
!enddo
!105 format(F8.4)

deallocate(tot_orb)
deallocate(nn_group)
deallocate(sl_group)
deallocate(nelimt)
deallocate(full_nn_group)
500 allocate(ssym(ncqs))
ssym = 0

  print*,'sourav_pi_2'
jj=1
loop4:do m19=1,ncqs
  if(m19.eq.1)ssym(1)=stsymsc(1)
    j=jj
    loop5:do i=1,j
      if(ssym(i).eq.stsymsc(m19)) then
        cycle loop4
      endif
    enddo loop5
    jj=jj+1
    ssym(i)=stsymsc(m19)
enddo loop4
nssym=jj

  print*,'sourav_pi_3'
allocate(order(ncqs))
order=0.0

!! arrange all structures in apropreate order
!! structures of same symmetry comes together

loop6:do k3=1,nssym
  loop7:do k4=1,nssym
    loop8:do k5=1,k3
      if(order(k5).eq.ssym(k4)) then
        cycle loop7
      endif
    enddo loop8
    if(ssym(k4).gt.order(k3))then
      order(k3)=ssym(k4)
    endif
  enddo loop7
enddo loop6

  print*,'sourav_pi_4'
deallocate(ssym)
!! each group get integer numbers serialy 
n=0
jj=0
do i=1,nssym
  jj=jj+1
  do j=1,ncqs
    if(order(i).eq.stsymsc(j))then
      symq(j)=jj
    endif
  enddo
enddo

  print*,'sourav_pi_5'
deallocate(order)

if(.not. allocated(symm_groups))then
   allocate(symm_groups(ncqs, ncqs))
   symm_groups = 0
endif
if(.not. allocated(symm_groups_count)) then
   allocate(symm_groups_count(ncqs))
   symm_groups_count = 0
endif
print*,'symm_groups',shape(symm_groups)
jj = 0
do i = 1, ncqs
  if (jj.lt.symq(i)) jj = symq(i)
enddo
symm_maxval = jj  !! symm_maxval = maximum number of symmetric groups

do i=1,ncqs
print*,'symq_pi',symq(i),symm_maxval
enddo


  print*,'sourav_pi_6'
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

  print*,'sourav_pi_7'
do i = 1, symm_maxval
   print*,'symm_groups_count_pi',symm_groups_count(i),(symm_groups(i, j), j = 1, symm_groups_count(i))
enddo

!do j=1,ncqs
!write(*,231)j,(str1(j,k),k=1,nae),symq(j)
!enddo

!231 format(20I3)

!CALL SYSTEM ("rm sympi.temp")
!deallocate(str1)

!deallocate(stsymsc)

return
end subroutine symmetry_cal_pi

end module mod_symmetry_cal_pi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
