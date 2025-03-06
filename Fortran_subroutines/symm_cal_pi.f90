!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! this subroutine designed to identify symmetric structure groups for the systems !
!! where only pi type active orbitals exists, If a system posses both sigma and    !
!! pi orbitals then need symmetry_cal_sig sunroutine.                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module symm_cal_pi
use commondat
use commondat1
use loops
implicit none
contains
subroutine symmetry_cal_pi(nl,str1,ncqs,stsymsc,symq,nssym)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!implicit none

!common/ats/atsymset,nsym,syn,at_sym
!common/loops/full_nn_group,fullgrp,natom,nelimt,sl_group,tot_orb,nn_group

integer::m19,m18,i,i1,i2,i3,i4,i5,nl,ncqs,k,k1,k2,n1,n2,n3,n4,n5,n6,&
nd,k3,k4,k5,n,nssym,j,jj,kkk,n11,n12,n17,n18,n9,b3,b4
integer::nn_count(50)
integer::nnscore,ipnum,connectivity_row(20,20),b1,b2,l1
double precision::tscore,strscore,score(100,2),new_score(100,2),lpna,&
ssym(15000),order(15000)
!integer::nn_group(50,10),nelimt(50),full_nn_group(1000),tot_orb(100),sl_group(50,10)
integer::symq(15000),numbond,loop_score,ncrow
real*8::coord_score,bndscore(20)
real::stsymsc1(15000)
!integer::atsymset(20,20),nsym,syn(50),at_sym(50),
integer::loop_score_row(20),coordination_mat(20,20)
integer, pointer:: str1(:,:)
double precision, pointer::stsymsc(:)

do i= 1, ncqs
print*,'str1',(str1(i, j),j=1,nae)
enddo
open(unit=15,file='sympi.temp',status='unknown')


!!! initialise the score array, it stores the scores of each atoms
do i=1,100
  do i1=1,2
    score(i,i1)=0.0
  enddo
enddo

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
        lpna=lpna+1.0/prime_num(1)
      endif
    enddo
    do i2=1,niach
      if(k.eq.iach(i2))then
        lpna=lpna+1.0/prime_num(3)
      endif
    enddo
  else

!! if any inactive biase (one or cluster of atoms) connected to any atom taken as a criteria
  lpna=biasval(i3)+dist_nnat(i3)
  endif

!! atomic number is a unique criteria of atoms
   score(i3,2)=lpna+at_num(i3)
   print*,i3,score(i3,2),lpna,at_num(i3)
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

print*,'active_atoms',(active_atoms(i),i=1,atom)
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

print*,'ncqs',ncqs,atom
do i=1,ncqs
  do i3=1,atom
   new_score(i3,2)=score(i3,2)
   print*,'score',i,i3,score(i3,2)
enddo

!! Scoring of the bonds of the structures started from here
!! the score of the bond is given by the score of the atoms involeved in the bond

strscore=0.0
ipnum=0
numbond=0

print*,'nl',nl,nae,nlast

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
    tscore=new_score(n1,2)+new_score(n2,2)+kkk/(new_score(n1,2)+new_score(n2,2))
    print*,'tscore1',tscore
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
  tscore=new_score(n1,2)+new_score(n2,2)+kkk/(new_score(n1,2)+new_score(n2,2))
  print*,'tscore2',tscore
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
  print*,'ttscore3',tscore
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
  
  
  do i2=1,natom
    if(b3.eq.tot_orb(i2))l1=i2
  enddo
  
  do i2=1,nelimt(l1)
    if(b4.eq.nn_group(l1,i2))then
      loop_score=1
      numbond=numbond+1
      loop_score_row(numbond)=loop_score+1
      coordination_mat(numbond,1)=b3
      coordination_mat(numbond,2)=b4
      bndscore(numbond)=tscore
      goto 555
    endif
  enddo
  
  call symm_loops(b3,b4,loop_score,connectivity_row,ncrow)
  
  !print*,'ncrow',ncrow
  !do i3=1,ncrow
  !print*,'loop_score 1',loop_score
  !print*,'connectivity_row',(connectivity_row(i3,i2),i2=1,loop_score)
  !enddo
  
  do i3=1,ncrow
    numbond=numbond+1
    loop_score_row(numbond)=loop_score+1
    do i2=1,loop_score+1
      coordination_mat(numbond,i2)=connectivity_row(i3,i2)
    enddo
    bndscore(numbond)=tscore
  enddo
  
  555 strscore=strscore+tscore
enddo

call coordination_val(coordination_mat,numbond,loop_score_row,bndscore,coord_score)

allocate(stsymsc(MaxStrOepo))

stsymsc(i)=1.0/strscore+coord_score
write(15,*)stsymsc(i)
enddo
rewind(15)

do i=1,ncqs
  read(15,105)stsymsc1(i)
enddo
105 format(F8.4)

jj=1
loop4:do m19=1,ncqs
  if(m19.eq.1)ssym(1)=stsymsc1(1)
    j=jj
    loop5:do i=1,j
      if(ssym(i).eq.stsymsc1(m19)) then
        cycle loop4
      endif
    enddo loop5
    jj=jj+1
    ssym(i)=stsymsc1(m19)
enddo loop4
nssym=jj

do i=1,10000
  order(i)=0.0
enddo

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

!! each group get integer numbers serialy 
n=0
jj=0
do i=1,nssym
  jj=jj+1
  do j=1,ncqs
    if(order(i).eq.stsymsc1(j))then
      symq(j)=jj
    endif
  enddo
enddo

!do j=1,ncqs
!write(*,231)j,(str1(j,k),k=1,nae),symq(j)
!enddo

!deallocate(stsymsc(MaxStrOepo))
!231 format(20I3)

!CALL SYSTEM ("rm sympi.temp")
!deallocate(str1)

deallocate(tot_orb)
deallocate(nn_group)
deallocate(sl_group)
deallocate(nelimt)
deallocate(full_nn_group)

return
end subroutine symmetry_cal_pi

end module symm_cal_pi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
