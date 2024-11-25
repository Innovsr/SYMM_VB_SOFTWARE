!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! this subroutine designed to identify symmetric structure groups for the systems !
!! where only pi type active orbitals exists, If a system posses both sigma and    !
!! pi orbitals then need symmetry_cal_sig sunroutine.                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine symmetry_cal_pi(nl,str1,ncqs,stsymsc,symq,nssym)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
use commondat1
implicit none

common/ats/atsymset,nsym,syn,at_sym
common/loops/full_nn_group,fullgrp,natom,nelimt,sl_group,tot_orb,nn_group

integer::m19,m18,i,i1,i2,i3,i4,i5,i6,i7,i8,ii,iii,iiii,nl,ncqs,k,k1,k2,n1,n2,n3,n4,n5,n6,l6,&
nd,fullgrp,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17,k18,k19,k20,k21,k22,k23,k24,&
k25,k26,k27,n,j1,j2,j3,j4,j5,j6,j7,j8,j9,j10,nssym,j,jj,nn1,nn2,kkk,n11,n12,n17,n18,n9,l,b3,b4
integer::nnat_bond_new(100,2),tot_orb(100),nn_count(50),sig_orb(100),sig_orb_1(100),nsig,nnatom,natom
integer::kk,nscore(100),nnscore,ipnum,connectivity_row(20,20),b1,b2,l1,l2,l3,l4,nnat_bond_1(100,2)
double precision::tscore,stsymsc(15000),strscore,bd_dist,score(100,2),new_score(100,2),lpna,&
ssym(15000),order(15000)
integer::nn(10),str1(15000,20),mbondq(15000),nn_group(50,10),nelimt(50),full_nn_group(1000),&
sl_group(50,10),symq(15000),numbond,loop_score,nrow,new_mat(10,10),ncrow

integer::atsymset(20,20),nsym,syn(50),at_sym(50),loop_score_row(20),coordination_mat(20,20)
real*8::atoset_symscr(1000),piscr,losc,lone_score(15000),symq1(15000),cot(1000),sig_sym_sc(15000)&
,sym_sets(15000),n7,n8,coord_score,bndscore(20)
real::stsymsc1(15000)

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


m18=0
m19=0
do i3=1,n5
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

natom=n5

m19=0
 do i3=1,n5
  do i4=1,nelimt(i3)
   m19=m19+1
   full_nn_group(m19)=nn_group(i3,i4)
  enddo
 enddo
fullgrp=m19


do i=1,ncqs
  do i3=1,atom
   new_score(i3,2)=score(i3,2)
enddo

!! Scoring of the bonds of the structures started from here
!! the score of the bond is given by the score of the atoms involeved in the bond

strscore=0.0
ipnum=0
numbond=0

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

231 format(20I3)

CALL SYSTEM ("rm sympi.temp")

200 return
end subroutine symmetry_cal_pi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
