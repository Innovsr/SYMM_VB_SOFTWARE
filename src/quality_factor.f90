!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! calling different quality calculations as opted and marging them into one score 
!! for each structue. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine quality_factor(nl,str1,ncqs,q_fac,str_quality_1,str_quality_2,bondq)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

common/ats/atsymset,nsym,syn,at_sym

integer::i,j,nl,qt,qt_cnt(10),ncqs,i9,m16,m17,m18,m19,qtg,str1(15000,20),equality,noqt&
,jj,nqul1,nqul2,i4,q_fac(15000),k6,k7,k8,str_quality_1(15000),str_quality_2(15000)&
,bondq(15000),qual_mat(6,15000),q,pref_radical(15000),qt_max,mbondq(15000)
integer::atsymset(20,20),nsym,syn(50),at_sym(50),fchrgq(15000)


print*,'itb,syb,nnb,radical,mnbond',itb,syb,nnb,radical,mnbond

if(qflg.eq.1)then 
  do j=1,ncqs
    q_fac(j)=1
  enddo
  return
endif

if(itb.ne.1.or.flg1.eq.1) call intra_bond_factor(nl,str1,ncqs,str_quality_1)
if(syb.ne.1.and.nsym.ne.1) call symm_break_factor(nl,str1,ncqs,str_quality_2)
if(input_flg.eq.1)then
  if(nnb.ne.1.or.flg1.eq.1) call nnat_bond_cal(nl,str1,ncqs,bondq)
endif
if(input_flg.eq.0)then
  if(nnb.ne.1.or.flg1.eq.1) call nnat_bond_cal_2(nl,str1,ncqs,bondq)
endif
if(radical.ne.1.and.prad.ne.0.and.nlast.ne.0.or.flg1.eq.1) call prio_rad_str(nl,str1,ncqs,pref_radical)
if(mnbond.ne.1.and.imbd.ne.0.or.flg1.eq.1) call main_bond_cal(nl,str1,ncqs,mbondq)

do j=1,ncqs
  do i=1,6
  qual_mat(i,j)=0
  enddo
enddo

do j=1,ncqs
i=itb
 qual_mat(i,j)=qual_mat(i,j)+str_quality_1(j)

i=syb
 qual_mat(i,j)=qual_mat(i,j)+str_quality_2(j)

i=nnb
 qual_mat(i,j)=qual_mat(i,j)+bondq(j)

i=radical
 qual_mat(i,j)=qual_mat(i,j)+pref_radical(j)

i=mnbond
 qual_mat(i,j)=qual_mat(i,j)+mbondq(j)
enddo

!213 format(10I5)

noqt=1
do i=2,6
  if(qual_mat(i,1).gt.0)then
  noqt=noqt+1
  endif
enddo

if(noqt.gt.2)then
  do i=2,noqt-1
  qt_max=1
    do j=1,ncqs
      if(qt_max.lt.qual_mat(i+1,j))qt_max=qual_mat(i+1,j)
    enddo
    do j=1,ncqs
      qual_mat(i+1,j)=qt_max*(qual_mat(i,j)-1)+qual_mat(i+1,j)
    enddo
  enddo
endif

do j=1,ncqs
  q_fac(j)=qual_mat(noqt,j)
enddo

!200 format(20I3)
 return
end subroutine quality_factor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
