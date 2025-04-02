!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! calling different quality calculations as opted and marging them into one score 
!! for each structue. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module mod_quality_factor
use commondat_mod
use mod_intra_bond_factor
use mod_symm_break_factor
use mod_main_bond_cal
use mod_prio_rad_str
use mod_nnat_bond_cal_2
use mod_nnat_bond_cal
implicit none

contains
subroutine quality_factor(nl, str1, ncqs, quality_fac, str_quality_1, str_quality_2, bondq, mbondq, pref_radical)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!common/ats/atsymset,nsym,syn,at_sym

integer::i,j,nl,ncqs,noqt, qt_max
!integer::atsymset(20,20),nsym,syn(50),at_sym(50)
integer, pointer :: str1(:,:), str_quality_1(:), str_quality_2(:)
integer, pointer :: bondq(:), mbondq(:), pref_radical(:), quality_fac(:)
integer, allocatable:: qual_mat(:,:)
!integer, allocatable:: quality_fac(:)


print*,'itb,syb,nnb,radical,mnbond',itb,syb,nnb,radical,mnbond
print*,'flg1`, nsym, input_flg, prad, imbd, nlast, ncqs',flg1, nsym, input_flg, prad, imbd, nlast, ncqs
do i=1,ncqs
print*,(str1(i,j),j=1,nae)
enddo

   allocate(str_quality_1(ncqs))
   allocate(str_quality_2(ncqs))
   allocate(bondq(ncqs))
   allocate(pref_radical(ncqs))
   allocate(mbondq(ncqs))
!   allocate(quality_fac(MaxStrOepo))
if (.not. associated(quality_fac)) then
!if (.not. allocated(quality_fac)) then
  allocate(quality_fac(ncqs))
  quality_fac = 0
endif
print*,'q_fac',(quality_fac(i),i=1, ncqs)

if(qflg.eq.1)then 
!  do j=1,ncqs
    quality_fac = 1
!  enddo
  return
endif


str_quality_1=0
print*,'i am here'
str_quality_2=0
bondq=0
pref_radical=0
mbondq=0

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

allocate(qual_mat(6, ncqs))

qual_mat = 0
!do j=1,ncqs
!  do i=1,6
!  qual_mat(i,j)=0
!  enddo
!enddo

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
  quality_fac(j)=qual_mat(noqt,j)
  write(*,200)j,(str1(j,i),i=1,nae), quality_fac(j)
enddo

deallocate(qual_mat)
200 format(20I3)
 return
end subroutine quality_factor
end module mod_quality_factor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
