!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module mod_nnat_bond_cal_2
use commondat_mod
implicit none

contains
subroutine nnat_bond_cal_2(nl,str1,ncqs,bondq)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This subroutine calculates the nearest naighbour score of the structures

integer::i, i1, i3, i4, i5, i7, iii, iiii, nl, ncqs
real*8::least
integer, allocatable:: sl(:)
real, allocatable:: bondq1_dist(:)
double precision, allocatable::bondq_dist(:)
double precision :: ii
integer, pointer::str1(:,:),bondq(:), nn(:)

print*,'enter nnat_bond_cal_2'

open(unit=13,file='nnbd.temp',status='unknown')

! initialisation
!allocate(bondq_dist(MaxStrOepo))
!allocate(bondq1_dist(MaxStrOepo))
!allocate(sl(MaxStrOepo))
allocate(bondq_dist(ncqs))
allocate(bondq1_dist(ncqs))
allocate(sl(ncqs))
allocate(nn(2))
bondq=0
sl = 0
bondq_dist = 0.0
!bondq1_dist = 0.0

!231 format (30I3)
i4=0
do i1=1,ncqs
  ii=0.0
  loop3:do i3=1+nl*2,(nae-nlast),2
    iii=0
    nn(1)=0
    nn(2)=0
    loop1:do i5=i3,i3+1
      do i4=1,atom
        do i7=1,atn(active_atoms(i4))
          if(str1(i1,i5).eq.atoset(active_atoms(i4),i7))then
            iii=iii+1
            if(mod(i5,2).eq.0) nn(2)=active_atoms(i4)
            if(mod(i5,2).eq.1) nn(1)=active_atoms(i4)
            cycle loop1
          endif
        enddo

      enddo
    enddo loop1
    if(iii.ne.2) cycle loop3
    if(nn(2).eq.nn(1))then

!!! changes done for EDEN
iab_length=1.00
!!!!!!!!!!!!!!!!!!!!!!!!!

      if (iab_length.eq.0.0)ii=ii+100.0
      if (iab_length.ne.0.0)ii=ii+iab_length
      cycle loop3
    endif
    ii=ii+dist_act_rel_mat(nn(1),nn(2))
    print*,'dist_act_rel_mat',dist_act_rel_mat(nn(1),nn(2))
  enddo loop3
  bondq_dist(i1)=ii
  write(13,*)bondq_dist(i1)
  write(*,*)bondq_dist(i1)
enddo
rewind(13)

do i1=1,ncqs
   read(13,105)bondq1_dist(i1)
enddo

105 format (F8.4)
iii=0
iiii=0
do
  least=10000.0
  loop4:do i=1,ncqs
    do i1=1,iii
      if(i.eq.sl(i1)) cycle loop4
    enddo
    if(bondq1_dist(i).le.least)least=bondq1_dist(i)
  enddo loop4
  if(least.ne.100.0.or.least.ne.0.0)iiii=iiii+1
  do i=1,ncqs
    if(bondq1_dist(i).eq.least)then
      iii=iii+1
      bondq(i)=iiii
      sl(iii)=i
    endif
  enddo
  if(iii.eq.ncqs) exit
enddo

do i=1,ncqs
write(*,231)i,(str1(i,i4),i4=1,nae),bondq(i)
!print*,'bondq(i)',i,bondq(i)
enddo

231 format(*(I0, 1x))

CALL SYSTEM ("rm nnbd.temp")
print*,'exit nnat_bond_cal_2'
!stop
deallocate(bondq_dist)
deallocate(bondq1_dist)
return
end subroutine nnat_bond_cal_2
end module mod_nnat_bond_cal_2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
