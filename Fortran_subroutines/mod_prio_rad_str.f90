!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module mod_prio_rad_str
use commondat_mod
implicit none

contains
subroutine prio_rad_str(nl,str1,ncqs,pref_radical)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! this subroutine calculate the scores for priority radical provided by the users
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer::i,j,i3,i4,i5,i7,ii,nl,ncqs,jj
integer, pointer:: str1(:,:), pref_radical(:)

!!!!! starting the arrangement of the structures according to priority radicles!!!


print*,'enter prio_rad_str'

! initialise the priority scores
pref_radical=nlast+1

if(nlpset.ne.0) then
  loop1:do i5=1,nlpset
    ii=0
    do i=1,nl*2,2
      do j=1,nl
        if(plpair(i5,j).eq.str1(1,i))then
          ii=ii+1
        endif
      enddo
    enddo
    if(ii.ne.nl) cycle loop1
    if(ii.eq.nl)then
      jj=i5
      if(plpair(i5,nl+1).ne.0.) goto 100
      if(plpair(i5,nl+1).eq.0.) goto 102
    endif
  enddo loop1
  return
endif


102 do i=1,ncqs
  do j=1,prad
    ii=0
    do i3=1,norad(j)
      do i4=nae-nlast+1,nae
        if(str1(i,i4).eq.prio_rad(j,i3))then
          ii=ii+1
        endif
      enddo
    enddo
    if(ii.eq.norad(j))then
      pref_radical(i)=pref_radical(i)-1
    endif
  enddo
enddo

if(nlpset.eq.0) then
   return
endif

100 if(plpair(jj,nl+1).eq.0) return


do i=1,ncqs
  do j=nl+1,lp(jj)
    i7=plpair(jj,j)
    ii=0
    do i3=1,norad(i7)
      do i4=nae-nlast+1,nae
        if(str1(i,i4).eq.prio_rad(i7,i3))then
          ii=ii+1
        endif
      enddo
    enddo
    if(ii.eq.norad(i7))then
      pref_radical(i)=pref_radical(i)-1
    endif
  enddo
enddo

!231 format(30I3)

print*,'exit prio_rad_str'

return
end subroutine prio_rad_str
end module mod_prio_rad_str
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
