!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! searching Rumer sets in the all possible VB set and mark the id
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module mod_Rumer_set_id
use commondat_mod
use mod_quality_factor
use quality_mod
implicit none

contains
subroutine Rumer_set_id(str2,i7,lonep,Rid,str)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer::i,ii,iii,iiii,i1,i2,i3,i4,i5,i6,i7,j,Rid,lonep,rumstrset(500,20)
integer::m119,m20, total_str
integer, pointer::str2(:,:),str(:,:)


print*,'enter Rumer_set_id', totrum,i7 
if (.not. allocated(set_num))then
  allocate(set_num(1000))
endif
!print*,'totrum',totrum,lonep,i7
j=0

Rid=0
rewind(31)
do i=1,totrum
  iiii=0
  read(31,*)
  do i1=1,i7
    read(31,*)(rumstrset(i1,i2),i2=1,nae)
!    print*,'rumstrset',(rumstrset(i1,i2),i2=1,nae)
  enddo
  loop3:do i1=1,i7
    do i2=1,i7
      iii=0
      loop2:do i3=lonep*2+1,nae-nlast,2
        do i4=lonep*2+1,nae-nlast,2
          ii=0
          loop1:do i5=i3,i3+1
            do i6=i4,i4+1
              !print*,'rumstrset(i2,i5).eq.str2(i1,i6)',rumstrset(i2,i5),str2(i1,i6)
              if(rumstrset(i2,i5).eq.str2(i1,i6))then
                ii=ii+1
                cycle loop1
              endif
            enddo
          enddo loop1
          if(ii.eq.2)then
            iii=iii+ii
            cycle loop2
          endif
        enddo
      enddo loop2
      if(iii.eq.nae-lonep*2-nlast)then
        iiii=iiii+1
        cycle loop3
      endif
    enddo
  enddo loop3
  !print*,'iiii,i7',iiii,i7
  if(iiii.eq.i7)then
    j=j+1
    Rid=1
    set_num(j)=i
  endif
enddo

nrs=j

!if(Rid.eq.1.and.Rumwrite.eq.1)then
!if(Rid.eq.0)then
!!  call quality_factor(lonep, str2, i7, qq, qq1, qq2, bondq4, mbondq, pref_radical)
!  total_str=MaxStrOepo
!!  call str_sl_check(lonep,str2,i7,str,total_str)
!  !if(symtype.eq.'check')then
!    !call sym_check(nl,str2,n,sym_str_sl,sym_str_num,nssym1)
!  !endif
!  rumset=rumset+1
!  do m119=1,i7
!  
!    if(niao.eq.0)then
!      write(10,900)str_quality_1(m119),bondq(m119),str_quality_2(m119),quality_fac(m119),'|',(str2(m119,m20),m20=1,nae)
!    endif
!    if(niao.gt.1)then
!      write(10,901)str_quality_1(m119),bondq(m119),str_quality_2(m119),quality_fac(m119),'|',1,':',niao,(str2(m119,m20),m20=1,nae)
!    endif
!    if(niao.eq.1)then
!      write(10,909)str_quality_1(m119),bondq(m119),str_quality_2(m119),quality_fac(m119),'|',1,1,(str2(m119,m20),m20=1,nae)
!    endif
!  enddo
!  write(10,*)'Set_number=',rumset
!  write(10,*)
!endif
!
!
!900 format(3x,I0,3x,I0,3x,I0,3x,I0,1x,a,1x,*(I0, 1x))
!901 format(3x,I0,3x,I0,3x,I0,3x,I0,1x,a,1x,I0,a,I0,1x,*(I0, 1x))
!909 format(3x,I0,3x,I0,3x,I0,3x,I0,1x,a,1x,I0,I0,1x,*(I0, 1x))

return
end subroutine Rumer_set_id

end module mod_Rumer_set_id
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
