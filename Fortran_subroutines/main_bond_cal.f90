!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! It counts how many preferred bond present in a structures and score them accordingly!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine main_bond_cal(nl,str1,ncqs,mbondq)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

integer::m19,m18,i1,i2,i3,i4,i5,i6,i7,i8,i9,ii,iii,iiii,nl,ncqs
integer::nn(10),str1(15000,20),mbondq(15000)

print*,'enter main_bond_cal'

!! initialise main or preferred bond quality score : total bond + 1
do i1=1,15000
  mbondq(i1)=1+(nae-nl*2-nlast)/2
enddo


do i1=1,ncqs
  i9=0
  do i3=1,100
    bond_count(i1,i3)=0
  enddo
  iiii=1+(nae-nl*2-nlast)/2
  do i3=1+nl*2,(nae-nlast),2
    i9=0
    nn(1)=0
    nn(2)=0
    loop3:do i4=1,nmbond*2,2
      i9=i9+1
      iii=0
      do i5=i3,i3+1
        do i7=i4,i4+1
          if(str1(i1,i5).eq.main_bond(i7))then
            iii=iii+1
            if(mod(i5,2).eq.0) nn(2)=i5
            if(mod(i5,2).eq.1) nn(1)=i5
          endif
        enddo
      enddo
      if(iii.ne.2) cycle loop3
      iiii=iiii-1
      bond_count(i1,i9)=bond_count(i1,i9)+1
    enddo loop3
  enddo
  mbondq(i1)=iiii
enddo

!do i1=1,ncqs
!write(*,231)(str1(i1,i2),i2=1,nae),mbondq(i1)
!write(*,231)(bond_count(i1,i2),i2=1,nmbond)
!enddo

231 format (30I3)

!print*,'exit main_bond_cal'
return
end subroutine main_bond_cal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
