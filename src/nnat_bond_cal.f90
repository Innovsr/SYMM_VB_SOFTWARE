!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine nnat_bond_cal(nl,str1,ncqs,bondq)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! This subroutine calculates the number of nearest neighbour bonds !!!!!

use commondat
implicit none

integer::m19,m18,i1,i2,i3,i4,i5,i6,i7,ii,iii,iiii,nl,ncqs
integer::nn(10),str1(15000,20),bondq(15000)

print*,'enter nnat_bond_cal'

! initialisation of scores
do i1=1,15000
  bondq(i1)=0
enddo

if(nnnatom.eq.0) then
  return
endif

231 format (30I3)
i4=0
do i1=1,ncqs
  iiii=1+(nae-nl*2-nlast)/2
  loop2:do i3=1+nl*2,(nae-nlast),2
    iii=0
    nn(1)=0
    nn(2)=0
    loop1:do i5=i3,i3+1
      do i4=1,atom
        do i7=1,atn(active_atoms(i4))
          if(str1(i1,i5).eq.atoset(active_atoms(i4),i7))then
            iii=iii+1
            if(mod(i5,2).eq.0) nn(2)=i4
            if(mod(i5,2).eq.1) nn(1)=i4
            cycle loop1
          endif
        enddo
      enddo
    enddo loop1
    if(nn(2).eq.nn(1)) cycle loop2
    if(iii.ne.2) cycle loop2
    do i6=1,nnnatom
      ii=0
      do i7=1,2
        do i4=1,2
          if(nn(i4).eq.nnat_bond(i6,i7))then
            ii=ii+1
          endif
        enddo
      enddo
      if(ii.eq.2)then
        iiii=iiii-1
        cycle loop2
      endif
    enddo

  enddo loop2
  bondq(i1)=iiii
enddo

print*,'exit nnat_bond_cal'
 return
end subroutine nnat_bond_cal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
