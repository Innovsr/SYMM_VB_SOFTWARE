!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! searching Rumer sets in the all possible VB set and mark the id
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Rumer_set_id(str2,i7,lonep,Rid,set_num)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

common/qfacRumid/qq1,bondq4,qq2,qq,rumset

integer::i,ii,iii,iiii,i1,i2,i3,i4,i5,i6,i7,j,Rid,orb1,orb2,sc,orbs(500),str2(2000,20),lonep,rumstrset(500,20),set_num(100)
integer::qq1(5000),qq2(5000),qq(5000),bondq4(5000),m119,m20,rumset
character(5)::a


!print*,'totrum',totrum,lonep,i7
j=0

Rid=0
rewind(31)
do i=1,totrum
  iiii=0
  read(31,*)
  do i1=1,i7
    read(31,*)(rumstrset(i1,i2),i2=1,nae)
  enddo
  loop3:do i1=1,i7
   do i2=1,i7
     iii=0
    loop2:do i3=lonep*2+1,nae-nlast,2
     do i4=lonep*2+1,nae-nlast,2
       ii=0
      loop1:do i5=i3,i3+1
        do i6=i4,i4+1
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
 if(iiii.eq.i7)then
   j=j+1
   Rid=1
   set_num(j)=i
 endif
enddo

nrs=j

if(Rid.eq.1.and.Rumwrite.eq.1)then
  rumset=rumset+1
  do m119=1,i7
  
    if(niao.eq.0)then
      write(23,900)qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',(str2(m119,m20),m20=1,nae)
    endif
    if(niao.gt.1)then
      write(23,901)qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,':',niao,(str2(m119,m20),m20=1,nae)
    endif
    if(niao.eq.1)then
      write(23,909)qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,1,(str2(m119,m20),m20=1,nae)
    endif
  enddo
  write(23,*)'Set_number=',rumset
  write(23,*)
endif


900 format(I3,2x,I3,2x,I3,2x,I3,x,a,x,25I4)
901 format(I3,2x,I3,2x,I3,2x,I3,x,a,x,I2,a,I2,x,25I4)
909 format(I3,2x,I3,2x,I3,2x,I3,x,a,x,I3,I3,x,25I4)

return
end subroutine Rumer_set_id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
