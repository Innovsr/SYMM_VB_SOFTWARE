!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine main_bond_str(nl,str1,ncqs,qual1,qual2,str2,q_fac)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

integer::jj,nl,str1(15000,15),str2(15000,15),ncqs,qual1(15000),qual2(15000),i,j,k,l &
,qul2(15000),qul1(15000),q_fac(15000),nqul,m19,k6,k7,k8,nqul1,nqul2,qt,qtg &
,qual3(15000),qual4(15000)

print*,'enter main_bond_str'
!do i=1,ncqs
!write(*,231)qual1(i),qual2(i),(str1(i,j),j=1,nae)
!enddo

jj=1
loop1:do m19=1,ncqs
  if(m19.eq.1)qul1(1)=qual1(1)
  j=jj
  do i=1,j
    if(qul1(i).eq.qual1(m19)) then
      cycle loop1
    endif
  enddo
  jj=jj+1
  qul1(i)=qual1(m19)
enddo loop1
nqul1=jj

jj=1
loop2:do m19=1,ncqs
  if(m19.eq.1)qul2(1)=qual2(1)
  j=jj
  do i=1,j
    if(qul2(i).eq.qual2(m19)) then
      cycle loop2
    endif
  enddo
  jj=jj+1
  qul2(i)=qual2(m19)
enddo loop2
nqul2=jj

do k6=1,15000
  qual3(k6)=0
  qual4(k6)=0
enddo

do k6=1,nqul1
  loop3:do k7=1,nqul1
    do k8=1,k6
      if(qual3(k8).eq.qul1(k7)) then
        cycle loop3
      endif
    enddo
    if(qul1(k7).gt.qual3(k6))then
      qual3(k6)=qul1(k7)
    endif
  enddo loop3
enddo

do k6=1,nqul2
  loop4:do k7=1,nqul2
    do k8=1,k6
      if(qual4(k8).eq.qul2(k7)) then
        cycle loop4
      endif
    enddo
    if(qul2(k7).gt.qual4(k6))then
      qual4(k6)=qul2(k7)
    endif
  enddo loop4
enddo

do k6=1,15000
  qul1(k6)=0
  qul2(k6)=0
enddo

k7=0
do k6=nqul1,1,-1
  k7=k7+1
  qul1(k7)=qual3(k6)
enddo

k7=0
do k6=nqul2,1,-1
  k7=k7+1
  qul2(k7)=qual4(k6)
enddo

k=0
qt=0
qtg=1
do m19=1,nqul1
  do j=1,nqul2
    if(qtg.ne.k)qt=qt+1
    qtg=k
    loop5:do i=1,ncqs
      if(qul1(m19).ne.qual1(i)) then
        cycle loop5
      endif
      if(qual2(i).eq.qul2(j))then
        k=k+1
        do l=1,nae
          str2(k,l)=str1(i,l)
        enddo
        q_fac(k)=qt
      endif
    enddo loop5
  enddo
enddo

!do i=1,nqs
!write(*,231)(str2(i,j),j=1,nae)
!enddo

231 format(20I3)
print*,'exit main_bond_str'

return
end subroutine main_bond_str
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
