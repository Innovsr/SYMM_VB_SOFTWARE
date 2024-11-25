!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! It creates the matrix of the vectorial form of the set of structures and sending them !
!! for the checking of linearly independency among them.                                 !
!! the matrix composed of the overlap of the structures. It is a symmetric overlap matrix!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mat_ind(nl,numstr,totstr,strno,Ifail,det_inv)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use commondat
use commondat1
implicit none

common/fail/faiil
integer::nl,str1(10000,15),Ifail,totstr,tndet,vec(10000,1000),numstr,strno(1000),detmnt1(10000,15)&
,i,j,k,l,m,n,i1,comdet(10000),mat(10000,10000),sdet(10000),mat_sign(10000,10000),S00(10000,10000)&
,det_sign1(10000),faiil,det_inv
integer::col(1000)
DOUBLE PRECISION :: SPL0,SPL,D1(100)
character(len=3)::ind

tndet=0

do i=1,numstr
  j=strno(i)
  loop1:do l=1,totstr*ndet
    if(strdet(l).ne.j) then
      cycle loop1
    endif
    tndet=tndet+1
    do m=1,nae
      detmnt1(tndet,m)=detmnt(l,m)
    enddo
    det_sign1(tndet)=det_sign(l)
  enddo loop1
enddo

l=0
n=0
loop2:do i=1,tndet
  do i1=1,l
    if(comdet(i1).eq.i) then
      cycle loop2
    endif
  enddo
  n=n+1
  m=0
  loop3:do j=i,tndet
    do k=1,nae
      if(detmnt1(i,k).ne.detmnt1(j,k)) then
        cycle loop3
      endif
    enddo
    l=l+1
    m=m+1
    mat(n,m)=strdet(j)
    mat_sign(n,m)=det_sign1(j)
    comdet(l)=j
  enddo loop3
sdet(n)=m
enddo loop2

102 format(50I5)
103 format(a,2x,50I5)

! initialisation of the matrix
do i=1,totstr
  do j=1,totstr
    ind_mat(i,j)=0.0
  enddo
enddo

!creating diagonal of the symmetric matrix which always be number of the determinent
! in each structures
do i=1,totstr
  ind_mat(i,i)=ndet
enddo

! creating upper half of the symmetric matrix
do i=1,numstr
  do j=i+1,numstr
    do k=1,n
      do l=1,sdet(k)
        if(mat(k,l).eq.i)then
          do m=1,sdet(k)
            if(mat(k,m).eq.j)then
              ind_mat(i,j)=ind_mat(i,j)+mat_sign(k,m)*mat_sign(k,l)
            endif
          enddo
        endif
      enddo
    enddo
  enddo
enddo

! creating lower halp of the matrix by copying upper half
Do i=1,numstr
 Do j=i+1,numstr
  ind_mat(j,i)=ind_mat(i,j)
 Enddo
Enddo


call Invmat(numstr,Ifail)
111 faiil=Ifail
return
end subroutine mat_ind
