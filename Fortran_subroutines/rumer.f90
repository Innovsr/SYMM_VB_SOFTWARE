!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rumer(permutation,n,j,setno,nl)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! All the Rumer sets are written in the file "Rumer_Sets.dat" varified
!! with the subroutine "Rumer_set_id" and then written in the file
!! "Rumer_Sets_all.dat" with full format of the output.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none
common/orb1/orbs1,rstr,nlonep

integer::permutation(100),n,i,i1,i2,j,orbs1(20),rstr(500,20),setno,m20,m19,rstr1(2000,20),nlonep
Double Precision::bond_rum(10000,100),bond_rum1(10000,100),k1,k2
integer::l,ll,jj,nl,rum_count(10000),set_num(100),Rid

print*,'enter rumer'
Rumwrite=0
do i=1,setno
  do i1=1,nlonep*2
    rstr1(i,i1)=rstr(i,i1)
  enddo
enddo

do i=1,setno
  do i1=1,nae
    do i2=1,n
      if(rstr(i,i1).eq.orbs1(i2))then
        rstr1(i,i1)=permutation(i2)
      endif
    enddo
  enddo
enddo


if(j.eq.1)then
  jj=0
  goto 103
endif

call Rumer_set_id(rstr1,setno,nlonep,Rid,set_num) 

if(Rid.eq.1) then
  return
endif



103 jj=jj+1

write(31,*)'set number',jj,(permutation(i),i=1,n),setno
write(23,*)'permutation',jj,'>',(permutation(i),i=1,n)
write(23,*)'*****************************************************'
write(23,*)
do m19=1,setno
write(31,914)(rstr1(m19,m20),m20=1,nae)
if(nfset.eq.5)write(10,915)1,':',niao,(rstr1(m19,m20),m20=1,nae)
if(niao.gt.1)then
endif
enddo
if(nfset.eq.5)write(10,*)'set number',jj

900 format(a,I3,10I5)

914 format(x,25I4)
915 format(x,I1,a,I3,x,25I4)
916 format(x,I3,I3,x,25I4)
print*,'exit rumer'
totrum=jj
!print*,'totrum',totrum
return
end subroutine rumer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
