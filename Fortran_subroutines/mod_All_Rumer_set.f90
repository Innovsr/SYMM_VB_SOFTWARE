!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! generating all possible rumer sets by permuting the orbital numbers and 
! rejecting the duplicate sets
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module mod_All_Rumer_set
use commondat_mod
use mod_print_rumer
use orb_mod
use rum_rad_mod
use final_str_mod
implicit none

contains
subroutine All_Rumer_set(nstr,str,nlonep)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! N! number of Runmer set generates here 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer::j,i,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,n1,k,k1
integer::setno,nlonep,nstr
integer, allocatable::permutation(:), orbs2(:)
integer, pointer::str(:,:),rstr(:,:)

print*,'enter All_Rumer_set',nstr
open(unit=31,file='Rumer_Sets.dat',status='unknown')
allocate(rumerstr(CovDim, nae))
rumerstr = 0
rumset = 0

i7 = 0
do i1=1,nstr
   if (rumer(i1)*rumer_rad(i1).eq.1)then
      i7=i7+1
      do i2=1,nae
         rumerstr(i7,i2)=str(i1,i2)
      enddo
   endif
enddo
setno = i7

nlpr=nlonep
allocate(rstr(MaxStrOepo, nae))
!print*,'i am here1'

do i=1,setno
 ! print*,(rumerstr(i,i1),i1=1,nae)
  do i1=1,nae
    rstr(i,i1)=rumerstr(i,i1)
  enddo
enddo
allocate(orbs2(nae))
!print*,'i am here2'

i=0
do j=nlonep*2+1,nae
  i=i+1
  orbs2(i)=rumerstr(1,j)
enddo
n1=i
!print*,'orbs2',(orbs2(i),i=1,n1)

allocate(permutation(nao))

i=100
do j=1,n1
  if(orbs2(j).lt.i)i=orbs2(j)
enddo
!print*,'orbs2',i

!deallocate(orbs2)
if(.not. allocated(orbs1))then
allocate(orbs1(nae))
orbs1 = 0
endif

orbs1(1)=i
k=1
do j=1,20
  i=i+1
  do i1=1,n1
    if(i.eq.orbs2(i1))then
      k=k+1
      orbs1(k)=i
    endif
  enddo
enddo

!print*,'orbs1',orbs1

j=0
do i1=1,n1
  permutation(1)=orbs1(i1)
  loop1:do i2=1,n1
    if(permutation(1).eq.orbs1(i2)) then
      cycle loop1
    endif
      permutation(2)=orbs1(i2)
    if(n1.eq.2)then
      j=j+1
!      print*,'jjjj',j
      call print_rumer(permutation,n1,j,setno,rstr,str)
      cycle loop1
    endif

    loop2:do i3=1,n1
      do i=1,2
        if(permutation(i).eq.orbs1(i3)) then
          cycle loop2
        endif
      enddo
      permutation(3)=orbs1(i3)
      if(n1.eq.3)then
        j=j+1
!        print*,'jjjj',j
        call print_rumer(permutation,n1,j,setno,rstr,str)
        cycle loop2
      endif

      loop3:do i4=1,n1
        do i=1,3
          if(permutation(i).eq.orbs1(i4)) then
            cycle loop3
          endif
        enddo
        permutation(4)=orbs1(i4)
        if(n1.eq.4)then
          j=j+1
!          print*,'jjjj',j
          call print_rumer(permutation,n1,j,setno,rstr,str)
          cycle loop3
        endif

        loop4:do i5=1,n1
          do i=1,4
            if(permutation(i).eq.orbs1(i5)) then
              cycle loop4
            endif
          enddo
          permutation(5)=orbs1(i5)
          if(n1.eq.5)then
            j=j+1
!            print*,'jjjj',j,(permutation(k1),k1=1,n1)
            call print_rumer(permutation,n1,j,setno,rstr,str)
            cycle loop4
          endif

          loop5:do i6=1,n1
            do i=1,5
              if(permutation(i).eq.orbs1(i6)) then
                cycle loop5
              endif
            enddo
            permutation(6)=orbs1(i6)
            if(n1.eq.6)then
              j=j+1
 !             print*,'jjjj',j
              call print_rumer(permutation,n1,j,setno,rstr,str)
              cycle loop5
            endif

            loop6:do i7=1,n1
              do i=1,6
                if(permutation(i).eq.orbs1(i7)) then
                  cycle loop6
                endif
              enddo
              permutation(7)=orbs1(i7)
              if(n1.eq.7)then
                j=j+1
!                print*,'jjjj',j
                call print_rumer(permutation,n1,j,setno,rstr,str)
                cycle loop6
              endif

              loop7:do i8=1,n1
                do i=1,7
                  if(permutation(i).eq.orbs1(i8)) then
                    cycle loop7
                  endif
                enddo
                  permutation(8)=orbs1(i8)
                if(n1.eq.8)then
                  j=j+1
!                  print*,'jjjj',j
                  call print_rumer(permutation,n1,j,setno,rstr,str)
                  cycle loop7
                endif

                loop8:do i9=1,n1
                  do i=1,8
                    if(permutation(i).eq.orbs1(i9)) then
                      cycle loop8
                    endif
                  enddo
                  permutation(9)=orbs1(i9)
                  if(n1.eq.9)then
                    j=j+1
                    call print_rumer(permutation,n1,j,setno,rstr,str)
                    cycle loop8
                  endif

                  loop9:do i10=1,n1
                    do i=1,9
                      if(permutation(i).eq.orbs1(i10)) then
                        cycle loop9
                      endif
                    enddo
                    permutation(10)=orbs1(i10)
                    if(n1.eq.10)then
                      j=j+1
                      call print_rumer(permutation,n1,j,setno,rstr,str)
                      cycle loop9
                    endif

                    loop10:do i11=1,n1
                      do i=1,10
                        if(permutation(i).eq.orbs1(i11)) then
                          cycle loop10
                        endif
                      enddo
                        permutation(11)=orbs1(i11)
                      if(n1.eq.11)then
                        j=j+1
                        call print_rumer(permutation,n1,j,setno,rstr,str)
                        cycle loop10
                      endif
                    enddo loop10
                  enddo loop9
                enddo loop8
              enddo loop7
            enddo loop6
          enddo loop5
        enddo loop4
      enddo loop3
    enddo loop2
  enddo loop1
enddo

print*,'exit All_Rumer_set'
deallocate(f_str)
deallocate(rumerstr)
close (31)
CALL SYSTEM ("rm Rumer_Sets.dat")
return
end subroutine All_Rumer_set
end module mod_All_Rumer_set
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
