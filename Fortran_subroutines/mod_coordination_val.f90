!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculationg coordination among bonds presents in a structure 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module mod_coordination_val
use commondat_mod
implicit none


contains
subroutine coordination_val(coordination_mat,numbond,bndscore,coord_score)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer::i,i1,i2,i3,numbond
integer, intent(in)::coordination_mat(:,:)
real*8::coord_score,c_score,score,val1, val2
real*8, intent(in)::bndscore(:)

print*,'enter coordination_val '!, (loop_score_row(i),i=1,numbond),coordination_mat

coord_score=0.0
do i=1,numbond
  if (coordination_mat(i,1).ne.0)then
    c_score=0.0
    !print*,'score_coordination_2',c_score
    do i1=i+1,numbond
      if (i.ne.i1)then
        if (coordination_mat(i1,1).ne.0)then
          score=0.0
          !print*,'score_coordination_2',score
          do i2=1,loop_score_row(i)
            do i3=1,loop_score_row(i1)
            !print*,'coordination_mat',coordination_mat(i,i2),coordination_mat(i1,i3),i, i1, i2, i3
              if (coordination_mat(i,i2).eq.coordination_mat(i1,i3))then
                call new_row(i2, i, val1)
                call new_row(i3, i1, val2)
                score=score+((bndscore(i)/val1) + (bndscore(i1)/val2))
!                print*,'score',bndscore(i),val1,bndscore(i1),val2,score
              endif
            enddo
          enddo
          if(score.ne.0.0)c_score=c_score+1.0/score
        endif
      endif
    enddo
    if(c_score.ne.0.0)coord_score=coord_score+c_score
  endif
enddo

print*,'exit coordination value',coord_score
return
end subroutine coordination_val
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This function replace the bond connectivity sequence to a palindromic 
! sequence for the calculation of connectivity score
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine new_row(l,n,val)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer, intent(in):: l, n
integer::i,i1,new_row_1(20), k
real*8::val

val = 0.0
k=int(loop_score_row(n)/2)+mod(loop_score_row(n),2)
!print*,'kkk',k
if (mod(loop_score_row(n),2).eq.0)then
  i1=0
  do i=1,loop_score_row(n)/2
    i1=i1+1
    new_row_1(i1)=k-ABS(k-i)
  enddo

  do i=loop_score_row(n)/2,1,-1
    i1=i1+1
    new_row_1(i1)=k-ABS(k-i)
  enddo

endif

if (mod(loop_score_row(n),2).eq.1)then
  do i=1,loop_score_row(n)
    new_row_1(i)=k-abs(k-i)
  enddo
endif

val = real(new_row_1(l))
!print*,'coordination_val', l, val

end subroutine new_row
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module mod_coordination_val
