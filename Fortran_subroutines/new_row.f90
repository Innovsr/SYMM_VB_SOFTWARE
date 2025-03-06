!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This function replace the bond connectivity sequence to a palindromic 
! sequence for the calculation of connectivity score
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function new_row(l,n) result(val)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer::n,l,i,i1,loop_score_row(20),new_row_1(20), val
common/lp/loop_score_row

k=int(loop_score_row(n)/2)+mod(loop_score_row(n),2)
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

val=new_row_1(l)
print*,'coordination_val', val

return
end function new_row
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
