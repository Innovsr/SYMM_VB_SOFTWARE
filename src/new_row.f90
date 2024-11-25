!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function new_row(l,n)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer::n,l,i,i1,i2,i3,loop_score_row(20),new_row_1(20),new_row
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

new_row=new_row_1(l)

return
end function new_row
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
