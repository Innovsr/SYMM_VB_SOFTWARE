!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine wigner(nnae,wig2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This subroutine calculates the number of permisible structures by wigner's
!!! theorem
use commondat
implicit none

integer::wig2,nm,s,i,j,k,l,b,m,nnae
real*8::a,factorial,ss


j=0
m=mod(nnae,2)
if(m.eq.1)then
ss=0.5
a=2.0*ss+1.0
do i=1,nnae,2
j=j+1
enddo
k=j+1
l=j-1
wig2=a*factorial(nnae)/(factorial(k)*factorial(l))
endif

if(m.eq.0)then
if(mult.eq.1)s=0
if(mult.eq.3)s=1
b=2*s+1
k=(nnae/2)+1+s
l=(nnae/2)-s

wig2=b*factorial(nnae)/(factorial(k)*factorial(l))
endif


return
end subroutine wigner
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
