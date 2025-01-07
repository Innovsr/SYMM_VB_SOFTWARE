!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine symm_break_factor(nl,str,tonstruc,str_quality_2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

integer::nl,i2,m8,l1,l2,l3,k1,k2,m13,m14,str(15000,20),str_quality_2(15000),tonstruc
integer::atsymset(20,20),nsym,syn(50),at_sym(50)

common/ats/atsymset,nsym,syn,at_sym

print*,'enter symm_break_factor'
do m8=1,15000
str_quality_2(m8)=1
enddo
do m8=1,tonstruc
!print*,'symm_break:sttr',(str(m8,i2),i2=1,nae)
l2=1
if(nsym.eq.1)goto 509
if(nsym.ne.1)l2=1+(nae-nlast-nl*2)/2
do k2=1+nl*2,nae-nlast,2
!print*,'k2',k2
do m13=1,nsym
l1=0
do m14=1,syn(m13)
if(syn(m13).eq.1)goto 507
do k1=k2,k2+1
if(str(m8,k1).eq.atsymset(m13,m14))then
l1=l1+1
endif
enddo
507 enddo
if(l1.eq.2) then
l2=l2-1 
goto 508
endif
enddo
508 enddo
str_quality_2(m8)=l2
509 enddo

!print*,'tonstruc',tonstruc
do m8=1,tonstruc
!!quality_fac(m8)=str_quality_1(m8)*str_quality_2(m8)-(str_quality_2(m8)-1)*(str_quality_1(m8)-1)
print*,'symm_break:sttr',(str(m8,i2),i2=1,nae)
print*,'qqqqqq',str_quality_2(m8)
enddo

print*,'exit symm_break_factor'
return
end subroutine symm_break_factor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
