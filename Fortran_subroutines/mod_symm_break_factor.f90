!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module mod_symm_break_factor
use commondat_mod
implicit none

contains
subroutine symm_break_factor(nl,str,tonstruc,str_quality_2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer::nl,i2,m8,l1,l2,k1,k2,m13,m14,tonstruc
!integer::atsymset(20,20),nsym,syn(50),at_sym(50)
integer, pointer::str(:,:),str_quality_2(:)

!common/ats/atsymset,nsym,syn,at_sym

print*,'enter symm_break_factor'

str_quality_2=1
loop1:do m8=1,tonstruc
   l2=1
   if(nsym.eq.1) cycle loop1
   if(nsym.ne.1)l2=1+(nae-nlast-nl*2)/2
   loop2:do k2=1+nl*2,nae-nlast,2
      do m13=1,nsym
         l1=0
         loop3:do m14=1,syn(m13)
            if(syn(m13).eq.1) cycle loop3
            do k1=k2,k2+1
               if(str(m8,k1).eq.atsymset(m13,m14))then
                  l1=l1+1
               endif
            enddo
         enddo loop3
         if(l1.eq.2) then
            l2=l2-1 
            cycle loop2
         endif
      enddo
   enddo loop2
   str_quality_2(m8)=l2
enddo loop1

do m8=1,tonstruc
   print*,'symm_break:sttr',(str(m8,i2),i2=1,nae)
   print*,'qqqqqq',str_quality_2(m8)
enddo

print*,'exit symm_break_factor'
return
end subroutine symm_break_factor
end module mod_symm_break_factor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
