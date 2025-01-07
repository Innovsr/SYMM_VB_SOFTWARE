!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine checking each bond of every stracture, finding if both of the orbitals associated with that 
!! bond are from the same atom and count how many bonds of that kind present in a perticular structure.
!! score each structure accordingly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine intra_bond_factor(nl,str,tonstruc,str_quality_1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

integer::nl,i2,m8,l1,l2,l3,k1,k2,m13,m14,str(15000,20),str_quality_1(15000),tonstruc

print*,'enter intra_bond_factor'

do m8=1,15000
  str_quality_1(m8)=1
enddo

l3=1
do m8=1,tonstruc
  l2=1
  loop2:do k2=1+nl*2,nae-nlast,2
    do m13=1,atom
      l1=0
        loop4:do m14=1,atn(active_atoms(m13))
          if(atn(active_atoms(m13)).eq.1) then
            cycle loop4
          endif
          do k1=k2,k2+1
            if(str(m8,k1).eq.atoset(active_atoms(m13),m14))then
              l1=l1+1
            endif
          enddo
        enddo loop4
       if(l1.eq.2) then
         l2=l2+l3
         cycle loop2
       endif
     enddo
   enddo loop2
   str_quality_1(m8)=l2
enddo


!print*,'exit intra_bond_factor'
return
end subroutine intra_bond_factor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
