!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_geometry_info(symat_array, coordx, coordy, coordz, symatno_array, array_size)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
use commondat1
implicit none

integer::i, array_size
real*8::coordx(array_size), coordy(array_size), coordz(array_size), symatno_array(array_size)
character(len=5)::symat_array(array_size)

symat=symat_array
symatno = symatno_array

do i = 1,array_size 
print*,'from_fortran:symat, coordx, coordy, coordz, symatno ',symat(i), coordx(i), coordy(i), coordz(i), symatno(i)
enddo

return
end subroutine get_geometry_info
