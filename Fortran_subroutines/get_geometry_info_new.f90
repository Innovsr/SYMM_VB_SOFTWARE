!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_geometry_info(symat_array, coordx, coordy, coordz, symatno_array)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat1
use commondat
implicit none

integer::i,n
real*8, allocatable::coordx(:), coordy(:), coordz(:),symatno_array(:)
character(len=5), allocatable::symat_array(:)

n = size(coordx)

symat = symat_array
symatno = symatno_array

do i = 1,n 
print*,'symat, coordx, coordy, coordz, symatno ',symat(i), coordx(i), coordy(i), coordz(i), symatno(i)
enddo

return
end subroutine get_geometry_info
