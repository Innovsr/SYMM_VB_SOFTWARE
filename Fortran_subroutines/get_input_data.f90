subroutine get_input_data(geometry_unit, nao_py, nae_py, nmul, output_file_name, &
                           atoset_row, atoset_col, norbsym_size, active_size, atn_size,&
                           orbsym_row, orbsym_col, atoset_array, norbsym_array,&
                           active_atom_array, atn_array, orbsym_array, symat_array, coordx,&
                           coordy, coordz, symatno_array, array_size, chinst)
                   !, symm_py, set_order_py, &
                    !       nset, mout_py, ovlp_py, itb_py, nnb_py, syb_py, mnbond_py, radical_py, nmbond_py)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
use commondat1
implicit none

! Input variables for geometry and control inputs
integer :: nao_py, nae_py, nmul
character(len=100) :: geometry_unit, output_file_name

! Input variables for orbitals information
integer :: atoset_row, atoset_col, norbsym_size, active_size, atn_size, orbsym_row, orbsym_col
integer :: atoset_array(atoset_row, atoset_col), norbsym_array(norbsym_size)
integer :: active_atom_array(active_size), atn_array(atn_size)
integer :: orbsym_array(orbsym_row, orbsym_col)

! Input variables for geometry information
integer :: i, array_size
real*8 :: coordx(array_size), coordy(array_size), coordz(array_size),symatno_array(array_size)
character(len=5) :: symat_array(array_size)

! Input variables for special keywords
integer :: chinst, symm_py, set_order_py, nset, mout_py, ovlp_py, itb_py
integer :: nnb_py, syb_py, mnbond_py, radical_py, nmbond_py

! Local variables
integer :: j

! Section 1: Handle control inputs
nao = nao_py
nae = nae_py
mult = nmul

print*, 'Control Inputs:'
print*, 'nao, nae, nmul:', nao, nae, mult
print*, 'geometry_unit:', geometry_unit
print*, 'file_name:', output_file_name

! Section 2: Handle orbitals information
atoset = atoset_array
norbsym = norbsym_array
active_atoms = active_atom_array
atn = atn_array
orbsym = orbsym_array

print*, 'Orbitals Information:'
do i = 1, atoset_row
  print*, 'atoset:', (atoset(i, j), j = 1, atoset_col)
end do

print*, 'norbsym_array:', (norbsym(i), i = 1, norbsym_size)
print*, 'active_atom_array:', (active_atoms(i), i = 1, active_size)
print*, 'atn_array:', (atn(i), i = 1, atn_size)

do i = 1, orbsym_row
  print*, 'orbsym:', (orbsym(i, j), j = 1, orbsym_col)
end do

! Section 3: Handle geometry information
symat = symat_array
symatno = symatno_array

print*, 'Geometry Information:'
do i = 1, array_size
  print*, 'symat, coordx, coordy, coordz, symatno:', symat(i), coordx(i), coordy(i), coordz(i), symatno(i)
end do

! Section 4: Handle special keywords
symm = symm_py
set_order = set_order_py
mset = mout_py
ovlp_int = ovlp_py
itb = itb_py
nnb = nnb_py
syb = syb_py
mnbond = mnbond_py
radical = radical_py
nmbond = nmbond_py
! symtype = symtype_py

print*, 'Special Keywords:'
print*, 'chinst, symm, set_order, nset, mout, ovlp, itb:', chinst, symm, set_order, nset, mset, ovlp_int, itb
print*, 'nnb, syb, mnbond, radical, nmbond:', nnb, syb, mnbond, radical, nmbond
!print*, 'symtype:', symcheck

return
end subroutine get_input_data
