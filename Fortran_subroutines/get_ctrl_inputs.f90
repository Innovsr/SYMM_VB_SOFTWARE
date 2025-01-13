!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_ctrl_inputs(geometry_unit, nao_py, nae_py, nmul, output_file_name, &
                chinst, symm_py, symtype_py, set_order_py, nset, mout_py, ovlp_py, itb_py, nnb_py,&
                syb_py, mnbond_py, radical_py, nmbond_py, symat_array, coordx, coordy, coordz, symatno_array,&
                atoset_array, norbsym_array, active_atom_array, atn_array, orbsym_array) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
use commondat1
implicit none

integer:: nao_py, nae_py, nmul
character(len = 100)::geometry_unit, output_file_name
integer::chinst, symm_py, set_order_py, nset, mout_py, ovlp_py, itb_py
integer::nnb_py, syb_py, mnbond_py, radical_py, nmbond_py!, main_bond_py
character(len = 6)::symtype_py
real*8::coordx(100), coordy(100), coordz(100), symatno_array(20)
character(len=5)::symat_array(20)
!integer::atoset_row, atoset_col,norbsym_size,atn_size,orbsym_row, orbsym_col,active_size
integer::atoset_array(200, 20), norbsym_array(50), atn_array(200)
integer::orbsym_array(20, 20), active_atom_array(30), i, j


nao = nao_py
nae = nae_py
mult = nmul
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
!main_bond = main_bond_py
symtype = symtype_py
symat=symat_array
symatno = symatno_array
atoset = atoset_array
norbsym = norbsym_array
active_atoms = active_atom_array
atn = atn_array
orbsym = orbsym_array

print*,'nao, nae, nmul',nao, nae, mult
print*,'geometry_unit',geometry_unit
print*,'file_name',output_file_name
print*,'chinst, symm, set_order, nset, mout, ovlp, itb',chinst, symm, set_order, nset, mset, ovlp_int, itb
print*,'nnb, syb, mnbond, radical, mnbond, main_bond',nnb, syb, mnbond, radical, mnbond!, main_bond
print*,'symtype',symtype
do i = 1,10 
print*,'from_fortran:symat, coordx, coordy, coordz, symatno ',symat(i), coordx(i), coordy(i), coordz(i), symatno(i)
enddo

!print*,'atoset_row, atoset_col,norbsym_size,atn_size,orbsym_row, orbsym_col,active_size',atoset_row, atoset_col,&
!       norbsym_size,atn_size,orbsym_row, orbsym_col,active_size

do i = 1,10 
print*,'from_Fortran: atoset ',(atoset(i,j),j=1,10)
enddo
print*,'from_fortran: norbsym_array',(norbsym(i),i=1,10)
print*,'from_fortran: active_atom_array',(active_atoms(i),i=1,10)
print*,'from_fortran: atn_array',(atn(i),i=1,10)

do i = 1,10 
print*,'from_Fortran: orbsym ',(orbsym(i,j),j=1,10)
enddo

stop
end subroutine get_ctrl_inputs
