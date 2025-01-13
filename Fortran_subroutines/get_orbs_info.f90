!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_orbs_info(atoset_array, norbsym_array, active_atom_array, atn_array,&
orbsym_array, atoset_row, atoset_col, norbsym_size, active_size, atn_size, orbsym_row, orbsym_col)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat1
use commondat
implicit none

integer::atoset_row, atoset_col,norbsym_size,atn_size,orbsym_row, orbsym_col,active_size
integer::atoset_array(atoset_row, atoset_col), norbsym_array(norbsym_size), atn_array(atn_size)
integer::orbsym_array(orbsym_row, orbsym_col), active_atom_array(active_size), i, j

atoset = atoset_array
norbsym = norbsym_array
active_atoms = active_atom_array
atn = atn_array
orbsym = orbsym_array

do i = 1,atoset_row 
print*,'from_Fortran: atoset ',(atoset(i,j),j=1,atoset_col)
enddo
print*,'from_fortran: norbsym_array',(norbsym(i),i=1,norbsym_size)
print*,'from_fortran: active_atom_array',(active_atoms(i),i=1,active_size)
print*,'from_fortran: atn_array',(atn(i),i=1,atn_size)

do i = 1,orbsym_row 
print*,'from_Fortran: orbsym ',(orbsym(i,j),j=1,orbsym_col)
enddo
return
end subroutine get_orbs_info
