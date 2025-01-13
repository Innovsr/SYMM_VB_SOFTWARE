!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_splkeywds(chinst, symm_py, symtype_py, set_order_py, nset, mout_py, ovlp_py, itb_py, nnb_py,&
                syb_py, mnbond_py, radical_py, nmbond_py)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
use commondat1
implicit none

integer::chinst, symm_py, set_order_py, nset, mout_py, ovlp_py, itb_py
integer::nnb_py, syb_py, mnbond_py, radical_py, nmbond_py!, main_bond_py
character(len = 6)::symtype_py

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

print*,'chinst, symm, set_order, nset, mout, ovlp, itb',chinst, symm, set_order, nset, mset, ovlp_int, itb
print*,'nnb, syb, mnbond, radical, mnbond, main_bond',nnb, syb, mnbond, radical, mnbond!, main_bond
print*,'symtype',symtype

return
end subroutine get_splkeywds
