module quality_mod
implicit none
integer, pointer::str_quality_1(:)=>null(), str_quality_2(:)=>null(), bondq(:)=>null()
integer, pointer :: pref_radical(:)=>null(), mbondq(:)=>null(), quality_fac(:)=>null()
integer::tqlty,bqlty,sqlty,tnqs,nssym
integer::qulsym(15000), sigsym(15000),tnqs_sig
integer, allocatable::symq(:)
end module quality_mod
