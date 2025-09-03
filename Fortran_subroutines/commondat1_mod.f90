!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module commondat1_mod
implicit none

integer::orbsym(20,20),norbsym(50),ovlp_int
character(len=5)::symat(100)
integer, allocatable::symatno(:),atm_nb_orbs1(:),atm_nb_orbs(:),ind_mat(:,:)
!double precision::ovlp_mat_norm(5000,100)

save
end module commondat1_mod
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
