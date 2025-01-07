!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module commondat1
implicit none

integer::orbs(50,500),orbsl,freezorb(100),frzn,orbn(50),sl1(10),orbsym(20,20),norbsym(50),frag_at(20,20)&
,frag_atn(20),atm_nb_orbs1(20),atm_nb_orbs(20),ovlp_int,nbasatom(20),nbasis,ovlp_bas_num(20,10),lnum&
,num_act_orb_typ(20),ind_mat(1000,1000)
character(len=5)::symat(20),all_atbas(500)
real*8::symatno(20)
double precision::orb_ovlp_mat(20,20),orb_ovlp_mat1(20,20),ovlp_mat_norm(5000,100)

save
end module commondat1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
