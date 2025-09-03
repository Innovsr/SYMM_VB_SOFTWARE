module commondat_mod
implicit none

integer::nao,nae,niao,mult,nlast,nlp,atom,ndet,symm,asymm,itb,syb,nnb,rad,cov_space
integer::flgst,flg_ion,flg_cov,nfset,noq0,noq1,noq2,noq3,nnnatom,prad,nlpset,prt_cnt
integer::flg1,mnbond,nmbond,totrum,vacorb,MaxStrOepo,nsym,set_number,CovDim,ion_space
integer::prio_rad(20,100),lp(100),plpair(100,10),norad(20),valence_state(54)
integer, allocatable::nnat_bond(:, :), loopsymsc(:), all_at_num(:), active_orbs(:)
integer::atoset(200,20),at_num(88),at_covrad(88),atn(200),active_atoms(30)
integer::nrs,mset,radical,ovopt,vpt,tot_atom,sig_sym_flg, nstrt,nactorb,qflg,u1
character(len=2),public ::at_list(88),at_list_bold(88),int_num(40)
integer, allocatable::qq11(:),qq12(:),qq10(:),bondq14(:),strdet(:),detmnt(:,:),det_sign(:)
real*8,allocatable::biasval(:),dist_nnat(:),dist_mat(:,:)
real*8::ovval, prime_num(142)
character(len=6)::symtype
double precision, allocatable::dist_act_rel_mat(:,:)
integer::Rumwrite,set_order,rum_jj,rumset,symm_rumset,symm_maxval
character(len = 100)::STDOUT
character(len = 300)::out_folder_path 
integer, pointer::atsymset(:,:),syn(:),at_sym(:)
integer, allocatable ::loop_score_row(:), set_num(:), rumerstr(:,:), symm_groups(:,:), symm_groups_count(:)

DATA valence_state/1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,&
4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5/

DATA at_list/'H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al',&
'Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni',&
'Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc',&
'Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce',&
'Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta',&
'W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra'/

DATA at_list_bold/'H','HE','LI','BE','B','C','N','O','F','NE','NA','MG','AL',&
'SL','P','S','CL','AR','K','CA','SC','TI','V','CR','MN','FE','CO','NR',&
'CU','ZN','GA','GE','AS','SE','BR','KR','RB','SR','Y','ZR','NB','MO','TC',&
'RU','RH','PD','AG','CD','IN','SN','SB','TE','I','XE','CS','BA','LA','CE',&
'PR','ND','PM','SM','EU','GD','TB','DY','HO','ER','TM','YB','LU','HF','TA',&
'W','RE','OS','IR','PT','AU','HG','TL','PB','BI','PO','AT','RN','FR','RA'/

DATA int_num/'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16',&
'17','18','19','20','21','22','23','24','25','26','27','28','29','30','31','32','33',&
'34','35','36','37','38','39','40'/

DATA prime_num/11.0, 13.0, 17.0, 19.0, 23.0, 29.0, 31.0, 37.0, 41.0, 43.0, 47.0, 53.0, 59.0, 61.0&
, 67.0, 71.0, 73.0, 79.0, 83.0, 89.0, 97.0, 101.0, 103.0, 107.0, 109.0, 113.0, 127.0, 131.0, 137.0, 139.0,&
 149.0, 151.0, 157.0, 163.0, 167.0, 173.0, 179.0, 181.0, 191.0, 193.0, 197.0, 199.0, 211.0, 223.0, 227.0, 229.0,&
 233.0, 239.0, 241.0, 251.0, 257.0, 263.0, 269.0, 271.0, 277.0, 281.0, 283.0, 293.0, 307.0, 311.0, 313.0, 317.0,&
 331.0, 337.0, 347.0, 349.0, 353.0, 359.0, 367.0, 373.0, 379.0, 383.0, 389.0, 397.0, 401.0, 409.0, 419.0, 421.0,&
 432.0, 433.0, 439.0, 443.0, 449.0, 457.0, 461.0, 463.0, 467.0, 479.0, 487.0, 491.0, 499.0, 503.0, 509.0, 521.0,&
 523.0, 541.0, 547.0, 557.0, 563.0, 569.0, 571.0, 577.0, 587.0, 593.0, 599.0, 601.0, 607.0, 613.0, 617.0, 619.0,&
 631.0, 641.0, 643.0, 647.0, 653.0, 659.0, 661.0, 673.0, 677.0, 683.0, 691.0, 701.0, 709.0, 719.0, 727.0, 733.0,&
 739.0, 743.0, 751.0, 757.0, 761.0, 769.0, 773.0, 787.0, 797.0, 809.0, 811.0, 821.0, 823.0, 827.0, 829.0, 839.0/

!!!! Covalent radii (in pm unit)  taken from 'Webelement' ... from the link
!"http://crystalmaker.com/support/tutorials/atomic-radii/index.html"

DATA at_covrad /37, 32, 134, 90, 82, 77, 75, 73, 71, 69, 154, 130, 118, 111,&
106, 102, 99, 97, 196, 174, 144, 136, 125, 127, 139, 125, 126, 121, 138, 131,&
126, 122, 119, 116, 114, 110, 211, 192, 162, 148, 137, 145, 156, 126, 135, 131,153,&
148, 144, 141, 138, 135, 133, 130, 225, 198, 169, 204, 203, 201, 199, 198, 198,&
196, 194, 192, 192, 189, 190, 187, 160, 150, 138, 146, 159, 128, 137,&
128, 144, 149, 148, 147, 146, 140, 150, 145, 260, 221/
save

end module commondat_mod

