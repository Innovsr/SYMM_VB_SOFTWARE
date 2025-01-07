module commondat
implicit none

integer::nao,nae,niao,mult,nlast,nlp,tncs,tnis,is,atom,nabsym,ndet,symm,repl,itb,syb,nnb,rad &
,flgst,noqult,rplc,flg_ion,flg_cov,imbd,nfset,noq0,noq1,noq2,noq3,noq,niabd,nialp,niach,serial &
,nrad,choice,pbond(1000),mbond(1000),nnnatom,nnat_bond(100,2),prio_rad(20,100),norad(20),prad &
,nlpset,lp(100),plpair(100,10),flg1,uoptstr,mnbond,nmbond,main_bond(100),iabd(50),ialp(50),iach(50)
integer::atoset(200,20),at_num(88),at_covrad(88),key_frag,frag_cntr(10)&
,at_ab_symset(20,20),tot_ndet,nrs,mset,loopsymsc(15000) &
,atn(50),absyn(50),qult(100),radical,actv_atom(30),num_frag_cntr,input_flg,ovopt,vpt,all_at_num(100)&
,strt_struc(100,15),nstrt,repl_struc(100,15),nnmat_act(100,2),nnmat_inact(100,2),nactorb,active_orbs(20),&
active_atoms(30),nnat_bond_inact(100,2),nnatominact,tot_atom,biasmat(100),sig_sym_flg,valence_state(54)&
,bond_count(1000,100),totrum,vacorb
character(len=2),public ::at_list(88),at_list_bold(88),int_num(40)
integer,public::atm_nb_sym(20),nalpha,nbeta,str_det_sec(15000,1000),num_norbsym1(100),qflg
!real,dimension(:,:),allocatable::ovlp_mat
real*8::ovval,dist_rel_mat(20,20),biasval(20)&
,dist_nnat(20)
real*8::prime_num(142),dist_mat(20,20),iab_length
character(len=6),public::symtype,norbsym1(100,20)
double precision::dist_act_rel_mat(20,20)
integer::strdet(10000),detmnt(10000,15),det_sign(10000),Rumwrite,set_order


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

end module commondat

