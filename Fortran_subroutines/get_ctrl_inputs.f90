!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_ctrl_inputs(geometry_unit, nao_py, nae_py, nmul, output_file_name, &
                chinst, symm_py, asymm_py, set_order_py, nset_py, mout_py, ovlp_py, itb_py, nnb_py,&
                syb_py, mnbond_py, radical_py, nmbond_py, symat_array, coordx_py, coordy_py, coordz_py, symatno_array,&
                atoset_array, norbsym_array, active_atom_array, atn_array, orbsym_array, flgst_py, total_atoms,&
                niao_py, active_atm_num, output_folder) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat_mod
use commondat1_mod
use mod_cov_struc
use coordinates_mod
implicit none

integer:: nao_py, nae_py, nmul, flgst_py, total_atoms, niao_py, active_atm_num
character(len = 100):: output_file_name,tempoutfile,all_struc_file
character(len = 300)::output_folder, outfile
integer::chinst, symm_py, asymm_py, set_order_py, nset_py, mout_py, ovlp_py, itb_py
integer::nnb_py, syb_py, mnbond_py, radical_py, nmbond_py, cov_nlp, cov_bond
real*8::symatno_array(100),coordx_py(100), coordy_py(100), coordz_py(100) 
character(len=5)::geometry_unit,symat_array(100)
integer::atoset_array(200, 20), norbsym_array(50), atn_array(200)
integer::orbsym_array(20, 20), active_atom_array(30), i, j, k, k5, k6
character(len=30)::lowercase

allocate(atsymset(nao_py, nao_py))
allocate(syn(nao_py))
allocate(at_sym(nao_py))
allocate(symatno(total_atoms))
allocate(all_at_num(total_atoms))

nfset = nset_py
niao = niao_py
tot_atom = total_atoms
flgst = flgst_py
nao = nao_py
nae = nae_py
mult = nmul
symm = symm_py
asymm = asymm_py
set_order = set_order_py
mset = mout_py
!ovlp_int = ovlp_py
itb = itb_py
nnb = nnb_py
syb = syb_py
mnbond = mnbond_py
radical = radical_py
nmbond = nmbond_py
!main_bond = main_bond_py
!symtype = symtype_py
symat=symat_array
symatno = symatno_array
atoset = atoset_array
norbsym = norbsym_array
active_atoms = active_atom_array
atn = atn_array
orbsym = orbsym_array
STDOUT = output_file_name
out_folder_path = output_folder
coordx = coordx_py
coordy = coordy_py
coordz = coordz_py

print*,'nao, nae, nmul',nao, nae, mult
print*,'geometry_unit',geometry_unit
print*,'file_name',output_file_name
print*,'chinst, symm, set_order, nset, mout, ovlp, itb',chinst, symm, set_order, nset_py, mset, ovlp_int, itb
print*,'nnb, syb, mnbond, radical, mnbond, main_bond',nnb, syb, mnbond, radical, mnbond!, main_bond
!print*,'symtype',symtype
do i = 1,10 
print*,'from_fortran:symat, coordx, coordy, coordz, symatno ',symat(i), coordx(i), coordy(i), coordz(i), symatno(i)
enddo
!print*,'atoset_row, atoset_col,norbsym_size,atn_size,orbsym_row, orbsym_col,active_size',atoset_row, atoset_col,&
!       norbsym_size,atn_size,orbsym_row, orbsym_col,active_size

do i = 1,tot_atom 
print*,'from_Fortran: atoset ',(atoset(i,j),j=1,10)
enddo
print*,'from_fortran: norbsym_array',(norbsym(i),i=1,10)
print*,'from_fortran: active_atom_array',(active_atoms(i),i=1,10)
print*,'from_fortran: atn_array',(atn(i),i=1,10)

do i = 1,10 
print*,'from_Fortran: orbsym ',(orbsym(i,j),j=1,10)
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    Output files   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! structures.dat:: stores all possible structures with identifications of rumer structurs                 !
! molecule_i.dat:: output file, i=1, 2, .... depeding on volume of sets; one file contain 75000 sets !
! out.temp : file to write sets for python_GUI.
! Rumer_Sets_all.dat:: unique Rumer sets for all possible permutations of orbitals                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
u1=1 !
outfile=trim(out_folder_path)//trim('/')//trim(STDOUT)//trim('_1.out')
open(unit=9+u1,file=outfile,status='unknown')
tempoutfile=trim(out_folder_path)//trim('/')//trim('out.temp')
open(unit=5,file=tempoutfile,status='unknown')
all_struc_file=trim(out_folder_path)//trim('/')//trim('structures.dat')
open(unit=7,file=all_struc_file,status='unknown')
open(unit=23,file='Rumer_Sets_all.dat',status='unknown')


noq0=10000
noq1=10000
noq2=10000
noq3=10000
if (chinst.eq.0) flg1 = 1
if (chinst.eq.1) flg1 = 0

if(itb.eq.1.and.syb.eq.1.and.nnb.eq.1.and.radical.eq.1.and.mnbond.eq.1.and.flg1.ne.1)then
  qflg=1
endif

vacorb = nao - nae
nlp = nae - nao
if (vacorb.gt.1) nlp = 0

nlast=(mult-1)

if(ovlp_py.eq.0)ovopt=0
if(ovlp_py.eq.1)ovopt=1
!if(ovlp.eq.1.and.nfset.eq.0)ovopt=0
!if(ovlp.eq.1.and.nfset.eq.1)ovopt=1
!if(ovlp.eq.1.and.nfset.eq.2)ovopt=1
!if(ovlp.eq.1.and.nfset.eq.3)ovopt=1
!!!!!!!!! "vpt" is the option for user specifying overlap value "ovval". It will work
!!!!!!!!!!!!!only for vpt=1. To lock it please put any other value
vpt=1
itb=itb+1
syb=syb+1
nnb=nnb+1
radical=radical+1
mnbond=mnbond+1


if (lowercase(geometry_unit).eq.'bohr') then
coordx=coordx*0.529177
coordy=coordy*0.529177
coordz=coordz*0.529177
endif
all_at_num=int(symatno)

k=0
loop1:do i=1,tot_atom
  do j=1,active_atm_num
    if(active_atoms(j).eq.i) then
      k=k+1
      at_num(k)=int(symatno(i))
      exit
    endif
  enddo
  if(k.eq.active_atm_num) exit
  print*,'iiiii',i
  !coordx(i)=coordx(i)*0.529177
  !coordy(i)=coordy(i)*0.529177
  !coordz(i)=coordz(i)*0.529177
enddo loop1
print*,'at_num',tot_atom,at_num,active_atm_num

k5=0
k6=0
do i=1,4
if(norbsym(i).ne.0)then
k5=k5+1
if(i.lt.4)then
k6=1
else
k6=k6+1
endif
at_sym(k5)=k6
syn(k5)=norbsym(i)
do j=1,norbsym(i)
atsymset(k5,j)=orbsym(i,j)
enddo
endif
enddo
nsym=k5

print*,'atsymset'
do i=1,nsym
print*,'atsymset',nsym,at_sym(i),(atsymset(i,j),j=1,syn(i))
enddo


call geocal
print*,'i am here',flgst
prt_cnt =0
!flgst: user input for 'str' keyword; 1: covalent (cov), 2: ionic (ion), 3: both (full)
if(flgst.eq.1.or.flgst.eq.2) then
  prt_cnt=prt_cnt+1
  flg_ion=0
  flg_cov=1
  cov_space=0
  call max_avail_str
  call cov_struc
endif

if(allocated(symq))then
  deallocate(symq)
endif

if(flgst.eq.1.or.flgst.eq.3) then
  prt_cnt=prt_cnt+1
  flg_ion=1
  flg_cov=0
  ion_space=0
  cov_nlp = nlp
  cov_bond = (nao - 2 * cov_nlp - nlast)/2
  print*,'cov_bond',cov_bond
  do i = 1, cov_bond
    nlp = cov_nlp + i
    print*,'number of lonepair',i
    call max_avail_str
    call cov_struc
  enddo
endif
close(7)
close(5)

end subroutine get_ctrl_inputs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine max_avail_str
!!!! calculate the maximum available structures

use commondat_mod
implicit none

integer::i,noeo, navo, spin, term1, term2, term3, product_term, denom, factorial, comb

!noeo = nao - nlp
noeo = nae - 2*nlp
navo = nao -nlp
spin = (mult - 1)/2
!term1 = comb(nao, nao - noeo)
term1 = comb(nao, nlp)
term2 = comb(noeo, int(2*spin))
term3 = comb(navo, noeo)

print*,'nao, nae, mult, nlp, noeo, navo',nao, nae, mult, nlp, noeo, navo
print*,'term1, term2, term3',term1, term2, term3
product_term = 1  ! Start with 1 for product computation

do i = 0, (noeo / 2) - spin - 1
   product_term = product_term * comb(noeo - 2 * spin - 2 * i, 2)
end do

denom = factorial(int((noeo / 2) - spin))
MaxStrOepo = int((term1 * term2 * term3 * product_term) / denom) !total structure one electron per orbital

call wigner(noeo, CovDim)
print*,'nnb',nnb

return
end subroutine max_avail_str
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
