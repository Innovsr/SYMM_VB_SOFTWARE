!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_ctrl_inputs(geometry_unit, nao_py, nae_py, nmul, output_file_name, &
                chinst, symm_py, symtype_py, set_order_py, nset_py, mout_py, ovlp_py, itb_py, nnb_py,&
                syb_py, mnbond_py, radical_py, nmbond_py, symat_array, coordx, coordy, coordz, symatno_array,&
                atoset_array, norbsym_array, active_atom_array, atn_array, orbsym_array, flgst_py, total_atoms,&
                niao_py, active_atm_num, output_folder) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
use commondat1
implicit none

common/ats/atsymset,nsym,syn,at_sym

integer:: nao_py, nae_py, nmul, flgst_py, total_atoms, niao_py, active_atm_num
character(len = 100)::geometry_unit, output_file_name
character(len = 300)::output_folder, outfile
integer::chinst, symm_py, set_order_py, nset_py, mout_py, ovlp_py, itb_py
integer::nnb_py, syb_py, mnbond_py, radical_py, nmbond_py!, main_bond_py
character(len = 6)::symtype_py
real*8::coordx(100), coordy(100), coordz(100), symatno_array(100)
character(len=5)::symat_array(100)
!integer::atoset_row, atoset_col,norbsym_size,atn_size,orbsym_row, orbsym_col,active_size
integer::atoset_array(200, 20), norbsym_array(50), atn_array(200)
integer::orbsym_array(20, 20), active_atom_array(30), i, j, k, k5, k6
integer::atsymset(20,20),nsym,syn(50),at_sym(50)

nfset = nset_py
niao = niao_py
tot_atom = total_atoms
flgst = flgst_py
nao = nao_py
nae = nae_py
mult = nmul
symm = symm_py
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
symtype = symtype_py
symat=symat_array
symatno = symatno_array
atoset = atoset_array
norbsym = norbsym_array
active_atoms = active_atom_array
atn = atn_array
orbsym = orbsym_array
STDOUT = output_file_name
out_folder_path = output_folder

print*,'nao, nae, nmul',nao, nae, mult
print*,'geometry_unit',geometry_unit
print*,'file_name',output_file_name
print*,'chinst, symm, set_order, nset, mout, ovlp, itb',chinst, symm, set_order, nset_py, mset, ovlp_int, itb
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    Output files   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! structures.dat:: stores all possible structures with identifications of rumer structurs                 !
! structure_set_i.dat:: output file, i=1, 2, .... depeding on volume of sets; one file contain 75000 sets !
! Rumer_Sets_all.dat:: unique Rumer sets for all possible permutations of orbitals                        !
! quality_str.dat:: all structures arranged according to their qualities                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
open(unit=7,file='structures.dat',status='unknown')
outfile=trim(out_folder_path)//trim('/')//trim(STDOUT)//trim('_1.out')
open(unit=10,file=outfile,status='unknown')
open(unit=23,file='Rumer_Sets_all.dat',status='unknown')
open(unit=35,file='quality_str.dat',status='unknown')


noq0=100
noq1=100
noq2=100
noq3=100

if(itb.eq.1.and.syb.eq.1.and.nnb.eq.1.and.radical.eq.1.and.mnbond.eq.1.and.flg1.ne.1)then
  qflg=1
endif

vacorb = nao - nae
nlp = nae - nao
if (vacorb.gt.1) nlp = 0

nlast=(mult-1)

if(ovlp_py.eq.0)ovopt=0
if(ovlp_py.eq.1)ovopt=1
#if(ovlp.eq.1.and.nfset.eq.0)ovopt=0
#if(ovlp.eq.1.and.nfset.eq.1)ovopt=1
#if(ovlp.eq.1.and.nfset.eq.2)ovopt=1
#if(ovlp.eq.1.and.nfset.eq.3)ovopt=1
!!!!!!!!! "vpt" is the option for user specifying overlap value "ovval". It will work
!!!!!!!!!!!!!only for vpt=1. To lock it please put any other value
vpt=1
itb=itb+1
syb=syb+1
nnb=nnb+1
radical=radical+1
mnbond=mnbond+1

write(7,*)'******************************************************************************************'
write(7,900)'*','active orbs =',nao,'active electrons =',nae,'multiplicity =',mult,'inactive orbs =',niao
write(7,*)'******************************************************************************************'
write(7,*)'******************************************************************************************'
900 format(a,1x,a,1x,I4,3x,a,1x,I4,3x,a,1x,I4,3x,a,1x,I4)

if (geometry_unit.eq.'Bohr') then
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

do i=1,nsym
print*,'atsymset',nsym,at_sym(i),(atsymset(i,j),j=1,syn(i))
enddo
!do i=1,atom
!act_at_num(i)=symatno(active_atoms(i))
!val_state_num(i)=valence_state(act_at_num(i))
!enddo

call geocal(coordx,coordy,coordz)

!flgst: user input for 'str' keyword; 1: covalent (cov), 2: ionic (ion), 3: both (full)
if(flgst.eq.1.or.flgst.eq.2) then
  call cov_struc
endif

!!! ionic str calculation has been stopped for this version

!if(flgst.eq.1.or.flgst.eq.3) then
!call ion_struc
!endif

stop
end subroutine get_ctrl_inputs
