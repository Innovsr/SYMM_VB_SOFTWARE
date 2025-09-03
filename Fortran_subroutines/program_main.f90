!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                            Chemical Insight Structure Sets                             !!
!!                                                                                        !!
!!                               Written by: Dr. Sourav Roy                               !!
!!                                                                                        !!
!!                    School Of Pharmacy, Hebrew University of Jerusalem                  !!
!!                                                                                        !!
!!                                    Jerusalem, Israel                                   !!
!!                                                                                        !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program main
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none

integer:: nao_py, nae_py, nmul, flgst_py, total_atoms, niao_py, active_atm_num,i, j
character(len = 100)::geometry_unit, output_file_name
character(len = 300)::output_folder
integer::chinst, symm_py, asymm_py, set_order_py, nset_py, mout_py, ovlp_py, itb_py
integer::nnb_py, syb_py, mnbond_py, radical_py, nmbond_py, flg!, symtype_py
!character(len = 6)::symtype_py
real*8, target::coordx(100), coordy(100), coordz(100), symatno_array(100)
character(len=5), target::symat_array(100)
integer, target::atoset_array(200, 20), norbsym_array(50), atn_array(200)
integer, target::orbsym_array(20, 20), active_atom_array(30)
!integer, allocatable, intent(out)::atsymset(:,:),nsym,syn(:),at_sym(:)
integer, pointer::sub_atoset_array(:), sub_norbsym_array(:), sub_atn_array(:)
integer, pointer::sub_orbsym_array(:), sub_active_atom_array(:)
integer, pointer::sub_atoset_array1(:,:),sub_orbsym_array1(:,:)
real*8, pointer::sub_coordx(:), sub_coordy(:), sub_coordz(:), sub_symatno_array(:)
character(len=5), pointer::sub_symat_array(:)

print*,'write your input choice'
read(*,*) flg

print*,'i am here1'
symat_array = ''
atoset_array = 0
print*,'i am here1'
norbsym_array = 0
atn_array = 0
orbsym_array = 0
active_atom_array = 0
coordx = 0.0 
coordy = 0.0 
coordz = 0.0 
symatno_array = 0.0
print*,'i am here2'

if(flg.eq.1) then
geometry_unit = 'Bohr'
nao_py = 5
nae_py = 5
nmul = 2
output_file_name = 'C5H5'
chinst = 1
symm_py = 1
asymm_py = 0
!symtype_py = 0
set_order_py = 0
nset_py = 0
mout_py = 1
ovlp_py = 0
itb_py = 0
nnb_py = 1
syb_py = 0
mnbond_py = 0
radical_py = 0
nmbond_py = 0
print*,'i am here3'

allocate(sub_symat_array(10))
allocate(sub_coordx(10)) 
allocate(sub_coordy(10)) 
allocate(sub_coordz(10)) 
allocate(sub_symatno_array(10))
!allocate(sub_atoset_array(5))
allocate(sub_norbsym_array(4))
allocate(sub_active_atom_array(5))
allocate(sub_atn_array(5))
allocate(sub_orbsym_array(5))

sub_symat_array => symat_array(1:10)
sub_coordx => coordx(1:10)
sub_coordy => coordy(1:10)
sub_coordz => coordz(1:10)
sub_symatno_array => symatno_array(1:10)
sub_atoset_array => atoset_array(1:5, 1)
sub_norbsym_array => norbsym_array(1:4)
sub_active_atom_array => active_atom_array(1:5)
sub_atn_array => atn_array(1:5)
sub_orbsym_array => orbsym_array(3 ,1:5)

sub_symat_array = (/'C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H'/)
sub_coordx = (/0.0, 2.2, 1.4, -1.4, -2.2, 0.0, 4.1, 2.5, -2.5, -4.1/)
sub_coordy = (/2.2, 0.6, -1.7, -1.7, 0.6, 4.2, 1.3, -3.4, -3.4, 1.3/)
sub_coordz = (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
sub_symatno_array = (/6, 6, 6, 6, 6, 1, 1, 1, 1, 1/)
sub_atoset_array = [16, 17, 18, 19, 20]
sub_norbsym_array = (/0, 0, 5, 0/)
sub_active_atom_array = (/1, 2, 3, 4, 5/)
sub_atn_array = (/1,1,1,1,1/)
sub_orbsym_array = [16, 17, 18, 19, 20]
flgst_py = 1
total_atoms = 10
niao_py = 15
active_atm_num = 5
output_folder = '/home/sourav/Desktop/SYMM_VB_SOFTWARE/Python_interface/C5H5_output'
endif

if(flg.eq.2)then
geometry_unit = 'Bohr'
nao_py = 8
nae_py = 10
nmul = 1
output_file_name = 'NCN'
chinst = 1
symm_py = 1
asymm_py = 0
!symtype_py = 0
set_order_py = 0
nset_py = 0
mout_py = 1
ovlp_py = 0
itb_py = 1
nnb_py = 2
syb_py = 3
mnbond_py = 0
radical_py = 0
nmbond_py = 0
print*,'i am here3'

allocate(sub_symat_array(3))
allocate(sub_coordx(3)) 
allocate(sub_coordy(3)) 
allocate(sub_coordz(3)) 
allocate(sub_symatno_array(3))
allocate(sub_atoset_array1(3, 3))
allocate(sub_norbsym_array(4))
allocate(sub_active_atom_array(3))
allocate(sub_atn_array(3))
allocate(sub_orbsym_array1(4, 4))

sub_symat_array => symat_array(1:3)
sub_coordx => coordx(1:3)
sub_coordy => coordy(1:3)
sub_coordz => coordz(1:3)
sub_symatno_array => symatno_array(1:3)
sub_atoset_array1 => atoset_array(1:3, 1:3)
sub_norbsym_array => norbsym_array(1:4)
sub_active_atom_array => active_atom_array(1:3)
sub_atn_array => atn_array(1:3)
sub_orbsym_array1 => orbsym_array(1:4 ,1:4)

sub_symat_array = (/'N', 'C', 'N'/)
sub_coordx = (/0.0, 0.0, 0.0/)
sub_coordy = (/0.0, 0.0, 0.0/)
sub_coordz = (/-2.3, 0.0, 2.3/)
sub_symatno_array = (/7, 6, 7/)
sub_atoset_array1 = transpose(reshape([6, 7, 8, 9, 10, 0, 11, 12, 13],[3,3]))
sub_norbsym_array = (/3, 3, 0, 2/)
sub_active_atom_array = (/1, 2, 3/)
sub_atn_array = (/3,2,3/)
sub_orbsym_array1 =  transpose(reshape([7, 9, 11, 0, 8, 10, 12, 0, 0, 0, 0, 0,  6, 13, 0, 0],[4,4]))
flgst_py = 1
total_atoms = 3
niao_py = 5
active_atm_num = 3
output_folder = '/home/sourav/Desktop/SYMM_VB_SOFTWARE/Python_interface/C5H5_output'

endif
if(flg.eq.3)then
geometry_unit = 'Bohr'
nao_py = 8
nae_py = 8
nmul = 1
output_file_name = 'C2'
chinst = 1
symm_py = 1
asymm_py = 0
!symtype_py = 0
set_order_py = 0
nset_py = 0
mout_py = 1
ovlp_py = 0
itb_py = 1
nnb_py = 0
syb_py = 2
mnbond_py = 0
radical_py = 0
nmbond_py = 0
print*,'i am here3'

allocate(sub_symat_array(2))
allocate(sub_coordx(2)) 
allocate(sub_coordy(2)) 
allocate(sub_coordz(2)) 
allocate(sub_symatno_array(2))
allocate(sub_atoset_array1(2, 4))
allocate(sub_norbsym_array(4))
allocate(sub_active_atom_array(2))
allocate(sub_atn_array(2))
allocate(sub_orbsym_array1(4, 4))

sub_symat_array => symat_array(1:2)
sub_coordx => coordx(1:2)
sub_coordy => coordy(1:2)
sub_coordz => coordz(1:2)
sub_symatno_array => symatno_array(1:2)
sub_atoset_array1 => atoset_array(1:2, 1:4)
sub_norbsym_array => norbsym_array(1:3)
sub_active_atom_array => active_atom_array(1:2)
sub_atn_array => atn_array(1:2)
sub_orbsym_array1 => orbsym_array(1:4 ,1:4)

sub_symat_array = (/'C', 'C'/)
sub_coordx = (/0.0, 0.0/)
sub_coordy = (/0.0, 0.0/)
sub_coordz = (/1.315, -1.315/)
sub_symatno_array = (/6, 6/)
sub_atoset_array1 = transpose(reshape([3, 5, 7, 9, 4, 6, 8, 10],[4,2]))
sub_norbsym_array = (/0, 2, 2, 4/)
sub_active_atom_array = (/1, 2/)
sub_atn_array = (/4, 4/)
!sub_orbsym_array1 =  transpose(reshape([0, 0, 0, 0, 5, 6, 7, 8, 0, 0, 0, 0, 3, 4, 9, 10],[4,4]))
sub_orbsym_array1 =  transpose(reshape([0, 0, 0, 0, 7, 8, 0, 0, 9, 10, 0, 0, 3, 4, 5, 6],[4,4]))
flgst_py = 1
total_atoms = 2
niao_py = 2
active_atm_num = 2
output_folder = '/home/sourav/Desktop/SYMM_VB_SOFTWARE/Python_interface/C5H5_output'

endif
if(flg.eq.4)then
geometry_unit = 'Bohr'
nao_py = 8
nae_py = 8
nmul = 1
output_file_name = 'SF4'
chinst = 1
symm_py = 1
asymm_py = 0
!symtype_py = 0
set_order_py = 0
nset_py = 0
mout_py = 1
ovlp_py = 0
itb_py = 1
nnb_py = 2
syb_py = 0
mnbond_py = 0
radical_py = 0
nmbond_py = 0
print*,'i am here3'

allocate(sub_symat_array(2))
allocate(sub_coordx(2)) 
allocate(sub_coordy(2)) 
allocate(sub_coordz(2)) 
allocate(sub_symatno_array(2))
allocate(sub_atoset_array1(2, 4))
allocate(sub_norbsym_array(4))
allocate(sub_active_atom_array(2))
allocate(sub_atn_array(2))
allocate(sub_orbsym_array1(4, 4))

sub_symat_array => symat_array(1:5)
sub_coordx => coordx(1:5)
sub_coordy => coordy(1:5)
sub_coordz => coordz(1:5)
sub_symatno_array => symatno_array(1:5)
sub_atoset_array1 => atoset_array(1:4, 1:5)
sub_norbsym_array => norbsym_array(1:4)
sub_active_atom_array => active_atom_array(1:5)
sub_atn_array => atn_array(1:5)
sub_orbsym_array1 => orbsym_array(1:8 ,1:4)

sub_symat_array = (/'S', 'F', 'F', 'F', 'F'/)
sub_coordx = (/0.0, 0.0, 0.0, -1.23, 1.23/)
sub_coordy = (/0.0, 1.68, -1.68, 0.0, 0.0/)
sub_coordz = (/0.39, 0.27, 0.27, -0.61, -0.61/)
sub_symatno_array = (/16, 9, 9, 9, 9/)
sub_atoset_array1 = transpose(reshape([23, 24, 25, 26, 27, 0, 0, 0, 28, 0, 0, 0, 29, 0, 0, 0, 30, 0, 0, 0],[4,5]))
sub_norbsym_array = (/0, 0, 0, 8/)
sub_active_atom_array = (/1, 2, 3, 4, 5/)
sub_atn_array = (/4, 1, 1, 1, 1/)
!sub_orbsym_array1 =  transpose(reshape([0, 0, 0, 0, 5, 6, 7, 8, 0, 0, 0, 0, 3, 4, 9, 10],[4,4]))
sub_orbsym_array1 = transpose(reshape([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,23,24,25,26,27,28,29,30],[8,4]))
flgst_py = 2
total_atoms = 5
niao_py = 22
active_atm_num = 5
output_folder = '/home/sourav/Desktop/SYMM_VB_SOFTWARE/Python_interface/SF4_output'

endif
if(flg.eq.5)then
geometry_unit = 'Bohr'
nao_py = 8
nae_py = 8
nmul = 1
output_file_name = 'CN+'
chinst = 1
symm_py = 1
asymm_py = 0
!symtype_py = 0
set_order_py = 0
nset_py = 0
mout_py = 1
ovlp_py = 0
itb_py = 1
nnb_py = 0
syb_py = 2
mnbond_py = 0
radical_py = 0
nmbond_py = 0
print*,'i am here3'

allocate(sub_symat_array(2))
allocate(sub_coordx(2)) 
allocate(sub_coordy(2)) 
allocate(sub_coordz(2)) 
allocate(sub_symatno_array(2))
allocate(sub_atoset_array1(2, 4))
allocate(sub_norbsym_array(4))
allocate(sub_active_atom_array(2))
allocate(sub_atn_array(2))
allocate(sub_orbsym_array1(4, 4))

sub_symat_array => symat_array(1:2)
sub_coordx => coordx(1:2)
sub_coordy => coordy(1:2)
sub_coordz => coordz(1:2)
sub_symatno_array => symatno_array(1:2)
sub_atoset_array1 => atoset_array(1:2, 1:4)
sub_norbsym_array => norbsym_array(1:3)
sub_active_atom_array => active_atom_array(1:2)
sub_atn_array => atn_array(1:2)
sub_orbsym_array1 => orbsym_array(1:4 ,1:4)

sub_symat_array = (/'C', 'N'/)
sub_coordx = (/0.0, 0.0/)
sub_coordy = (/0.0, 0.0/)
sub_coordz = (/-1.268, 1.087/)
sub_symatno_array = (/6, 7/)
sub_atoset_array1 = transpose(reshape([3, 5, 7, 9, 4, 6, 8, 10],[4,2]))
sub_norbsym_array = (/0, 2, 2, 4/)
sub_active_atom_array = (/1, 2/)
sub_atn_array = (/4, 4/)
!sub_orbsym_array1 =  transpose(reshape([0, 0, 0, 0, 5, 6, 7, 8, 0, 0, 0, 0, 3, 4, 9, 10],[4,4]))
sub_orbsym_array1 =  transpose(reshape([0, 0, 0, 0, 7, 8, 0, 0, 9, 10, 0, 0, 3, 4, 5, 6],[4,4]))
flgst_py = 2
total_atoms = 2
niao_py = 2
active_atm_num = 2
output_folder = '/home/sourav/Desktop/SYMM_VB_SOFTWARE/Python_interface/CN+_output'

endif
print*,'i am here4'

print*,'symat_array',symat_array
print*,'coordx',coordx
print*,'coordy',coordy
print*,'coordz',coordz
print*,'symatno_array',symatno_array
do i=1, 10
print*,'atoset_array',(atoset_array(i, j), j=1,10)
enddo
print*,'norbsym_array',norbsym_array
print*,'active_atom_array',active_atom_array
print*,'atn_array',atn_array
do i=1, 10
print*,'orbsym_array',(orbsym_array(i, j), j= 1,10)
enddo
print*,'i am here5'
call get_ctrl_inputs(geometry_unit, nao_py, nae_py, nmul, output_file_name, &
                chinst, symm_py, asymm_py, set_order_py, nset_py, mout_py, ovlp_py, itb_py, nnb_py,&
                syb_py, mnbond_py, radical_py, nmbond_py, symat_array, coordx, coordy, coordz, symatno_array,&
                atoset_array, norbsym_array, active_atom_array, atn_array, orbsym_array, flgst_py, total_atoms,&
                niao_py, active_atm_num, output_folder)

stop
end program main
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
