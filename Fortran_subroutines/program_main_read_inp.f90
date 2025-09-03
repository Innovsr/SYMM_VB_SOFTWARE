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

integer:: nao, nae, nmul, flgst, total_atoms, niao, active_atm_num,i, j, pos
character(len = 200)::output_file_name, inputfilename
character(len = 300)::output_folder,line, cwd
character(len = 5)::lname, geometry_unit
character(len = 3)::str
integer::chinst, symm, asymm, set_order, nset, mout, ovlp, itb
integer::nnb, syb, mnbond, radical, nmbond, flg!, symtype_py
character(len = 30)::orb_typ(100),orb,lowercase
real*8, target::coordx(100), coordy(100), coordz(100), symatno_array(100)
character(len=5), target::symat_array(100)
integer, target::atoset_array(200, 20), norbsym_array(50), atn_array(200)
integer, target::orbsym_array(20, 20), active_atom_array(30)
integer::lineno,natoms, mol_info_line, mol_geo_line, act_orbs_line, keywd_line,at_index(100)
integer::minval, max_atom, l, ios, read_flg, k,kk,ierr
integer::line_num, geo_line, geo_line1, geo_line2, info_flg, geo_flg, orb_flg, key_flg
logical::fileexists


!print*,'i am here1'
symat_array = ''
atoset_array = 0
!print*,'i am here1'
norbsym_array = 0
atn_array = 0
orbsym_array = 0
active_atom_array = 0
coordx = 0.0 
coordy = 0.0 
coordz = 0.0 
symatno_array = 0.0
read_flg = 0
!print*,'i am here2'

call getarg(1,inputfilename)
l=len(TRIM(inputfilename))
lname=inputfilename(l-3:l)

if(lname.ne.'.cis') then
print*,'please put the .xmi file as input'
stop
endif
INQUIRE(FILE=TRIM(inputfilename),EXIST=fileexists)
IF (fileexists) THEN
open(unit=21,file=TRIM(inputfilename),status='old')
ELSE
PRINT*,'SORRY This input file does not exist or you may not provide the &
filename at all'
stop
ENDIF

lineno = 0
mol_info_line = 0
mol_geo_line = 0
act_orbs_line = 0
keywd_line = 0
do
    read(21, '(A)', iostat=ios) line
    if (ios /= 0) exit
    lineno = lineno + 1
    if (INDEX(line, '$info') > 0) then
            mol_info_line = lineno
    end if
    if (INDEX(line, '$geo') > 0) then
            mol_geo_line = lineno
    end if
    if (INDEX(line, '$orbs') > 0) then
            act_orbs_line = lineno
    end if
    if (INDEX(line, '$keywd') > 0) then
            keywd_line = lineno
    end if
end do
rewind(21)

if(mol_info_line == 0) then
   write(*,*)"molecular info is not given or you forget to mark the section as'$info'"
   stop
endif
if(mol_geo_line == 0) then
   write(*,*)"molecular geometry is not given or you forget to mark the section as '$geo'"
   stop
endif
if(act_orbs_line == 0) then
   write(*,*)"Active orbitals info are not given or you forget to mark the section as '$orbs'"
   stop
endif
if(keywd_line == 0) then
   write(*,*)"molecular geometry is not given or you forget to mark the section as '$keywd'"
   stop
endif

   
chinst = 0 !chinst = 0 gives Rumer sets (Default) chinst = 1 provides Chem. Insight sets
symm = 0    !symm = 1 will provide symmetric set symm = 0 will provide asymmetric set
asymm = 0
set_order = 0 ! set_order decide order of symmetric groups 2=smallest set to largest
              ! 3=largest set to smallest
nset = 0 ! number of set to display 0=single set, 1=All same quality sets, 2=all possible sets
mout =1 !if nset=2 the total number of sets can be huge in one single file maximum 75000 sets could printed.
        ! mout provide how many this type of output file will be generated. 1=only one file will 75000 sets
ovlp=0 !ovlp=1 calculate overlap of the structures inside a set and prited bydefault it is zero for now.
itb=1 !intra-bond first priority default
nnb=2 !neighboring atom bonds second priority default
syb=3 !Symmetry breaking bond third priority default
mnbond=0 !User defined bond not included default
radical=0 !User defined radical not included default
nmbond=0
str='cov' !flgst defines covalent (cov)=2, ionic (ion)=3 or all (both)=1 type of structures 'cov' is default

nao=0
nae=0
niao=0
nmul=0
output_file_name=''
read_flg = 0
line_num = 0
geo_line = 0
geo_line1 = 0
geo_line2 = 0
info_flg = 0
geo_flg = 0
orb_flg = 0
key_flg = 0
total_atoms = 0
geometry_unit = ''
do 
    read(21, '(A)', iostat=ios) line
    if (ios /= 0) exit
    line_num = line_num+1

  if (line_num > mol_info_line.and.info_flg == 0) then
    pos=index(line,'nao=')
    if (pos>0) then
      read(line(pos+4:),*) nao
    endif

    pos=index(line,'nae=')
    if (pos>0) then
      read(line(pos+4:),*) nae
    endif

    pos=index(line,'niao=')
    if (pos>0) then
      read(line(pos+5:),*) niao
    endif

    pos=index(line,'nmul=')
    if (pos>0) then
      read(line(pos+5:),*) nmul
    endif

    pos=index(line,'name=')
    if (pos>0) then
      read(line(pos+5:),*) output_file_name
    endif
 endif
    if(nao /= 0 .and. nae /= 0 .and. niao /= 0 .and. nmul /= 0 .and. output_file_name /= '')then
            info_flg = 1
    endif
 if(line_num > mol_geo_line .and. geo_flg == 0) then
    pos=index(line,'natoms=')
    if (pos>0) then
      read(line(pos+7:),*) total_atoms
      geo_line1 = line_num
    endif
    pos=index(line,'unit=')
    if (pos>0) then
      read(line(pos+5:),*) geometry_unit
      geo_line2 = line_num
    endif
 endif
 if(total_atoms /= 0 .and. geometry_unit /= '') then
       geo_line=max(geo_line1, geo_line2)
       geo_flg = 1
 endif
 if(line_num > keywd_line) then

    pos=index(line,'chinst=')
    if (pos>0) then
      read(line(pos+7:),*) chinst
    endif

    pos=index(line,'symm')
    if (pos>0) then
      symm = 1
    endif

    pos=index(line,'nset=')
    if (pos>0) then
      read(line(pos+5:),*) nset
    endif

    pos=index(line,'mout=')
    if (pos>0) then
      read(line(pos+5:),*) mout
    endif

    pos=index(line,'set_order=')
    if (pos>0) then
      read(line(pos+10:),*) set_order
    endif

    pos=index(line,'itb=')
    if (pos>0) then
      read(line(pos+4:),*) itb
    endif

    pos=index(line,'nab=')
    if (pos>0) then
      read(line(pos+4:),*) nnb
    endif

    pos=index(line,'sbb=')
    if (pos>0) then
      read(line(pos+4:),*) syb
    endif

    pos=index(line,'udb=')
    if (pos>0) then
      read(line(pos+4:),*) mnbond
    endif

    pos=index(line,'udr=')
    if (pos>0) then
      read(line(pos+4:),*) radical
    endif

    pos=index(line,'str=')
    if (pos>0) then
      read(line(pos+4:),*) str
    endif

 endif
enddo
rewind(21)

if(nao == 0)then
  write(*,*)"number of active orbitals 'nao' are not given"
  stop
endif
if(nae == 0) then
  write(*,*)"number of active electrons 'nae' are not given"
  stop
endif
if(niao == 0) then
  write(*,*)"number of in-active orbitals 'niao' are not given"
  stop
endif
if(nmul == 0) then
  write(*,*)"multiplicity 'nmul' is not given"
  stop
endif
if(output_file_name == '') then
  write(*,*)"name of the molecule is not given"
  stop
endif
if (total_atoms == 0) then
   write(*,*),' "total_atoms" in the $geo section before geometry is not provided'
   stop
endif

do i=1,geo_line
   read(21,*)
enddo
i=0
do
  read(21,'(A)') line
  if(len(trim(line))==0)cycle
  i=i+1
  read(line,*)symat_array(i),symatno_array(i),coordx(i),coordy(i),coordz(i)
!  print*,'coord',symat_array(i),symatno_array(i),coordx(i),coordy(i),coordz(i)
  if(i==total_atoms)exit
enddo
rewind(21)

do i=1,act_orbs_line
   read(21,*)
enddo
i=0
do 
  read(21,'(A)') line
  if(len(trim(line))==0)cycle
  i=i+1
  read(line,*) at_index(i),orb_typ(i)
!  print*,'at_index, orb_type',at_index(i),orb_typ(i)
  if(i==nao)exit
enddo
 max_atom = maxval(at_index)
 kk = 0
 loop1:do i =1, max_atom
      k = 0
      do j = 1, nao
        if (i == at_index(j))then
           k = k+1
           atoset_array(i,k) = j+niao
        endif
      enddo
      atn_array(i)=k
      do j = 1, nao
        if(i == at_index(j))then
          kk=kk+1
          active_atom_array(kk)=i
          cycle loop1
        endif
      enddo
 enddo loop1

active_atm_num = kk

do i = 1, 4
   if (i==1) orb = 'px'
   if (i==2) orb = 'py'
   if (i==3) orb = 'pz'
   if (i==4) orb = 's'
   k = 0
   do j = 1, nao
     pos = index(lowercase(orb_typ(j)),orb)
     if (pos > 0)then
        k = k+1
        orbsym_array(i,k) = j+niao
     endif
   enddo
   norbsym_array(i)=k
enddo
rewind(21)

!do i=1, nao
!print*,'atoset_array',(atoset_array(i, j), j=1,10)
!enddo
!do i=1, 4
!print*,'orbsym_array',(orbsym_array(i, j), j= 1,10)
!enddo
!
!print*,'atn_array',atn_array
!print*,'active_atom_array',active_atom_array
!print*,'norbsym_array',norbsym_array
!
!print*,'keywd_line',keywd_line

if (symm == 0) asymm = 1
if (symm == 1) asymm = 0
if(str=='cov') flgst=2
if(str=='ion') flgst=3
if(str=='all') flgst=1

call getcwd(cwd, ierr)

if (ierr == 0) then
   output_folder = trim(cwd)
else
   print *, "Error: could not get current working directory."
end if

print*,'total_atoms',total_atoms,geometry_unit
print*,'chinst',chinst
print*,'symm',symm
print*,'nset',nset
print*,'mout',mout
print*,'intra atomic bond',itb
print*,'neighbour atomic bond',nnb
print*,'symmetry breaking bond',syb
print*,'user defined bond',mnbond
print*,'user defined radical',radical
print *, "Current working directory: ", trim(cwd)
print*,'symm,asymm',symm,asymm
print*,'structure type',str,flgst
print*,'symat_array',symat_array
print*,'coordx',coordx
print*,'coordy',coordy
print*,'coordz',coordz
print*,'symatno_array',symatno_array
print*,'i am here5'

call get_ctrl_inputs(geometry_unit, nao, nae, nmul, output_file_name, &
                chinst, symm, asymm, set_order, nset, mout, ovlp, itb, nnb,&
                syb, mnbond, radical, nmbond, symat_array, coordx, coordy, coordz, symatno_array,&
                atoset_array, norbsym_array, active_atom_array, atn_array, orbsym_array, flgst, total_atoms,&
                niao, active_atm_num, output_folder)

stop
end program main
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
