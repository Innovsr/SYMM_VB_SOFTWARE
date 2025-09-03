!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine geocal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat_mod
use commondat1_mod
use coordinates_mod
implicit none

integer:: i,j,ij,j1,k,ii1,i1,i2,i3,i4,k1,k2,l,l1,l2,l3,nland,nncatm,nnatominact,ncatm,a,factorial
real*8::least,bond_dist,ddist,kk
integer, allocatable::active_orb2(:),active_orb1(:),nnat_bond_inact(:,:),biasmat(:)
integer, allocatable::atm_nb_sym(:),nisland(:),islands(:,:),nnmat_act(:,:),nnmat_inact(:,:)
!real*8::coordx(100),coordy(100),coordz(100)
real*8, allocatable::island_mat(:,:),dist_rel_mat(:,:)
integer, allocatable::active_atoms_sort(:),nnact_1(:),nnact(:)
!character(len=5)::line

print*,'enter geocal', nactorb,atom,tot_atom
a = (tot_atom*tot_atom-1)/2
allocate(nnmat_act(a,2))
allocate(nnmat_inact(a,2))
allocate(nnat_bond_inact(tot_atom,2))
allocate(active_orbs(nao))
allocate(atm_nb_orbs1(nao))
allocate(atm_nb_orbs(nao))
allocate(dist_nnat(tot_atom))
sig_sym_flg=0

print*,'geocal1'
l3=0
l2=0
do k1=1,tot_atom
  l1=0
  if(atoset(k1,1).eq.0) cycle
    l2=l2+1
  do k2=1,20
    if(atoset(k1,k2).ne.0)then
      l1=l1+1
      l3=l3+1
      active_orbs(l3)=atoset(k1,k2)
      atm_nb_orbs(l3)=k1
    endif
  enddo
  atn(k1)=l1
enddo
atom=l2
nactorb=l3

!do i=1,nactorb
!  active_orb2(i)=active_orbs(i)
!  atm_nb_orbs1(i)=atm_nb_orbs(i)
!enddo

print*,'atm_nb_orbs',nao,(atm_nb_orbs(i),i=1,nao)
allocate(active_orb2(nactorb))
active_orb2=active_orbs
atm_nb_orbs1=atm_nb_orbs

print*,'geocal2'
!print*,'active_orbs',(active_orbs(i),i=1,nactorb)
!print*,'atm_nb_orbs',(atm_nb_orbs(i),i=1,nactorb)
!print*,'active_orb2',(active_orb2(i),i=1,nactorb)
!print*,'atm_nb_orbs1',(atm_nb_orbs1(i),i=1,nactorb)
!print*, atom, nactorb
!stop

allocate(active_orb1(nactorb))
k1=0
k=0
do
  least=100
  loop2:do i=1,nactorb
    do j=1,k1
      if(active_orb1(j).eq.i) then
        cycle loop2
      endif
    enddo
    if(active_orbs(i).lt.least)then
      least=active_orbs(i)
      k=i
    endif
  enddo loop2
  
  k1=k1+1
  active_orb1(k1)=k
  if(k1.eq.nactorb) exit
enddo

print*,'geocal3'
do i=1,nactorb
  active_orbs(i)=0
  atm_nb_orbs(i)=0
enddo
do i=1,nactorb
  active_orbs(i)=active_orb2(active_orb1(i))
  atm_nb_orbs(i)=atm_nb_orbs1(active_orb1(i))
enddo


allocate(atm_nb_sym(nactorb))
!! atm_nb_sym() store the symmetry property of each active orbitals
!! mentioned by num 1, 2, 3, 4 represents px, py, pz, s respectively 

do i=1,nactorb
  do j=1,4
    do k=1,norbsym(j)
      if(norbsym(j).ne.0)then
        if(active_orbs(i).eq.orbsym(j,k))then
          atm_nb_sym(i)=j
        endif
      endif
    enddo
  enddo
enddo

print*,'atm_nb_sym',(atm_nb_sym(i),i=1,nactorb)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! if we take all 'pi' orbitals have non-zero overlapping!!!!

print*,'geocal4'
do i=1,nactorb
  if(atm_nb_sym(i).lt.4)atm_nb_sym(i)=1
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!990 format(a,x,I2,x,a,x,I2,x,a)

print*,'active_atoms',(active_atoms(i),i=1,atom)
do i = 1, tot_atom
print*,'coordx:',coordx(i),'coordy:',coordy(i),'coordz:',coordz(i)
enddo

allocate(dist_rel_mat(tot_atom, tot_atom))
allocate(dist_mat(tot_atom, tot_atom))
allocate(dist_act_rel_mat(tot_atom, tot_atom))

j=0
do ij=1,tot_atom-1
  do ii1=ij+1,tot_atom
    ddist=sqrt((coordx(ij)-coordx(ii1))**2.0+(coordy(ij)-coordy(ii1))**2.0+&
    (coordz(ij)-coordz(ii1))**2.0)
    dist_mat(ij,ii1)=ddist
    bond_dist=(at_covrad(all_at_num(ij))+at_covrad(all_at_num(ii1)))/100.0
    if(ddist.le.bond_dist)then
      k1=0
      do i2=1,atom
        if(ij.eq.active_atoms(i2))k1=k1+1
        if(ii1.eq.active_atoms(i2))k1=k1+1
      enddo
      if(k1.ne.2)then
        dist_act_rel_mat(ij,ii1)=0.0
      else
        dist_act_rel_mat(ij,ii1)=1.0
      endif
      dist_rel_mat(ij,ii1)=1.0
    else
      k1=0
      do i2=1,atom
        if(ij.eq.active_atoms(i2))k1=k1+1
        if(ii1.eq.active_atoms(i2))k1=k1+1
      enddo
      if(k1.ne.2)then
        dist_act_rel_mat(ij,ii1)=0.0
      else
        dist_act_rel_mat(ij,ii1)=ddist/bond_dist
      endif
        dist_rel_mat(ij,ii1)=ddist/bond_dist
    endif
  enddo
enddo


do i=1,tot_atom
  dist_act_rel_mat(i,i)=0.0
  dist_rel_mat(i,i)=0.0
  dist_mat(i,i)=0.0
enddo

do i=1,tot_atom-1
  do i1=i+1,tot_atom
    dist_act_rel_mat(i1,i)=dist_act_rel_mat(i,i1)
    dist_rel_mat(i1,i)=dist_rel_mat(i,i1)
    dist_mat(i1,i)=dist_mat(i,i1)
  enddo
enddo

do i = 1, tot_atom
print*,'dist_act_rel_mat',(dist_act_rel_mat(i, j), j=1, tot_atom)
enddo
do i = 1, tot_atom
print*,'dist_rel_mat',(dist_rel_mat(i, j), j=1, tot_atom)
enddo
do i = 1, tot_atom
print*,'dist_mat',(dist_mat(i, j), j=1, tot_atom)
enddo


!!!!!!!!!!!! below active_atoms() moved to the regular order as
!given in INFO ... for dist matrix formation it was (may be) necessary
! to change it in the order of active orbitals.
! 'atom'= active atoms

allocate(active_atoms_sort(atom))
j=0

print*,'sourav1'
!413 k=100
do
  k=100
  loop4:do i=1,atom
    do i1=1,j
      if(active_atoms_sort(i1).eq.active_atoms(i)) then
        cycle loop4
      endif
    enddo
    if(k.gt.active_atoms(i))k=active_atoms(i)
  enddo loop4
  j=j+1
  active_atoms_sort(j)=k
  if(j.eq.atom) exit
enddo

print*,'sourav2'
do i=1,atom
  active_atoms(i)=active_atoms_sort(i)
enddo
print*,'active_atoms',(active_atoms(i),i=1,atom)
deallocate(active_atoms_sort)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

k=0
do i=1,tot_atom
  least=1000.0
  do i1=1,tot_atom
    if(dist_act_rel_mat(i,i1).lt.least.and.dist_act_rel_mat(i,i1).ne.0.0)least=dist_act_rel_mat(i,i1)
  enddo
  do i1=1,tot_atom
    k1=0
    do i2=1,atom
      if(i.eq.active_atoms(i2))k1=k1+1
      if(i1.eq.active_atoms(i2))k1=k1+1
    enddo
    if(k1.eq.2.and.dist_act_rel_mat(i,i1).eq.least)then
      k=k+1
      nnmat_act(k,1)=i
      nnmat_act(k,2)=i1
      print*,'nnmat_act',k,nnmat_act(k,1),nnmat_act(k,2)
    endif
  enddo
enddo

print*,'sourav3'
allocate(nnact(k))
nnact = 0

allocate(nnact_1(k))
nnact_1 = 0
!print*,'kkk',k, atom
!do i = 1, k
!print*,'nnmat_act',(nnmat_act(i,j),j=1,2)
!enddo

do i=1,k
   nnact(i)=nnmat_act(i,1)**2 + nnmat_act(i,2)**2
enddo


print*,'sourav4'
j1=0
loop5:do i=1,k
  do i1=1,j1
    if(nnact(i).eq.nnact(nnact_1(i1))) then
      cycle loop5
    endif
  enddo
  j1=j1+1
  nnact_1(j1)=i
  print*,'nnact_1',nnact_1(j1)
enddo loop5
nnnatom=j1

print*,'nnnatom',nnnatom
allocate(nnat_bond(nnnatom, 2))
do i=1,nnnatom
  do j=1,2
    nnat_bond(i,j)=nnmat_act(nnact_1(i),j)
  enddo
enddo


print*,'sourav5'
k2=0
do i=1,tot_atom
  least=1000.0
  do i1=1,tot_atom
    if(dist_rel_mat(i,i1).lt.least.and.dist_rel_mat(i,i1).ne.0.0)least=dist_rel_mat(i,i1)
  enddo
  do i1=1,tot_atom
  k1=0
    do i2=1,atom
      if(i.eq.active_atoms(i2))k1=k1+1
      if(i1.eq.active_atoms(i2))k1=k1+1
    enddo

    if(k1.ne.2.and.dist_rel_mat(i,i1).eq.least)then
      k2=k2+1
      nnmat_inact(k2,1)=i
      nnmat_inact(k2,2)=i1
    endif
  enddo
enddo

deallocate(dist_rel_mat)
nnact = 0
nnact_1 = 0

do i=1,k2
  nnact(i)=nnmat_inact(i,1)**2+nnmat_inact(i,2)**2
enddo

print*,'sourav6'
j1=0
loop13:do i=1,k2
  do i1=1,j1
    if(nnact(i).eq.nnact(nnact_1(i1))) cycle loop13
  enddo
  j1=j1+1
  nnact_1(j1)=i
enddo loop13
nnatominact=j1
do i=1,nnatominact
  do j=1,2
    nnat_bond_inact(i,j)=nnmat_inact(nnact_1(i),j)
  enddo
enddo

!do i=1,nnnatom
!print*,'read_info:nnat_bond',nnnatom,(nnat_bond(i,i1),i1=1,2)
!enddo
!print*,'nnatominact',nnatominact
!do i=1,nnatominact
!print*,'read_info:nnat_bond_inact',nnatominact,(nnat_bond_inact(i,i1),i1=1,2)
!enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Routine to produce the the inactive bias to active atooms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate(biasmat(tot_atom))
allocate(biasval(tot_atom))
biasval=0.0

print*,'sourav7'
do i=1,atom
  biasmat=0
  k=0
  k=k+1
  biasmat(k)=active_atoms(i)
  loop10:do i1=1,nnatominact
    do j=1,2
      if(active_atoms(i).eq.nnat_bond_inact(i1,j)) then
        exit
      else
        cycle loop10
      endif
    enddo
    if(j.eq.1)then
      k=k+1
      biasmat(k)=nnat_bond_inact(i1,2)
      biasval(i)=biasval(i)+all_at_num(nnat_bond_inact(i1,2))+dist_mat(active_atoms(i),nnat_bond_inact(i1,2))
    else
      k=k+1
      biasmat(k)=nnat_bond_inact(i1,1)
      biasval(i)=biasval(i)+all_at_num(nnat_bond_inact(i1,1))+dist_mat(active_atoms(i),nnat_bond_inact(i1,1))
    endif
  enddo loop10
  ncatm=k
  nncatm=0
  do
    do i2=1+nncatm,ncatm
      loop9:do i1=1,nnatominact
        do j=1,2
          if(biasmat(i2).eq.nnat_bond_inact(i1,j))then
            exit
          else
            cycle loop9
          endif
        enddo
        if(j.eq.1)then
          do i3=1,k
            if(biasmat(i3).eq.nnat_bond_inact(i1,2)) cycle loop9
          enddo
          k=k+1
          biasmat(k)=nnat_bond_inact(i1,2)
          biasval(i)=biasval(i)+all_at_num(nnat_bond_inact(i1,2))+dist_mat(active_atoms(i),nnat_bond_inact(i1,2))
        else
          do i3=1,k
            if(biasmat(i3).eq.nnat_bond_inact(i1,1)) cycle loop9
          enddo
          k=k+1
          biasmat(k)=nnat_bond_inact(i1,1)
          biasval(i)=biasval(i)+all_at_num(nnat_bond_inact(i1,1))+dist_mat(active_atoms(i),nnat_bond_inact(i1,1))
        endif
      enddo loop9
    enddo
    nncatm=ncatm
    ncatm=k
    if(nncatm.eq.ncatm) exit
  enddo
enddo

print*,'sourav8'
deallocate(biasmat)
deallocate(all_at_num)
!do i=1,atom
!print*,'biasval_active',biasval(i)
!enddo
!do i=1,nnnatom
!print*,'nnat_bond_ild',(nnat_bond(i,j),j=1,2)
!enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Routine to take care of the Islands of active atooms !!!!
!!!!!!!!! Islands means the cluster of atoms !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate(nisland(tot_atom))
allocate(islands(nnnatom, tot_atom))

!do i1=1,50
!  do i=1,20

islands=0
nisland=0

!  enddo
!enddo

i4=0
loop6:do i=1,atom
  if(i.ne.1) then
    do i1=1,i4
      do j=1,nisland(i1)
        if(active_atoms(i).eq.islands(i1,j)) cycle loop6
      enddo
    enddo
  endif
  k=0
  k=k+1
  i4=i4+1
  islands(i4,k)=active_atoms(i)
  loop7:do i1=1,nnnatom
    do j=1,2
      if(active_atoms(i).eq.nnat_bond(i1,j)) then
        exit
      else
        cycle loop7
      endif
    enddo
    if(j.eq.1)then
      k=k+1
      islands(i4,k)=nnat_bond(i1,2)
    else
      k=k+1
      islands(i4,k)=nnat_bond(i1,1)
    endif
  enddo loop7
  ncatm=k
  nncatm=0
  do
    do i2=1+nncatm,ncatm
      loop8:do i1=1,nnnatom
        do j=1,2
          if(islands(i4,i2).eq.nnat_bond(i1,j)) then
            exit
          else
            cycle loop8
          endif
        enddo
        if(j.eq.1)then
          do i3=1,k
            if(islands(i4,i3).eq.nnat_bond(i1,2)) cycle loop8
          enddo
          k=k+1
          islands(i4,k)=nnat_bond(i1,2)
        else
          do i3=1,k
            if(islands(i4,i3).eq.nnat_bond(i1,1)) cycle loop8
          enddo
          k=k+1
          islands(i4,k)=nnat_bond(i1,1)
        endif
      enddo loop8
    enddo
    nncatm=ncatm
    ncatm=k
    if(nncatm.eq.ncatm) exit
  enddo
  nisland(i4)=k
enddo loop6
nland=i4

print*,'sourav9'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!lowest distance between two islands
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(island_mat(nland, nland))

if(nland.gt.1)then
  j1=nnnatom
  do i=1,nland
    island_mat(i,i)=0.0
  enddo
  do i=1,nland-1
    do j=i+1,nland
      least=100.0
      do i1=1,nisland(i)
        do i2=1,nisland(j)
          if(dist_act_rel_mat(islands(i,i1),islands(j,i2)).lt.least.and.dist_act_rel_mat&
            (islands(i,i1),islands(j,i2)).ne.0.0)least=dist_act_rel_mat(islands(i,i1),islands(j,i2))
        enddo
      enddo
      island_mat(i,j)=least
      island_mat(j,i)=least
    enddo
  enddo

!do i=1,nland
!print*,'island_mat(i,j)',(island_mat(i,j),j=1,nland)
!enddo

print*,'sourav10'
  do i=1,nland
    least=100.0
    do j=1,nland
      if(island_mat(i,j).lt.least.and.island_mat(i,j).ne.0.0)then
        least=island_mat(i,j)
        k1=i
        k2=j
        island_mat(j,i)=0.0
      endif
    enddo

    dO i1=1,nisland(k1)
      do i2=1,nisland(k2)
        if(dist_act_rel_mat(islands(k1,i1),islands(k2,i2)).eq.least)then
          j1=j1+1
          nnat_bond(j1,1)=islands(k1,i1)
          nnat_bond(j1,2)=islands(k2,i2)
        endif
      enddo
    enddo
  enddo
  nnnatom=j1
endif

!do i=1,nnnatom
!print*,'nnat_bond',(nnat_bond(i,j),j=1,2)
!enddo

print*,'atm_nb_orbs',nactorb,(atm_nb_orbs(i),i=1,nactorb)
print*,'atm_nb_sym',nactorb,(atm_nb_sym(i),i=1,nactorb)
print*,'sourav10.5',tot_atom, atom

print*,'sourav11'
do i=1,atom
  kk=0.0
  do j=1,nnnatom
    do k=1,2
      if(active_atoms(i).eq.nnat_bond(j,k))then
        if(k.eq.1)l=2
        if(k.eq.2)l=1
        kk=kk+dist_mat(active_atoms(i),nnat_bond(j,l))
      endif
    enddo
  enddo
  dist_nnat(i)=kk
  print*,'dist_nnat',i,kk
enddo

print*,'dist_nnat',(dist_nnat(i),i=1,atom),nactorb
do i=1,nactorb-1
  do j=i+1,nactorb
    if(atm_nb_orbs(i).eq.atm_nb_orbs(j))then
      if(atm_nb_sym(i).eq.4.and.atm_nb_sym(j).eq.4)then
        sig_sym_flg=1
        return
      endif
    endif
  enddo
enddo

print*,'sourav12'
deallocate(nnmat_act)
deallocate(nnmat_inact)
deallocate(atm_nb_sym)
deallocate(atm_nb_orbs1)
print*,'exit geocal'
return
end subroutine geocal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
