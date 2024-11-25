!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine geocal(coordx,coordy,coordz)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
use commondat1
implicit none

integer:: i,j,ij,j1,j2,j3,k,ii1,i1,i2,i3,i4,k1,k2,k3,k4,k9,k10,nland,nisland(50),islands(20,50),nncatm,ncatm,&
nnmat(100,100),nnact(100),nnact_1(100),active_orb2(20),active_orb1(20),l,l1,l2,l3,l4,val_st(20)&
,ind(100),act_at_num(20),val_state_num(100),n,n1,n2,n3,n4&
,val_state_num1(100),linenumi,linenumj,active_atoms_sort(10000)
real*8::least,bond_dist,ddist,kk,island_mat(20,20)
real*8::coordx(100),coordy(100),coordz(100),oval,a(4)
character(len=5)::act_orb_typ(20,10),line

print*,'enter geocal,key_frag',key_frag
sig_sym_flg=0

if(key_frag.eq.1)then
  do i=1,200
    do j=1,20
      atoset(i,j)=0
    enddo
  enddo
  k3=0
  do i=1,tot_atom
    ind(i)=1
  enddo
  loop1:do k=orbsl-nao+1,orbsl
    if(frag_atn(orbs(k,1)).gt.1)then
      write(*,990)'The orbital',k-1,'is a fragment orbital and we took the atom number'&
      ,frag_at(orbs(k,1)+1,1),'as the central atom.'
      write(*,*)'If it is not right please put the central atom number at first place &
      in the row in the $frag section and run again.'
    endif
    if(atoset(frag_at(orbs(k,1)+1,1),ind(frag_at(orbs(k,1)+1,1))).ne.0)&
    ind(frag_at(orbs(k,1)+1,1))=ind(frag_at(orbs(k,1)+1,1))+1
    if(num_frag_cntr.ne.0)then
      do j=1,frag_atn(orbs(k,1)+1)
        do i=1,num_frag_cntr
          if(frag_cntr(i).eq.frag_at(orbs(k,1)+1,j))then
            frag_at(orbs(k,1)+1,1)=frag_cntr(i)
          endif
        enddo
      enddo
    endif

    atoset(frag_at(orbs(k,1)+1,1),ind(frag_at(orbs(k,1)+1,1)))=k-1

    do i=1,k3
      if(active_atoms(i).eq.frag_at(orbs(k,1)+1,1)) cycle loop1
    enddo
    k3=k3+1
    active_atoms(k3)=frag_at(orbs(k,1)+1,1)

  enddo loop1
  do i=1,k3
  atn(active_atoms(i))=ind(active_atoms(i))
  act_at_num(i)=symatno(active_atoms(i))
  val_state_num(i)=valence_state(act_at_num(i))
  val_state_num1(i)=val_state_num(i)
  enddo
  atom=k3
endif

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

do i=1,nactorb
  active_orb2(i)=active_orbs(i)
  atm_nb_orbs1(i)=atm_nb_orbs(i)
enddo

k1=0
k=0
436 least=100
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
if(k1.lt.nactorb)goto 436

do i=1,nactorb
  active_orbs(i)=0
  atm_nb_orbs(i)=0
enddo
do i=1,nactorb
  active_orbs(i)=active_orb2(active_orb1(i))
  atm_nb_orbs(i)=atm_nb_orbs1(active_orb1(i))
enddo


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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! if we take all 'pi' orbitals have non-zero overlapping!!!!

do i=1,nactorb
  if(atm_nb_sym(i).lt.4)atm_nb_sym(i)=1
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
990 format(a,x,I2,x,a,x,I2,x,a)

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

!!!!!!!!!!!! below active_atoms() moved to the regular order as
!given in INFO ... for dist matrix formation it was (may be) necessary
! to change it in the order of active orbitals.

j=0

413 k=100
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
if(j.ne.atom)goto 413

do i=1,atom
  active_atoms(i)=active_atoms_sort(i)
enddo


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
    endif
  enddo
enddo

do i=1,k
   nnact(i)=nnmat_act(i,1)**2.0+nnmat_act(i,2)**2.0
enddo


j1=0
loop5:do i=1,k
  do i1=1,j1
    if(nnact(i).eq.nnact(nnact_1(i1))) then
      cycle loop5
    endif
  enddo
  j1=j1+1
  nnact_1(j1)=i
enddo loop5
nnnatom=j1
do i=1,nnnatom
  do j=1,2
    nnat_bond(i,j)=nnmat_act(nnact_1(i),j)
  enddo
enddo


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
do i=1,100
nnact(i)=0
nnact_1(i)=0
enddo

do i=1,k2
nnact(i)=nnmat_inact(i,1)**2.0+nnmat_inact(i,2)**2.0
enddo

j1=0
do i=1,k2
do i1=1,j1
if(nnact(i).eq.nnact(nnact_1(i1)))goto 234
enddo
j1=j1+1
nnact_1(j1)=i
234 enddo
nnatominact=j1
do i=1,nnatominact
do j=1,2
nnat_bond_inact(i,j)=nnmat_inact(nnact_1(i),j)
enddo
enddo

do i=1,nnnatom
print*,'read_info:nnat_bond',nnnatom,(nnat_bond(i,i1),i1=1,2)
enddo
print*,'nnatominact',nnatominact
do i=1,nnatominact
print*,'read_info:nnat_bond_inact',nnatominact,(nnat_bond_inact(i,i1),i1=1,2)
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
n=0
n1=0
do i=1,tot_atom
if(atoset(i,1).eq.0)goto 189
n1=n1+1
do i3=1,88
if(symat(i).eq.at_list(i3))then
n=n+1
actv_atom(n)=i
at_num(n)=i3
if(i3.le.54)val_st(n)=valence_state(i3)
endif
enddo
189 enddo
if(n1.ne.n)then
n=0
do i=1,tot_atom
if(atoset(i,1).eq.0)goto 190
do i3=1,88
if(symat(i).eq.at_list_bold(i3))then
n=n+1
actv_atom(n)=i
at_num(n)=i3
if(i3.le.54)val_st(n)=valence_state(i3)
endif
enddo
190 enddo

endif

do i=1,n
print*,'at_num:geo',at_num(i),actv_atom(i),val_st(i)
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Routine to produce the the inactive bias to active atooms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do i1=1,20
biasval(i1)=0.0
enddo

do i=1,atom
  do i1=1,100
    biasmat(i1)=0
  enddo
k=0
k=k+1
biasmat(k)=active_atoms(i)
  do i1=1,nnatominact
    do j=1,2
      if(active_atoms(i).eq.nnat_bond_inact(i1,j))goto 235
    enddo
  goto 239
235 if(j.eq.1)then
     k=k+1
     biasmat(k)=nnat_bond_inact(i1,2)
     biasval(i)=biasval(i)+all_at_num(nnat_bond_inact(i1,2))+dist_mat(active_atoms(i),nnat_bond_inact(i1,2))
    else
     k=k+1
     biasmat(k)=nnat_bond_inact(i1,1)
     biasval(i)=biasval(i)+all_at_num(nnat_bond_inact(i1,1))+dist_mat(active_atoms(i),nnat_bond_inact(i1,1))
    endif
239 enddo
ncatm=k
nncatm=0
237 do i2=1+nncatm,ncatm
      do i1=1,nnatominact
        do j=1,2
          if(biasmat(i2).eq.nnat_bond_inact(i1,j))goto 236
        enddo
      goto 238
236 if(j.eq.1)then
        do i3=1,k
          if(biasmat(i3).eq.nnat_bond_inact(i1,2))goto 238
        enddo
        k=k+1
        biasmat(k)=nnat_bond_inact(i1,2)
     biasval(i)=biasval(i)+all_at_num(nnat_bond_inact(i1,2))+dist_mat(active_atoms(i),nnat_bond_inact(i1,2))
     else
         do i3=1,k
           if(biasmat(i3).eq.nnat_bond_inact(i1,1))goto 238
         enddo
         k=k+1
         biasmat(k)=nnat_bond_inact(i1,1)
     biasval(i)=biasval(i)+all_at_num(nnat_bond_inact(i1,1))+dist_mat(active_atoms(i),nnat_bond_inact(i1,1))
     endif
238 enddo
  enddo
nncatm=ncatm
ncatm=k
if(nncatm.ne.ncatm)goto 237
enddo

do i=1,atom
print*,'biasval_active',biasval(i)
enddo
do i=1,nnnatom
print*,'nnat_bond_ild',(nnat_bond(i,j),j=1,2)
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Routine to take care of the Islands of active atooms !!!!
!!!!!!!!! Islands means the cluster of atoms !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do i1=1,50
do i=1,20
islands(i,i1)=0
enddo
enddo

i4=0
do i=1,atom
if(i.eq.1)goto 246
do i1=1,i4
do j=1,nisland(i1)
if(active_atoms(i).eq.islands(i1,j))goto 245
enddo
enddo
246 k=0
k=k+1
i4=i4+1
islands(i4,k)=active_atoms(i)
do i1=1,nnnatom
do j=1,2
if(active_atoms(i).eq.nnat_bond(i1,j))goto 240
enddo
goto 241
240 if(j.eq.1)then
k=k+1
islands(i4,k)=nnat_bond(i1,2)
else
k=k+1
islands(i4,k)=nnat_bond(i1,1)
endif
241 enddo
ncatm=k
nncatm=0
242 do i2=1+nncatm,ncatm
do i1=1,nnnatom
do j=1,2
if(islands(i4,i2).eq.nnat_bond(i1,j))goto 243
enddo
goto 244
243 if(j.eq.1)then
do i3=1,k
if(islands(i4,i3).eq.nnat_bond(i1,2))goto 244
enddo
k=k+1
islands(i4,k)=nnat_bond(i1,2)
else
do i3=1,k
if(islands(i4,i3).eq.nnat_bond(i1,1))goto 244
enddo
k=k+1
islands(i4,k)=nnat_bond(i1,1)
endif
244 enddo
enddo
nncatm=ncatm
ncatm=k
if(nncatm.ne.ncatm)goto 242
nisland(i4)=k
245 enddo
nland=i4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!lowest distance between two islands
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

do i=1,nland
print*,'island_mat(i,j)',(island_mat(i,j),j=1,nland)
enddo

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

do i=1,nnnatom
print*,'nnat_bond',(nnat_bond(i,j),j=1,2)
enddo

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

do i=1,nactorb
do j=1,nactorb
if(atm_nb_orbs(i).eq.atm_nb_orbs(j))then
if(atm_nb_sym(i).eq.4.and.atm_nb_sym(j).eq.4)then
sig_sym_flg=1
goto 550
endif
endif
enddo
enddo

print*,'exit geocal'
550 return
end subroutine geocal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
