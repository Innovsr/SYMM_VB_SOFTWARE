!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! reading INFO file (geometry and AO informations) and calculating distance matrix, bias value etc.!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
use commondat1
implicit none

integer::i,i1,i2,i3,k,k1,k2,k3,l,j,unkwn&
,nbasatom1(20),cnt(20),ufrzorbn,ovlp_basis_num(20,10),MNP,iO
integer::fragorb(30),frag_cent(50),frag_cent1(50),k4,k5,nfragcent,n,k9,k10,linenumi,linenumj
integer::bfiexists,n1,n2,n3,n4,l1,l2,j1,i4,sig(20),l3,l4,act_at_num(20),val_state_num(20) &
,pi1(20),pi2(20),k6,k7,nsig,npi1,npi2,k8,ovlp_at_bas_num(20,10),nbas_slab(100)
real*8::oval,a(4)
!biasval(20),least
character(len=5)::atbas(500),act_orb_typ(20,10),atbas_frz(50),line
logical :: fileexists
integer::atsymset(20,20),nsym,syn(50),at_sym(50)
real*8::coordx(100),coordy(100),coordz(100)

common/ats/atsymset,nsym,syn,at_sym
common/coordinate/coordx,coordy,coordz

print*,'enter read_info'

! reading one body integral file to check for frozen orbitals
ovlp_int=0
INQUIRE(FILE='x1e.int',EXIST=fileexists)
IF (fileexists) THEN
   open (unit=29,file='x1e.int',status='old')
   read(29,*)ufrzorbn
   ovlp_int=0
   MNP=0
   do
     read(29,'(a)',iostat=io)

     if(io.ne.0)exit
     MNP=MNP+1
     enddo
     rewind(29)


     do i=1,MNP
       read(29,*)line
       if(trim(line(1:5)).eq.'SSF0')lnum=i
     enddo

ELSE
PRINT*,' x1e.int file does not exist as a result I could not able &
to check if your system has frozen orbitals or not hope you have previded &
correct information otherwise result will be wrong'
ufrzorbn=0
ENDIF

! reading INFO file 
INQUIRE(FILE='INFO',EXIST=fileexists)
IF (fileexists) THEN
   open (unit=39,file='INFO',status='old')
ELSE
   PRINT*,'SORRY INFO file does not exist'
   stop
ENDIF
 
read(39,*)tot_atom,nbasis,unkwn
k=0
do i=1,tot_atom
  read(39,*)nbasatom(i)
  k=k+nbasatom(i)
  nbas_slab(i)=k
enddo

bfiexists=0
do i=1,10
  if(sl1(i).eq.4)bfiexists=1
enddo

if(ufrzorbn.ne.0)then
  if(ufrzorbn.lt.k.and.bfiexists.eq.0)then
    write(*,*)'it seems that you have frozen orbitals but you forget to put the $bfi section'
    stop
  else
    k=0
    do i=1,tot_atom
      l=0
      k=k+nbasatom(i)
      cnt(i)=k
      do j=1,frzn
        if(freezorb(j).le.cnt(i))then
          if(i.gt.1)then
            if(freezorb(j).le.cnt(i-1)) then
              exit
            endif
          endif
          l=l+1
        endif
      enddo
      nbasatom(i)=nbasatom(i)-l
    enddo
  endif
endif

if(key_frag.eq.0)then
  do i=1,200
    do j=1,20
      atoset(i,j)=0
    enddo
  enddo

  k5=0
  j=0
  k1=0
  k3=0
  k4=0
  do i1=1,tot_atom
    k=0
    j=j+nbasatom(i1)
    do i=orbsl-nao+1,orbsl
      k2=0
      do l=1,orbn(i)
        if(orbs(i,l).gt.j.or.orbs(i,l).le.j-nbasatom(i1))goto 200
        k2=k2+1
      enddo
      k5=k5+1
      k=k+1

      do i4=1,k3
        if(i1.eq.active_atoms(i4))goto 207
      enddo
      k3=k3+1
      active_atoms(k3)=i1

207 atoset(i1,k)=i-1
      atn(i1)=k 
      active_orbs(k5)=i-1
      atm_nb_orbs(k5)=i1
      goto 201
200 if(k2.ne.0.and.num_frag_cntr.eq.0)then
        k1=k1+1
        fragorb(k1)=i-1
        if(orbs(i,1).le.j.and.orbs(i,1).gt.j-nbasatom(i1))then
          k5=k5+1
          k=k+1
          k4=k4+1

do i4=1,k3
if(i1.eq.active_atoms(i4))goto 206
enddo
k3=k3+1
active_atoms(k3)=i1

206 frag_cent(k4)=i1
atoset(i1,k)=i-1
atn(i1)=k 
active_orbs(k5)=i-1
atm_nb_orbs(k5)=i1
endif
endif
if(k2.ne.0.and.num_frag_cntr.ne.0)then
do i3=1,num_frag_cntr
j1=0
do i2=1,frag_cntr(i3)
j1=j1+nbasatom(i2)
enddo
k5=k5+1
k=k+1
k4=k4+1
do i4=1,k3
if(frag_cntr(i3).eq.active_atoms(i4))goto 205
enddo
k3=k3+1
active_atoms(k3)=frag_cntr(i3)

205 frag_cent(k4)=frag_cntr(i3)
atoset(frag_cntr(i3),k)=i-1
active_orbs(k5)=i-1
atm_nb_orbs(k5)=frag_cntr(i3)
atn(frag_cntr(i3))=k

enddo
endif
201 enddo
enddo

atom=k3
nactorb=k5

k5=1
frag_cent1(k5)=frag_cent(1)
do i=2,k4
do j=1,k5
if(frag_cent(i).eq.frag_cent1(j))goto 320
enddo
k5=k5+1
frag_cent1(k5)=frag_cent(i-1)
320 enddo
nfragcent=k5

if(k1.ne.0.and.num_frag_cntr.eq.0)then
if(k1.eq.1)print*,'You have fragment orbital ',fragorb(1),' . Please mention the &  
central atoms of the fragmants or start the orbitals with basis function of the central atom'
if(k1.gt.1)print*,'You have fragment orbitals ',(fragorb(k),k=1,k1),' . Please mention the &  
central atoms of the fragmants or start the orbitals with basis function of the central atom'
endif

l2=0
do k1=1,tot_atom
l1=0
if(atoset(k1,1).eq.0)goto 321
l2=l2+1
do k2=1,20
if(atoset(k1,k2).ne.0)then
l1=l1+1
endif
enddo
atn(k1)=l1
321 enddo
atom=l2

endif

do i=1,tot_atom
read(39,*)symat(i),symatno(i),coordx(i),coordy(i),coordz(i)
coordx(i)=coordx(i)*0.529177
coordy(i)=coordy(i)*0.529177
coordz(i)=coordz(i)*0.529177
all_at_num(i)=int(symatno(i))
enddo
if(key_frag.eq.0)then
do i=1,atom
act_at_num(i)=symatno(active_atoms(i))
val_state_num(i)=valence_state(act_at_num(i))
enddo
endif


n=0
do i=1,nbasis
if(frzn.ne.0)then
do j=1,frzn
if(freezorb(j).eq.i)then
read(39,*)atbas_frz(i)
all_atbas(i)=atbas_frz(i)
goto 197
endif
enddo
endif
n=n+1
read(39,*)atbas(n)
all_atbas(i)=atbas(n)
197 enddo

if(key_frag.eq.0)then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! determination of the type of the orbitals 'pi or sigma' !!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
norbsym(4)=0
norbsym(1)=0
norbsym(2)=0
norbsym(3)=0

k=0
k5=0
k8=0
k4=0
k6=0
k1=0
k3=0
do i=orbsl-nao+1,orbsl
k7=0
n=0
do i2=1,atom
do i3=1,atn(active_atoms(i2))
print*,'ovlp_act_atm,atoset',active_atoms(i2),atoset(active_atoms(i2),i3)
if(atoset(active_atoms(i2),i3).eq.i-1)then
print*,'active_atoms(i2)',active_atoms(i2),nbasatom(3)
do j=1,active_atoms(i2)
k7=k7+1
n=n+nbasatom(j)
enddo
endif
enddo
enddo

do i2=1,orbn(i)
if(orbs(i,i2).le.n.and.orbs(i,i2).gt.n-nbasatom(k7))then
if(atbas(orbs(i,i2)).eq.'S')then
k=k+1
orbsym(4,k)=i-1
goto 208
endif
endif
enddo
do i2=1,orbn(i)
if(orbs(i,i2).le.n.and.orbs(i,i2).gt.n-nbasatom(k7))then
if(atbas(orbs(i,i2)).eq.'X')then
k3=k3+1
orbsym(1,k3)=i-1
goto 208
endif
endif
enddo
do i2=1,orbn(i)
if(orbs(i,i2).le.n.and.orbs(i,i2).gt.n-nbasatom(k7))then
if(atbas(orbs(i,i2)).eq.'Y')then
k4=k4+1
orbsym(2,k4)=i-1
goto 208
endif
endif
enddo
do i2=1,orbn(i)
if(orbs(i,i2).le.n.and.orbs(i,i2).gt.n-nbasatom(k7))then
if(atbas(orbs(i,i2)).eq.'Z')then
k8=k8+1
orbsym(3,k8)=i-1
goto 208
endif
endif
enddo

208 enddo
norbsym(4)=k
norbsym(1)=k3
norbsym(2)=k4
norbsym(3)=k8


if(ovlp_int.eq.1.and.ovopt.eq.1)then
n2=0
do i=1,atom
n=0
n1=0
n4=0
do i2=1,active_atoms(i)-1
n=n+nbasatom(i2)
enddo
do i2=1,active_atoms(i)
n1=n1+nbasatom(i2)
enddo
do i2=1,frzn
if(freezorb(i2).ge.n.and.freezorb(i2).le.n1)then
n4=n4+1
endif
enddo
n2=n2+n4
do i3=1,atn(active_atoms(i))
l1=0
l2=0
l3=0
l4=0
k9=0
if(n4.eq.0)then
do i2=1,orbn(atoset(active_atoms(i),i3)+1)
print*,'orbn(i)',i2,orbn(atoset(active_atoms(i),i3)+1),orbs(atoset(active_atoms(i),i3)+1,i2)+n2
n3=orbs(atoset(active_atoms(i),i3)+1,i2)+n2
if(n3.lt.n.or.n3.gt.n1)goto 531
if(all_atbas(n3).eq.'S')then
l1=l1+1
if(l1.eq.val_state_num(i))goto 530
endif
if(all_atbas(n3).eq.'X')then
l2=l2+1
if(l2.eq.val_state_num(i))goto 530
endif
if(all_atbas(n3).eq.'Y')then
l3=l3+1
if(l3.eq.val_state_num(i)-1)goto 530
endif
if(all_atbas(n3).eq.'Z')then
l4=l4+1
if(l4.eq.val_state_num(i)-1)goto 530
endif
goto 531
530 k9=k9+1
act_orb_typ(atoset(active_atoms(i),i3),k9)=all_atbas(n3)
ovlp_bas_num(atoset(active_atoms(i),i3),k9)=n3
531 enddo
endif
if(n4.ne.0)then
do i2=1,orbn(atoset(active_atoms(i),i3)+1)
n3=orbs(atoset(active_atoms(i),i3)+1,i2)+n2
if(n3.lt.n.or.n3.gt.n1)goto 631
if(all_atbas(n3).eq.'S'.and.l1.eq.0)then
l1=1
goto 630
endif
if(all_atbas(n3).eq.'X'.and.l2.eq.0)then
l2=1
goto 630
endif
if(all_atbas(n3).eq.'Y'.and.l3.eq.0)then
l3=1
goto 630
endif
if(all_atbas(n3).eq.'Z'.and.l4.eq.0)then
l4=1
goto 630
endif
goto 631
630 k9=k9+1
act_orb_typ(atoset(active_atoms(i),i3),k9)=all_atbas(n3)
ovlp_bas_num(atoset(active_atoms(i),i3),k9)=n3
631 enddo
endif
num_act_orb_typ(atoset(active_atoms(i),i3))=k9

enddo

enddo

endif

endif


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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

print*,'exit read_info'
call geocal(coordx,coordy,coordz)
return
end subroutine read_info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
