!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! It identifies the nearest naighbour connections among the sigma orbitals presents  !
!! in a atom and updated the connectivity matrix accordingly                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine nnat_bond_sig(nnat_bond_new,n,sig_orb,nsig,nnat_bond_new_1,j)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

integer::i,i1,i2,i3,i4,i5,i6,j,k,l,nnat_bond_new(100,2),nnat_bond_new_1(100,2),n
integer::atoset_sig(20,10),atn_sig(20),orbs_sig(10),sig_orb(100),nsig,k1,k2,c1,c2,s_orbs(0)
integer::k3,k4
real*8::coordx(100),coordy(100),coordz(100),a,b,c,d,sig_coord(50,3),norm,kk,sig_dist_mat(40,40)

common/coordinate/coordx,coordy,coordz

! finding the sigma orbitals associated with an atom
k=0
do i=1,atom
  if(atn(i).gt.2)then
    j=0
    do i1=1,atn(i)
      do i2=1,nsig
!print*,'atoset',atoset(active_atoms(i),i1),sig_orb(i2)
        if(atoset(active_atoms(i),i1).eq.sig_orb(i2))then
          j=j+1
          orbs_sig(j)=sig_orb(i2)
        endif
      enddo
    enddo
    if(j.gt.2)then
      k=k+1
      do i4=1,j
        atoset_sig(k,i4)=orbs_sig(i4)
      enddo
      atn_sig(k)=j
    endif
  endif
enddo


do i=1,k
print*,'atoset_sig',(atoset_sig(i,j),j=1,atn_sig(i))
enddo

j=0
loop1:do i2=1,n
  do i3=1,2
    do i=1,k
      do i1=1,atn_sig(i)
        if(atoset_sig(i,i1).eq.nnat_bond_new(i2,i3))then
          k1=nnat_bond_new(i2,i3)
          if(i3.eq.1)k2=nnat_bond_new(i2,2)
          if(i3.eq.2)k2=nnat_bond_new(i2,1)
          do i4=1,atn_sig(i)
          if(k2.eq.atoset_sig(i,i4)) then
            cycle loop1
          endif
          enddo
          endif
        enddo
      enddo
    enddo
    j=j+1
nnat_bond_new_1(j,1)=nnat_bond_new(i2,1)
nnat_bond_new_1(j,2)=nnat_bond_new(i2,2)
enddo loop1



l=0
k1=0
do i=1,k
  do i1=1,atn_sig(i)
    loop3:do i2=1,n
      do i3=1,2
        if(atoset_sig(i,i1).eq.nnat_bond_new(i2,i3))then
          k1=nnat_bond_new(i2,i3)
          if(i3.eq.1)k2=nnat_bond_new(i2,2)
          if(i3.eq.2)k2=nnat_bond_new(i2,1)
          if(k2.gt.nao+niao)k2=k2-nao-niao
          do i4=1,atn_sig(i)
            if(k2.eq.atoset_sig(i,i4)) cycle loop3
          enddo

          c1=0
          c2=0
          do i5=1,atom
            do i6=1,atn(i5)
              if(k1.eq.atoset(active_atoms(i5),i6))c1=i5
              if(k2.eq.atoset(active_atoms(i5),i6))c2=i5
            enddo
          enddo

          if(c2.eq.0)c2=k2
          if(c1.eq.c2) cycle loop3
          l=l+1
          a=coordx(c2)-coordx(c1)
          b=coordy(c2)-coordy(c1)
          c=coordz(c2)-coordz(c1)
          norm=1.0/sqrt(a**2.0 + b**2.0 + c**2.0)
          d=0.3
          s_orbs(l)=k1
          sig_coord(l,1)=coordx(c1)+a*d*norm
          sig_coord(l,2)=coordy(c1)+b*d*norm
          sig_coord(l,3)=coordz(c1)+c*d*norm

!print*,'norm,a,b,c',norm,a,b,c
!print*,'a,b,c',a,b,c
!print*,'c1 coord',coordx(c1),coordy(c1),coordz(c1)
!print*,'c2 coord',coordx(c2),coordy(c2),coordz(c2)
!print*,'sig_coord',l,sig_coord(l,1),sig_coord(l,2),sig_coord(l,3)
!

        endif
      enddo
    enddo loop3
  enddo


  do i1=1,l
    do i2=1,l
      sig_dist_mat(i1,i2)=0.0
    enddo
  enddo
  
  do i1=1,l
    do i2=1,l
      sig_dist_mat(i1,i2)=0.0
      sig_dist_mat(i1,i2)=sqrt((sig_coord(i1,1)-sig_coord(i2,1))**2.0 + &
      (sig_coord(i1,2)-sig_coord(i2,2))**2.0 + &
      (sig_coord(i1,3)-sig_coord(i2,3))**2.0)
    enddo
  enddo
  
  !do i1=1,l
  !print*,'sig_dist_mat',(sig_dist_mat(i1,i2),i2=1,l)
  !enddo
  
  do i1=1,l
    kk=1000.0
    do i2=1,l
      if(i1.eq.i2) cycle
      if(sig_dist_mat(i1,i2).lt.kk)kk=sig_dist_mat(i1,i2)
    enddo
    loop4:do i2=1,l
      if(kk.eq.sig_dist_mat(i1,i2))then
        do i3=1,j
          do i4=1,2
            if(s_orbs(i1).eq.nnat_bond_new_1(i3,i4))then
              k3=nnat_bond_new_1(i3,i4)
              if(i4.eq.1)k4=nnat_bond_new_1(i3,2)
              if(i4.eq.2)k4=nnat_bond_new_1(i3,1)
            endif
          enddo
          if(s_orbs(i2).eq.k4) cycle loop4
        enddo
        j=j+1
        nnat_bond_new_1(j,1)=s_orbs(i1)
        nnat_bond_new_1(j,2)=s_orbs(i2)
      endif
    enddo loop4
  enddo

enddo


return
end subroutine nnat_bond_sig
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
