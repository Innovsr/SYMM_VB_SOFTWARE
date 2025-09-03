!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module mod_write_symm_xmi_1
use commondat_mod
use quality_mod
use infosymm_mod
use mod_mat_ind
use str_module
!use mod_MatLDR 
implicit none

contains
subroutine write_symm_xmi_1(i7,m21,m3,sf1,sf2,str3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer::nl,i,i1,i2,i3,i7,i8,m119,m19,m20,m21,count,j,jj,flg,m3,x,Ifail
integer::ttqlty,sf1,sf2,j1,j2,grp(30)
!real*8::ovlp
!Double Precision, allocatable::D(:)
!Double Precision::D(10000),ovlp
character(10)::a
integer, pointer::str3(:,:)
!integer, allocatable::strno(:),str12(:,:),col9(:)

print*,'enter write_symm_xmi_1',strnn,MaxStrOepo,CovDim

do m19=i7+1,MaxStrOepo
  strno(m19)=0
  do i8=1,nae
    str12(m19,i8)=0
  enddo
enddo

totstr=i7

jj=0
sf1=0
sf2=0
flg=0

if (symm.eq.1) then
  if(nqset(m3).gt.strnn)then
    sf1=1
    write(*,558)'number of structures of this group is =',nqset(m3),'so the group is rejected'
    return
  endif
endif
558 format(a,x,I2,x,a)

if (asymm.eq.1) then
  totstr=totstr+1
  strno(totstr)=m3
endif

if (symm.eq.1) then
  if(m3.eq.1)then
    do i=1,strset(m3)
      totstr=totstr+1
      strno(totstr)=i
    enddo
  else
    do i=strset(m3-1)+1,strset(m3)
      totstr=totstr+1
      strno(totstr)=i
    enddo
  endif
endif

if (symm.eq.1) then
  j=0
  do i=1,totstr
    if(i.ne.1)then
      j1=group_num(strno(i-1))
      j2=group_num(strno(i))
      if(j1.eq.j2) then
        cycle
      endif
    endif
    j=j+1
    grp(j)=group_num(strno(i))
  enddo

  write(*,556)'groups involved = ',(grp(i),i=1,j)
556 format(a,3x,20I3)
  write(*,*)
  
  if(totstr.gt.strnn)then
     sf1=1
     write(*,557)"includig this group the total number of structures becomes",totstr,&
             "which is biggere than",strnn,"so group has been rejected"
     return
  endif
557 format(a,x,I2,x,a,x,I2,x,a)
  
  write(*,559)'grps','str_num','structures'
  do i=1,totstr
     write(*,191)group_num(strno(i)),strno(i),(str3(strno(i),j),j=1,nae)
  enddo
endif
191 format(I3,3x,3x,I3,3X,20I3)
559 format(a,5x,a,6x,a)

call mat_ind(totstr,ncqss,strno,Ifail)

if(Ifail.eq.0)print*,'structures are INDEPENDENT'
if(Ifail.eq.1)print*,'structures are DEPENDENT'

if(Ifail.eq.1) then
   sf1=1
   return
endif

if(totstr.gt.mincmplt)then
   mincmplt=totstr
   do i=1,mincmplt
      mincmplt_set(i)=strno(i)
   enddo
endif

!if(totstr.gt.1.and.ovopt.eq.vpt) then
! !allocate(D(totstr))
! call MatLDR('str',strno,totstr,D)
!
! ovlp=1.0
! do i=1,totstr
!    ovlp=ovlp*D(i)
! enddo
!
! !deallocate(D)
!
! if(1.0-ovlp.gt.ovval)then 
! sf1=1
! return
! endif
!endif

if(symm.eq.1) m21=m21+nqset(m3)
if(asymm.eq.1) m21=m21+1

 if(m21.gt.strnn)then
    sf1=1
    return
 endif
! if (.not. allocated(col9))then
!   allocate(col9(CovDim))
!   col9 = 0 
! endif

if (asymm.eq.1)then
  i7=i7+1
  do i3=1,nae
     str12(i7,i3)=str3(m3,i3)
  enddo
  col9(i7)=m3
  qq10(i7)=quality_fac(m3)
  qq11(i7)=str_quality_1(m3)
  qq12(i7)=str_quality_2(m3)
  bondq14(i7)=bondq(m3)
endif

if (symm.eq.1) then
 if(m3.eq.1)then
    loop1:do m19=1,strset(m3)
       if(quality_fac(m19).ne.qul(m3)) cycle loop1
       i7=i7+1
        do i3=1,nae
           str12(i7,i3)=str3(m19,i3)
        enddo
        col9(i7)=m19
        qq10(i7)=quality_fac(m19)
        qq11(i7)=str_quality_1(m19)
        qq12(i7)=str_quality_2(m19)
        bondq14(i7)=bondq(m19)
        flg=1
        count=3
    enddo loop1
 else
    loop2:do m19=strset(m3-1)+1,strset(m3)
       if(quality_fac(m19).ne.qul(m3)) cycle loop2
       i7=i7+1
        do i3=1,nae
           str12(i7,i3)=str3(m19,i3)
        enddo
        col9(i7)=m19
        qq10(i7)=quality_fac(m19)
        qq11(i7)=str_quality_1(m19)
        qq12(i7)=str_quality_2(m19)
        bondq14(i7)=bondq(m19)
        flg=1
        count=3
    enddo loop2
 endif
endif

if(i7.eq.strnn) then
  call write_output(sf1, sf2, i7, str12, col9)
endif


return
end subroutine write_symm_xmi_1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_output(sf1, sf2, i7, str12, col9)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none

integer::ttqlty, sf1, sf2, m119, m19, m20, i7, i
integer, allocatable:: str12(:,:),col9(:)
character(10)::a
character(len=300)::outfile

print*,'enter write_output',i7,ttqlty, ttqlty0
incmplt=0
ttqlty=0
tqlty=0
bqlty=0
sqlty=0
do m119=1,i7
  ttqlty=ttqlty+qq10(m119)
enddo
if(ttqlty.le.ttqlty0)then
  set_number=set_number+1
  mns=mns+1
  if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
  if(nfset.eq.1.and.set_number.gt.1)then
    if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
  endif

  do m19=1,i7
    print*,'enter loop',m19
    if(niao.eq.0)then
      write(9+u1,900)qq11(m19),bondq14(m19),qq12(m19),qq10(m19),'|',(str12(m19,m20),m20=1,nae)
    endif
    if(niao.gt.1)then
      write(9+u1,901)qq11(m19),bondq14(m19),qq12(m19),qq10(m19),'|',1,':',niao,(str12(m19,m20),m20=1,nae)
    endif
    if(niao.eq.1)then
      write(9+u1,909)qq11(m19),bondq14(m19),qq12(m19),qq10(m19),'|',1,1,(str12(m19,m20),m20=1,nae)
    endif
  enddo

!   if(ovopt.eq.vpt) then
!    !allocate(D(i7))
!    call MatLDR('str',col9,i7,D)
!
!    ovlp=1.0
!    do i=1,i7
!     ovlp=ovlp*D(i)
!    enddo
!    write(9+u1,912)'Overlap of this set of the structures =',1.0-ovlp
!    !deallocate(D)
!   endif

!    if(nfset.eq.3)then
!    Rumwrite=1
!    call Rumer_set_id(str12,i7,nl,Rid)
!    endif

  if(Rid.eq.0)write(9+u1,910)'quality value',ttqlty
  if(Rid.eq.1)write(9+u1,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
  bqlty=0
  write(9+u1,*)'Set_number=',set_number
  write(9+u1,*)'    '
  if(flg_cov.eq.1) write(5,913)'Set_number_c',set_number,(col9(m19),m19=1,i7)
  if(flg_ion.eq.1) write(5,913)'Set_number_i',set_number,(col9(m19),m19=1,i7)
  print*,'col9',(col9(m19),m19=1,i7)
endif

if(nfset.eq.0.and.ttqlty.le.ttqlty0)then 
  sf2=1
  sf1=0
  return
endif

if(mns.eq.max_set)then
  if(u1.eq.mset-1)then
    sf2=0
    sf1=1
    return
  endif
  close(9+u1)
  mns=0
  u1=u1+1
  write(a,'(I0)')u1
!  if(u1.le.9)then
!   write(a,'(I1)')u1
!  endif
!  if(u1.gt.9.and.u1.lt.100)then
!   write(a,'(I2)')u1
!  endif
!  if(u1.gt.99.and.u1.lt.1000)then
!   write(a,'(I3)')u1
!  endif
  if(u1.eq.1000)then
    write(*,*)'Maximum output file has been set to 1000, it seems that you have more sets. Please increase the limit'
    sf2=1
    return
  endif
  outfile=trim(out_folder_path)//trim('/')//trim(STDOUT)//trim('_')//trim(a)//trim('.out')
  open(unit=9+u1,file=trim(outfile),status='unknown')
endif

900 format(2x,I0,3x,I0,3x,I0,3x,I0,8x,a,1x,*(I0, 1x))
901 format(2x,I0,3x,I0,3x,I0,3x,I0,8x,a,1x,I0,a,I0,1x,*(I0, 1x))
909 format(2x,I0,3x,I0,3x,I0,3x,I0,8x,a,1x,I0,I0,1x,*(I0, 1x))
910 format(a,I0)
920 format(a,2x,I5,2x,a,2x,*(I0, 1x))
912 format(a,3x,F10.3)
913 format(a,2x,I0,4x,*(I0, 1x))

!deallocate(str12)
end subroutine write_output

end module mod_write_symm_xmi_1
