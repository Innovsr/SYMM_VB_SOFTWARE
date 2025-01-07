!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_symm_xmi_1(i7,m21,m3,sf1,sf2,str3,q_fac2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

integer::nl,strn,ncqs,tostr,initstr,i,i1,i2,i3,i4,i5,i6,i7,i8,i9,m119,m18,m19,m20,m21,m23,m24,count&
,qul(100),nqul,j,jj,jjj,fg,flg,ii5,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,x,y,strs
integer::i5up16,i16i7,m16m21,i5up15,i15i7,m15m21,i5up14,i14i7,m14m21,i5up13,i13i7,m13m21,i5up12,i12i7,m12m21,&
i5up11,i11i7,m11m21,i5up10,i10i7,m10m21,i5up9,i9i7,m9m21,i5up8,i8i7,m8m21,i5up7,i7i7,m7m21,i5up6,i6i7,m6m21,&
i5up5,i5i7,m5m21,i5up4,i4i7,m4m21,i5up3,i3i7,m3m21,i5up2,i2i7,m2m21,i5up1,i1i7,m1m21,i5up17,i17i7,m17m21
integer::str3(15000,20),q_fac2(15000),finalvec(15000),strset(1000),col(1000),sigsym(15000),tnqs_sig,&
ffvec2(15000,1000),bondq(15000),bondq4(15000),nqset(15000),str5(2000,20),nstr7,&
tndet,totstr,Ifail,indpnt,strno(1000),str_quality_1(15000),str_quality_2(15000),ttqlty0,ttqlty&
,tqlty,bqlty,sqlty,hqlty,tnqs,nssym,qulsym(15000),symq(15000),set_number,ttqlty1,det_inv,ttqlty2,ttqlty3
integer::rumer(15000),rumer_rad(15000),quality_fac(15000),rumset,u1,max_set,mns,mm,group_num(15000)
!integer,dimension(:),allocatable::qq1,qq2,qq
integer::qq1(5000),qq2(5000),qq(5000),str2(2000,20),Rid,set_num(100),sf1,sf2,incmplt,mincmplt,mincmplt_set(500)
integer::j1,j2,grp(30)
real*8::ovlp
Double Precision::D(1000)
character(10)::dd,a
character(len=100)::outfile

common/quality/str_quality_1,str_quality_2,bondq,tqlty,bqlty,sqlty,tnqs,nssym,qulsym,symq,&
sigsym,tnqs_sig
common/str/str5,nstr7
common/infosymm/set_number,hqlty,ttqlty0,ttqlty1,ttqlty2,ttqlty3,ncqs,qul,Rid&
,mns,u1,max_set,rumset,nqset,strset,strn,totstr,strno,incmplt,mincmplt,mincmplt_set,group_num

print*,'*****************************'
print*,'enter write_symm_xmi_1'
print*,''
print*,''


!do i=i4+1,i5
!do i=i5+1,15000
!finalvec(i)=0
!enddo
!do m19=1,20
!print*,(str2(m19,i8),i8=1,8)
!enddo

!print*,'i7,m21',i7,m21
!stop
do m19=i7+1,1000
strno(m19)=0
 do i8=1,15
 str2(m19,i8)=0
 enddo
enddo

!i7=i3
totstr=i7
!m21=mm

jj=0
sf1=0
sf2=0
flg=0

!print*,'nqset,strn',nqset(m3),strn
if(nqset(m3).gt.strn)then
sf1=1
!print*,'rejected'
write(*,558)'number of structures of this group is =',nqset(m3),'so the group is rejected'
goto 100
endif
558 format(a,x,I2,x,a)

do i=strset(m3-1)+1,strset(m3)
totstr=totstr+1
strno(totstr)=i
!print*,'strno',strno(totstr)
enddo

j=0
do i=1,totstr
if(i.ne.1)then
j1=group_num(strno(i-1))
j2=group_num(strno(i))
if(j1.eq.j2)goto 555
endif
j=j+1
grp(j)=group_num(strno(i))
555 enddo
write(*,556)'groups involved = ',(grp(i),i=1,j)
556 format(a,3x,20I3)
write(*,*)

if(totstr.gt.strn)then
sf1=1
!print*,'rejected'
write(*,557)'includig this group the total number of structures &
becomes',totstr,'which is biggere than',strn,'so group has been rejected'
goto 100
endif
557 format(a,x,I2,x,a,x,I2,x,a)

write(*,559)'grps','str_num','structures'
do i=1,totstr
write(*,191)group_num(strno(i)),strno(i),(str3(strno(i),j),j=1,nae)
enddo
191 format(I3,3x,3x,I3,3X,20I3)
559 format(a,5x,a,6x,a)

call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)

if(Ifail.eq.0)print*,'structures are INDEPENDENT'
if(Ifail.eq.1)print*,'structures are DEPENDENT'

if(Ifail.eq.1) then
sf1=1
goto 100
endif

if(totstr.gt.mincmplt)then
mincmplt=totstr
do i=1,mincmplt
mincmplt_set(i)=strno(i)
enddo
endif

if(totstr.gt.1.and.ovopt.eq.vpt) then
 call MatLDR('str',strno,totstr,D)

 ovlp=1.0
 do i=1,totstr
 ovlp=ovlp*D(i)
 enddo

 if(1.0-ovlp.gt.ovval)then 
 sf1=1
 goto 100
 endif
endif

 m21=m21+nqset(m3)
 if(m21.gt.strn)then
 sf1=1
 goto 100
 endif

 do m19=strset(m3-1)+1,strset(m3)
 if(q_fac2(m19).ne.qul(m3))goto 603
 i7=i7+1
! print*,'new struc',i7,m19,q_fac2(m19),str_quality_1(m19),str_quality_2(m19),bondq(m19)
  do i3=1,nae
  str2(i7,i3)=str3(m19,i3)
  enddo
 col(i7)=m19
 qq(i7)=q_fac2(m19)
 qq1(i7)=str_quality_1(m19)
 qq2(i7)=str_quality_2(m19)
 bondq4(i7)=bondq(m19)
 flg=1
 count=3
 603 enddo

 if(i7.eq.strn) then
  ttqlty=0
  tqlty=0
  bqlty=0
  sqlty=0
  do m119=1,i7
   ttqlty=ttqlty+qq(m119)
  enddo
!   print*,'ttqlty,ttqlty0',ttqlty,ttqlty0
  if(ttqlty.le.ttqlty0)then
    set_number=set_number+1
    mns=mns+1
    if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
    if(nfset.eq.1.and.set_number.gt.1)then
     if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
!     write(9,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
    endif

print*,'i7***********',i7,m21,nfset
   do m19=1,i7
incmplt=0
    if(niao.eq.0)then
     write(9+u1,900)qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
    endif
    if(niao.gt.1)then
     write(9+u1,901)qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
    endif
    if(niao.eq.1)then
     write(9+u1,909)qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
    endif
   enddo

   if(ovopt.eq.vpt) then
    call MatLDR('str',col,i7,D)

    ovlp=1.0
    do i=1,i7
     ovlp=ovlp*D(i)
    enddo
    write(9+u1,912)'Overlap of this set of the structures =',1.0-ovlp
   endif

    if(nfset.eq.3)then
    Rumwrite=1
    call Rumer_set_id(str2,i7,nl,Rid,set_num)
    endif

   if(Rid.eq.0)write(9+u1,910)'quality value',ttqlty
   if(Rid.eq.1)write(9+u1,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
   bqlty=0
   write(9+u1,*)'Set_number=',set_number
   write(9+u1,*)'    '
  endif

  if(nfset.eq.0.and.ttqlty.le.ttqlty0)then 
   sf2=1
   sf1=0
   goto 100
  endif
! endif

print*,'i am here',mns,max_set,u1,mset,mset-1
 if(mns.eq.max_set)then
  if(u1.eq.mset-1)then
print*,'i am here',mns,max_set,u1,mset,mset-1
   sf2=0
   sf1=1
   goto 100
  endif
  close(9+u1)
  mns=0
  u1=u1+1
  if(u1.le.9)then
   write(a,'(I1)')u1
  endif
  if(u1.gt.9.and.u1.lt.100)then
   write(a,'(I2)')u1
  endif
  if(u1.gt.99.and.u1.lt.1000)then
   write(a,'(I3)')u1
  endif
  if(u1.eq.1000)then
   write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
   sets. Please increase the limit'
   sf2=1
   goto 100
  endif
  outfile='structure_set_'//trim(a)//trim('.dat')
  open(unit=9+u1,file=outfile,status='unknown')
 endif

!sf2=1
goto 100
endif

900 format(I3,x,I3,x,I3,x,I3,x,a,x,25I4)
901 format(I3,x,I3,x,I3,x,I3,x,a,x,I1,a,I3,x,25I4)
909 format(I3,x,I3,x,I3,x,I3,x,a,x,I3,I3,x,25I4)
910 format(a,I3)
920 format(a,2x,I5,2x,a,2x,100I5)
911 format(15x,I3,7x,I3,7x,I3)
912 format(a,3x,F10.3)

100 return
end subroutine write_symm_xmi_1
