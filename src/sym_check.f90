!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine sym_check(nl,str2,n,sym_str_sl,sym_str_num,numsymset)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

common/quality/str_quality_1,str_quality_2,bondq,tqlty,bqlty,sqlty,tnqs,nssym,qulsym,symq,&
sigsym,tnqs_sig

logical :: fileexists
character(len=35)::inputfilename
integer::i,i1,i2,i3,i4,i5,i6,i7,i8,l,l1,m,bnd,n,nl,numsymset,nsymset,nsymstr,str6(100,20),str2(15000,20)
integer::str_quality_1(15000),str_quality_2(15000),bondq(15000),tqlty,bqlty,sqlty,tnqs,nssym,sym_check_flg,&
qulsym(15000),symq(15000),sigsym(15000),tnqs_sig,sym_str_sl1(1000,100),sym_str_sl(1000,100),sym_str_num(1000),&
sym_str_num1(1000),val,stc

print*, 'enter sym_check'
inputfilename='sym_str_set.dat'
INQUIRE(FILE=TRIM(inputfilename),EXIST=fileexists)
IF (fileexists) THEN
open(unit=43,file='sym_str_set.dat',status='unknown')
ELSE
PRINT*,'SORRY sym_str_set.dat file does not exist'
print*,'please write the file: at the top need to mention the number of sets and then number of structures in each sets'
print*,'please leave a empty line after each segment'
stop
ENDIF



read(43,*)nsymset,nsymstr
stc=0
bnd=(nae-nl*2-nlast)/2

101 do i=1,nsymset
read(43,*)
do i2=1,numsymset
do i4=1,sym_str_num(i2)
sym_str_sl1(i2,i4)=sym_str_sl(i2,i4)
enddo
enddo
do i1=1,nsymstr
read(43,*)(str6(i1,i2),i2=1,nae)
do i3=1,n
m=0
do i4=nl*2+1,nl*2+bnd*2,2
do i7=nl*2+1,nl*2+bnd*2,2
l=0
do i5=i4,i4+1
do i8=i7,i7+1
if(str6(i1,i5).eq.str2(i3,i8))then
l=l+1
goto 109
endif
enddo
109 enddo
if(l.eq.2)then
m=m+1
goto 111
endif
enddo
goto 112
111 enddo
if(m.eq.bnd)then
do i2=1,numsymset
do i4=1,sym_str_num(i2)
if(sym_str_sl1(i2,i4).eq.i3)then
sym_str_sl1(i2,i4)=0
endif
enddo
enddo
endif
112 enddo
enddo

do i2=1,1000
sym_str_num1(i2)=0
enddo

do i2=1,numsymset
do i4=1,sym_str_num(i2)
if(sym_str_sl1(i2,i4).ne.0)val=1
if(sym_str_sl1(i2,i4).eq.0)val=0
sym_str_num1(i2)=sym_str_num1(i2)+val
enddo
enddo

!do i2=1,numsymset
!print*,sym_str_num1(i2),sym_str_num(i2)
!enddo

!do i3=1,numsymset
!do i2=1,sym_str_num(i)
!write(*,200),i,')',(str2(sym_str_sl(i3,i2),i4),i4=1,nae)
!enddo
!enddo

sym_check_flg=1
do i2=1,numsymset
if(sym_str_num(i2)-sym_str_num1(i2).ne.0.and.sym_str_num(i2)-sym_str_num1(i2).ne.sym_str_num(i2))then
sym_check_flg=0
goto 114
endif
enddo
114 do i3=1,nsymstr
write(10,202)(str6(i3,i2),i2=1,nae)
enddo
if(sym_check_flg.eq.1)then
stc=stc+1
write(10,201)'set',i,'is symmetric'
endif
if(sym_check_flg.eq.0)write(10,201)'set',i,'is not symmetric'
write(10,*)

enddo
write(10,*)'Total number of symmetric sets =',stc

200 format (I2,a,2x,20I3)
201 format (a,I5,2x,a)
202 format (20I3)

!do i2=1,numsymset
!print*,(sym_str_sl(i2,i3),i3=1,sym_str_num(i2))
!enddo
!print*,'******************'
!do i2=1,numsymset
!do i3=1,sym_str_num(i2)
!print*,i2,sym_str_sl(i2,i3),(str2((sym_str_sl(i2,i3)),i4),i4=1,nae)
!enddo
!enddo

print*, 'exit sym_check'
stop
end subroutine sym_check
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
