!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_rumer_xmi(nl,str,nstr,rumer,rumer_rad,q_fac)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none
common/quality/str_quality_1,str_quality_2,bondq,tqlty,bqlty,sqlty,tnqs,nssym,qulsym,symq,&
sigsym,tnqs_sig

integer::i,i7,nstr,m19,m20,nl,rumer(15000),rumer_rad(15000),q_fac(15000),str_quality_1(15000),bondq(15000),allrum,&
str(15000,20),str_quality_2(15000),tqlty,bqlty,sqlty,tnqs,nssym,qulsym(15000),symq(15000),col(1000),sigsym(15000),tnqs_sig&
,rumstr(500,20)
real*8::ovlp
Double Precision::D(1000)

print*,'enter write_rumer_xmi',nstr
allrum=1
tqlty=0
bqlty=0
sqlty=0
i7=0
do m19=1,nstr
if(rumer(m19)*rumer_rad(m19).eq.1)then
i7=i7+1
col(i7)=m19
do m20=1,nae
rumstr(i7,m20)=str(m19,m20)
enddo
if(nfset.eq.3.or.nfset.eq.5)goto 200
if(niao.eq.0)then
write(10,914)str_quality_1(m19),bondq(m19),str_quality_2(m19),'|',(str(m19,m20),m20=1,nae)
endif
if(niao.gt.1)then
write(10,915)str_quality_1(m19),bondq(m19),str_quality_2(m19),'|',1,':',niao,(str(m19,m20),m20=1,nae)
endif
if(niao.eq.1)then
write(10,916)str_quality_1(m19),bondq(m19),str_quality_2(m19),'|',1,1,(str(m19,m20),m20=1,nae)
endif
tqlty=tqlty+str_quality_1(m19)
bqlty=bqlty+bondq(m19)
sqlty=sqlty+str_quality_2(m19)
endif
200 enddo
open(unit=31,file='Rumer_Sets.dat',status='unknown')
if(allrum.eq.1)call All_Rumer_set(rumstr,i7,nl)
if(nfset.eq.5)stop
if(nfset.eq.3)goto 201
if(ovopt.eq.1)then
call MatLDR('str',col,i7,D)

ovlp=1.0
do i=1,i7
ovlp=ovlp*D(i)
enddo
write(10,912)'Overlap of this set of the structures =',1.0-ovlp
endif
write(10,910)'qualities:',' intra_bond =',tqlty,'nn_bond =',bqlty,'sym_break=',sqlty

914 format(I3,x,I3,x,I3,x,a,x,25I4)
915 format(I3,x,I3,x,I3,x,a,x,I1,a,I1,x,25I4)
916 format(I3,x,I3,x,I3,x,a,x,I3,I3,x,25I4)
910 format(a,a,I3,x,a,I3,x,a,I3)
912 format(a,3x,F10.3)
print*,'exit write_rumer_xmi'
    
201 return
end subroutine write_rumer_xmi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
