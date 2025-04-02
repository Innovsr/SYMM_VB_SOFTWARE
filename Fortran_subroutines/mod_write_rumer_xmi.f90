!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module mod_write_rumer_xmi
use commondat_mod
use quality_mod
!use mod_MatLDR 
use rum_rad_mod
implicit none

contains
subroutine write_rumer_xmi(nl,str,nstr)!,rumer,rumer_rad)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer::i,i7,nstr,m19,m20,nl,set_number
integer, allocatable::col(:)
integer, pointer::str(:,:)

print*,'enter write_rumer_xmi',nstr, CovDim

allocate(col(CovDim))
tqlty = 0
bqlty = 0
sqlty = 0
i7=0
loop1:do m19=1,nstr
   if (rumer(m19)*rumer_rad(m19).eq.1)then
      i7=i7+1
      col(i7)=m19
      !do m20=1,nae
      !   rumerstr(i7,m20)=str(m19,m20)
      !enddo
      if(nfset.eq.3.or.nfset.eq.5) cycle loop1
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
enddo loop1
set_number = 1
write(5,913)'Set_number',set_number,(col(m19),m19=1,i7)
!open(unit=31,file='Rumer_Sets.dat',status='unknown')
if(nfset.eq.5) stop
if(nfset.eq.3) then
   return
endif

!if(ovopt.eq.1)then
!   !allocate(D(i7))
!   call MatLDR('str',col,i7,D)
!   ovlp=1.0
!   do i=1,i7
!      ovlp=ovlp*D(i)
!   enddo
!   write(10,912)'Overlap of this set of the structures =',1.0-ovlp
!   !deallocate(D)
!endif
!write(10,910)'qualities:',' intra_bond =',tqlty,'nn_bond =',bqlty,'sym_break=',sqlty
!if(allrum.eq.1) call All_Rumer_set(i7,nl)

913 format(a,2x,I0,4x,*(I0, 1x))
914 format(2x,I0,3x,I0,3x,I0,8x,a,1x,*(I0, 1x))
915 format(2x,I0,3x,I0,3x,I0,8x,a,1x,I0,a,I0,1x,*(I0, x))
916 format(2x,I3,3x,I3,3x,I3,8x,a,1x,I3,I3,1x,*(I0, x))
910 format(2x,a,a,I0,1x,a,I0,1x,a,I0)
912 format(a,3x,F10.3)
print*,'exit write_rumer_xmi'
    
return
end subroutine write_rumer_xmi
end module mod_write_rumer_xmi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
