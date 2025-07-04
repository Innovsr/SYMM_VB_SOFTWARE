!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! All possible sets containg various linearly independent symmetric structures groups generated here.!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module mod_write_symm_xmi_new
use commondat_mod
use quality_mod
use mod_rumer_structures
!use write_rumer
use mod_write_symm_xmi_1
use infosymm_mod
use rum_rad_mod
implicit none

contains
subroutine write_symm_xmi_new(nl,strn,str3,ncqs)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer::nl,strn,ncqs,i,i3,i7,m19,m20,m21,indpnt,nqul,j,jj,jjj,flg
integer::m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,x
integer::i16i7,m16m21,i15i7,m15m21,i14i7,m14m21,i13i7,m13m21,i12i7
integer::m12m21,i11i7,m11m21,i10i7,m10m21,i9i7,m9m21,i8i7,m8m21,i7i7,sf1,sf2
integer::m7m21,i6i7,m6m21,i5i7,m5m21,i4i7,m4m21,i3i7,m3m21,i2i7,m2m21,i1i7,m1m21
integer, allocatable::qq1(:),qq2(:),qq(:),bondq4(:)
character(10)::a
integer, pointer::str3(:,:)!q_fac2(:), 
!integer, allocatable::rumer(:),rumer_rad(:)

print*,'enter_symm_xmi_new',strn, ncqs, nfset, flg_cov, flg_ion

if(flg_cov.eq.1)then
  cov_space=cov_space+1
  write(5,913)'Cov_Space',cov_space
  write(9+u1,914)'================= Cov Space num:',cov_space,'================'
endif
if(flg_ion.eq.1)then
  ion_space=ion_space+1
  write(5,913)'Ion_Space',ion_space
  write(9+u1,914)'================= Ion Space num:',ion_space,'================'
endif
913 format(a,2x,I0)
914 format(a,2x,I0,2x,a)

if (.not. allocated(nqset))then
   allocate(nqset(MaxStrOepo))
   nqset = 0
endif
if (.not. allocated(strset))then
   allocate(strset(MaxStrOepo))
   strset = 0
endif
if (.not. allocated(group_num)) then
   allocate(group_num(MaxStrOepo))
   group_num = 0
endif
if (.not. allocated(qul)) then
  allocate(qul(MaxStrOepo))
  qul = 0
endif
if (.not. allocated(mincmplt_set)) then
allocate(mincmplt_set(MaxStrOepo))
endif

Rid=0 ! Rumer id; to varify the set is Rumer of not. Rid=1 call the Rumer sets

! if the loops could not able to find out a set it will provide the biggest incomplete set
! which is counted and stored with 'mincmplt' and mincmplt_set array.
mincmplt=0
mincmplt_set=0

ncqss=ncqs ! total number of structures
strnn=strn ! total number of permiissible structures in a set

incmplt=1 ! running variable to count maximum number of independent structures in each set
mns=0

max_set=75000 ! maximum number of set will be written in one file


set_number=0
bqlty=0
tqlty=0
sqlty=0
hqlty=0
indpnt=2
if(noq0.gt.strnn)then
   ttqlty0=noq0
else
   ttqlty0=strnn+noq0
endif
ttqlty1=strnn+noq1
ttqlty2=strnn+noq2
ttqlty3=strnn+noq3


if(symm.eq.1) write(*,*)'sl  structures           group_numbers'
if(asymm.eq.1) write(*,*)'sl  structures           quality_factors'
do i=1,ncqss
   write(*,231)i,(str3(i,j),j=1,nae),quality_fac(i)
   group_num(i)=quality_fac(i)
enddo

231 format(30I3)

if(symm.eq.1) then
jj=1
loop1:do m19=1,ncqss
  if(m19.eq.1)qul(1)=quality_fac(1)
    j=jj
    loop2:do i=1,j
      if(qul(i).eq.quality_fac(m19)) then
        cycle loop1
      endif
    enddo loop2
    jj=jj+1
    qul(i)=quality_fac(m19)
enddo loop1
nqul=jj

! counting of groups and the number of structures in each group
do i=1,nqul
   jjj=0
   jj=0
   do m19=1,ncqss
      if(qul(i).eq.quality_fac(m19))then
         jjj=jjj+1
         jj=m19
      endif
   enddo
   nqset(i)=jjj
   strset(i)=jj
enddo
endif

if (asymm.eq.1) nqul=ncqs

flg=0
totstr=0
i7=0
m21=0
i1i7=0
m1m21=0


i1i7=i7
m1m21=m21

! nested loops to find group combinations starts here
do m1=1,nqul
   i7=0
   m21=0
   print*,'m1m1',m1
   print*,'loop1'

   call write_symm_xmi_1(i7,m21,m1,sf1,sf2,str3)
   if (sf1.eq.1) cycle
   if (sf2.eq.1) return

   i2i7=i7
   m2m21=m21

   do m2=m1+1,nqul
      i7=i2i7
      m21=m2m21
      print*,'m2m2',m2
      print*,'loop2'

      call write_symm_xmi_1(i7,m21,m2,sf1,sf2,str3)
      if (sf1.eq.1) cycle
      if (sf2.eq.1) return
      
      i3i7=i7
      m3m21=m21

      do m3=m2+1,nqul
         i7=i3i7
         m21=m3m21
         print*,'m3m3',m3
         print*,'loop3'
         
         call write_symm_xmi_1(i7,m21,m3,sf1,sf2,str3)
         if (sf1.eq.1) cycle
         if (sf2.eq.1) return
         
         i4i7=i7
         m4m21=m21

         do m4=m3+1,nqul
            i7=i4i7
            m21=m4m21
            print*,'m4m4',m4
            print*,'loop4'
            
            call write_symm_xmi_1(i7,m21,m4,sf1,sf2,str3)
            if (sf1.eq.1) cycle
            if (sf2.eq.1) return
            
            i5i7=i7
            m5m21=m21

            do m5=m4+1,nqul
               i7=i5i7
               m21=m5m21
               print*,'m5m5',m5
               print*,'loop5'
               
               call write_symm_xmi_1(i7,m21,m5,sf1,sf2,str3)
               if (sf1.eq.1) cycle
               if (sf2.eq.1) return
               
               i6i7=i7
               m6m21=m21

               do m6=m5+1,nqul
                  i7=i6i7
                  m21=m6m21
                  
                  call write_symm_xmi_1(i7,m21,m6,sf1,sf2,str3)
                  if (sf1.eq.1) cycle
                  if (sf2.eq.1) return
                  
                  i7i7=i7
                  m7m21=m21

                  do m7=m6+1,nqul
                     i7=i7i7
                     m21=m7m21
                     
                     !print*,'loop7'
                     call write_symm_xmi_1(i7,m21,m7,sf1,sf2,str3)
                     if (sf1.eq.1) cycle
                     if (sf2.eq.1) return
                     
                     i8i7=i7
                     m8m21=m21

                     do m8=m7+1,nqul
                        i7=i8i7
                        m21=m8m21
                        
                        call write_symm_xmi_1(i7,m21,m8,sf1,sf2,str3)
                        if (sf1.eq.1) cycle
                        if (sf2.eq.1) return
                        
                        i9i7=i7
                        m9m21=m21

                        do m9=m8+1,nqul
                           i7=i9i7
                           m21=m9m21
                           
                           !print*,'loop9'
                           call write_symm_xmi_1(i7,m21,m9,sf1,sf2,str3)
                           if (sf1.eq.1) cycle
                           if (sf2.eq.1) return
                           
                           i10i7=i7
                           m10m21=m21

                           do m10=m9+1,nqul
                              i7=i10i7
                              m21=m10m21
                              
                              !print*,'loop10'
                              call write_symm_xmi_1(i7,m21,m10,sf1,sf2,str3)
                              if (sf1.eq.1) cycle
                              if (sf2.eq.1) return
                              
                              i11i7=i7
                              m11m21=m21

                              do m11=m10+1,nqul
                                 i7=i11i7
                                 m21=m11m21
                                 
                                 !print*,'loop11'
                                 call write_symm_xmi_1(i7,m21,m11,sf1,sf2,str3)
                                 if (sf1.eq.1) cycle
                                 if (sf2.eq.1) return
                                 
                                 i12i7=i7
                                 m12m21=m21

                                 do m12=m11+1,nqul
                                    i7=i12i7
                                    m21=m12m21
                                    
                                    print*,'loop12'
                                    call write_symm_xmi_1(i7,m21,m12,sf1,sf2,str3)
                                    if (sf1.eq.1) cycle
                                    if (sf2.eq.1) return
                                    
                                    i13i7=i7
                                    m13m21=m21

                                    do m13=m12+1,nqul
                                       i7=i13i7
                                       m21=m13m21
                                       
                                       print*,'loop13'
                                       call write_symm_xmi_1(i7,m21,m13,sf1,sf2,str3)
                                       if (sf1.eq.1) cycle
                                       if (sf2.eq.1) return
                                       
                                       i14i7=i7
                                       m14m21=m21

                                       do m14=m13+1,nqul
                                          i7=i14i7
                                          m21=m14m21
                                          
                                          call write_symm_xmi_1(i7,m21,m14,sf1,sf2,str3)
                                          if (sf1.eq.1) cycle
                                          if (sf2.eq.1) return
                                          
                                          i15i7=i7
                                          m15m21=m21

                                          do m15=m14+1,nqul
                                             i7=i15i7
                                             m21=m15m21
                                             
                                             print*,'loop15'
                                             call write_symm_xmi_1(i7,m21,m15,sf1,sf2,str3)
                                             if (sf1.eq.1) cycle
                                             if (sf2.eq.1) return
                                             
                                             i16i7=i7
                                             m16m21=m21

                                             do m16=m15+1,nqul
                                                i7=i16i7
                                                m21=m16m21
                                                
                                                print*,'loop16'
                                                call write_symm_xmi_1(i7,m21,m16,sf1,sf2,str3)
                                                if (sf1.eq.1) cycle
                                                if (sf2.eq.1) return
                                                
                                                i16i7=i7
                                                m16m21=m21

                                              enddo

                                           enddo
 
                                        enddo
   
                                     enddo
 
                                  enddo
 
                               enddo

                            enddo

                         enddo

                      enddo

                   enddo

                enddo

             enddo

          enddo

       enddo

    enddo

enddo


deallocate(nqset)
deallocate(strset)
deallocate(qul)
deallocate(group_num)
close(21)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! if the above loop become unsuccessfull to provide a full dimention set the below part wil then!
!! print the maximum generated set in the output                                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(incmplt.eq.1)then
if (.not. allocated(qq)) then
  allocate(qq(MaxStrOepo))
  qq = 0
endif
if (.not. allocated(qq1)) then
  allocate(qq1(MaxStrOepo))
  qq1 = 0
endif
if (.not. allocated(qq2)) then
  allocate(qq2(MaxStrOepo))
  qq2 = 0
endif
if (.not. allocated(bondq4)) then
  allocate(bondq4(MaxStrOepo))
  bondq4 = 0
endif
write(10,*)'----------- incomplete set -----------'
do i=1,mincmplt
 qq(i)=quality_fac(mincmplt_set(i))
 qq1(i)=str_quality_1(mincmplt_set(i))
 qq2(i)=str_quality_2(mincmplt_set(i))
 bondq4(i)=bondq(mincmplt_set(i))
    if(niao.eq.0)then
     write(10,900)qq1(i),bondq4(i),qq2(i),qq(i),'|',(str3(mincmplt_set(i),m20),m20=1,nae)
    endif
    if(niao.gt.1)then
     write(10,901)qq1(i),bondq4(i),qq2(i),qq(i),'|',1,':',niao,(str3(mincmplt_set(i),m20),m20=1,nae)
    endif
    if(niao.eq.1)then
     write(10,909)qq1(i),bondq4(i),qq2(i),qq(i),'|',1,1,(str3(mincmplt_set(i),m20),m20=1,nae)
    endif
enddo
endif

900 format(I3,x,I3,x,I3,x,I3,x,a,x,25I4)
901 format(I3,x,I3,x,I3,x,I3,x,a,x,I1,a,I3,x,25I4)
909 format(I3,x,I3,x,I3,x,I3,x,a,x,I3,I3,x,25I4)

! removing temporary file
!CALL SYSTEM ("rm Rumer_Sets.dat")


return
end subroutine write_symm_xmi_new
end module mod_write_symm_xmi_new
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
