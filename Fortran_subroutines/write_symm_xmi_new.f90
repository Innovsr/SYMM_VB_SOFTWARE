!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! All possible sets containg various linearly independent symmetric structures groups generated here.!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_symm_xmi_new(nl,strn,str3,ncqs,q_fac2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

common/quality/str_quality_1,str_quality_2,bondq,tqlty,bqlty,sqlty,tnqs,nssym,qulsym,symq,&
sigsym,tnqs_sig
common/infosymm/set_number,hqlty,ttqlty0,ttqlty1,ttqlty2,ttqlty3,ncqss,qul,Rid&
,mns,u1,max_set,rumset,nqset,strset,strnn,totstr,strno,incmplt,mincmplt,mincmplt_set,group_num
common/str/str5,nstr7

integer::nl,strn,strnn,ncqs,ncqss,tostr,initstr,i,i1,i2,i3,i4,i5,i6,i7,i8,i9,m119,m18,m19,m20,m21,m23,m24,count&
,qul(100),nqul,j,jj,jjj,fg,flg,ii5,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,x,y,strs
integer::i5up16,i16i7,m16m21,i5up15,i15i7,m15m21,i5up14,i14i7,m14m21,i5up13,i13i7,m13m21,i5up12,i12i7,m12m21,&
i5up11,i11i7,m11m21,i5up10,i10i7,m10m21,i5up9,i9i7,m9m21,i5up8,i8i7,m8m21,i5up7,i7i7,m7m21,i5up6,i6i7,m6m21,&
i5up5,i5i7,m5m21,i5up4,i4i7,m4m21,i5up3,i3i7,m3m21,i5up2,i2i7,m2m21,i5up1,i1i7,m1m21,i5up17,i17i7,m17m21
integer::str3(15000,20),q_fac2(15000),finalvec(15000),strset(1000),col(1000),sigsym(15000),tnqs_sig,&
ffvec2(15000,1000),bondq(15000),bondq4(15000),nqset(15000),str5(2000,20),nstr7,qfac3(15000),group_num(15000)&
,tndet,totstr,Ifail,indpnt,strno(1000),str_quality_1(15000),str_quality_2(15000),ttqlty0,ttqlty&
,tqlty,bqlty,sqlty,hqlty,tnqs,nssym,qulsym(15000),symq(15000),set_number,ttqlty1,det_inv,ttqlty2,ttqlty3
integer::rumer(15000),rumer_rad(15000),quality_fac(15000),rumset,u1,max_set,mns,mincmplt,mincmplt_set(500)
!integer,dimension(:),allocatable::qq1,qq2,qq
integer::qq1(5000),qq2(5000),qq(5000),str2(2000,20),Rid,set_num(100),sf1,sf2,incmplt
real*8::ovlp
Double Precision::D(1000)
character(10)::dd,a
character(len=100)::outfile


print*,'enter_symm_xmi_new'

Rid=0 ! Rumer id; to varify the set is Rumer of not. Rid=1 call the Rumer sets

! if the loops could not able to find out a set it will provide the biggest incomplete set
! which is counted and stored with 'mincmplt' and mincmplt_set array.
mincmplt=0
do i=1,500
mincmplt_set(i)=0
enddo

ncqss=ncqs ! total number of permiissible structures in a set
strnn=strn ! total number of structures

incmplt=1 ! running variable to count maximum number of independent structures in each set
mns=0
u1=1

max_set=75000 ! maximum number of set will be written in one file

if(nfset.eq.3.or.nfset.eq.5)then
rumset=0
call rumer_structures(nl,str3,ncqss,rumer,rumer_rad)
call write_rumer_xmi(nl,str3,ncqss,rumer,rumer_rad,quality_fac)
endif

!various keywds
set_number=0
bqlty=0
tqlty=0
sqlty=0
hqlty=0
indpnt=2
!ovlpval=1.0
if(noq0.gt.strnn)then
ttqlty0=noq0
else
ttqlty0=strnn+noq0
endif
ttqlty1=strnn+noq1
ttqlty2=strnn+noq2
ttqlty3=strnn+noq3

write(*,*)'sl  structures           group_numbers'
do i=1,ncqss
write(*,231)i,(str3(i,j),j=1,nae),q_fac2(i)
!write(*,231),i,(str3(i,j),j=1,nae),q_fac2(i),str_quality_1(i),str_quality_2(i),bondq(i)
group_num(i)=q_fac2(i)
enddo
231 format(30I3)

jj=1
loop1:do m19=1,ncqss
  if(m19.eq.1)qul(1)=q_fac2(1)
    j=jj
    loop2:do i=1,j
      if(qul(i).eq.q_fac2(m19)) then
        cycle loop1
      endif
    enddo loop2
    jj=jj+1
    qul(i)=q_fac2(m19)
enddo loop1
nqul=jj
!print*,'nqul',nqul

! counting of groups and the number of structures in each group
do i=1,nqul
jjj=0
jj=0
do m19=1,ncqss
!print*,qul(i),q_fac2(m19)
if(qul(i).eq.q_fac2(m19))then
jjj=jjj+1
jj=m19
endif
enddo
nqset(i)=jjj
strset(i)=jj
!print*,'i,nqset(i),strset(i)',i,nqset(i),strset(i)
enddo

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

!p rint*,'loop1'
   call write_symm_xmi_1(i7,m21,m1,sf1,sf2,str3,q_fac2)
   if (sf1.eq.1) cycle
   if (sf2.eq.1) return

   i2i7=i7
   m2m21=m21

   do m2=m1+1,nqul
      i7=i2i7
      m21=m2m21

      call write_symm_xmi_1(i7,m21,m2,sf1,sf2,str3,q_fac2)
      if (sf1.eq.1) cycle
      if (sf2.eq.1) return
      
      i3i7=i7
      m3m21=m21

      do m3=m2+1,nqul
         i7=i3i7
         m21=m3m21
         
         call write_symm_xmi_1(i7,m21,m3,sf1,sf2,str3,q_fac2)
         if (sf1.eq.1) cycle
         if (sf2.eq.1) return
         
         i4i7=i7
         m4m21=m21

         do m4=m3+1,nqul
            i7=i4i7
            m21=m4m21
            
            call write_symm_xmi_1(i7,m21,m4,sf1,sf2,str3,q_fac2)
            if (sf1.eq.1) cycle
            if (sf2.eq.1) return
            
            i5i7=i7
            m5m21=m21

            do m5=m4+1,nqul
               i7=i5i7
               m21=m5m21
               
               call write_symm_xmi_1(i7,m21,m5,sf1,sf2,str3,q_fac2)
               if (sf1.eq.1) cycle
               if (sf2.eq.1) return
               
               i6i7=i7
               m6m21=m21

               do m6=m5+1,nqul
                  i7=i6i7
                  m21=m6m21
                  
                  call write_symm_xmi_1(i7,m21,m6,sf1,sf2,str3,q_fac2)
                  if (sf1.eq.1) cycle
                  if (sf2.eq.1) return
                  
                  i7i7=i7
                  m7m21=m21

                  do m7=m6+1,nqul
                     i7=i7i7
                     m21=m7m21
                     
                     !print*,'loop7'
                     call write_symm_xmi_1(i7,m21,m7,sf1,sf2,str3,q_fac2)
                     if (sf1.eq.1) cycle
                     if (sf2.eq.1) return
                     
                     i8i7=i7
                     m8m21=m21

                     do m8=m7+1,nqul
                        i7=i8i7
                        m21=m8m21
                        
                        call write_symm_xmi_1(i7,m21,m8,sf1,sf2,str3,q_fac2)
                        if (sf1.eq.1) cycle
                        if (sf2.eq.1) return
                        
                        i9i7=i7
                        m9m21=m21

                        do m9=m8+1,nqul
                           i7=i9i7
                           m21=m9m21
                           
                           !print*,'loop9'
                           call write_symm_xmi_1(i7,m21,m9,sf1,sf2,str3,q_fac2)
                           if (sf1.eq.1) cycle
                           if (sf2.eq.1) return
                           
                           i10i7=i7
                           m10m21=m21

                           do m10=m9+1,nqul
                              i7=i10i7
                              m21=m10m21
                              
                              !print*,'loop10'
                              call write_symm_xmi_1(i7,m21,m10,sf1,sf2,str3,q_fac2)
                              if (sf1.eq.1) cycle
                              if (sf2.eq.1) return
                              
                              i11i7=i7
                              m11m21=m21

                              do m11=m10+1,nqul
                                 i7=i11i7
                                 m21=m11m21
                                 
                                 !print*,'loop11'
                                 call write_symm_xmi_1(i7,m21,m11,sf1,sf2,str3,q_fac2)
                                 if (sf1.eq.1) cycle
                                 if (sf2.eq.1) return
                                 
                                 i12i7=i7
                                 m12m21=m21

                                 do m12=m11+1,nqul
                                    i7=i12i7
                                    m21=m12m21
                                    
                                    print*,'loop12'
                                    call write_symm_xmi_1(i7,m21,m12,sf1,sf2,str3,q_fac2)
                                    if (sf1.eq.1) cycle
                                    if (sf2.eq.1) return
                                    
                                    i13i7=i7
                                    m13m21=m21

                                    do m13=m12+1,nqul
                                       i7=i13i7
                                       m21=m13m21
                                       
                                       print*,'loop13'
                                       call write_symm_xmi_1(i7,m21,m13,sf1,sf2,str3,q_fac2)
                                       if (sf1.eq.1) cycle
                                       if (sf2.eq.1) return
                                       
                                       i14i7=i7
                                       m14m21=m21

                                       do m14=m13+1,nqul
                                          i7=i14i7
                                          m21=m14m21
                                          
                                          call write_symm_xmi_1(i7,m21,m14,sf1,sf2,str3,q_fac2)
                                          if (sf1.eq.1) cycle
                                          if (sf2.eq.1) return
                                          
                                          i15i7=i7
                                          m15m21=m21

                                          do m15=m14+1,nqul
                                             i7=i15i7
                                             m21=m15m21
                                             
                                             print*,'loop15'
                                             call write_symm_xmi_1(i7,m21,m15,sf1,sf2,str3,q_fac2)
                                             if (sf1.eq.1) cycle
                                             if (sf2.eq.1) return
                                             
                                             i16i7=i7
                                             m16m21=m21

                                             do m16=m15+1,nqul
                                                i7=i16i7
                                                m21=m16m21
                                                
                                                print*,'loop16'
                                                call write_symm_xmi_1(i7,m21,m16,sf1,sf2,str3,q_fac2)
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


close(21)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! if the above loop become unsuccessfull to provide a full dimention set the below part wil then!
!! print the maximum generated set in the output                                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(incmplt.eq.1)then
write(10,*)'----------- incomplete set -----------'
do i=1,mincmplt
 qq(i)=q_fac2(mincmplt_set(i))
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
CALL SYSTEM ("rm Rumer_Sets.dat")


100 return
end subroutine write_symm_xmi_new
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
