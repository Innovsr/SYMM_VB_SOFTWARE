!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! All possible Covalent structures are generated here it has three parts 
!!! lone-paires, covalent-bonds, and unpaired or radicals, each part generated 
!!! separately 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module covalent_str
use commondat
use quality
use symm_cal_pi
use symm_cal_sig
use str_select
implicit none

contains
subroutine cov_struc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!common/quality/str_quality_1,str_quality_2,bondq,tqlty,bqlty,sqlty,tnqs,nssym,qulsym,symq,&
!sigsym,tnqs_sig

integer::ijk,ij,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,perm_nstr,totset,&
k1,k2,k4,k5,a,d,e,f,nnn,elporb,x,alstr,y
integer::factorial,i,j,k,l,j1,m,j2,ii,j3,j4,wig2,tnqs,lonep,bonds,allp, c
!integer::symq(15000),sigsym(15000),tnqs_sig,nssym&
!,qulsym(15000),str_quality_1(15000),str_quality_2(15000),bondq(15000),tqlty,bqlty,sqlty
!integer,dimension(:,:),allocatable::strc,fstr,strct,num
double precision, pointer::symsc(:)
integer, pointer::fullcovstr(:, :)
integer,allocatable::n(:), nn(:), num(:,:), strc(:,:), strct(:,:)
character(len=10)::alstr_string,allowed_string
!!! memory allocation !!!!!!!
print*,'enter cov_str'
x=100000
y=100

allocate(strct(x,y))
allocate(num(x,y))

allocate(n(MaxStrOepo))
allocate(nn(MaxStrOepo))

!!! initialization of arrays !!!
strct = 0
n = 0
nn=0
num=0

!!! ionic and covalent flags
flg_ion=0
flg_cov=1

!!!! new structures generation part starts !!!!!
vacorb=nao-nae
lonep=nae-nao
if(vacorb.gt.0)lonep=0
bonds=(nae-lonep*2-nlast)/2
d=nao
e=2
f=d-e
c=factorial(d)/(factorial(e)*factorial(f))

allocate(strc(c, 2))
strc = 0
!!!! production of the set of bonded orbitals strats !!!!!
i4=0
do i1=1,nao-1
  do i2=i1+1,nao
    i4=i4+1
    do i3=1,2
    i5=i2
      if(i3.eq.1)i5=i1
    strc(i4,i3)=i5
    n(i4)=i1
    enddo
  enddo
enddo
totset=i4

print*,'totset',totset
!!!! production of the set of bonded orbitals ends !!!!!
!!!! production of the covalent bonding part start !!!!!
!!!! nested loops generates all combinations of bons !!!
!!!! for the covalent past of the structures; 7 nested !
!!!! loops can generate upto seven bond structures !!!!!

i=0
j=0
! loop starts
do i1=1,totset
  j=1
  nn(1)=i1
  if(j.eq.bonds)then
  i=i+1
  k4=0
    do k2=1,bonds
      do k5=1,2
      k4=k4+1
      strct(i,k4)=strc(nn(k2),k5)
      enddo
    enddo
  cycle 
  endif

  ! loop 2
  do i2=i1+1,totset
    j=2
    l=0
  
    do k=1,2
      if(strc(nn(1),k).eq.n(i2))then
      l=l+1
      endif
    enddo
  
    if(l.eq.0) then
      l=0
      do k=1,2
        do k1=1,2
          if(strc(i2,k).eq.strc(nn(1),k1))then
            l=l+1
          endif
        enddo
      enddo
  
      if(l.eq.0) then
        nn(2)=i2
        if(j.eq.bonds)then
          i=i+1
          k4=0
          do k2=1,bonds
            do k5=1,2
              k4=k4+1
              strct(i,k4)=strc(nn(k2),k5)
            enddo
          enddo
          cycle
        endif
         
      else
        cycle
      endif
  
    else
      cycle
    endif
  
    ! loop 3
    do i3=i2+1,totset
      j=3
      l=0
      do k1=1,j-1
        do k=1,2
          if(strc(nn(k1),k).eq.n(i3))then
            l=l+1
          endif
        enddo
      enddo
      if(l.eq.0) then
        l=0
        do k2=1,j-1
          do k=1,2
            do k1=1,2
              if(strc(i3,k).eq.strc(nn(k2),k1))then
                l=l+1
              endif
            enddo
          enddo
        enddo
        if(l.eq.0) then
          nn(3)=i3
          if(j.eq.bonds)then
            i=i+1
            k4=0
            do k2=1,bonds
              do k5=1,2
                k4=k4+1
                strct(i,k4)=strc(nn(k2),k5)
              enddo
            enddo
            cycle
          endif
    
          else
            cycle
          endif
    
        else
          cycle
        endif
    
      ! loop 4
      do i4=i3+1,totset
        j=4
        l=0
        do k1=1,j-1
          do k=1,2
            if(strc(nn(k1),k).eq.n(i4))then
              l=l+1
            endif
          enddo
        enddo
        if(l.eq.0) then
          l=0
          do k2=1,j-1
            do k=1,2
              do k1=1,2
                if(strc(i4,k).eq.strc(nn(k2),k1))then
                  l=l+1
                endif
              enddo
            enddo
          enddo
          if(l.eq.0) then
            nn(4)=i4
            if(j.eq.bonds)then
              i=i+1
              k4=0
              do k2=1,bonds
                do k5=1,2
                  k4=k4+1
                  strct(i,k4)=strc(nn(k2),k5)
                enddo
              enddo
              cycle
            endif
            else
              cycle
            endif
      
          else
            cycle
          endif
      
        !loop 5
        do i5=i4+1,totset
          j=5
          l=l+1
          do k1=1,j-1
            do k=1,2
              if(strc(nn(k1),k).eq.n(i5))then
                l=l+1
              endif
            enddo
          enddo
          if(l.eq.0) then
            l=0
            do k2=1,j-1
              do k=1,2
                do k1=1,2
                  if(strc(i5,k).eq.strc(nn(k2),k1))then
                    l=l+1
                  endif
                enddo
              enddo
            enddo
            if(l.eq.0)then
              nn(5)=i5
              if(j.eq.bonds)then
                i=i+1
                k4=0
                do k2=1,bonds
                  do k5=1,2
                    k4=k4+1
                    strct(i,k4)=strc(nn(k2),k5)
                  enddo
                enddo
                cycle
              endif
              else
                cycle
              endif
        
            else
              cycle
            endif
        
          !loop 6
          do i6=i5+1,totset
            j=6
            l=0
            do k1=1,j-1
              do k=1,2
                if(strc(nn(k1),k).eq.n(i6))then
                  l=l+1
                endif
              enddo
            enddo
            if(l.eq.0) then
              l=0
              do k2=1,j-1
                do k=1,2
                  do k1=1,2
                    if(strc(i6,k).eq.strc(nn(k2),k1))then
                      l=l+1
                    endif
                  enddo
                enddo
              enddo
              if(l.eq.0) then
                nn(6)=i6
                if(j.eq.bonds)then
                  i=i+1
                  k4=0
                  do k2=1,bonds
                    do k5=1,2
                      k4=k4+1
                      strct(i,k4)=strc(nn(k2),k5)
                    enddo
                  enddo
                  cycle
                endif
          
                else
                  cycle
                endif
          
              else
                cycle
              endif
          
            ! loop7
            do i7=i6+1,totset
              j=7
              l=0
              do k1=1,j-1
                do k=1,2
                  if(strc(nn(k1),k).eq.n(i7))then
                    l=l+1
                  endif
                enddo
              enddo
              if(l.eq.0) then
                l=0
                do k2=1,j-1
                  do k=1,2
                    do k1=1,2
                      if(strc(i7,k).eq.strc(nn(k2),k1))then
                        l=l+1
                      endif
                    enddo
                  enddo
                enddo
                if(l.eq.0) then
                  nn(7)=i7
                  if(j.eq.bonds)then
                    i=i+1
                    k4=0
                    do k2=1,bonds
                      do k5=1,2
                        k4=k4+1
                        strct(i,k4)=strc(nn(k2),k5)
                      enddo
                    enddo
                    cycle
                  endif
                  else
                    cycle
                  endif
            
                else
                  cycle
                endif


            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
enddo
!!!! production of the covalent bonding part of the structures ends !!!!!
!!!! bond part of the structures stored in 'strct' matrix

deallocate(strc)
alstr=i ! total number of possible structures revealed

allocate(strc(alstr, nae))
strc=0
n=0
nn=0


do i1=1,alstr
  do k4=1,bonds*2
    strc(i1,k4)=strct(i1,k4)
  enddo
enddo

!!!!! production of the lone pairs and radical part of the structures starts !!!!!!
if(nlast.ne.0.or.lonep.ne.0)then
  i=0
  do ii=1,alstr
    k5=0
    i2=0
     loop19:do i1=1,nao
       do k1=1,bonds*2
         if(i1.eq.strct(ii,k1)) cycle loop19
       enddo
       i2=i2+1
       num(ii,i2)=i1
     enddo loop19
     allp=i2
!!!!! production of the lone pair part of the structures starts !!!!!!
!!!!! 10 nested loops can creat 10 lone paires in each structures !!!!
!!!!! due to combinations of the lone paires the number of structures!
!!!!! may increase.

    if(lonep.ne.0)then
    j=0
    m=0
    
    ! loop 1
      loop20:do i1=1,allp
      j=1
      n(1)=i1
        if (j.eq.lonep) then
        i=i+1
        j4=0
          do j1=lonep,1,-1
            do j2=1,2
            j4=j4+1
            strc(i,j4)=num(ii,n(j1))
            enddo
          enddo
        j2=0
          do j3=1,bonds*2
          j2=j3+j4
          strc(i,j2)=strct(ii,j3)
          enddo
        cycle loop20
        endif
      
        ! loop 2
        loop21:do i2=i1+1,allp
          j=2
          if(i2.eq.n(1)) then
            cycle loop21
          endif
            n(2)=i2
            if (j.eq.lonep) then
              i=i+1
              j4=0
              do j1=1,lonep
                do j2=1,2
                j4=j4+1
                strc(i,j4)=num(ii,n(j1))
                enddo
              enddo
              j2=0
              do j3=1,bonds*2
                j2=j3+j4
                strc(i,j2)=strct(ii,j3)
              enddo
              cycle loop21
            endif
        
          ! loop 3
          loop22:do i3=i2+1,allp
            j=3
            do k=1,2
              if(i3.eq.n(k)) then
                cycle loop22
              endif
            enddo
            n(3)=i3
            if (j.eq.lonep) then
              i=i+1
              j4=0
              do j1=1,lonep
                do j2=1,2
                  j4=j4+1
                  strc(i,j4)=num(ii,n(j1))
                enddo
              enddo
              j2=0
              do j3=1,bonds*2
                j2=j3+j4
                strc(i,j2)=strct(ii,j3)
              enddo
              cycle loop22
            endif
          
            ! loop 4
            loop23:do i4=i3+1,allp
              j=4
              do k=1,3
                if(i4.eq.n(k)) then 
                  cycle loop23
                endif
              enddo
              n(4)=i4
              if (j.eq.lonep) then
                i=i+1
                j4=0
                do j1=1,lonep
                  do j2=1,2
                    j4=j4+1
                    strc(i,j4)=num(ii,n(j1))
                  enddo
                enddo
                j2=0
                do j3=1,bonds*2
                  j2=j3+j4
                  strc(i,j2)=strct(ii,j3)
                enddo
                cycle loop23
              endif
            
              !loop 5
              loop24:do i5=i4+1,allp
                j=5
                do k=1,4
                  if(i5.eq.n(k)) then
                    cycle loop24
                  endif
                enddo
                n(5)=i5
                if (j.eq.lonep) then
                  i=i+1
                  j4=0
                  do j1=1,lonep
                    do j2=1,2
                      j4=j4+1
                      strc(i,j4)=num(ii,n(j1))
                    enddo
                  enddo
                j2=0
                  do j3=1,bonds*2
                    j2=j3+j4
                    strc(i,j2)=strct(ii,j3)
                  enddo
                  cycle loop24
                endif
              
                ! loop 6
                loop25:do i6=i5+1,allp
                  j=6
                  do k=1,5
                    if(i6.eq.n(k)) then
                      cycle loop25
                    endif
                  enddo
                  n(6)=i6
                  if (j.eq.lonep) then
                    i=i+1
                    j4=0
                    do j1=1,lonep
                      do j2=1,2
                        j4=j4+1
                        strc(i,j4)=num(ii,n(j1))
                      enddo
                    enddo
                    j2=0
                    do j3=1,bonds*2
                      j2=j3+j4
                      strc(i,j2)=strct(ii,j3)
                    enddo
                    cycle loop25
                  endif
                
                  !loop 7
                  loop26:do i7=i6+1,allp
                    j=7
                    do k=1,6
                      if(i7.eq.n(k)) then
                        cycle loop26
                      endif
                    enddo
                    n(7)=i7
                    if (j.eq.lonep) then
                      i=i+1
                      j4=0
                      do j1=1,lonep
                        do j2=1,2
                        j4=j4+1
                        strc(i,j4)=num(ii,n(j1))
                        enddo
                      enddo
                      j2=0
                      do j3=1,bonds*2
                        j2=j3+j4
                        strc(i,j2)=strct(ii,j3)
                      enddo
                      cycle loop26
                    endif
                  
                    loop27:do i8=i7+1,allp
                      j=8
                      do k=1,7
                        if(i8.eq.n(k)) then
                          cycle loop27
                        endif
                      enddo
                      n(8)=i8
                      if (j.eq.lonep) then
                        i=i+1
                        j4=0
                        do j1=1,lonep
                          do j2=1,2
                            j4=j4+1
                            strc(i,j4)=num(ii,n(j1))
                          enddo
                        enddo
                        j2=0
                        do j3=1,bonds*2
                          j2=j3+j4
                          strc(i,j2)=strct(ii,j3)
                        enddo
                        cycle loop27
                      endif
                    
                      ! loop 9
                      loop28:do i9=i8+1,allp
                        j=9
                        do k=1,8
                          if(i9.eq.n(k)) then
                            cycle loop28
                          endif
                        enddo
                        n(9)=i9
                        if (j.eq.lonep) then
                          i=i+1
                          j4=0
                          do j1=1,lonep
                            do j2=1,2
                              j4=j4+1
                              strc(i,j4)=num(ii,n(j1))
                            enddo
                          enddo
                          j2=0
                          do j3=1,bonds*2
                            j2=j3+j4
                            strc(i,j2)=strct(ii,j3)
                          enddo
                          cycle loop28
                        endif
                      
                        ! loop 10
                        loop29:do i10=i9+1,allp
                          j=10
                          do k=1,9
                            if(i10.eq.n(k)) then
                              cycle loop29
                            endif
                          enddo
                          n(10)=i10
                          if (j.eq.lonep) then
                            i=i+1
                            j4=0
                            do j1=1,lonep
                              do j2=1,2
                                j4=j4+1
                                strc(i,j4)=num(ii,n(j1))
                              enddo
                            enddo
                            j2=0
                            do j3=1,bonds*2
                              j2=j3+j4
                              strc(i,j2)=strct(ii,j3)
                            enddo
                            cycle loop29
                          endif
                        
                        enddo loop29
                      enddo loop28
                    enddo loop27
                  enddo loop26
                enddo loop25
              enddo loop24
            enddo loop23
          enddo loop22
        enddo loop21
      enddo loop20
    endif
  enddo
  if (lonep.ne.0) alstr=i ! new number of structures revealed
  !!!!! production of the radical part of the structures starts !!!!!!
  if(nlast.ne.0)then
    do i=1,alstr
      j4=nae-nlast
      loop31:do i1=1,nae
        do i2=1,nae-nlast
          if(strc(i,i2).eq.i1) cycle loop31
        enddo
        j4=j4+1
        strc(i,j4)=i1
      enddo loop31
    enddo
  endif
endif
!!!!! production of the lone pairs and radical part of the structures ends !!!!!!
deallocate(num)
!990 format(30I3)
deallocate(strct)
!!!!!!! new structures generation part ends !!!!!!!!!!!!!!!!!!!!!!!
d=nao
e=nlp
f=nao-nlp
elporb=nao-nlp
c=factorial(d)/(factorial(e)*factorial(f))
call wigner(elporb,wig2)
deallocate(n)
deallocate(nn)

!!! header of the 'structures.dat' file is written here
write(7,*)'**********************************************************'
write(7,*)'*      All possible structures with their qualities      *'
write(7,*)'**********************************************************'
write(7,*)'  '
write(7,*)'----------------------------------------------------------------------------------------------'
write(7,309)'*','active orbs =',nao,'active electrons =',nae,'multiplicity =',mult,'inactive orbs =',niao
write(7,*)'----------------------------------------------------------------------------------------------'
write(7,*)

WRITE(alstr_string, '(I0)') alstr
!    string_value = ADJUSTL(atstr_string)

WRITE(allowed_string, '(I0)') wig2*c
!    string_value = ADJUSTL(atstr_string)

write(7,307)'Total number of covalent structure of the system = ',trim(alstr_string)
write(7,307)'Total number of allowed covalent structres are = ',trim(allowed_string)
write(7,*)' '
write(7,308)'Sl No.','various qualities','Overall qualities','Rumer/Non Rumer', 'Structures'
write(7,306)'IAB','NAB','SBB','PDB','PDR'
307 format(10x,a,a)
308 format(2x,a,5x,a,5x,a,3x,a,10x,a)
306 format(8x,a,1x,a,1x,a,1x,a,1x,a)
309 format(a,1x,a,1x,I4,3x,a,1x,I4,3x,a,1x,I4,3x,a,1x,I4)

allocate (fullcovstr(alstr,nae))

do i=1,alstr
  do i1=1,nae
    fullcovstr(i,i1)=strc(i,i1)+niao
  enddo
enddo

if(symm.eq.1) then
  !allocate(symsc(MaxStrOepo))
  if(sig_sym_flg.eq.1)call symmetry_cal_sig(nlp,fullcovstr,alstr,symsc,symq,nssym)
  if(sig_sym_flg.ne.1)call symmetry_cal_pi(nlp,fullcovstr,alstr,symsc,symq,nssym)

  nnn=0
  do ij=1,tnqs
    do ijk=nssym,1,-1
      do j2=1,alstr
        if(symq(j2).eq.ijk)then
          nnn=nnn+1
          if(niao.eq.0)then
            write(7,900)'cov structure',nnn,')',qulsym(j2),(fullcovstr(j2,j1),j1=1,nao-nlp+nlp*2)
          endif
          if(niao.ne.0)then
            write(7,901)'cov structure',nnn,')',qulsym(j2),1,':',niao,(fullcovstr(j2,j1),j1=1,nao-nlp+nlp*2)
          endif
        endif
      enddo
      write(7,*)'******************************************************'
    enddo
  enddo
else
  do j2=1,alstr
    if(niao.eq.0)then
      write(7,902)'cov structure',j2,')',(fullcovstr(j2,j1),j1=1,nao-nlp+nlp*2)
    endif
    if(niao.ne.0)then
      write(7,903)'cov structure',j2,')',1,':',niao,(fullcovstr(j2,j1),j1=1,nao-nlp+nlp*2)
    endif
  enddo
  write(7,*)'******************************************************'
endif
!print*,alstr, nae
!do i=1,alstr
!print*,'fullcovstr',(fullcovstr(i,j),j=1,nae)
!enddo
900 format(a,I5,a,x,I3,3x,30(I5))
901 format(a,I5,a,x,I3,3x,I3,x,a,I3,30(I5))
902 format(a,I5,a,x,30(I5))
903 format(a,I5,a,x,I3,x,a,I3,30(I5))


if(nfset.eq.3)write(23,*)'********** RUMER STRUCTURES ******************'
if(nfset.eq.5)write(10,*)'**********ALL POSSIBLE RUMER STRUCTURES ******************'
if(flg1.eq.1)then
write(10,*)'********** RUMER STRUCTURES ******************'
endif
if(flg1.eq.0.and.nfset.ne.5)then
write(10,*)'******* CHEM. QUAL. STRUCTURES **************'
endif

perm_nstr=wig2*c
write(10,*)perm_nstr,' covalent structures' 
write(23,*)perm_nstr,' covalent structures' 
if(nfset.ne.5)write(10,913)'IAB','NNB','SBB','fqual'
write(23,913)'IAB','NNB','SBB','fqual'
913 format(x,a,x,a,x,a,x,a)
deallocate(strc)
call str_selection(fullcovstr,nlp,alstr,perm_nstr)

!deallocate(fullcovstr)

print*,'exit cove_struc'
return
end subroutine cov_struc
end module covalent_str
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


