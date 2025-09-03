!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! All possible Covalent structures are generated here it has three parts 
!!! lone-paires, covalent-bonds, and unpaired or radicals, each part generated 
!!! separately 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module mod_cov_struc
use commondat_mod
use quality_mod
use str_module
!use mod_symmetry_cal_sig
use mod_str_selection
implicit none

contains
subroutine cov_struc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer::ijk,ij,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,perm_nstr,totset
integer::k1,k2,k4,k5,a,d,e,f,nnn,elporb,x,alstr,y,alstr_new, comb, max_val
integer::factorial,i,j,k,l,j1,m,j2,ii,j3,j4,wig2,tnqs,lonep,bonds,allp, c, c1
double precision, pointer::symsc(:)
integer, pointer::fullcovstr(:, :)
integer, allocatable::n(:),  strc(:,:)
!integer::nn(1500), num(1500,10), strct(1500,10)
integer, allocatable::nn(:), num(:,:), strct(:,:)
integer, allocatable::lonep_mat(:,:),lonep_score(:), strc_new(:,:)
character(len=10)::alstr_string,allowed_string
print*,'enter cov_str',MaxStrOepo, nae

!print*,'sourav1'
allocate(num(MaxStrOepo, nae))
!print*,'sourav2'
!
allocate(nn(MaxStrOepo))
!print*,'sourav3'
allocate(strct(MaxStrOepo, nae))

!!! initialization of arrays !!!
strct = 0
nn=0
num=0

!!!! new structures generation part starts !!!!!
vacorb=nao-nae
lonep=nlp
if(vacorb.gt.0)lonep=0
bonds=(nae-lonep*2-nlast)/2
d=nao
e=2
f=d-e
c=factorial(d)/(factorial(e)*factorial(f))
print*,'bonds',bonds

allocate(strc(c, 2))
allocate(n(c))
n = 0
strc = 0
!!!! production of the set of bonded orbitals strats !!!!!
if (bonds.ne.0) then
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
do i=1,i4
  print*,'strc_bonded_orbs',(strc(i,i1),i1=1,2)
enddo

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
print*,'alstr',alstr
!allocate(strc(alstr+alstr*(nae-bonds*2)*(nlast+lonep), nae))
allocate(strc(MaxStrOepo, nae))
endif

if(bonds.eq.0)then
  deallocate(strc)
  !alstr = comb(nae,lonep)
  alstr = 1
  !!allocate(strc(comb(nae,lonep)*(nlast+lonep), nae))
  allocate(strc(MaxStrOepo, nae))
endif

strc=0
n=0
nn=0

if(bonds.ne.0) then
  print*,'shape',shape(strc)
do i1=1,alstr
  do k4=1,bonds*2
    strc(i1,k4)=strct(i1,k4)
  enddo
enddo
endif

do i1=1,alstr
  print*,'bonded_parts',(strc(i1,k4),k4=1,bonds*2)
enddo

!!!!! production of the lone pairs and radical part of the structures starts !!!!!!
if(nlast.ne.0.or.lonep.ne.0)then
  i=0
  do ii=1,alstr
    print*,'alstr',ii
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
          print*,'strc..',strc(i,j2)
          enddo
        cycle loop20
        endif
      
        ! loop 2
        print*,'loop2'
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
                  print*,'strc',i,strc(i,j4)
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

  !do i=1,alstr
  !print*,'strc',(strc(i,i2),i2=1,nae)
  !enddo

  !! sorting of the structures according to it's lone-pairs
  !! using s.Roy's prioratise integer merging formula.
  if (lonep.gt.0) then
    allocate(lonep_mat(lonep+1, alstr))
    lonep_mat=0 ! matrix to store lonepairs for all structures along row

    do j=1,alstr
        lonep_mat(1,j) = lonep_mat(1,j) + strc(j,1) ! 1st row for 1st lone paire 
      do i=1,lonep
        lonep_mat(i+1,j) = lonep_mat(i+1,j) + strc(j,i+(i-1)) !1st lone pair comes towice
      enddo                                                   !1st and 2nd rows, it generalise
    enddo                                                     !for all system:: 1st lp merges 
                                                              !with it's own.
    !do i=1,lonep+1
    !print*,'lonep_mat',(lonep_mat(i,j),j=1,alstr)
    !enddo

    !! lower lonepair has the priority. if there are three lonepairs the sorting 
    !! process will put more priority to the 1st lone-pair, then the 2nd lone-pair
    !! and then the third lone-pair 
    do i=1,lonep
      max_val=1
      do j=1,alstr
        if(max_val.lt.lonep_mat(i+1,j))max_val=lonep_mat(i+1,j) ! finding highest value of less  
      enddo                                                     ! priority lone-pair's
      do j=1,alstr
        lonep_mat(i+1,j)=max_val*(lonep_mat(i,j)-1)+lonep_mat(i+1,j)
      enddo
    enddo
    allocate(lonep_score(alstr))

    do i=1,alstr
      lonep_score(i)=lonep_mat(lonep,i)
    enddo
    deallocate(lonep_mat)

    max_val=1
    do i=1,alstr
      if(max_val.lt.lonep_score(i))max_val=lonep_score(i)
    enddo

    allocate(strc_new(alstr,nae))

    k=0
    do i=1,max_val
      do j=1,alstr
        if(i.eq.lonep_score(j))then
          k=k+1
          do i1=1,nae
            strc_new(k,i1)=strc(j,i1)
          enddo
        endif
      enddo
    enddo
    do i=1,alstr
      do j=1,nae
        strc(i,j)=strc_new(i,j)
      enddo
    enddo

    deallocate(strc_new)
    deallocate(lonep_score)

  endif


  do i=1,alstr
  print*,'strc**',i,(strc(i,i2),i2=1,nae)
  enddo
  deallocate(strct)

  !!!!! construction of the radical part of the structures starts !!!!!!
  if(nlast.ne.0)then
    allocate(strct(alstr*(nae-bonds*2)*(nlast+lonep), nae))
    strct = strc
    strc = 0
    alstr_new = 0
    do i=1,alstr
      j4=nae-nlast
      loop31:do i1=1,nae
        do i2=1,nae-nlast
          if(strct(i,i2).eq.i1) cycle loop31
        enddo
        j4=j4+1
        alstr_new = alstr_new + 1
        do i2 = 1, nae
           strc(alstr_new,i2) = strct(i,i2)
        enddo
        strc(alstr_new,j4)=i1
        if (j4.eq.nae) j4 = nae-nlast
        print*,'sets',alstr_new,(strc(alstr_new,i2),i2=1,nae)
      enddo loop31

    enddo
    alstr = alstr_new
  endif
endif
!!!!! production of the lone pairs and radical part of the structures ends !!!!!!
deallocate(num)
!990 format(30I3)
!!!!!!! new structures generation part ends !!!!!!!!!!!!!!!!!!!!!!!
!c=1
d=nao
e=nao-nlp
elporb=nae-nlp*2
print*,'elporb',elporb,nlp,d,e
!c=factorial(d)/(factorial(e)*factorial(f))
c = comb(d, e)

print*,'cccc',c
!!! for ionic structure increasing number of nlp, this combination
!!! multiplied for extra radicals 
d=nao-nlp
e=nae-2*nlp
c1=comb(d, e)
call wigner(elporb,wig2)
print*,'wig',c,c1,wig2
deallocate(n)
deallocate(nn)
print*,'sourav1'

!!! header of the 'structures.dat' file is written here
if(prt_cnt.eq.1) then
  write(7,*)'**********************************************************'
  write(7,*)'*      All possible structures with their qualities      *'
  write(7,*)'**********************************************************'
  write(7,*)'  '
  write(7,*)'----------------------------------------------------------------------------------------------'
  write(7,309)'*','active orbs =',nao,'active electrons =',nae,'multiplicity =',mult,'inactive orbs =',niao
  write(7,*)'----------------------------------------------------------------------------------------------'
  write(7,*)
endif

print*,'sourav2'
WRITE(alstr_string, '(I0)') alstr
!    string_value = ADJUSTL(atstr_string)

WRITE(allowed_string, '(I0)') wig2*c*c1
!    string_value = ADJUSTL(atstr_string)

if (flg_cov.eq.1) then
  write(7,307)'Total number of covalent structure of the system = ',trim(alstr_string)
  write(7,307)'Total number of allowed covalent structres are = ',trim(allowed_string)
endif
if (flg_ion.eq.1) then
  write(7,310)'Total number of ionic structure of the system = ',trim(alstr_string),&
          'with number of lone-pair = ',nlp
  write(7,307)'Total number of allowed ionic structres are = ',trim(allowed_string)
endif
write(7,*)' '
write(7,308)'Sl No.','various qualities','Overall qualities','Rumer/Non Rumer', 'Structures'
write(7,306)'IAB','NAB','SBB','PDB','PDR'
307 format(10x,a,a)
308 format(2x,a,5x,a,5x,a,3x,a,10x,a)
306 format(8x,a,1x,a,1x,a,1x,a,1x,a)
309 format(a,1x,a,1x,I4,3x,a,1x,I4,3x,a,1x,I4,3x,a,1x,I4)
310 format(10x,a,a,2x,a,I0)

print*,'sourav3',alstr
allocate (fullcovstr(alstr,nae))

do i=1,alstr
  do i1=1,nae
    fullcovstr(i,i1)=strc(i,i1)+niao
  enddo
  print*,'fullcovstr',i,(fullcovstr(i,i1),i1=1,nae)
enddo

print*,'sourav4'
!if(symm.eq.1) then
!  !allocate(symsc(MaxStrOepo))
!  if(sig_sym_flg.eq.1)call symmetry_cal_sig(nlp,fullcovstr,alstr,symsc,symq,nssym)
!  if(sig_sym_flg.ne.1)call symmetry_cal_pi(nlp,fullcovstr,alstr,symsc,symq,nssym)
!
!  nnn=0
!  do ij=1,tnqs
!    do ijk=nssym,1,-1
!      do j2=1,alstr
!        if(symq(j2).eq.ijk)then
!          nnn=nnn+1
!          if(niao.eq.0)then
!            write(7,900)'cov structure',nnn,')',qulsym(j2),(fullcovstr(j2,j1),j1=1,nao-nlp+nlp*2)
!          endif
!          if(niao.ne.0)then
!            write(7,901)'cov structure',nnn,')',qulsym(j2),1,':',niao,(fullcovstr(j2,j1),j1=1,nao-nlp+nlp*2)
!          endif
!        endif
!      enddo
!      write(7,*)'******************************************************'
!    enddo
!  enddo
!else
!  do j2=1,alstr
!    if(niao.eq.0)then
!      write(7,902)'cov structure',j2,')',(fullcovstr(j2,j1),j1=1,nao-nlp+nlp*2)
!    endif
!    if(niao.ne.0)then
!      write(7,903)'cov structure',j2,')',1,':',niao,(fullcovstr(j2,j1),j1=1,nao-nlp+nlp*2)
!    endif
!  enddo
!  write(7,*)'******************************************************'
!endif
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

print*,'sourav5'
perm_nstr=wig2*c*c1
if(flg_cov.eq.1) write(9+u1,*)perm_nstr,' covalent structures' 
if(flg_ion.eq.1) write(9+u1,*)perm_nstr,' ionic structures' 
write(23,*)perm_nstr,' covalent structures' 
if(nfset.ne.5)write(10,913)'IAB','NNB','SBB','fqual'
write(23,913)'IAB','NNB','SBB','fqual'
913 format(x,a,x,a,x,a,x,a)
deallocate(strc)
print*,'sourav6'
if (.not. allocated(col9)) then
allocate(col9(perm_nstr))
col9 = 0
endif
if (.not. allocated(strno)) then
allocate(strno(alstr))
strno = 0
endif
if (.not. allocated(str12)) then
allocate(str12(alstr, nae))
str12 = 0
endif
call str_selection(fullcovstr,nlp,alstr,perm_nstr)

deallocate(col9)
deallocate(strno)
deallocate(str12)

print*,'exit cove_struc'
return
end subroutine cov_struc
end module mod_cov_struc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


