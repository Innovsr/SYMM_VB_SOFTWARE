!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine symm_loops(k1,k2,loop_score,connectivity_row,nrow)

use commondat
use loops
implicit none

!common/loops/full_nn_group,fullgrp,natom,nelimt,sl_group,tot_orb,nn_group

integer::k1,k2,i,i1,i2,i4,i5,i6,j,loop_score1,nrow, count
!integer::nn_group(50,10),nelimt(50),full_nn_group(1000),sl_group(50,10),tot_orb(100)
integer::nloop(100),orbs1(200),orbs2(200),orbs3(40),n,nsl,nsl1,ii,iii,flg,loop_score
integer::ii0,ii1,ii2,ii3,ii4,ii5,ii6,ii7,ii8,ii9,ii10,ii11,ii12,j1,j2,&
j3,j4,j5,j6,j7,j8,j9,j10,j11,j12,j13,l,maxlp
integer::ls1,ls2,ls3,ls4,ls5,ls6,ls7,ls8,ls9,ls10,ls11,ls12,ls13,ls14,connectivity(20,20),&
connectivity_row(20,20),nn,nnn(10)

print*,'enter_symm_loops',k1,k2
nsl=0
nsl1=0
ii=0
iii=0
i4=0
i5=0
flg=0
maxlp=15
loop_score1=0
l=0

!print*,'natom',natom
do i=1,natom
  if(k1.eq.tot_orb(i)) exit
enddo

!print*,'iii',i,nelimt(i)
loop1:do j=1,nelimt(i)
  do i6=1,maxlp
    orbs2(i6)=0
  enddo
  ls1=0
  if(k2.eq.nn_group(i,j))then
    iii=iii+1
    ls1=ls1+1
    connectivity(iii,1)=k1
    connectivity(iii,2)=k2
    nloop(iii)=ls1+1
    cycle loop1
  else
    if(k1.ne.nn_group(i,j))then
      l=1
      ii=ii+1
      orbs3(ii)=nn_group(i,j)
    endif
  endif
enddo loop1

!print*,'iiiii',ii


ii0=ii
! loop 1
loop32:do j1=1,ii0
  ii=1
  do i6=1,maxlp
    orbs2(i6)=0
  enddo
  orbs2(1)=orbs3(j1)
  l=0
  ls2=ls1
  count = 0
  do i=1,natom
    if(orbs3(j1).eq.tot_orb(i)) then
      exit
    else
      count = count +1
      if (count.eq.natom) cycle loop32
    endif
  enddo
  do i1=1,nelimt(i)
    if(k2.eq.nn_group(i,i1))then
      iii=iii+1
      ls2=ls2+1
      connectivity(iii,1)=k1
      connectivity(iii,2)=orbs2(1)
      connectivity(iii,3)=k2
      nloop(iii)=ls2+1
      cycle loop32
    else
      if(k1.ne.nn_group(i,i1))then
        l=1
        ii=ii+1
        orbs1(ii)=nn_group(i,i1)
      endif
    endif
  enddo


  if(l.eq.1)ls2=ls2+l
  ii1=ii




  !!! loop 2
  loop31:do j2=2,ii1
    ii=ii1
    do i6=2,maxlp
      orbs2(i6)=0
    enddo
    do i6=1,1
      if(orbs2(i6).eq.orbs1(j2)) cycle loop31
    enddo
    orbs2(2)=orbs1(j2)
    l=0
    ls3=ls2
    count = 0
    do i=1,natom
      if(orbs1(j2).eq.tot_orb(i)) then
        exit
      else
      count = count +1
      if (count.eq.natom) cycle loop32
      !  cycle loop32
      endif
    enddo
    do i1=1,nelimt(i)
      if(k2.eq.nn_group(i,i1))then
        iii=iii+1
        ls3=ls3+1
        connectivity(iii,1)=k1
        do i6=1,2
          connectivity(iii,i6+1)=orbs2(i6)
        enddo
        connectivity(iii,4)=k2
        nloop(iii)=ls3+1
        cycle loop31
      else
        if(k1.ne.nn_group(i,i1))then
          ii=ii+1
          l=1
          orbs1(ii)=nn_group(i,i1)
        endif
      endif
    enddo
    if(l.eq.1)ls3=ls3+l
    ii2=ii




    !!! loop 3
    loop30:do j3=ii1+1,ii2
      ii=ii2
      do i6=3,maxlp
        orbs2(i6)=0
      enddo
      do i6=1,2
        if(orbs2(i6).eq.orbs1(j3)) cycle loop30
      enddo
      orbs2(3)=orbs1(j3)
      l=0
      ls4=ls3
      count = 0
      do i=1,natom
        if(orbs1(j3).eq.tot_orb(i)) then
          exit
        else
          count = count +1
          if (count.eq.natom) cycle loop31
         ! cycle loop31
        endif
      enddo
      do i1=1,nelimt(i)
        if(k2.eq.nn_group(i,i1))then
          iii=iii+1
          ls4=ls4+1
          connectivity(iii,1)=k1
          do i6=1,3
            connectivity(iii,i6+1)=orbs2(i6)
          enddo
          connectivity(iii,5)=k2
          nloop(iii)=ls4+1
          cycle loop30
        else
          if(k1.ne.nn_group(i,i1))then
            ii=ii+1
            l=1
            orbs1(ii)=nn_group(i,i1)
          endif
        endif
      enddo
    
      ii3=ii
      if(l.eq.1)ls4=ls4+l
    



      !!! loop 4
      loop29:do j4=ii2+1,ii3
        ii=ii3
        do i6=4,maxlp
          orbs2(i6)=0
        enddo
        do i6=1,3
          if(orbs2(i6).eq.orbs1(j4)) cycle loop29
        enddo
        orbs2(4)=orbs1(j4)
        l=0
        ls5=ls4
        count = 0
        do i=1,natom
          if(orbs1(j4).eq.tot_orb(i)) then
            exit
          else
            count = count +1
            if (count.eq.natom) cycle loop30
           !cycle loop30
          endif
        enddo
        do i1=1,nelimt(i)
          if(k2.eq.nn_group(i,i1))then
            iii=iii+1
            ls5=ls5+1
            connectivity(iii,1)=k1
            do i6=1,4
              connectivity(iii,i6+1)=orbs2(i6)
            enddo
            connectivity(iii,6)=k2
            nloop(iii)=ls5+1
            cycle loop29
          else
            if(k1.ne.nn_group(i,i1))then
              ii=ii+1
              l=1
              orbs1(ii)=nn_group(i,i1)
            endif
          endif
        enddo
        if(l.eq.1)ls5=ls5+l
        ii4=ii



      
        !!! loop 5
        loop28:do j5=ii3+1,ii4
          ii=ii4
          do i6=5,maxlp
            orbs2(i6)=0
          enddo
          do i6=1,4
            if(orbs2(i6).eq.orbs1(j5)) cycle loop28
          enddo
          orbs2(5)=orbs1(j5)
          l=0
          ls6=ls5
          count = 0
          do i=1,natom
            if(orbs1(j5).eq.tot_orb(i)) then
              exit
            else
              count = count +1
              if (count.eq.natom) cycle loop29
             !cycle loop29
            endif
          enddo
          do i1=1,nelimt(i)
            if(k2.eq.nn_group(i,i1))then
              iii=iii+1
              ls6=ls6+1
              connectivity(iii,1)=k1
              do i6=1,5
                connectivity(iii,i6+1)=orbs2(i6)
              enddo
              connectivity(iii,7)=k2
              nloop(iii)=ls6+1
              cycle loop28
            else
              if(k1.ne.nn_group(i,i1))then
                ii=ii+1
                l=1
                orbs1(ii)=nn_group(i,i1)
              endif
            endif
          enddo
          if(l.eq.1)ls6=ls6+l
          ii5=ii



        
          !!! loop 6
          loop27:do j6=ii4+1,ii5
            ii=ii5
            do i6=6,maxlp
              orbs2(i6)=0
            enddo
            do i6=1,5
              if(orbs2(i6).eq.orbs1(j6)) cycle loop27
            enddo
            orbs2(6)=orbs1(j6)
            l=0
            ls7=ls6
            count = 0
            do i=1,natom
              if(orbs1(j6).eq.tot_orb(i)) then
                exit
              else
                count = count +1
                if (count.eq.natom) cycle loop28
               !cycle loop28
              endif
            enddo
            do i1=1,nelimt(i)
              if(k2.eq.nn_group(i,i1))then
                iii=iii+1
                ls7=ls7+1
                connectivity(iii,1)=k1
                do i6=1,6
                  connectivity(iii,i6+1)=orbs2(i6)
                enddo
                connectivity(iii,8)=k2
                nloop(iii)=ls7+1
                cycle loop27
              else
                if(k1.ne.nn_group(i,i1))then
                  ii=ii+1
                  l=1
                  orbs1(ii)=nn_group(i,i1)
                endif
              endif
            enddo
            if(l.eq.1)ls7=ls7+l
            ii6=ii



          
            !!! loop 7
            loop26:do j7=ii5+1,ii6
              ii=ii6
              do i6=7,maxlp
                orbs2(i6)=0
              enddo
              do i6=1,6
                if(orbs2(i6).eq.orbs1(j7)) cycle loop26
              enddo
              orbs2(7)=orbs1(j7)
              l=0
              ls8=ls7
              count = 0
              do i=1,natom
                if(orbs1(j7).eq.tot_orb(i)) then
                  exit
                else
                  count = count +1
                  if (count.eq.natom) cycle loop27
               !  cycle loop27
                endif
              enddo
              do i1=1,nelimt(i)
                if(k2.eq.nn_group(i,i1))then
                  iii=iii+1
                  ls8=ls8+1
                  connectivity(iii,1)=k1
                  do i6=1,7
                    connectivity(iii,i6+1)=orbs2(i6)
                  enddo
                  connectivity(iii,9)=k2
                  nloop(iii)=ls8+1
                  cycle loop26
                else
                  if(k1.ne.nn_group(i,i1))then
                    ii=ii+1
                    l=1
                    orbs1(ii)=nn_group(i,i1)
                  endif
                endif
              enddo
              if(l.eq.1)ls8=ls8+l
              ii7=ii
            



              !!! loop 8
              loop25:do j8=ii6+1,ii7
                ii=ii7
                do i6=8,maxlp
                  orbs2(i6)=0
                enddo
                do i6=1,7
                  if(orbs2(i6).eq.orbs1(j8)) cycle loop25
                enddo
                orbs2(8)=orbs1(j8)
                l=0
                ls9=ls8
                count = 0
                do i=1,natom
                  if(orbs1(j8).eq.tot_orb(i)) then
                    exit
                  else
                    count = count +1
                    if (count.eq.natom) cycle loop26
                   !cycle loop26
                  endif
                enddo
                do i1=1,nelimt(i)
                  if(k2.eq.nn_group(i,i1))then
                    iii=iii+1
                    ls9=ls9+1
                    connectivity(iii,1)=k1
                    do i6=1,8
                      connectivity(iii,i6+1)=orbs2(i6)
                    enddo
                    connectivity(iii,10)=k2
                    nloop(iii)=ls9+1
                    cycle loop25
                  else
                    if(k1.ne.nn_group(i,i1))then
                      ii=ii+1
                      l=1
                      orbs1(ii)=nn_group(i,i1)
                    endif
                  endif
                enddo
                if(l.eq.1)ls9=ls9+l
                ii8=ii
              



                !!! loop 9
                loop24:do j9=ii7+1,ii8
                  ii=ii8
                  do i6=9,maxlp
                    orbs2(i6)=0
                  enddo
                  do i6=1,8
                    if(orbs2(i6).eq.orbs1(j9)) cycle loop24
                  enddo
                  orbs2(9)=orbs1(j9)
                  l=0
                  ls10=ls9
                  count = 0
                  do i=1,natom
                    if(orbs1(j9).eq.tot_orb(i)) then
                      exit
                    else
                      count = count +1
                      if (count.eq.natom) cycle loop25
                     !cycle loop25
                    endif
                  enddo
                  do i1=1,nelimt(i)
                    if(k2.eq.nn_group(i,i1))then
                      iii=iii+1
                      ls10=ls10+1
                      connectivity(iii,1)=k1
                      do i6=1,9
                        connectivity(iii,i6+1)=orbs2(i6)
                      enddo
                      connectivity(iii,11)=k2
                      nloop(iii)=ls10+1
                      cycle loop24
                    else
                      if(k1.ne.nn_group(i,i1))then
                        ii=ii+1
                        l=1
                        orbs1(ii)=nn_group(i,i1)
                      endif
                    endif
                  enddo
                  if(l.eq.1)ls10=ls10+l
                  ii9=ii


                
                  !!! loop 10
                  loop23:do j10=ii8+1,ii9
                    ii=ii9
                    do i6=10,maxlp
                      orbs2(i6)=0
                    enddo
                    do i6=1,9
                      if(orbs2(i6).eq.orbs1(j10)) cycle loop23
                    enddo
                    orbs2(10)=orbs1(j10)
                    l=0
                    ls11=ls10
                    count = 0
                    do i=1,natom
                      if(orbs1(j10).eq.tot_orb(i)) then
                        exit
                      else
                        count = count +1
                        if (count.eq.natom) cycle loop24
                       !cycle loop24
                      endif
                    enddo
                    do i1=1,nelimt(i)
                      if(k2.eq.nn_group(i,i1))then
                        iii=iii+1
                        ls11=ls11+1
                        connectivity(iii,1)=k1
                        do i6=1,10
                          connectivity(iii,i6+1)=orbs2(i6)
                        enddo
                        connectivity(iii,12)=k2
                        nloop(iii)=ls11+1
                        cycle loop23
                      else
                        if(k1.ne.nn_group(i,i1))then
                          ii=ii+1
                          l=1
                          orbs1(ii)=nn_group(i,i1)
                        endif
                      endif
                    enddo
                    if(l.eq.1)ls11=ls11+l
                    ii10=ii



                    !!! loop 11
                    loop22:do j11=ii9+1,ii10
                      ii=ii10
                      do i6=11,maxlp
                        orbs2(i6)=0
                      enddo
                      do i6=1,10
                        if(orbs2(i6).eq.orbs1(j11)) cycle loop22
                      enddo
                      orbs2(11)=orbs1(j11)
                      l=0
                      ls12=ls11
                      count = 0
                      do i=1,natom
                        if(orbs1(j11).eq.tot_orb(i)) then
                          exit
                        else
                          count = count +1
                          if (count.eq.natom) cycle loop22
                         !cycle loop23
                        endif
                      enddo
                      do i1=1,nelimt(i)
                        if(k2.eq.nn_group(i,i1))then
                          iii=iii+1
                          ls12=ls12+1
                          connectivity(iii,1)=k1
                          do i6=1,11
                            connectivity(iii,i6+1)=orbs2(i6)
                          enddo
                          connectivity(iii,13)=k2
                          nloop(iii)=ls12+1
                          cycle loop22
                        else
                          if(k1.ne.nn_group(i,i1))then
                            ii=ii+1
                            l=1
                            orbs1(ii)=nn_group(i,i1)
                          endif
                        endif
                      enddo
                      if(l.eq.1)ls12=ls12+l
                      ii11=ii
                    



                      !!! loop 12
                      loop21:do j12=ii10+1,ii11
                        ii=ii11
                        do i6=12,maxlp
                          orbs2(i6)=0
                        enddo
                        do i6=1,11
                          if(orbs2(i6).eq.orbs1(j12)) cycle loop21
                        enddo
                        orbs2(12)=orbs1(j12)
                        l=0
                        ls13=ls12
                        count = 0
                        do i=1,natom
                          if(orbs1(j12).eq.tot_orb(i)) then 
                            exit 
                          else 
                            count = count +1
                            if (count.eq.natom) cycle loop22
                           !cycle loop22
                          endif
                        enddo
                        do i1=1,nelimt(i)
                          if(k2.eq.nn_group(i,i1))then
                            iii=iii+1
                            ls13=ls13+1
                            connectivity(iii,1)=k1
                            do i6=1,12
                              connectivity(iii,i6+1)=orbs2(i6)
                            enddo
                            connectivity(iii,14)=k2
                            nloop(iii)=ls13+1
                            cycle loop21
                          else
                            if(k1.ne.nn_group(i,i1))then
                              ii=ii+1
                              l=1
                              orbs1(ii)=nn_group(i,i1)
                            endif
                          endif
                        enddo
                        if(l.eq.1)ls13=ls13+l
                        ii12=ii
                      


                        !print*,'loop13'
                        loop20:do j13=ii11+1,ii12
                          ii=ii12
                          do i6=13,maxlp
                            orbs2(i6)=0
                          enddo
                          do i6=1,12
                            if(orbs2(i6).eq.orbs1(j13)) cycle loop20
                          enddo
                          orbs2(13)=orbs1(j13)
                          l=0
                          ls14=ls13
                          count = 0
                          do i=1,natom
                            if(orbs1(j13).eq.tot_orb(i)) then
                              exit
                            else
                              count = count +1
                              if (count.eq.natom) cycle loop21
                             !cycle loop21
                            endif
                          enddo
                          do i1=1,nelimt(i)
                            if(k2.eq.nn_group(i,i1))then
                              iii=iii+1
                              ls14=ls14+1
                              connectivity(iii,1)=k1
                              do i6=1,13
                                connectivity(iii,i6+1)=orbs2(i6)
                              enddo
                              connectivity(iii,15)=k2
                              nloop(iii)=ls14+1
                              cycle loop20
                            else
                              if(k1.ne.nn_group(i,i1))then
                                ii=ii+1
                                orbs1(ii)=nn_group(i,i1)
                              endif
                            endif
                          enddo
                        enddo loop20
                      
                      enddo loop21
                    
                    enddo loop22
                  
                  enddo loop23
                
                enddo loop24
              
              enddo loop25
            
            enddo loop26
          
          enddo loop27
        
        enddo loop28
      
      enddo loop29
    
    enddo loop30

  enddo loop31

enddo loop32


do i1=1,iii
  print*,'nloopi_new',nloop(i1)
enddo

n=1000
do i1=1,iii
  if(n.gt.nloop(i1))n=nloop(i1)
enddo

loop_score=n
nn=0
l=0
do i1=1,iii
  if(loop_score.eq.nloop(i1))then
    l=l+1
    nnn(l)=i1
    nn=nn+1
  endif
enddo

print*,'nnn**',(nnn(i),i=1,l)

do i2=1,20
  do i1=1,20
    connectivity_row(i1,i2)=0
  enddo
enddo
!print*,'',nnn,nloop(nnn)
!if (nn.eq.1)then
do i2=1,l
  do i1=1,nloop(nnn(i2))+1
    connectivity_row(i2,i1)=connectivity(nnn(i2),i1)
!print*,'connectivity_row',connectivity_row(i1)
  enddo
enddo

nrow=l
!endif

print*,'nloop',(nloop(i1),i1=1,iii)
do i1=1,iii
print*,'connectivity',(connectivity(i1,i2),i2=1,nloop(i1)+1)
enddo
do i2=1,l
print*,'connectivity_row',(connectivity_row(i2,i1),i1=1,nloop(nnn(i2))+1)
enddo

print*,'shortestgest loop',n
!if(k1.eq.24.and.k2.eq.29)stop

return
end subroutine symm_loops
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
