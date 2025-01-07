!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! here thee structures are converted into vector form to be used in the gramian matrix  !
! to veryfy their linear independency. In this process structures are first converted   !
! into their determinants and the vectorial forms are actually indecating the presence  !
! of the determinants in the all possible determinants. Determinants has two part Alpha !
! and Beta means spin up and spin down.                                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine vector_rep(nl,str1,totstr,vec)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! vector representation of the structures for independency calculations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
use commondat1
implicit none 

common/ats/atsymset,nsym,syn,at_sym
integer::j,j1,j2,i,i1,i2,i3,i4,i5,i6,i7,i8,k,k6,k7,k8,k11,k12,k13,k14,k15,l,nend,c1,c2,c3,c4,c5,c6,c,d,e,f,totstr,&
cof,ntdet,det,det1,nl,permtn1,permtn2,permtn,tndet,kk,kkk,x,y,z,detno,n1,n2,n3,n4,n5,n6,n7,n8,i9,i10,i11,&
c7,c8,c9,c10
integer::rum1(100),rum2(100),m(500),bonds(100),coef(500),vec(15000,1000),str1(15000,20)
integer::score(40),score_sym(40),col(20),detcount,col1(2,20)
real*8::factorial,norm_const,detn
Double Precision::D1(1000),aldet_new,aldet(1000)
integer::atsymset(20,20),nsym,syn(50),at_sym(50)


integer,dimension(:,:),allocatable::str2,beta,beta1,alfa,vector,all_det,sign_det,all_det_new

print*,'enter vector_rep'

! memory allocations for local variables
x=100000
y=10000
z=1000

allocate(sign_det(x,z))
allocate(str2(x,y))
allocate(beta(x,y))
allocate(beta1(x,y))
allocate(alfa(x,y))
allocate(vector(x,y))
allocate(all_det(x,y))
allocate(all_det_new(x,y))

str2(x,y)=0
beta(x,y)=0
beta1(x,y)=0
alfa(x,y)=0
vector(x,y)=0
all_det(x,y)=0
str_det_sec(x,z)=0

! initialisation of the structure determinant array strdet
do i=1,10000
  strdet(i)=0
enddo

!!!!!!! "rum1()" is the orbital numbers except lone paires !!!!!
!!!!!!! "rum2()" is the orbital numbers except lone paires arranged in increasing order !!!!!


do k6=1,100
  rum1(k6)=0
enddo

do k6=1,nae-nl*2
  loop1:do k7=nl*2+1,nae
    do k8=1,k6
      if(rum1(k8).eq.str1(1,k7)) then
        cycle loop1
      endif
    enddo
    if(str1(1,k7).gt.rum1(k6))then
      rum1(k6)=str1(1,k7)
    endif
  enddo loop1
enddo

do k6=1,100
  rum2(k6)=0
enddo

k7=0
do k6=nae-nl*2,1,-1
  k7=k7+1
  rum2(k7)=rum1(k6)
enddo
nend=k7+mod(k7,2)


!!!!!!! "all_det()" are the all possible determinants of the set. Ex. for 5
!electrons in 5 orbitals case : there must be 3 alfa and two beta in each
!determinant so the possible number of determinants are= !5/(!3*!2)=10 so all_det()
!is the 10 defferent sets of determinant. 

k11=0
k12=0

! nested loop starts here
do i1=1,k7
  k11=1
  m(1)=i1
  if (k11.eq.nend/2) then
    k12=k12+1
    do k15=1,nend/2
      all_det(k12,k15)=rum2(m(k15))
    enddo
    cycle
  endif

! loop 2
  do i2=i1,k7
    k11=2
    if(i2.eq.m(1)) then
      cycle
    endif
      m(2)=i2
    if (k11.eq.nend/2) then
      k12=k12+1
      do k15=1,nend/2
        all_det(k12,k15)=rum2(m(k15))
      enddo
      cycle
    endif

    ! loop 3
    do i3=i2,k7
      k11=3
      do k13=1,2
        if(i3.eq.m(k13)) then
          cycle
        endif
      enddo
      m(3)=i3
      if (k11.eq.nend/2) then
        k12=k12+1
        do k15=1,nend/2
          all_det(k12,k15)=rum2(m(k15))
        enddo
        cycle
      endif

      ! loop 4
      do i4=i3,k7
        k11=4
        do k14=1,3
          if(i4.eq.m(k14)) then
            cycle
          endif
        enddo
        m(4)=i4
        if (k11.eq.nend/2) then
          k12=k12+1
          do k15=1,nend/2
            all_det(k12,k15)=rum2(m(k15))
          enddo
          cycle
        endif

        ! loop 5
        do i5=i4,k7
          k11=5
          do k14=1,4
            if(i5.eq.m(k14)) then
              cycle
            endif
          enddo
          m(5)=i5
          if (k11.eq.nend/2) then
            k12=k12+1
            do k15=1,nend/2
              all_det(k12,k15)=rum2(m(k15))
            enddo
            cycle
          endif

          ! loop 6
          do i6=i5,k7
            k11=6
            do k14=1,5
              if(i6.eq.m(k14)) then
                cycle
              endif
            enddo
            m(6)=i6
            if (k11.eq.nend/2) then
              k12=k12+1
              do k15=1,nend/2
                all_det(k12,k15)=rum2(m(k15))
              enddo
              cycle
            endif

            ! loop 7
            do i7=i6,k7
              k11=7
              do k14=1,6
                if(i7.eq.m(k14)) then
                  cycle
                endif
              enddo
              m(7)=i7
              if (k11.eq.nend/2) then
                k12=k12+1
                do k15=1,nend/2
                  all_det(k12,k15)=rum2(m(k15))
                enddo
                cycle
              endif

              ! loop 8
              do i8=i7,k7
                k11=8
                do k14=1,7
                  if(i8.eq.m(k14)) then
                    cycle
                  endif
                enddo
                m(8)=i8
                if (k11.eq.nend/2) then
                  k12=k12+1
                  do k15=1,nend/2
                    all_det(k12,k15)=rum2(m(k15))
                  enddo
                  cycle
                endif

                ! loop 9
                do i9=i8,k7
                  k11=9
                  do k14=1,8
                    if(i9.eq.m(k14)) then
                      cycle
                    endif
                  enddo
                  m(9)=i9
                  if (k11.eq.nend/2) then
                    k12=k12+1
                    do k15=1,nend/2
                      all_det(k12,k15)=rum2(m(k15))
                    enddo
                    cycle
                  endif

                  ! cycle 10
                  do i10=i9,k7
                    k11=10
                    do k14=1,9
                      if(i10.eq.m(k14)) then
                        cycle
                      endif
                    enddo
                    m(10)=i10
                    if (k11.eq.nend/2) then
                      k12=k12+1
                      do k15=1,nend/2
                        all_det(k12,k15)=rum2(m(k15))
                      enddo
                      cycle
                    endif

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

ntdet=k12  ! number of determinants found

d=nae-nl*2
e=nend/2
f=d-e
c=factorial(d)/(factorial(e)*factorial(f))
do i1=1,c
  i4=nend/2
  loop2:do i3=1,nae-nl*2
    do i2=1,nend/2
      if(all_det(i1,i2).eq.rum2(i3)) then
        cycle loop2
      endif
    enddo
    i4=i4+1
    all_det(i1,i4)=rum2(i3)
  enddo loop2
  nbeta=(nae-(nl*2)-(mult-1))/2
  nalpha=nbeta+(mult-1)

  ! if lone pairs are present in the system will be added in the determinant below
  if (nl.ne.0)then
    i2=0
    do i3=1,nalpha
      i2=i2+1
      all_det_new(i1,i2)=all_det(i1,i3)
    enddo
  
    do i3=1,nl,2
      i2=i2+1
      all_det_new(i1,i2)=str1(1,i3)
    enddo
    do i3=nalpha+1,nalpha+nbeta
      i2=i2+1
      all_det_new(i1,i2)=all_det(i1,i3)
    enddo
    do i3=1,nl,2
      i2=i2+1
      all_det_new(i1,i2)=str1(1,i3)
    enddo
  endif
  
  if(nl.eq.0)then
    do i3=1,nae
      all_det_new(i1,i3)=all_det(i1,i3)
    enddo
  endif
  ! finally the determinants are stored in all_det_new
enddo

tot_ndet=i1 ! total number of determinants are revealed 

! now it starting calculate vectors below
do i1=1,15000
  do i2=1,1000
    vector(i1,i2)=0
  enddo
enddo


!!!!!!  'struc_albe () ' is the determinants of each structures presenting only alfas !!!!!!
!!!!!!  for two bonds and one unpaired electrons Ex. (1-2)(3-4)5, two bonds
!!!!!!  between orb 1,2 and orbs 3,4 and unpaired electron in in 5. then the
!!!!!!  determinants are (a1b2-a2b1)(a3b4-a4b3)a5 taking unpaired as alfa. so the
!!!!!!  determinants are a1b2a3b4a5-a1b2a4b3a5-a2b1a3b4a5+a2b1a4b3a5 to delet the a,b
!!!!!!  notation or to make it in only number sets we put all alfas first and beta
!!!!!!  second, therefore the number set representing the above determinants are =
!!!!!!  135/24 - 145/23 - 235/14 + 245/13, 'cof()' is the sign of the determinants.  !!!!!!!
i4=0
tndet=0
do i7=1,totstr
  j=0
  i2=0
  detcount=0
  do i4=1,100
    bonds(100)=0
  enddo
  do i1=nl*2+1,nae-nlast
    i2=i2+1
    bonds(i2)=str1(i7,i1)
  enddo

  do i1=1,500
    m(i1)=0
  enddo

  k11=0
  k12=0
  cof=1
  do i1=1,2
  cof=cof*(-1)**(i1-1)  ! cof is the +ve or -ve sign if the determinants
    c1=cof
    k11=1
    m(1)=i1
    if (k11.eq.(nae-nl*2-nlast)/2) then
      k12=k12+1
      do k15=1,(nae-nl*2-nlast)/2
        str2(k12,k15)=bonds(m(k15))
        coef(k12)=cof
      enddo
      cof=c1
      cycle
    endif

    do i2=3,4
      cof=cof*(-1)**(i2-1)
      c2=cof
      k11=2
      if(i2.eq.m(1))then
        cof=c2
        cycle
      endif
      m(2)=i2
      if (k11.eq.(nae-nl*2-nlast)/2) then
        k12=k12+1
        do k15=1,(nae-nl*2-nlast)/2
          str2(k12,k15)=bonds(m(k15))
          coef(k12)=cof
        enddo
        cof=c2
        cycle
      endif

      do i3=5,6
        cof=cof*(-1)**(i3-1)
        c3=cof
        k11=3
        do k13=1,2
          if(i3.eq.m(k13)) then
            cof=c3
            cycle
          endif
        enddo
        m(3)=i3
        if (k11.eq.(nae-nl*2-nlast)/2) then
          k12=k12+1
          do k15=1,(nae-nl*2-nlast)/2
            str2(k12,k15)=bonds(m(k15))
            coef(k12)=cof
          enddo
          cof=c3
          cycle
        endif

        do i4=7,8
          cof=cof*(-1)**(i4-1)
          c4=cof
          k11=4
          do k14=1,3
            if(i4.eq.m(k14)) then
              cof=c4
              cycle
            endif
          enddo
          m(4)=i4
          if (k11.eq.(nae-nl*2-nlast)/2) then
            k12=k12+1
            do k15=1,(nae-nl*2-nlast)/2
              str2(k12,k15)=bonds(m(k15))
              coef(k12)=cof
            enddo
            cof=c4
            cycle
          endif

          do i5=9,10
            cof=cof*(-1)**(i5-1)
            c5=cof
            k11=5
            do k14=1,4
              if(i5.eq.m(k14)) then
                cof=c5
                cycle
              endif
            enddo
            m(5)=i5
            if (k11.eq.(nae-nl*2-nlast)/2) then
              k12=k12+1
              do k15=1,(nae-nl*2-nlast)/2
                str2(k12,k15)=bonds(m(k15))
                coef(k12)=cof
              enddo
              cof=c5
              cycle
            endif

              do i6=11,12
                cof=cof*(-1)**(i6-1)
                c6=cof
                k11=6
                do k14=1,5
                  if(i6.eq.m(k14)) then
                    cof=c6
                    cycle
                  endif
                enddo
                m(6)=i6
                if (k11.eq.(nae-nl*2-nlast)/2) then
                  k12=k12+1
                  do k15=1,(nae-nl*2-nlast)/2
                    str2(k12,k15)=bonds(m(k15))
                    coef(k12)=cof
                  enddo
                  cof=c6
                  cycle
                endif

                do i8=13,14
                  cof=cof*(-1)**(i8-1)
                  c7=cof
                  k11=7
                  do k14=1,6
                    if(i8.eq.m(k14)) then
                      cof=c7
                      cycle
                    endif
                  enddo
                  m(7)=i8
                  if (k11.eq.(nae-nl*2-nlast)/2) then
                    k12=k12+1
                    do k15=1,(nae-nl*2-nlast)/2
                      str2(k12,k15)=bonds(m(k15))
                      coef(k12)=cof
                    enddo
                    cof=c7
                    cycle
                  endif

                  do i9=15,16
                    cof=cof*(-1)**(i9-1)
                    c8=cof
                    k11=8
                    do k14=1,7
                      if(i9.eq.m(k14)) then
                        cof=c8
                        cycle
                      endif
                    enddo
                    m(8)=i9
                    if (k11.eq.(nae-nl*2-nlast)/2) then
                      k12=k12+1
                      do k15=1,(nae-nl*2-nlast)/2
                        str2(k12,k15)=bonds(m(k15))
                        coef(k12)=cof
                      enddo
                      cof=c8
                      cycle
                    endif

                    do i10=17,18
                      cof=cof*(-1)**(i10-1)
                      c9=cof
                      k11=9
                      do k14=1,8
                        if(i10.eq.m(k14)) then
                          cof=c9
                          cycle
                        endif
                      enddo
                      m(9)=i10
                      if (k11.eq.(nae-nl*2-nlast)/2) then
                        k12=k12+1
                        do k15=1,(nae-nl*2-nlast)/2
                          str2(k12,k15)=bonds(m(k15))
                          coef(k12)=cof
                        enddo
                        cof=c9
                        cycle
                      endif

                      do i11=19,20
                        cof=cof*(-1)**(i11-1)
                        c10=cof
                        k11=10
                        do k14=1,9
                          if(i11.eq.m(k14)) then
                            cof=c10
                            cycle
                          endif
                        enddo
                        m(10)=i11
                        if (k11.eq.(nae-nl*2-nlast)/2) then
                          k12=k12+1
                          do k15=1,(nae-nl*2-nlast)/2
                            str2(k12,k15)=bonds(m(k15))
                            coef(k12)=cof
                          enddo
                          cof=c10
                          cycle
                        endif

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

ndet=k12  ! number of determinants of each structures revealed 

! alpha and beta finally achieving here
do k14=1,k12
  tndet=tndet+1
  i8=(nae-nl*2-nlast)/2
  do i6=nae-nlast+1,nae
    i8=i8+1
    str2(k14,i8)=str1(i7,i6)
  enddo
  i5=0
  do i1=1,i8
    do i2=nl*2+1,nae-nlast
      if(str2(k14,i1).eq.str1(i7,i2))then
        if(mod(i2,2).eq.1)i3=i2+1
        if(mod(i2,2).eq.0)i3=i2-1
      endif
    enddo  
    beta1(k14,i1)=str1(i7,i3)
  enddo

  k7=0
  if(nl.ne.0)then
    do k6=1,nl*2,2
      k7=k7+1
      alfa(k14,k7)=str1(i7,k6)
      beta(k14,k7)=str1(i7,k6)
    enddo
  endif

  do k6=1,(nae-nl*2-nlast)/2
    k7=k7+1
    alfa(k14,k7)=str2(k14,k6)
    beta(k14,k7)=beta1(k14,k6)
  enddo


! in the systen has radical by convention it will be added to the alfa
! nlast is radical

if(nlast.ne.0)then
  do k6=nae-nlast+1,nae
    k7=k7+1
    alfa(k14,k7)=str1(i7,k6)
  enddo
endif

do k6=1,100
  rum1(k6)=0
enddo

k7=0
i5=0
kkk=0
permtn1=0
770 k=100
kk=0

loop4: do k6=1,(nae-nl*2-nlast)/2+nl+nlast
  do k8=1,k7
    if(rum1(k8).eq.alfa(k14,k6)) then
      cycle loop4
    endif
  enddo
  kk=kk+1
  if(k.gt.alfa(k14,k6))then
    kkk=kk
    k=alfa(k14,k6)
  endif
enddo loop4
k7=k7+1
i5=i5+1
rum1(k7)=k
permtn1=permtn1+(kkk-1)
detmnt(tndet,i5)=k
strdet(tndet)=i7


if(k7.lt.((nae-nl*2-nlast)/2)+nl+nlast)goto 770


do k6=1,100
  rum1(k6)=0
enddo

k7=0
permtn2=0
kkk=0
780 k=1000
kk=0
loop5:do k6=1,(nae-nl*2-nlast)/2+nl
  do k8=1,k7
    if(rum1(k8).eq.beta(k14,k6)) then 
      cycle loop5
    endif
  enddo
  kk=kk+1
  if(k.gt.beta(k14,k6))then
    kkk=kk
    k=beta(k14,k6)
  endif
enddo loop5
k7=k7+1
i5=i5+1
rum1(k7)=k
permtn2=permtn2+(kkk-1)
detmnt(tndet,i5)=k
strdet(tndet)=i7

if(k7.lt.((nae-nl*2-nlast)/2)+nl+nlast)goto 780

permtn=permtn1+permtn2

det_sign(tndet)=(-1)**permtn
j=j+1
sign_det(i7,j)=det_sign(tndet)
detcount=detcount+1
str_det_sec(i7,detcount)=tndet

enddo
300 format (a,20I3)


do i1=1,k12
  i4=nend/2
  loop6:do i3=1,nae-nl*2
    do i2=1,nend/2
      if(str2(i1,i2).eq.rum2(i3)) then
        cycle loop6
      endif
    enddo
    i4=i4+1
    str2(i1,i4)=rum2(i3)
  enddo loop6
enddo

! vector form of each structures now formed
i3=0
do i1=1,ndet
  do i2=1,ntdet
    det=0
    det1=0
    do i3=1,nend/2
      do i4=1,nend/2
        if(str2(i1,i3).eq.all_det(i2,i4))then
          det=det+1
        endif
      enddo
    enddo
    if(det.eq.nend/2)then
      do i3=1+nend/2,nae-nl*2
        do i4=1+nend/2,nae-nl*2
          if(str2(i1,i3).eq.all_det(i2,i4))then
            det1=det1+1
          endif
        enddo
      enddo
      if(det1.eq.(nae-nl*2)-(nend/2))then
        i3=i3+1
        vector(i7,i2)=1*coef(i1)
      endif
    endif
  enddo
enddo

i3=0
do i2=1,ntdet
  if(vector(i7,i2).ne.0)then
    i3=i3+1
    vec(i7,i3)=i2
  endif
enddo

detno=i3

detn= real(detno,4)
norm_const=1/sqrt(detn)

enddo


!do l=1,totstr*ndet
!enddo
123 format(a,2x,50I5)
321 format(I5,a,500I3)
!!!end of independency test

339 print*,'exit vector_rep'
deallocate(sign_det)
deallocate(str2)
deallocate(beta)
deallocate(beta1)
deallocate(alfa)
deallocate(vector)
deallocate(all_det)

return
end subroutine vector_rep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
