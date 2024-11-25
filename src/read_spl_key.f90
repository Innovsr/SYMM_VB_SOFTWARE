!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_spl_key()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use commondat
use commondat1
implicit none

integer::i,MDP,io,a,b,i1,j,j1,i2,i3,d,i5,i6,i4,i7,st_num1(100),st_num2(100),num(500),m5,m6,k,kk,p,q,&
l,ll,stn,k1,k2,k3,k4,qult1(100),nbond,nel,n,nset(10),pent_set(5),info_val(4),argnum,num_2,ovlp,runn
character(len=4)::abc,at_orb(100),im_bnd(100),ato,str,at_bas_sym(100),strc,rumr,info_set(4)
character(len=200)::charst(1000)
character(len=2)::prio
character(len=200)::line34,line35,line5,lowercase
character(len=5)::pent,stn1,stn2,stn3
character(len=200)::sttr,Seq(20),aa
real*8::num_1(500)
logical :: fileexists
character(len=35)::inputfilename
integer::atsymset(20,20),nsym,syn(50),at_sym(50),MDP1

common/ats/atsymset,nsym,syn,at_sym
print*,'enter input_1'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! default values !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

symtype='loose'
set_order=0
iab_length=0.0
noq0=1000
noq1=1000
noq2=1000
noq3=1000
symm=0
ovval=1.000
ovopt=0
nlpset=0
nfset=0
runn=0
mset=2
qflg=0
ovlp_int=0
if(input_flg.eq.1)flg1=1
if(input_flg.eq.1)nsym=1
if(input_flg.eq.1)nstrt=0
if(input_flg.eq.1)flgst=1
if(input_flg.eq.1)niach=0
if(input_flg.eq.1)niabd=0
if(input_flg.eq.1)nialp=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MDP=0
i6=0
i7=0
n=0


! reading .xmi file starts here
if(input_flg.eq.0)then

  call getarg(1,inputfilename)
  INQUIRE(FILE=TRIM(inputfilename),EXIST=fileexists)
  IF (fileexists) THEN
    open(unit=21,file=TRIM(inputfilename),status='old')
  ELSE
    write(*,*)'SORRY This input file does not exist or you may not provide the &
    filename at all'
    stop
  ENDIF
  
  ! counting number of lines in the .xmi file
  do
    read(21,'(a)',iostat=io)
    if(io.ne.0)exit
    MDP=MDP+1
  enddo
  rewind(21)
  
  ! reading line from .xni and storing in charst for analyse
  do i=1,MDP
    read(21,'(a)')charst(i)
    if(index(lowercase(charst(i)),'$end').ne.0)then
      MDP1=i
    endif
  enddo
  rewind(21)
endif




do i=MDP1+1,MDP

  kk=0
  q=0
  do j=1,100
    st_num1(j)=0
    st_num2(j)=0
  enddo
  do j=1,500
    num(j)=0
  enddo
  
  
  if(index(charst(i),'#').eq.1) then
    cycle
  endif
    if(index(charst(i),'#').gt.1)then
      ll=index(charst(i),'#')
      do k=1,ll
        line5=charst(i)(k:k)
        if(line5.eq.'')then
          if(trim(charst(i)(1:k-1)).ne.'')then
            charst(i)=charst(i)(1:ll-1)
            exit
          endif
        endif
      enddo
    endif
    charst(i)=lowercase(charst(i))
  
  call Keycheck(charst(i))
  
  if (index(charst(i),'sym').ne.0)then
    symm=1
    if(index(charst(i),'loose').ne.0)then
      symtype='loose'
    endif
    if(index(charst(i),'tight').ne.0)then
      symtype='tight'
    endif
    if(index(charst(i),'check').ne.0)then
      symtype='check'
    endif
    if(index(charst(i),'qual').ne.0)then
      set_order=0
    endif
    if(index(charst(i),'stob').ne.0)then
      set_order=1
    endif
    if(index(charst(i),'btos').ne.0)then
      set_order=2
    endif
  
  rewind(21)
  endif
  
  
  if(index(charst(i),'mout').ne.0)then
    do j=1,i-1
      read(21,*)
    enddo
    Read(21,*)aa,mset
    if(mset.eq.0)mset=1
    mset=mset+1
    rewind(21)
  endif
  
  
  
  if(index(charst(i),'nset').ne.0)then
    do j=1,i-1
      read(21,*)
    enddo
    read(21,'(a)')line34
    line35=line34
    ll=len(trim(line34))
    do k=1,ll+1
      line5=line34(k:k)
      if(line5.eq.'')then
        q=q+1
        if(k.eq.1)p=1
        if(p.ne.1) then
          if(q.eq.1)then
            kk=1
            st_num1(1)=1
          endif
          kk=kk+1
          st_num1(kk)=k+1
          st_num2(kk)=k
        else
          kk=kk+1
          st_num1(kk)=k+1
          st_num2(kk)=k
        endif
        cycle
      endif
    enddo
    stn=0
    do k=1,500
      num(k)=0
    enddo
    m5=0
    do k=1,kk
      a=0
      b=0
      if(st_num2(k+1)-st_num2(k).ne.1.and.st_num2(k).le.ll)then
        a=st_num1(k)
        b=st_num2(k+1)
        if(line34(a:b).ne.'nset')then
          if(line34(a:b).ne.'')then
            stn=stn+1
            read(line34(a:b),'(I10)')num(stn)
            m5=m5+1
            nset(m5)=num(stn)
            if(m5.eq.1)nfset=nset(1)
            if(nset(1).ne.0)then
              if(m5.eq.2)noq0=nset(2)
              if(m5.eq.3)noq1=nset(3)
              if(m5.eq.4)noq2=nset(4)
              if(m5.eq.5)noq3=nset(5)
            endif
          endif
        endif
      endif
    enddo
    rewind(21)
  endif 
  
  
  
  if (index(charst(i),'iab').ne.0.or.index(charst(i),'sbb').ne.0&
  .or.index(charst(i),'nnb').ne.0.or.index(charst(i),'udr').ne.0&
  .or.index(charst(i),'udb').ne.0)then
    runn=runn+1
    if(runn.eq.1)then
      itb=0
      syb=0
      nnb=0
      radical=0
      mnbond=0
    endif
  
    do j=1,i-1
      read(21,*)
    enddo
    read(21,'(a)')line34
    line35=line34
    ll=len(trim(line34))
    q=0
    kk=0
    p=0
    do k=1,ll+1
      line5=line34(k:k)
      if(line5.eq.'')then
        q=q+1
        if(k.eq.1) p=1
        if(p.ne.1) then
          if(q.eq.1) then
            kk=1
            st_num1(1)=1
          endif
          kk=kk+1
          st_num1(kk)=k+1
          st_num2(kk)=k
        else
          kk=kk+1
          st_num1(kk)=k+1
          st_num2(kk)=k
        endif
        cycle
      endif
    enddo
    stn=0
    do k=1,500
      num(k)=0
    enddo
    m5=0
    do k=1,kk
      a=0
      b=0
      if(st_num2(k+1)-st_num2(k).ne.1.and.st_num2(k).le.ll)then
        a=st_num1(k)
        b=st_num2(k+1)
        if(line34(a:b).eq.'iab')then
          pent=line34(a:b)
          stn=stn+1
          cycle
        endif
        if(line34(a:b).eq.'sbb')then
          pent=line34(a:b)
          stn=stn+1
          cycle
        endif
        if(line34(a:b).eq.'nnb')then
          pent=line34(a:b)
          stn=stn+1
          cycle
        endif
        if(line34(a:b).eq.'udr')then
          pent=line34(a:b)
          stn=stn+1
          cycle
        endif
        if(line34(a:b).eq.'udb')then
          pent=line34(a:b)
          stn=stn+1
          cycle
        endif
        read(line34(a:b),'(I10)')num(stn)
      endif
    enddo
  
    if(pent.eq.'iab')itb=num(stn)
    if(pent.eq.'sbb')syb=num(stn)
    if(pent.eq.'nnb')nnb=num(stn)
    if(pent.eq.'udr')radical=num(stn)
    if(pent.eq.'udb')mnbond=num(stn)
  endif
  rewind(21)
  
  if (index(charst(i),'iab_len').ne.0)then
    do j=1,i-1
      read(21,*)
    enddo
    read(21,*)line34,iab_length
    rewind(21)
  endif
enddo



vacorb=nao-nae
nlp=nae-nao
if(vacorb.gt.1)nlp=0
if(qult(1).eq.0)then
  noqult=(1+(nae-nlast-nlp*2)/2)*(1+(nae-nlast-nlp*2)/2)
  do j=1,noqult
    qult(j)=j
  enddo
endif

do i=1,5
  pent_set(i)=0
enddo

if(prad.eq.0)radical=0
if(imbd.eq.0)mnbond=0
pent_set(1)=itb
pent_set(2)=syb
pent_set(3)=nnb
pent_set(4)=radical
pent_set(5)=mnbond

l=1
do j=1,5
  if(l.lt.pent_set(j))l=pent_set(j)
enddo

k=0
do j=1,l
  do i=1,5
    if(pent_set(i).eq.j)then
      k=k+1
      exit
    endif
  enddo
  do i=1,5
    if(pent_set(i).eq.j)then
      pent_set(i)=k
    endif
  enddo
enddo 

itb=pent_set(1)
syb=pent_set(2)
nnb=pent_set(3)
radical=pent_set(4)
mnbond =pent_set(5)
       
itb=itb+1
syb=syb+1
nnb=nnb+1
radical=radical+1
mnbond=mnbond+1

if(itb.eq.1.and.syb.eq.1.and.nnb.eq.1.and.radical.eq.1.and.mnbond.eq.1.and.flg1.ne.1)then
  qflg=1
endif
nlast=(mult-1)

if(ovlp.eq.0)ovopt=0
if(ovlp.eq.2)ovopt=1
if(ovlp.eq.1.and.nfset.eq.0)ovopt=0
if(ovlp.eq.1.and.nfset.eq.1)ovopt=1
if(ovlp.eq.1.and.nfset.eq.2)ovopt=1
if(ovlp.eq.1.and.nfset.eq.3)ovopt=1
!!!!!!!!! "vpt" is the option for user specifying overlap value "ovval". It will work
!!!!!!!!!!!!!only for vpt=1. To lock it please put any other value
vpt=1

write(7,*)'******************************************************************************************'
write(7,900)'*','active orbs =',nao,'active electrons =',nae,'multiplicity =',mult,'inactive orbs =',niao
write(7,*)'******************************************************************************************'
write(7,*)'******************************************************************************************'
900 format(a,1x,a,1x,I4,3x,a,1x,I4,3x,a,1x,I4,3x,a,1x,I4)

return
end subroutine read_spl_key
