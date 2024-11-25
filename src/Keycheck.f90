!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Its checking any misspelled keywords provided by the user !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Keycheck(line)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none
integer::i,j,k,l,nkey
character(len=200)::line,line1
character(len=10)::Keywd(20),lower(26)
data Keywd/'nset','mout','loose','tight','check','sym','chinst','ovlp','btos',&
'stob','qual','iab','nnb','sbb','udr','udb','prad','lpst','imbd','imps'/
data lower/'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q',&
'r','s','t','u','v','w','x','y','z'/

line1=''
nkey=20
l=len(trim(line))
k=0
loop1:do j=1,l+1
  do i=1,26
    if(line(j:j).eq.lower(i))then
      k=k+1
      line1(k:k)=lower(i)
      cycle loop1
    endif
  enddo
  if(len(trim(line1)).gt.0)then
    do i=1,nkey
      if(trim(line1).eq.Keywd(i))then
        line1=''
        k=0
        cycle loop1
      endif
    enddo
    write(*,111)'UNKNOWN KEYWORD:',trim(line1)
    write(*,112)'Please choose the keywords from ',(keywd(i),i=1,20)
    stop
  else
  cycle loop1
  endif
enddo loop1

111 format(a,2x,a)
112 format(a,2x,20a)


return
end subroutine Keycheck
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
