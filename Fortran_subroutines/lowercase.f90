!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This function changes keyword into lower case before cheking them so that 
!! user can put case-free letters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function lowercase(line)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none
integer::i,j,k,l
character(len=200)::line,lowercase
character(len=2)::line1,lower(26),upper(26)
data upper/'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q'&
,'R','S','T','U','V','W','X','Y','Z'/
data lower/'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q',&
'r','s','t','u','v','w','x','y','z'/

l=len(line)

loop1:do j=1,l
  line1=line(j:j)
  do i=1,26
    if(line1.eq.upper(i))then
      line(j:j)=lower(i)
      cycle loop1
    endif
  enddo
enddo loop1

lowercase=line
return
end function lowercase
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
