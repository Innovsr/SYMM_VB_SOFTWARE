!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module mod_rumer_structures
use commondat_mod
use rum_rad_mod
implicit none

contains
subroutine rumer_structures(nl,str,nstr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer::kk,k2,k3,k4,k5,k6,k7,k8,k9,nl,nstr,rum,rum1,rum2,rum3,rum4
integer, allocatable::rumer1(:,:),rumer2(:,:),rrad(:)
integer, pointer :: str(:,:)

!allocate(rumer1(MaxStrOepo, nae))
!allocate(rumer2(MaxStrOepo, nae))
allocate(rumer1(nstr, nae))
allocate(rumer2(nstr, nae))
allocate(rrad(nae))
allocate(rumer(nstr))
allocate(rumer_rad(nstr))

print*,'enter rumer_structures',nstr

rumer=0
rumer_rad=0
rumer1=0

do k3=1,nstr
  do k6=1,nae-nl*2-nlast
    loop1:do k7=nl*2+1,nae-nlast
      do k8=1,k6
        if(rumer1(k3,k8).eq.str(k3,k7)) cycle loop1
      enddo
      if(str(k3,k7).gt.rumer1(k3,k6))then
        rumer1(k3,k6)=str(k3,k7)
        !print*,'**',rumer1(k3,k6)
      endif
    enddo loop1
  enddo
enddo

rumer2=0

!231 format(20I3)
do k3=1,nstr
  k7=0
  do k6=nae-nl*2-nlast,1,-1
    k7=k7+1
    rumer2(k3,k7)=rumer1(k3,k6)
  enddo
enddo
!do k3=1,nstr
!print*,'**',(rumer2(k3,k6),k6=1,nae-nl*2-nlast)
!enddo

rumer=0
do k9=1,nstr
  rum=0
  do k6=nl*2+1,nae-nlast,2
    do k2=1,k7
      k8=k6
      if(str(k9,k8).eq.rumer2(k9,k2))rum1=k2
    enddo
    do k2=1,k7
      k8=k6+1
      if(str(k9,k8).eq.rumer2(k9,k2))rum2=k2
    enddo
    do k5=k6+2,nae-nlast,2
      do k2=1,k7
        k8=k5
        if(str(k9,k8).eq.rumer2(k9,k2))rum3=k2
      enddo
      do k2=1,k7
        k8=k5+1
        if(str(k9,k8).eq.rumer2(k9,k2))rum4=k2
      enddo
      if((rum1-rum3)*(rum1-rum4)*(rum2-rum3)*(rum2-rum4).gt.0)then
        rum=rum+1
      endif
    enddo
  enddo
  kk=0
  do k6=1,((nae-nl*2-nlast)/2)-1
    kk=kk+k6
  enddo

  if(rum.eq.kk)then
    rumer(k9) = 1
  endif
enddo

rumer_rad = 1
if(nlast.ne.0)then
  rumer2=0
  rumer1=0
  do k3=1,nstr
    do k6=1,nae-nl*2
      loop2:do k7=nl*2+1,nae
        do k8=1,k6
          if(rumer1(k3,k8).eq.str(k3,k7)) cycle loop2
        enddo
        if(str(k3,k7).gt.rumer1(k3,k6))then
          rumer1(k3,k6)=str(k3,k7)
        endif
      enddo loop2
    enddo
  enddo
  rumer2=0

!  do k3=1,MaxStrOepo
  do k3=1,nstr
    k7=0
    do k6=nae-nl*2,1,-1
      k7=k7+1
      rumer2(k3,k7)=rumer1(k3,k6)
    enddo
  enddo
  
  do k9=1,nstr
    rum=0
    do k6=nae-nlast+1,nae
      do k2=1,k7
        if(str(k9,k6).eq.rumer2(k9,k2))rrad(k6)=k2
      enddo
    enddo
  
    do k6=nl*2+1,nae-nlast,2
      do k2=1,k7
        k8=k6
        if(str(k9,k8).eq.rumer2(k9,k2))rum1=k2
      enddo
      do k2=1,k7
        k8=k6+1
        if(str(k9,k8).eq.rumer2(k9,k2))rum2=k2
      enddo
      do k4=nae-nlast+1,nae
        if((rrad(k4)-rum1)*(rrad(k4)-rum2).gt.0)then
          rum=rum+1
        endif
      enddo
    enddo
    
    if(rum.lt.((nae-nlast-nl*2)/2)*nlast)then
      rumer_rad(k9)=0
    endif
  enddo


endif

print*,'exit rumer_structures'
!deallocate(rumer)
!deallocate(rumer_rad)
deallocate(rumer1)
deallocate(rumer2)
deallocate(rrad)
return
end subroutine rumer_structures
end module mod_rumer_structures
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
