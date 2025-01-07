!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rumer_structures(nl,str,nstr,rumer,rumer_rad)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

integer::kk,k2,k3,k4,k5,k6,k7,k8,k9,nl,nstr,rum,rum1,rum2,rum3,rum4,rrad(100),rrad1,x,y
integer::str(15000,20),rumer(15000),rumer_rad(15000),rumer1(15000,100),rumer2(15000,100)

print*,'enter rumer_structures'

do k6=1,15000
rumer(k6)=0
rumer_rad(k6)=0
enddo
do k3=1,15000
do k6=1,100
rumer1(k3,k6)=0
enddo
enddo
do k3=1,nstr
do k6=1,nae-nl*2-nlast
do k7=nl*2+1,nae-nlast
do k8=1,k6
if(rumer1(k3,k8).eq.str(k3,k7))goto 389
enddo
if(str(k3,k7).gt.rumer1(k3,k6))then
rumer1(k3,k6)=str(k3,k7)
!print*,'**',rumer1(k3,k6)
endif
389 enddo
enddo
enddo
do k3=1,15000
do k6=1,100
rumer2(k3,k6)=0
enddo
enddo

231 format(20I3)
do k3=1,nstr
k7=0
do k6=nae-nl*2-nlast,1,-1
k7=k7+1
rumer2(k3,k7)=rumer1(k3,k6)

!print*,rumer2(k3,k7)
enddo
enddo
do k3=1,nstr
print*,'**',(rumer2(k3,k6),k6=1,nae-nl*2-nlast)
enddo
do k6=1,15000
rumer(k6)=0
enddo
!print*,'k7',k7
do k9=1,nstr
!write(*,231),(str(k9,k6),k6=1,nae)
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

!print*,'rum1,rum2,rum',rum1,rum2,rum
if(rum.eq.kk)then
rumer(k9)=1
!print*,'rumer(k9)',rumer(k9)
endif
enddo

do k6=1,15000
rumer_rad(k6)=1
enddo
if(nlast.ne.0)then
do k3=1,15000
do k6=1,100
rumer2(k3,k6)=0
rumer1(k3,k6)=0
enddo
enddo
do k3=1,nstr
do k6=1,nae-nl*2
do k7=nl*2+1,nae
do k8=1,k6
if(rumer1(k3,k8).eq.str(k3,k7))goto 390
enddo
if(str(k3,k7).gt.rumer1(k3,k6))then
rumer1(k3,k6)=str(k3,k7)
!print*,rumer1(k6)
endif
390 enddo
enddo
enddo
do k3=1,15000
do k6=1,100
rumer2(k3,k6)=0
enddo
enddo
do k3=1,15000
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
!print*,str(k9,k6),rumer2(k9,k2)
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
!print*,'rrad,rum1,rum2',rrad(k4),rum1,rum2
if((rrad(k4)-rum1)*(rrad(k4)-rum2).gt.0)then
rum=rum+1
endif
enddo
enddo
if(rum.lt.((nae-nlast-nl*2)/2)*nlast)then
rumer_rad(k9)=0
!print*,'rumer_rad(k9)',rumer_rad(k9)
endif
enddo


endif

print*,'exit rumer_structures'
return
end subroutine rumer_structures
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
