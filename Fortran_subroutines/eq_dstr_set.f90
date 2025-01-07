!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! creating equally distributed structure sets : the sets contains the structures !
!! bearing bonds wich cover space of bods symmetrycally.                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine eq_dstr_set(tns,lnp,npstr,str2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

integer::i,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18,i19,&
l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14,l15,l16,l17,l18,l19,l20,j
integer::tns,lnp,npstr,fail,strsln(1000),str2(15000,20),num_str,numpstr
common/nst/num_str,numpstr

numpstr=npstr
num_str=0

print*,'enter eq_dstr_set'

! loops started here
do i1=1,tns
   do i=1,200
      strsln(i)=0
   enddo
   j=0
   j=j+1
   strsln(j)=i1

   l1=j

   ! loop 2
   do i2=i1+1,tns
      do i=l1+1,j
         strsln(i)=0
      enddo
      j=l1
      j=j+1
      strsln(j)=i2
      call eq_dst_check(j,strsln,fail,lnp,str2)
      if (fail.eq.1) then
         cycle
      endif

      l2=j

      ! loop 3
      do i3=i2+1,tns
         do i=l2+1,j
            strsln(i)=0
         enddo
         j=l2
         j=j+1
         strsln(j)=i3
         call eq_dst_check(j,strsln,fail,lnp,str2)
         if (fail.eq.1) then
            cycle
         endif

         l3=j

         ! loop 4
         do i4=i3+1,tns
            do i=l3+1,j
               strsln(i)=0
            enddo
            j=l3
            j=j+1
            strsln(j)=i4
            call eq_dst_check(j,strsln,fail,lnp,str2)
            if (fail.eq.1) then
               cycle
            endif

            l4=j

            ! loop 5
            do i5=i4+1,tns
               do i=l4+1,j
                  strsln(i)=0
               enddo
               j=l4
               j=j+1
               strsln(j)=i5 
               call eq_dst_check(j,strsln,fail,lnp,str2)
               if (fail.eq.1) then
                  cycle
               endif

               l5=j

               ! loop 6
               do i6=i5+1,tns
                  do i=l5+1,j
                     strsln(i)=0
                  enddo
                  j=l5
                  j=j+1
                  strsln(j)=i6
                  call eq_dst_check(j,strsln,fail,lnp,str2)
                  if (fail.eq.1) then
                     cycle
                  endif

                  l6=j

                  ! loop 7
                  do i7=i6+1,tns
                     do i=l6+1,j
                        strsln(i)=0
                     enddo
                     j=l6
                     j=j+1
                     strsln(j)=i7
                     call eq_dst_check(j,strsln,fail,lnp,str2)
                     if (fail.eq.1) then
                        cycle
                     endif

                     l7=j

                     ! loop 8
                     do i8=i7+1,tns
                        do i=l7+1,j
                           strsln(i)=0
                        enddo
                        j=l7
                        j=j+1
                        strsln(j)=i8
                        call eq_dst_check(j,strsln,fail,lnp,str2)
                        if (fail.eq.1) then
                           cycle
                        endif

                        l8=j

                        ! loop 9
                        do i9=i8+1,tns
                           do i=l8+1,j
                              strsln(i)=0
                           enddo
                           j=l8
                           j=j+1
                           strsln(j)=i9
                           call eq_dst_check(j,strsln,fail,lnp,str2)
                           if (fail.eq.1) then
                              cycle
                           endif

                           l9=j

                           do i10=i9+1,tns
                              do i=l9+1,j
                                 strsln(i)=0
                              enddo
                              j=l9
                              j=j+1
                              strsln(j)=i10
                              call eq_dst_check(j,strsln,fail,lnp,str2)
                              if (fail.eq.1) then
                                 cycle
                              endif

                              l10=j

                              do i11=i10+1,tns
                                 do i=l10+1,j
                                    strsln(i)=0
                                 enddo
                                 j=l10
                                 j=j+1
                                 strsln(j)=i11
                                 call eq_dst_check(j,strsln,fail,lnp,str2)
                                 if (fail.eq.1) then
                                    cycle
                                 endif

                                 l11=j

                                 do i12=i11+1,tns
                                    do i=l11+1,j
                                       strsln(i)=0
                                    enddo
                                    j=l11
                                    j=j+1
                                    strsln(j)=i12
                                    call eq_dst_check(j,strsln,fail,lnp,str2)
                                    if (fail.eq.1) then
                                       cycle
                                    endif

                                    l12=j

                                    do i13=i12+1,tns
                                       do i=l12+1,j
                                          strsln(i)=0
                                       enddo
                                       j=l12
                                       j=j+1
                                       strsln(j)=i13
                                       call eq_dst_check(j,strsln,fail,lnp,str2)
                                       if (fail.eq.1) then
                                          cycle
                                       endif

                                       l13=j

                                       do i14=i13+1,tns
                                          do i=l13+1,j
                                             strsln(i)=0
                                          enddo
                                          j=l13
                                          j=j+1
                                          strsln(j)=i14
                                          call eq_dst_check(j,strsln,fail,lnp,str2)
                                          if (fail.eq.1) then
                                             cycle
                                          endif

                                          l14=j

                                          do i15=i14+1,tns
                                             do i=l14+1,j
                                                strsln(i)=0
                                             enddo
                                             j=l14
                                             j=j+1
                                             strsln(j)=i15
                                             call eq_dst_check(j,strsln,fail,lnp,str2)
                                             if (fail.eq.1) then
                                                cycle
                                             endif

                                             l15=j

                                             do i16=i15+1,tns
                                                do i=l15+1,j
                                                   strsln(i)=0
                                                enddo
                                                j=l15
                                                j=j+1
                                                strsln(j)=i16
                                                call eq_dst_check(j,strsln,fail,lnp,str2)
                                                if (fail.eq.1) then
                                                   cycle
                                                endif

                                                l16=j

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


print*,'exit eq_dstr_set'
return
end subroutine eq_dstr_set
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
