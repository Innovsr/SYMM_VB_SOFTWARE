!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module mod_print_rumer
use commondat_mod
use mod_Rumer_set_id
use orb_mod
use quality_mod
use final_str_mod
implicit none

contains
subroutine print_rumer(permutation,n,j,setno, rstr, str)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! All the Rumer sets are written in the file "Rumer_Sets.dat" varified
!! with the subroutine "Rumer_set_id" and then written in the file
!! "Rumer_Sets_all.dat" with full format of the output.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer::n, i, i1, i2, j, k, setno, m, m1, m20, m19, nlonep
integer::jj,Rid, count, count1
integer, allocatable, intent(in)::permutation(:)
integer, pointer::rstr1(:,:),rstr(:,:),str(:,:)
integer, allocatable::col(:)!, symq(:)
logical::symmetry_group
!double precision, pointer::symsc(:)

print*,'enter print_rumer'
!do i=1,setno
!print*,'symm_groups_count',symm_groups_count(i)
!enddo

nlonep = nlpr
allocate(rstr1(MaxStrOepo, nae))
rstr1 = 0

!print*,'enter rumer',(permutation(i),i =1, n),j
Rumwrite=1
do i=1,setno
  do i1=1,nlonep*2
    rstr1(i,i1)=rstr(i,i1)
  enddo
enddo

do i=1,setno
  do i1=1,nae
    do i2=1,n
      if(rstr(i,i1).eq.orbs1(i2))then
        rstr1(i,i1)=permutation(i2)
        
      endif
    enddo
  enddo
enddo

!do i =1,setno
!print*,'rumer_set_rstr1',(rstr1(i,j),j=1,nae)
!enddo

if(j.eq.1)then
  rum_jj=0
  goto 103
endif
call Rumer_set_id(rstr1,setno,nlonep,Rid, str) 

if(Rid.eq.1) then
  return
endif



103 rum_jj=rum_jj+1

write(31,*)'set number',rum_jj,(permutation(i),i=1,n),setno
write(10,*)'permutation',rum_jj,'>',(permutation(i),i=1,n)
write(10,*)'*****************************************************'
write(23,*)
do m19=1,setno
write(31,914)(rstr1(m19,m20),m20=1,nae)
if(nfset.eq.5)write(10,915)1,':',niao,(rstr1(m19,m20),m20=1,nae)
if(niao.gt.1)then
endif
enddo
if(nfset.eq.5)write(10,*)'set number',rum_jj

!total_str=MaxStrOepo
rumset=rumset+1   ! counter of Rumer sets
call str_sl_check(nlonep, rstr1, setno, rumset, col)

!print*,'i am sourav'
!if(symtype.eq.'check')then
!  if(sig_sym_flg.eq.1)call symmetry_cal_sig(nlonep, f_str, strplp, symsc, symq, nssym)
!  if(sig_sym_flg.ne.1)call symmetry_cal_pi(nlonep, f_str, strplp, symsc, symq, nssym)
!endif
!print*,'symq',sig_sym_flg,'|',symq
!stop

if (symm.eq.1) then
  print*,'i am sourav1', col, symm_maxval, setno
  count1 = 0
  symmetry_group = .false.
  do j = 1, symm_maxval
     print*,'jj',j,symm_maxval
     count = 0
     loop2:do i = 1, setno
        print*,'setno',i,symm_groups_count(j)
        do k = 1, symm_groups_count(j)
           print*,'kkk',k,symm_groups_count(j),col(i),symm_groups(j,k)
           if(col(i).eq.symm_groups(j,k)) then
              count = count + 1
              print*,'i am sourav2'
              cycle loop2
           endif
        enddo
     enddo loop2
     print*,'count',count, symm_groups_count(j)
     if (count.eq.symm_groups_count(j)) then
        count1=count1+count
        print*,symmetry_group
     endif
     if (count1.eq.setno) then
        symmetry_group = .true.
        exit
     endif
  enddo 

  if( symmetry_group ) then
    symm_rumset = symm_rumset + 1 !counter of symmetric Rumer sets
    write(5,913)'Set_number',symm_rumset,(col(i),i=1,setno)
  endif
endif
913 format(a,2x,I0,4x,*(I0, 1x))

print*,'i am sourav3'
!deallocate(symm_groups) 
do m=1,setno

  if(niao.eq.0)then
    write(10,900)str_quality_1(col(m)),bondq(col(m)),str_quality_2(col(m)),quality_fac(col(m)),'|',(f_str(col(m),m1),m1=1,nae)
  endif
  if(niao.gt.1)then
    write(10,901)str_quality_1(col(m)),bondq(col(m)),str_quality_2(col(m)),quality_fac(col(m)), &
            '|',1,':',niao,(f_str(col(m),m1),m1=1,nae)
  endif
  if(niao.eq.1)then
    write(10,909)str_quality_1(m),bondq(m),str_quality_2(m),quality_fac(m),'|',1,1,(rstr1(m,m1),m1=1,nae)
  endif
enddo
if (symm.eq.1) then
   if( symmetry_group ) then
      write(10,*)'it is a symmetric set'
      write(*,*)'it is a symmetric set'
   else
      write(10,*)'it is a non-symmetric set'
      write(*,*)'it is a non-symmetric set'
   endif
endif        


write(10,*)'Set_number=',rumset
write(10,*)


900 format(3x,I0,3x,I0,3x,I0,3x,I0,1x,a,1x,*(I0, 1x))
901 format(3x,I0,3x,I0,3x,I0,3x,I0,1x,a,1x,I0,a,I0,1x,*(I0, 1x))
909 format(3x,I0,3x,I0,3x,I0,3x,I0,1x,a,1x,I0,I0,1x,*(I0, 1x))

!:900 format(a,I3,10I5)

914 format(x,25I4)
915 format(x,I1,a,I3,x,25I4)
!916 format(x,I3,I3,x,25I4)
!deallocate(symm_groups_count)
!deallocate(symm_groups)
print*,'exit print_rumer'
totrum=rum_jj
!print*,'rumer totrum',totrum
return
end subroutine print_rumer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine str_sl_check(nlonep, str2, nstr, rumset, col)
integer:: i, j, k, l, nlonep, nstr, match_count, rumset
integer, pointer::str2(:,:)
logical :: row_matched
integer, allocatable::col(:)

print*,'enter str_sl_check'!, size(f_str), strplp

!do i=1,nstr
!print*,'str2',i,(str2(i,j),j=1,nae)
!enddo

if (.not. allocated(col))then
  allocate(col(CovDim))
endif
col = 0

match_count = 0

loop1:do i = 1, nstr
!print*,'str2',i,(str2(i,j),j=1,nae)
    !row_matched = .false.

    loop2:do j = 1, strplp
        row_matched = .true.
        loop3:do k = nlonep*2+1, nae-nlast, 2  ! Step by 2 (pairwise check)
            do l = nlonep*2+1, nae-nlast, 2
               if (str2(i, k) == f_str(j, l) .and. str2(i, k+1) == f_str(j, l+1)) then
                      cycle loop3
               end if
               if (str2(i, k) == f_str(j, l+1) .and. str2(i, k+1) == f_str(j, l)) then
                      cycle loop3
               end if
            end do 
            row_matched = .false.
            cycle loop2
        end do loop3

        if (row_matched) then
            match_count = match_count + 1
            col(match_count) = j
            cycle loop1
        end if
    end do loop2

    !if (.not. row_matched) then
    !    write(5,*), "Row ", j, " in array1 has no pairwise match in array2!"
    !    !stop
    !end if
end do loop1

if (symm.ne.1) write(5,913)'Set_number',rumset,(col(i),i=1,match_count)
print*,'Set_number',rumset
913 format(a,2x,I0,4x,*(I0, 1x))

print*,'exit str_sl_check'
!if (rumset.eq.2)stop
return
end subroutine str_sl_check

end module mod_print_rumer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
