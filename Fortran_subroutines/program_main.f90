!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                            Chemical Insight Structure Sets                             !!
!!                                                                                        !!
!!                               Written by: Dr. Sourav Roy                               !!
!!                                                                                        !!
!!                    School Of Pharmacy, Hebrew University of Jerusalem                  !!
!!                                                                                        !!
!!                                    Jerusalem, Israel                                   !!
!!                                                                                        !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program main
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
use commondat1
implicit none

integer :: ierr, rank, size 
integer::wig1,nnae,STDOUT,argnum,nstr7,str5(2000,20)
logical :: fileexists
character(len=35)::inputfilename
common/str/str5,nstr7


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    Output files   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! structures.dat:: stores all possible structures with identifications of rumer structurs                 !
! structure_set_i.dat:: output file, i=1, 2, .... depeding on volume of sets; one file contain 75000 sets !
! Rumer_Sets_all.dat:: unique Rumer sets for all possible permutations of orbitals                        !
! quality_str.dat:: all structures arranged according to their qualities                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
open(unit=7,file='structures.dat',status='unknown')
open(unit=10,file='structure_set_1.dat',status='unknown')
open(unit=23,file='Rumer_Sets_all.dat',status='unknown')
open(unit=35,file='quality_str.dat',status='unknown')

serial=0
input_flg=0
nstr7=0
Rumwrite=0

!!! code takes inputfile (.xmi) as command line arguments
argnum=iargc()
if(argnum.eq.0)then
  write(*,*)'probably you forget to put .xmi file put it in commandline and rerun symmstr code'
  stop
else
  call read_xmi
endif

!flgst: user input for 'str' keyword; 1: covalent (cov), 2: ionic (ion), 3: both (full)
if(flgst.eq.1.or.flgst.eq.2) then
  call cov_struc
endif

!!! ionic str calculation has been stopped for this version

!if(flgst.eq.1.or.flgst.eq.3) then
!call ion_struc
!endif


call close_file

stop
end program main
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
