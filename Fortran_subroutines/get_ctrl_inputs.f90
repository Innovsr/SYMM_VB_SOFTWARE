!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_ctrl_inputs(geometry_unit, nao, nae, nmul, output_file_name)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!use commondat
!use commondat1
implicit none

integer:: nao, nae, nmul
character(len = 100)::geometry_unit, output_file_name

print*,'nao, nae, nmul',nao, nae, nmul
print*,'geometry_unit',geometry_unit
print*,'file_name',output_file_name

stop
end subroutine get_ctrl_inputs
