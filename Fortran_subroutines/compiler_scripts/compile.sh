#!/bin/bash

# Exit on any error
set -e

# Define the source files
SOURCE_FILES="
commondat1.f90
commondat.f90
quality.f90
factorial.f90
comb.f90
new_row.f90
close_file.f90
wigner.f90
sym_check.f90
TwoSideJacobi.f90
MatTran.f90
MatLDR.f90
vector_rep.f90
rumer_structures.f90
Rumer_set_id.f90
main_bond.f90
check_str_bond.f90
symm_loops.f90
coordination_val.f90
nnat_bond_sig.f90
geocal.f90
lowercase.f90
main_bond_str.f90
Invmat.f90
mat_ind.f90
eq_dst_check.f90
eq_dstr_set.f90
rumer.f90
All_Rumer_set.f90
write_rumer_xmi.f90
write_symm_xmi_1.f90
write_symm_xmi_new.f90
qult_str_arrange.f90
intra_bond_fac.f90
symm_break_fac.f90
nnat_bond.f90
nnat_bond_2.f90
prio_rad.f90
quality_fac.f90
symm_cal_sig.f90
symm_cal_pi.f90
str_select.f90
cov_struc.f90
get_ctrl_inputs.f90
program_main.f90
"

# Define the output executable name
EXECUTABLE="symm_vb_exec"

echo "Compiling Fortran source files into an executable..."

# Compile using gfortran
#gfortran -o "$EXECUTABLE" $SOURCE_FILES
gfortran -g -Wall -fcheck=all -fbounds-check -fbacktrace $SOURCE_FILES -o $EXECUTABLE

# Check if the executable was created successfully
if [[ ! -f "$EXECUTABLE" ]]; then
  echo "Error: Compilation failed. Executable file not created."
  exit 1
fi
