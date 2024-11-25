#!/bin/bash

# Define the source files
SOURCE_FILES="
commondat1.f90 
commondat.f90 
factorial.f90 
new_row.f90
close_file.f90 
wigner.f90 
sym_check.f90 
TwoSideJacobi.f90 
MatTran.f90 
Keycheck.f90
MatLDR.f90 
vector_rep.f90 
rumer_structures.f90 
Rumer_set_id.f90
main_bond_cal.f90 
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
intra_bond_factor.f90
symm_break_factor.f90
nnat_bond_cal.f90
nnat_bond_cal_2.f90
prio_rad_str.f90
quality_factor.f90
symmetry_cal_sig.f90
symmetry_cal_pi.f90
str_selection.f90
read_spl_key.f90
ol.f90
cov_struc.f90
read_info.f90
program_main.f90
"

# Create an object file list from the source files
OBJECT_FILES=$(echo $SOURCE_FILES | sed 's/\.f90/\.o/g')

# Compile the source files into object files
#gfortran -ffpe-trap=underflow,denormal -g -fbacktrace -c $SOURCE_FILES
gfortran -ffpe-summary=none -c $SOURCE_FILES
#gfortran -g -fcheck=all -c $SOURCE_FILES
#mpif90.openmpi -g -fcheck=all -c $SOURCE_FILES

# Link the object files into an executable
#gfortran -ffpe-trap=underflow,denormal -g -fbacktrace -o symm_str $OBJECT_FILES
gfortran -ffpe-summary=none -o symm_str $OBJECT_FILES
#gfortran -g -fcheck=all -o symm_str $OBJECT_FILES
#mpif90.openmpi -g -fcheck=all -o symm_str $OBJECT_FILES

# Clean up object files if no longer needed
rm -f $OBJECT_FILES

echo "Compilation complete. Executable: symm_str"

