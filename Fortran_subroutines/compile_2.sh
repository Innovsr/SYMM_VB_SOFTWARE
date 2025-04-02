#!/bin/bash

# Exit on any error
set -e

# Define the source files as a space-separated list
SOURCE_FILES=(
  commondat1_mod.f90
  commondat_mod.f90
  final_str_mod.f90
  orb_mod.f90
  rum_rad_mod.f90
  check_bd_mod.f90
  quality_mod.f90
  loops_mod.f90
  infosymm_mod.f90
  factorial.f90
  comb.f90
  #close_file.f90
  wigner.f90
  #sym_check.f90
  #TwoSideJacobi.f90
  #MatTran.f90
  #mod_MatLDR.f90
  mod_vector_rep.f90
  mod_rumer_structures.f90
  mod_intra_bond_factor.f90
  mod_symm_break_factor.f90
  mod_nnat_bond_cal.f90
  mod_nnat_bond_cal_2.f90
  mod_prio_rad_str.f90
  mod_main_bond_cal.f90
  mod_quality_factor.f90
  mod_Rumer_set_id.f90
  mod_check_str_bond.f90
  mod_symm_loops.f90
  mod_coordination_val.f90
  mod_nnat_bond_sig.f90
  geocal.f90
  #lowercase.f90
  #main_bond_str.f90
  Invmat.f90
  mod_mat_ind.f90
  mod_eq_dst_check.f90
  mod_eq_dstr_set.f90
  mod_print_rumer.f90
  mod_All_Rumer_set.f90
  mod_write_rumer_xmi.f90
  mod_write_symm_xmi_1.f90
  mod_write_symm_xmi_new.f90
  mod_qult_str_arrange.f90
  #intra_bond_fac.f90
  #symm_break_fac.f90
  #nnat_bd_cal.f90
  #nnat_bond_2.f90
  #prio_radical.f90
  #mod_quality_factor.f90
  mod_symmetry_cal_sig.f90
  mod_symmetry_cal_pi.f90
  mod_str_selection.f90
  mod_cov_struc.f90
  get_ctrl_inputs.f90
  program_main.f90
)

# Define the output executable name
EXECUTABLE="symm_vb_exec"

echo "Compiling Fortran source files into an executable..."

# Compile using gfortran with debugging and bounds checking
#gfortran -o "$EXECUTABLE" "${SOURCE_FILES[@]}" 
#gfortran -g -Wall -fcheck=all -fbounds-check -fbacktrace "${SOURCE_FILES[@]}" -o "$EXECUTABLE"
#gfortran -g -Wall -fcheck=all -fbounds-check -fbacktrace -O0  "${SOURCE_FILES[@]}" -o "$EXECUTABLE"
gfortran -g -Wall -fcheck=all -fbounds-check -fbacktrace -O0 -finit-real=nan -finit-integer=-9999 -ffpe-trap=zero,overflow,invalid "${SOURCE_FILES[@]}" -o "$EXECUTABLE"

echo "Compilation successful. Executable: $EXECUTABLE"

