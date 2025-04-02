#!/bin/bash

# Exit on any error
set -e

# Define the source files as a space-separated list
SOURCE_FILES=(
  commondat1.f90
  commondat.f90
  orb.f90
  mod_rum_rad.f90
  check_mod.f90
  quality.f90
  final_str.f90
  loops.f90
  infosymm.f90
  factorial.f90
  comb.f90
  close_file.f90
  wigner.f90
  sym_check.f90
  TwoSideJacobi.f90
  MatTran.f90
  vector_rep.f90
  rumer_structures.f90
  main_bond.f90
  check_str_bond.f90
  mod_symm_loops.f90
  mod_coordination_val.f90
  nnat_bond_sig.f90
  geocal.f90
  lowercase.f90
  main_bond_str.f90
  Invmat.f90
  ind_matrix.f90
  eq_dst_check.f90
  eq_dstr_set.f90
  write_symm_xmi_1.f90
  write_symm_xmi_new.f90
  qult_str_arrange.f90
  intra_bond_fac.f90
  symm_break_fac.f90
  nnat_bond.f90
  nnat_bond_2.f90
  prio_rad.f90
  mod_quality_factor.f90
  Rum_set_id.f90
  mod_print_rumer.f90
  All_Rumer_set.f90
  write_rumer_xmi.f90
  symm_cal_sig.f90
  symm_cal_pi.f90
  str_select.f90
  cov_struc.f90
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

